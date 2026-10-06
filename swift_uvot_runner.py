"""
swift_uvot_runner.py

The one way the Swift UVOT pipeline runs a HEASoft tool. Every step uses it.

Each call gets:

  * Its own working folder (a fresh temporary directory). uvotsource and
    its sub-tasks write temporary files into the current folder
    (this.<pid>.<n>.fits, src.<pid>.<n>.reg, ...), and a crashed call
    leaves them behind. Separate folders keep parallel calls apart and
    keep a failed call's leftovers in one place.
  * Short names for its inputs. Each input file is linked into the working
    folder under a short name, and the tool is given that name. HEASoft
    tools have fixed-size path buffers (in the XRT pipeline xselect and
    grppha failed beyond ~90-150 characters); with links, the paths a tool
    sees are short wherever the data live.
  * A private parameter-file folder ($PFILES). HEASoft tasks keep their
    parameters in .par files, and parallel calls sharing $HOME/pfiles can
    collide (the XRT pipeline's parallel-extraction "Bug B"). uvotsource
    wrote no .par files in our tests, but its sub-tasks and other tools do.
  * No terminal prompts: HEADASNOQUERY=1 and HEADASPROMPT=/dev/null.
    Without a controlling terminal (nohup, cron, a dropped SSH session, a
    script started by another program) a HEASoft task otherwise stops with

        Task uvotsource 0.0 terminating with status 6
        Unable to redirect prompts to the /dev/tty (at headas_stdio.c:152)

    and, because uvotsource is a Perl script, still exits with status 0.
  * stdin from /dev/null, a wall-clock timeout, and a new process group,
    so a stuck call is killed together with everything it started.
  * A log, if asked for: the command, working folder, exit status,
    elapsed time and the complete output.

A call succeeds only if the tool exits 0, prints no "terminating with
status N" line with N other than 0, and leaves every expected output file,
non-empty. Its outputs are then moved to their destinations and the
working folder is removed. A failed call keeps its working folder (named in
the result and the log) unless keep_failed=False.

Usage:

    from swift_uvot_runner import run_tool, run_many

    result = run_tool(
        'uvotsource',
        {'image': 'sky.img.gz[3]', 'srcreg': 'src.reg', 'bkgreg': 'bkg.reg',
         'expfile': 'ex.img.gz[3]', 'sigma': 3, 'apercorr': 'NONE',
         'history': 'no', 'outfile': 'phot.fits', 'clobber': 'yes'},
        inputs={'sky.img.gz': sky_path, 'ex.img.gz': ex_path,
                'src.reg': src_path, 'bkg.reg': bkg_path},
        outputs={'phot.fits': dest_path},
        log=log_path)
    if not result.ok:
        print(result.reason)
        print(result.tail)

    results = run_many([dict(tool=..., params=..., ...), ...], nproc=16)

Positional arguments (e.g. for quzcif) can be given as a list of strings,
which are passed to the tool verbatim.
"""

import os
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

# Environment variable naming the folder for working folders (default: the
# system temporary directory). Keep it short; see the module docstring.
WORK_ROOT_ENV = 'SWIFT_UVOT_TMP'

DEFAULT_TIMEOUT = 300

# A HEASoft task announces its own failure with this line, sometimes while
# its wrapper still exits 0.
_TERMINATING = re.compile(r'terminating with status (-?\d+)')

# Lines worth showing when a call fails.
_ERRORISH = re.compile(
    r'error|fail|terminating|unable|cannot|not found|no such|overflow|'
    r'abort|warning|invalid', re.I)


@dataclass
class CallResult:
    """What happened in one HEASoft call."""
    tool: str
    ok: bool
    reason: str = ''          # why it failed ('' when ok)
    returncode: int = None
    timed_out: bool = False
    elapsed: float = 0.0
    command: str = ''
    stdout: str = ''
    stderr: str = ''
    log: str = None           # log file, if one was written
    workdir: str = None       # kept working folder (failed calls only)
    outputs: dict = field(default_factory=dict)   # name -> destination

    @property
    def tail(self):
        """The informative last lines of the output, for summaries."""
        return error_tail(self.stdout + '\n' + self.stderr)


def error_tail(text, n=6):
    """
    The last n lines of text that look like errors or warnings, or the
    last n non-blank lines if none do.
    """
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    errs = [ln for ln in lines if _ERRORISH.search(ln)]
    return '\n'.join((errs or lines)[-n:])


def _format_params(params):
    """argv strings for a tool: dict -> key=value, list -> verbatim."""
    if params is None:
        return []
    if isinstance(params, dict):
        argv = []
        for key, value in params.items():
            if isinstance(value, bool):
                value = 'yes' if value else 'no'
            argv.append('%s=%s' % (key, value))
        return argv
    return [str(p) for p in params]


def _signal_name(returncode):
    try:
        return signal.Signals(-returncode).name
    except (ValueError, TypeError):
        return 'signal %d' % -returncode


def _write_log(path, result, workdir):
    with open(path, 'w') as f:
        f.write('# %s\n' % result.tool)
        f.write('# Command: %s\n' % result.command)
        f.write('# Working folder: %s%s\n' % (
            workdir, '' if result.workdir else ' (removed)'))
        f.write('# Result: %s\n' % ('OK' if result.ok
                                    else 'FAILED (%s)' % result.reason))
        f.write('# Return code: %s\n' % result.returncode)
        f.write('# Elapsed: %.2f s\n\n' % result.elapsed)
        if result.stdout:
            f.write('=== STDOUT ===\n%s\n' % result.stdout.rstrip())
        if result.stderr:
            f.write('\n=== STDERR ===\n%s\n' % result.stderr.rstrip())


def run_tool(tool, params=None, inputs=None, outputs=None, dest_dir='.',
             log=None, timeout=DEFAULT_TIMEOUT, work_root=None,
             keep_failed=True, env_extra=None):
    """
    Run one HEASoft tool in isolation; see the module docstring.

    tool        task name, e.g. 'uvotsource'
    params      dict (key=value) or list (verbatim) of arguments; refer to
                inputs by their short names
    inputs      {short name: path} of existing files to link into the
                working folder
    outputs     files the tool must create in its working folder: a list
                of names (moved to dest_dir/<name>) or {name: destination}
    log         path for the call's log (optional)
    timeout     seconds before the call and everything it started are
                killed
    work_root   folder for the working folder (default: $SWIFT_UVOT_TMP or
                the system temporary directory)
    keep_failed keep a failed call's working folder for inspection
    env_extra   extra environment variables for this call
    """
    inputs = inputs or {}
    if isinstance(outputs, dict):
        out_map = dict(outputs)
    else:
        out_map = {name: os.path.join(dest_dir, name)
                   for name in (outputs or [])}
    argv = [tool] + _format_params(params)
    result = CallResult(tool=tool, ok=False, command=' '.join(argv), log=log)

    missing = [p for p in inputs.values() if not os.path.exists(p)]
    if missing:
        result.reason = 'missing input %s' % ', '.join(missing)
        if log:
            _write_log(log, result, '(not created)')
        return result
    if shutil.which(tool) is None:
        result.reason = '%s not on PATH (is HEASoft set up?)' % tool
        if log:
            _write_log(log, result, '(not created)')
        return result

    root = work_root or os.environ.get(WORK_ROOT_ENV) or tempfile.gettempdir()
    os.makedirs(root, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix='swuvot_%s_' % tool, dir=root)
    pfiles = os.path.join(workdir, 'pfiles')
    os.mkdir(pfiles)
    for name, path in inputs.items():
        os.symlink(os.path.abspath(path), os.path.join(workdir, name))

    env = os.environ.copy()
    headas = env.get('HEADAS', '')
    env['PFILES'] = ('%s;%s/syspfiles' % (pfiles, headas)) if headas else pfiles
    env['HEADASNOQUERY'] = '1'
    env['HEADASPROMPT'] = '/dev/null'
    if env_extra:
        env.update(env_extra)

    start = time.monotonic()
    proc = subprocess.Popen(argv, cwd=workdir, env=env,
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            errors='replace', start_new_session=True)
    try:
        result.stdout, result.stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        result.timed_out = True
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        result.stdout, result.stderr = proc.communicate()
    result.elapsed = time.monotonic() - start
    result.returncode = proc.returncode

    output = result.stdout + '\n' + result.stderr
    statuses = [int(s) for s in _TERMINATING.findall(output)]
    if result.timed_out:
        result.reason = 'timeout after %d s' % timeout
    elif proc.returncode < 0:
        result.reason = 'killed by %s' % _signal_name(proc.returncode)
    elif proc.returncode != 0:
        result.reason = 'exit status %d' % proc.returncode
    elif any(s != 0 for s in statuses):
        result.reason = ('task reported "terminating with status %d" '
                         '(exit status was 0)'
                         % next(s for s in statuses if s != 0))
    else:
        absent = [name for name in out_map
                  if not os.path.isfile(os.path.join(workdir, name))
                  or os.path.getsize(os.path.join(workdir, name)) == 0]
        if absent:
            result.reason = 'no output %s' % ', '.join(absent)
    if 'buffer overflow detected' in output:
        note = ('"buffer overflow detected": a path the tool wrote was too '
                'long (is $CALDB short?)')
        result.reason = (result.reason + ' -- ' + note) if result.reason \
            else note

    if not result.reason:
        for name, dest in out_map.items():
            dest_parent = os.path.dirname(os.path.abspath(dest))
            os.makedirs(dest_parent, exist_ok=True)
            shutil.move(os.path.join(workdir, name), dest)
            result.outputs[name] = dest
        result.ok = True
        shutil.rmtree(workdir, ignore_errors=True)
    elif keep_failed:
        result.workdir = workdir
    else:
        shutil.rmtree(workdir, ignore_errors=True)

    if log:
        os.makedirs(os.path.dirname(os.path.abspath(log)), exist_ok=True)
        _write_log(log, result, workdir)
    return result


def run_many(calls, nproc=1, progress=None):
    """
    Run many calls, nproc at a time. calls is a list of dicts of run_tool
    keyword arguments (including 'tool'). Returns the CallResults in the
    same order. Each call is isolated as in run_tool, so the order in which
    they finish does not matter. progress, if given, is called as
    progress(n_done, n_total) after each call.
    """
    done = [0]
    lock = threading.Lock()

    def one(call):
        res = run_tool(**call)
        if progress:
            with lock:
                done[0] += 1
                progress(done[0], len(calls))
        return res

    if nproc <= 1:
        return [one(call) for call in calls]
    with ThreadPoolExecutor(max_workers=nproc) as pool:
        return list(pool.map(one, calls))
