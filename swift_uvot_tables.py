"""
swift_uvot_tables.py

The plain-text tables the Swift UVOT pipeline steps write and read
(inventory, photometry, master table). Same style as the XRT pipeline's
master tables: '#' comment lines, one header line, then whitespace-aligned
columns, so they are easy to read in a terminal and to diff. A value
containing spaces (a reason or comment) is written in double quotes.
Missing values are written as '-'.

    write_table(path, rows, columns, comments=[...])
    rows = read_table(path)          # list of dicts of strings
"""

import os
import shlex

MISSING = '-'


def _cell(value):
    if value is None or value == '':
        return MISSING
    text = str(value)
    if any(c.isspace() for c in text) or text.startswith('"'):
        return '"%s"' % text.replace('"', "'")
    return text


def write_table(path, rows, columns, comments=()):
    """
    Write rows (dicts) as an aligned text table with the given columns,
    preceded by '# ' comment lines. Written to a temporary name first and
    renamed, so a reader never sees half a table.
    """
    cells = [[_cell(row.get(col)) for col in columns] for row in rows]
    widths = [max([len(col)] + [len(r[i]) for r in cells])
              for i, col in enumerate(columns)]
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        for line in comments:
            fh.write(('# ' + line).rstrip() + '\n')
        fh.write('  '.join(col.ljust(w) for col, w in zip(columns, widths))
                 .rstrip() + '\n')
        for r in cells:
            fh.write('  '.join(c.ljust(w) for c, w in zip(r, widths))
                     .rstrip() + '\n')
    os.replace(tmp, path)


def read_table(path):
    """Rows of a table written by write_table, as dicts of strings."""
    rows, columns = [], None
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = shlex.split(line, posix=True)
            if columns is None:
                columns = fields
                continue
            if len(fields) != len(columns):
                raise ValueError('%s:%d: %d values for %d columns'
                                 % (path, lineno, len(fields), len(columns)))
            rows.append(dict(zip(columns, fields)))
    return rows
