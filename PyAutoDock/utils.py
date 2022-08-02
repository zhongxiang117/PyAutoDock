import os


def file_gen_new(filename,noseq=None):
    """return a non-exist filename to avoid file overwriting"""
    if noseq:
        if not os.path.isfile(filename):
            return filename

    base,ext = os.path.splitext(filename)
    i = 1
    while True:
        new = '{:}-{:}{:}'.format(base,1,ext)
        if not os.path.isfile(new):
            break
        i += 1
    return new


