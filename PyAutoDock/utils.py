import os
import math


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


def calc_ddd_Mehler_Solmajer(distance,tol=None):
    """Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910."""
    tol = tol if tol else 0.00000001
    epsilon = 1.0
    mlambda = 0.003627
    epsilon0 = 78.4
    A = -8.5525
    B = epsilon0 - A
    rk = 7.7839
    lambada_B = -mlambda * B
    epsilon = A + B/(1.0+rk*math.exp(lambada_B*distance))
    if epsilon < tol: epsilon = 1.0
    return epsilon


def read_map_data(file,exit_on_error=None):
    """read map data in format: List(x(y(z(num_maps))))"""
    xmax,ymax,zmax = 0,0,0
    data = []
    with open(file,'rt') as f:
        for line in f:
            line = line.strip()
            if not len(line) or line[0] == '#': continue
            ltmp = line.replace(',',' ').split()
            try:
                x = int(ltmp[0])
                y = int(ltmp[1])
                z = int(ltmp[2])
                d = list(map(float,ltmp[3:]))
            except (ValueError,IndexError):
                if exit_on_error: break
                continue
            else:
                xmax = max(x,xmax)
                ymax = max(y,ymax)
                zmax = max(z,zmax)
                data.append(d)
    if (xmax+1)*(ymax+1)*(zmax+1) != len(data):
        return []
    dx = (ymax+1) * (zmax+1)
    dy = zmax + 1
    return [
        [[data[i*dx+j*dy+k] for k in range(zmax+1)] for j in range(ymax+1)]
        for i in range(xmax+1)
    ]


def write_map_data(maps,file=None):
    """write map data in format: `x,y,z,d1,d2, ...' """
    if not file: file = 'maps.txt'
    with open(file,'wt') as f:
        for x in range(len(maps[0])):
            for y in range(len(maps[1])):
                for z in range(len(maps[2])):
                    f.write('{:},{:},{:},{:}\n'.format(x,y,z,
                        ','.join([str(t) for t in maps[x][y][z]]))
                    )


def write_table_vdW_hbond_data(table,file=None):
    """write vdW hydrogen bond table data in format: `point,receptor_atomtype,ligand_atomtypes, ...' """
    if not file: file = 'vdW-hbond-table.txt'
    with open(file,'wt') as f:
        for i in range(len(table)):
            for j in range(len(table[i])):
                    f.write('{:},{:},{:}\n'.format(i,j,
                        ','.join(['{:}'.format(t) for t in table[i][j]]))
                    )


def read_table_vdW_hbond_data(file,exit_on_error=None):
    """read map data in format: List(num_points(num_receptor_atomtypes(num_ligand_atomtypes)))"""
    num_points,num_receptor_atomtypes = 0,0
    data = []
    with open(file,'rt') as f:
        for line in f:
            line = line.strip()
            if not len(line) or line[0] == '#': continue
            ltmp = line.replace(',',' ').split()
            try:
                x = int(ltmp[0])
                y = int(ltmp[1])
                d = list(map(float,ltmp[2:]))
            except (ValueError,IndexError):
                if exit_on_error: break
                continue
            else:
                num_points = max(x,num_points)
                num_receptor_atomtypes = max(y,num_receptor_atomtypes)
                data.append(d)
    if (num_points+1)*(num_receptor_atomtypes+1) < len(data):
        return []
    dx = num_receptor_atomtypes + 1
    return [[data[i*dx+j] for j in range(num_receptor_atomtypes+1)] for i in range(num_points+1)]



