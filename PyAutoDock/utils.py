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


