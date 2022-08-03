

def allclose(a,b,tol=None):
    """emulate numpy.allclose: recursively compare a and b on given tolerance"""
    tol = tol if tol else 0.000001
    if isinstance(a, (int,float)) and isinstance(b, (int,float)):
        if abs(a-b) > tol:
            return False
        else:
            return True
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        if not a:
            return True
        for i,vi in enumerate(a):
            vj = b[i]
            if not allclose(vi,vj,tol):
                return False
        return True
    else:
        return False



