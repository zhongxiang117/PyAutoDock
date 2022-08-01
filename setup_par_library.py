import logger

import collections
import os

logger = logger.mylogger()

AD_MAP_COEFFS = collections.namedtuple(
    'AD_MAP_COFEES',(
        'FE_coeff_vdW FE_coeff_hbond FE_coeff_estat FE_coeff_desolv FE_coeff_tors'
    ),
    defaults=[-1.0 for i in range(5)]
)

AD_MAP_ATOM = collections.namedtuple(
    'AD_MAP_ATOM',(
        'atomtype Rii epsii vol solpar Rij_hb epsij_hb hbond '  # space is important
        'rec_index map_index bond_index'
    ),
    defaults=[-1.0 for i in range(11)]
)

# for hydrogen bonds
NON = 0     # none
DS = 1      # spherical donor
D1 = 2      # directional donor
AS = 3      # spherical acceptor
A1 = 4      # acceptor of 1 directional hbond
A2 = 5      # acceptor of 2 directional hbond


class SetupParLibrary:
    """Setup parameter library, onetime setup, forever using

    Args:
        filename(str):  input parameter file

    Attributes:
        e:  coefficients

    Method:
        get_atom_par:   return nametuple of atom parameters
    """
    def __init__(self,filename=None,*args,**kwargs):
        self.filename = filename
        self.atoms = []
        self._read(self.filename)
        self._atomtypes = [i.atomtype.lower() for i in self.atoms]

    def get_atom_par(self,atomtype):
        if atomtype.lower() in self._atomtypes:
            return self.atoms[self._atomtypes.index(atomtype.lower())]
        logger.warning('not defined: atomtype: {:}'.format(atomtype))
        return None

    def _read(self,filename=None):
        if not filename: filename = self.filename
        if not os.path.isfile(filename): return
        fileds_coeff = [i.lower() for i in AD_MAP_COEFFS._fields]
        coeff_dict = {}
        pars = []
        with open(filename,'rt') as f:
            for line in f:
                new = self._strip_comments(line)
                if not len(new): continue
                ltmp = new.split()
                if len(ltmp) == 2:
                    if ltmp[0].lower() in fileds_coeff:
                        key = AD_MAP_COEFFS._fields[fileds_coeff.index(ltmp[0].lower())]
                        coeff_dict[key] = float(ltmp[1])
                elif len(ltmp) == 12 and ltmp[0].lower() == 'atom_par':
                    pars.append([ltmp[1],*list(map(float,ltmp[2:8])),*list(map(int,ltmp[8:]))])
                else:
                    logger.warning('ignoring: line: {:}'.format(line))
        self.e = AD_MAP_COEFFS(**coeff_dict)

        # check
        for i in self.e._fields:
            if getattr(self.e,i) < 0.0:
                logger.critical('not found: {:}'.format(i))

        # process
        self.atoms = []
        for a in pars:
            a[2] *= self.e.FE_coeff_vdW
            a[6] *= self.e.FE_coeff_hbond
            self.atoms.append(AD_MAP_ATOM(*a))

    def _strip_comments(self,line):
        ndx = line.find('#')
        if ndx == -1:
            return line.strip()
        return line[:ndx].strip()


