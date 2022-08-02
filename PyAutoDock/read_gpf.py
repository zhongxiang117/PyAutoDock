from PyAutoDock import logger
from PyAutoDock.utils import file_gen_new

import collections
import os

logger = logger.mylogger()

# Future keys can be added into list, order matters, case-insensitive
AD_GPF_KEYS = [
    {'key':'npts',              'default':[0, 0, 0],},
    {'key':'gridfld',           'default':None,     },
    {'key':'spacing',           'default':0.375,    },
    {'key':'receptor_types',    'default':[],       },  # can be guessed from file
    {'key':'ligand_types',      'default':[],       },  # can be guessed from file
    {'key':'receptor',          'default':None,     },
    {'key':'gridcenter',        'default':None,     },  # means `auto`
    {'key':'smooth',            'default':0.5,      },
    {'key':'map',               'default':[],       },
    {'key':'elecmap',           'default':None,     },
    {'key':'dsolvmap',          'default':None,     },
    {'key':'dielectric',        'default':-0.1465,  },
]


class ReadGPF:
    """read Grid Parameter File

    keys will be smartly picked, case-insensitive, orders does not matter
    """
    def __init__(self,filename=None,*args,**kwargs):
        self.filename = filename
        self.gpf = collections.OrderedDict()
        for k in AD_GPF_KEYS:
            self.gpf[k['key'].lower()] = k['default']
        self._read()

    def _read(self):
        if not os.path.isfile(self.filename):
            logger.critical('not the file: {:}'.format(self.filename))
            return
        with open(self.filename,'rt') as f:
            for line in f:
                new = self._strip_comments(line)
                if not new: continue
                ltmp = new.split()
                key = ltmp[0].lower()
                if key == 'npts':
                    if len(ltmp) == 4:
                        try:
                            t = list(map(int,ltmp[1:]))
                        except ValueError:
                            logger.critical('npts: value error: {:}'.format(line))
                        else:
                            if False in ([i%2 == 0 for i in t]):
                                logger.critical('npts: not the even number: {:}'.format(line))
                            else:
                                self.gpf[key] = t
                    else:
                        logger.critical('npts: number of args should be 3: {:}'.format(line))
                elif key == 'gridfld':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            logger.warning('gridfld: file exit: {:}'.format(ltmp[1]))
                            new = file_gen_new(ltmp[1])
                            logger.info(
                                'gridfld: changing filename from {:} to {:} to'
                                'avoid overwritting'.format(ltmp[1],new)
                            )
                            self.gpf[key] = new
                        else:
                            self.gpf[key] = ltmp[1]
                    else:
                        logger.critical('gridfld: number of args should be 2: {:}'.format(line))
                elif key == 'spacing':
                    if len(ltmp) == 2:
                        try:
                            self.gpf[key] = float(ltmp[1])
                        except ValueError:
                            logger.critical('spacing: not a number: {:}'.format(line))
                    else:
                        logger.critical('spacing: number of args should be 2: {:}'.format(line))
                elif key == 'receptor_types':
                    self.gpf[key].extend(ltmp[1:])
                elif key == 'ligand_types':
                    self.gpf[key].extend(ltmp[1:])
                elif key == 'receptor':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            self.gpf[key] = ltmp[1]
                        else:
                            logger.critical('not a file: {:}'.format(line))
                    else:
                        logger.critical('receptor: number of args should be 2: {:}'.format(line))
                elif key == 'gridcenter':
                    if len(ltmp) == 2:
                        if ltmp[1].lower() not in ['a','auto']:
                            logger.critical('gridcenter: wrong setting: {:}'.format(line))
                    elif len(ltmp) == 4:
                        try:
                            self.gpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('gridcenter: value error: {:}'.format(line))
                    else:
                        logger.critical('gridcenter: number of args: {:}'.format(line))
                elif key == 'smooth':
                    if len(ltmp) == 2:
                        try:
                            self.gpf[key] = float(ltmp[1])
                        except ValueError:
                            logger.critical('smooth: not a number: {:}'.format(line))
                    else:
                        logger.critical('smooth: number of args should be 2: {:}'.format(line))
                elif key == 'map':
                    for i in ltmp[1:]:
                        if os.path.isfile(i):
                            logger.warning('elecmap: file exit: {:}'.format(i))
                            new = file_gen_new(i)
                            logger.info(
                                'elecmap: changing filename from {:} to {:} to'
                                'avoid overwritting'.format(i,new)
                            )
                            self.gpf[key] = new
                        else:
                            self.gpf[key] = i
                elif key == 'elecmap':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            logger.warning('elecmap: file exit: {:}'.format(ltmp[1]))
                            new = file_gen_new(ltmp[1])
                            logger.info(
                                'elecmap: changing filename from {:} to {:} to'
                                'avoid overwritting'.format(ltmp[1],new)
                            )
                            self.gpf[key] = new
                        else:
                            self.gpf[key] = ltmp[1]
                    else:
                        logger.critical('elecmap: number of args should be 2: {:}'.format(line))
                elif key == 'dsolvmap':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            logger.warning('dsolvmap: file exit: {:}'.format(ltmp[1]))
                            new = file_gen_new(ltmp[1])
                            logger.info(
                                'dsolvmap: changing filename from {:} to {:} to'
                                'avoid overwritting'.format(ltmp[1],new)
                            )
                            self.gpf[key] = new
                        else:
                            self.gpf[key] = ltmp[1]
                    else:
                        logger.critical('dsolvmap: number of args should be 2: {:}'.format(line))
                elif key == 'dielectric':
                    if len(ltmp) == 2:
                        try:
                            self.gpf[key] = float(ltmp[1])
                        except ValueError:
                            logger.critical('dielectric: not a number: {:}'.format(line))
                    else:
                        logger.critical('dielectric: number of args should be 2: {:}'.format(line))
                else:
                    logger.warning('currently not built: ignoring: {:}'.format(line))

    def _strip_comments(self,line):
        ndx = line.find('#')
        if ndx == -1:
            return line.strip()
        return line[:ndx].strip()


