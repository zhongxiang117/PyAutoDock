from PyAutoDock import logger
from PyAutoDock.utils import file_gen_new

import collections
import os

logger = logger.mylogger()

# Future keys can be added into list, order matters, case-insensitive
AD_DPF_KEYS = [
    #{'key':'autodock_parameter_version',    'default':'4.2',},
    {'key':'outlev',                'default':None, },
    {'key':'parameter_file',        'default':None, },  # library
    {'key':'intelec',               'default':True, },
    {'key':'seed',                  'default':[],   },
    {'key':'smooth',                'default':0.5,  },
    {'key':'ligand_types',          'default':[],   },  # can be guessed from file
    {'key':'fld',                   'default':None, },
    {'key':'map',                   'default':[],   },
    {'key':'elecmap',               'default':None, },
    {'key':'desolvmap',             'default':None, },
    {'key':'move',                  'default':None, },
    {'key':'about',                 'default':[],   },
    {'key':'trans0',                'default':[],   },
    {'key':'quaternion0',           'default':[],   },
    {'key':'axisangle0',            'default':[],   },
    {'key':'quat0',                 'default':[],   },  # synonym for `axisangle0`
    {'key':'dihe0',                 'default':[],   },
    {'key':'torsdof',               'default':None, },
    {'key':'rmstol',                'default':None, },
    {'key':'extnrg',                'default':None, },
    {'key':'e0max',                 'default':[],   },
    {'key':'ga_pop_size',           'default':None, },
    {'key':'ga_num_evals',          'default':None, },
    {'key':'ga_num_generations',    'default':None, },
    {'key':'ga_elitism',            'default':None, },
    {'key':'ga_mutation_rate',      'default':None, },
    {'key':'ga_crossover_rate',     'default':None, },
    {'key':'ga_window_size',        'default':None, },
    {'key':'ga_cauchy_alpha',       'default':None, },
    {'key':'ga_cauchy_beta',        'default':None, },
    {'key':'set_ga',                'default':True, },
    {'key':'unbound_model',         'default':None, },
    {'key':'analysis',              'default':True, },
]


class ReadDPF:
    """read Grid Parameter File

    keys will be smartly picked, case-insensitive, orders does not matter
    """
    def __init__(self,filename=None,loglevel=None,*args,**kwargs):
        level = logger.level
        if loglevel: logger.setLevel(loglevel)
        self.filename = filename
        self.dpf = collections.OrderedDict()
        for k in AD_DPF_KEYS:
            self.dpf[k['key'].lower()] = k['default']
        self._read()
        logger.setLevel(level)      # restore

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
                if key == 'outlev':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            logger.info('>DPF: outlev: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('outlev: not a number: {:}'.format(line))
                    else:
                        logger.critical('outlev: number of args should be 2: {:}'.format(line))
                elif key == 'smooth':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: smooth: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('smooth: not a number: {:}'.format(line))
                    else:
                        logger.critical('smooth: number of args should be 2: {:}'.format(line))
                elif key == 'intelec':
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2:
                        if ltmp[1].lower() in ['on','true','t','o']:
                            self.dpf[key] = True
                        else:
                            self.dpf[key] = False
                    else:
                        logger.critical('intelec: wrong line: {:}'.format(line))
                elif key == 'seed':
                    if len(ltmp) == 3:
                        try:
                            self.dpf[key] = list(map(int,ltmp[1:]))
                            logger.info('>DPF: seed: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('seed: not a number: {:}'.format(line))
                    else:
                        logger.critical('seed: number of args should be 3: {:}'.format(line))
                elif key == 'ligand_types':
                    self.dpf[key].extend(ltmp[1:])
                    logger.info('>DPF: ligand_types: {:}'.format(self.dpf[key]))
                elif key == 'fld':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            self.dpf[key] = ltmp[1]
                            logger.info('>DPF: fld: {:}'.format(self.dpf[key]))
                        else:
                            logger.critical('not a file: {:}'.format(line))
                    else:
                        logger.critical('fld: number of args should be 2: {:}'.format(line))
                elif key == 'map':
                    for i in ltmp[1:]:
                        if os.path.isfile(i):
                            self.dpf[key].append(i)
                            logger.info('>DPF: map: {:}'.format(self.dpf[key]))
                        else:
                            logger.critical('map: file not exit: {:}'.format(i))
                elif key == 'elecmap':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            self.dpf[key] = ltmp[1]
                            logger.info('>DPF: elecmap: {:}'.format(self.gpf[key]))
                        else:
                            logger.critical('elecmap: file not exit: {:}'.format(ltmp[1]))
                    else:
                        logger.critical('elecmap: number of args should be 2: {:}'.format(line))
                elif key == 'desolvmap':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            self.dpf[key] = ltmp[1]
                            logger.info('>DPF: desolvmap: {:}'.format(self.gpf[key]))
                        else:
                            logger.critical('desolvmap: file not exit: {:}'.format(ltmp[1]))
                    else:
                        logger.critical('desolvmap: number of args should be 2: {:}'.format(line))
                elif key == 'move':
                    if len(ltmp) == 2:
                        if os.path.isfile(ltmp[1]):
                            self.dpf[key] = ltmp[1]
                            logger.info('>DPF: move: {:}'.format(self.gpf[key]))
                        else:
                            logger.critical('move: file not exit: {:}'.format(ltmp[1]))
                    else:
                        logger.critical('move: number of args should be 2: {:}'.format(line))
                elif key == 'about':
                    if len(ltmp) == 4:
                        try:
                            self.dpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('about: value error: {:}'.format(line))
                    else:
                        logger.critical('about: number of args should be 3: {:}'.format(line))
                elif key == 'tran0':
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2 and ltmp[1].lower() == 'random':
                        pass
                    elif len(ltmp) == 4:
                        try:
                            self.dpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('tran0: value error: {:}'.format(line))
                    else:
                        logger.critical('tran0: number of args should be 4: {:}'.format(line))
                elif key == 'quaternion0':
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2 and ltmp[1].lower() == 'random':
                        pass
                    if len(ltmp) == 5:
                        try:
                            self.dpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('quaternion0: value error: {:}'.format(line))
                    else:
                        logger.critical('quaternion0: number of args should be 5: {:}'.format(line))
                elif key in ['axisangle', 'quat0']:
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2 and ltmp[1].lower() == 'random':
                        pass
                    if len(ltmp) == 5:
                        try:
                            self.dpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('quaternion0: value error: {:}'.format(line))
                    else:
                        logger.critical('quaternion0: number of args should be 5: {:}'.format(line))
                elif key == 'dihe0':
                    if 'random' in [i.lower() for i in ltmp[1:]]:
                        pass
                    else:
                        try:
                            self.dpf[key] = list(map(float,ltmp[1:]))
                        except ValueError:
                            logger.critical('dihe0: value error: {:}'.format(line))
                elif key == 'torsdof':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            logger.info('>DPF: torsdof: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('torsdof: not a number: {:}'.format(line))
                    else:
                        logger.critical('torsdof: number of args should be 2: {:}'.format(line))
                elif key == 'rmstol':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: rmstol: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('rmstol: not a number: {:}'.format(line))
                    else:
                        logger.critical('rmstol: number of args should be 2: {:}'.format(line))
                elif key == 'extnrg':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: extnrg: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('extnrg: not a number: {:}'.format(line))
                    else:
                        logger.critical('extnrg: number of args should be 2: {:}'.format(line))
                elif key == 'e0max':
                    if len(ltmp) == 3:
                        try:
                            self.dpf[key].append(float(ltmp[1]))
                            self.dpf[key].append(int(ltmp[2]))
                            if self.dpf[key][1] <= 0: raise ValueError
                        except ValueError:
                            logger.critical('e0max: value error: {:}'.format(line))
                    else:
                        logger.critical('e0max: number of args should be 2: {:}'.format(line))
                elif key == 'ga_pop_size':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            if self.dpf[key] <= 0: raise ValueError
                            logger.info('>DPF: ga_pop_size: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_pop_size: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_pop_size: number of args should be 2: {:}'.format(line))
                elif key == 'ga_num_evals':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            if self.dpf[key] <= 0: raise ValueError
                            logger.info('>DPF: ga_num_evals: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_num_evals: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_num_evals: number of args should be 2: {:}'.format(line))
                elif key == 'ga_num_generations':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            if self.dpf[key] <= 0: raise ValueError
                            logger.info('>DPF: ga_num_generations: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_num_generations: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_num_generations: number of args should be 2: {:}'.format(line))
                elif key == 'ga_elitism':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            logger.info('>DPF: ga_elitism: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_elitism: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_elitism: number of args should be 2: {:}'.format(line))
                elif key == 'ga_mutation_rate':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: ga_mutation_rate: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_mutation_rate: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_mutation_rate: number of args should be 2: {:}'.format(line))
                elif key == 'ga_crossover_rate':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: ga_crossover_rate: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_crossover_rate: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_crossover_rate: number of args should be 2: {:}'.format(line))
                elif key == 'ga_window_size':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            logger.info('>DPF: ga_window_size: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_window_size: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_window_size: number of args should be 2: {:}'.format(line))
                elif key == 'ga_cauchy_alpha':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: ga_cauchy_alpha: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_cauchy_alpha: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_cauchy_alpha: number of args should be 2: {:}'.format(line))
                elif key == 'ga_cauchy_beta':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = float(ltmp[1])
                            logger.info('>DPF: ga_cauchy_beta: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('ga_cauchy_beta: not a number: {:}'.format(line))
                    else:
                        logger.critical('ga_cauchy_beta: number of args should be 2: {:}'.format(line))
                elif key == 'set_ga':
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2:
                        if ltmp[1].lower() in ['on','true','t','o']:
                            self.dpf[key] = True
                        else:
                            self.dpf[key] = False
                    else:
                        logger.critical('set_ga: wrong line: {:}'.format(line))
                elif key == 'unbound_model':
                    if len(ltmp) == 2:
                        self.dpf[key] = ltmp[1].lower()
                    else:
                        logger.critical('unbound_model: wrong line: {:}'.format(line))
                elif key == 'do_global_only':
                    if len(ltmp) == 2:
                        try:
                            self.dpf[key] = int(ltmp[1])
                            logger.info('>DPF: do_global_only: {:}'.format(self.dpf[key]))
                        except ValueError:
                            logger.critical('do_global_only: not a number: {:}'.format(line))
                    else:
                        logger.critical('do_global_only: number of args should be 2: {:}'.format(line))
                elif key == 'analysis':
                    if len(ltmp) == 1:
                        pass
                    elif len(ltmp) == 2:
                        if ltmp[1].lower() in ['on','true','t','o']:
                            self.dpf[key] = True
                        else:
                            self.dpf[key] = False
                    else:
                        logger.critical('analysis: wrong line: {:}'.format(line))


    def _strip_comments(self,line):
        ndx = line.find('#')
        if ndx == -1:
            return line.strip()
        return line[:ndx].strip()


