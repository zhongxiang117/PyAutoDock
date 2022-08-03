from PyAutoDock import logger

import os

logger = logger.mylogger()


class ReadPDBQT:
    def __init__(self,filename=None,*args,**kwargs):
        self.filename = filename
        self.atoms = []
        self.xmin = 0.0
        self.xmax = 0.0
        self.ymin = 0.0
        self.ymax = 0.0
        self.zmin = 0.0
        self.zmax = 0.0
        self.qmin = 0.0
        self.qmax = 0.0
        self.total_charges = 0.0
        self.center = [0.0, 0.0, 0.0]
        self.counter = {'undefined':0, }
        self._read()

    def _read(self):
        if not os.path.isfile(self.filename): return
        with open(self.filename,'rt') as f:
            for line in f:
                if len(line) < 76: continue
                if line[:4].lower() not in ['atom','heta','char']: continue

                atomname = line[12:16].strip()          # strip is important
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    charge = float(line[70:76])
                except ValueError:
                    logger.warning('wrong setting: ignoring: line: {:}'.format(line))
                else:
                    bo = False
                    if len(line) >= 78:
                        atomtype = line[77:79].strip()  # strip is important
                        if not atomtype:
                            bo = True
                    else:
                        bo = True
                    if bo:
                        # guess atomtype, striping numeric numbers
                        atomtype = atomname.strip('0123456789').strip()
                    if len(atomtype) > 2: atomtype = None
                    self.atoms.append([atomtype,atomname,x,y,z,charge])
        if not len(self.atoms): return

        atomtypelist = [i[0] for i in self.atoms]
        self.atomset = set(atomtypelist)
        for i in self.atomset:
            key = i if i else 'undefined'
            self.counter[key] = atomtypelist.count(i)

        self._process()
        self.info()

    def _process(self):
        xtot = self.atoms[0][2]
        ytot = self.atoms[0][3]
        ztot = self.atoms[0][4]
        self.xmin, self.xmax = xtot, xtot
        self.zmin, self.ymax = ytot, ytot
        self.zmin, self.zmax = ztot, ztot
        for v in self.atoms:
            self.xmin = min(self.xmin, v[2])
            self.xmax = max(self.xmax, v[2])
            xtot += v[2]
            self.ymin = min(self.ymin, v[3])
            self.ymax = max(self.ymax, v[3])
            ytot += v[3]
            self.zmin = min(self.zmin, v[4])
            self.zmax = max(self.zmax, v[4])
            ztot += v[4]
            self.qmin = min(self.qmin, v[5])
            self.qmax = max(self.qmax, v[5])
            self.total_charges += v[5]
        n = len(self.atoms)
        self.center = [xtot/n, ytot/n, ztot/n]

    def info(self):
        logger.info('\n\nFor input Receptor file: {:}'.format(self.filename))
        logger.info('Receptor coordinates fit within the following volume:\n')
        logger.info('                   _______({:.6f} {:.6f} {:.6f})'.format(self.xmax,self.ymax,self.zmax))
        logger.info('                  /|     /|')
        logger.info('                 / |    / |')
        logger.info('                /______/  |')
        logger.info('                |  |___|__| Center = ({:.6f} {:.6f} {:.6f})'.format(*self.center))
        logger.info('                |  /   |  /')
        logger.info('                | /    | /')
        logger.info('                |/_____|/')
        logger.info('({:.6f} {:.6f} {:.6f})\n'.format(self.xmin,self.ymin,self.zmin))
        logger.info('Minimum Coordinates: ({:.6f} {:.6f} {:.6f})'.format(self.xmin,self.ymin,self.zmin))
        logger.info('Maximum Coordinates: ({:.6f} {:.6f} {:.6f})\n'.format(self.xmax,self.ymax,self.zmax))
        logger.info('  Atomtype    Counters in Receptor')
        for k,v in self.counter.items():
            logger.info(' {:^10}         {:}'.format(k,v))
        logger.info('\n')
        logger.info(f'Total number of atoms: {len(self.atoms)}')
        logger.info(f'Total number of atom types: {len(self.atomset)}')
        logger.info('Minimum charge: {:.6f} elementary charge'.format(self.qmin))
        logger.info('Maximum charge: {:.6f} elementary charge'.format(self.qmax))
        logger.info('Total charges : {:.6f} elementary charge'.format(self.total_charges))

    def centroid(self,x=None,y=None,z=None):
        if x:
            dx = self.center[0] - x
            for i in range(len(self.atoms)):
                self.atoms[i][2] -= dx
        if y:
            dy = self.center[1] - y
            for i in range(len(self.atoms)):
                self.atoms[i][3] -= dy
        if z:
            dz = self.center[2] - z
            for i in range(len(self.atoms)):
                self.atoms[i][4] -= dz
        if x or y or z:
            self._process()
            self.info()


