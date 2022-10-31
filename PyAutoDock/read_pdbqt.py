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
        """centroid molecule at [x,y,z], if any one of them exists.
           otherwise, centroid at the center.
        """
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
        else:
            self.translate(*self.center)

    def translate(self,x=None,y=None,z=None):
        """translate molecule to x, y, or z, if any one of them exists"""
        if x:
            for i in range(len(self.atoms)):
                self.atoms[i][2] -= x
        if y:
            for i in range(len(self.atoms)):
                self.atoms[i][3] -= y
        if z:
            for i in range(len(self.atoms)):
                self.atoms[i][4] -= z
        if x or y or z:
            self._process()
            self.info()


class ReadPDBQTLigand:
    #TODO add default atomtypes, found in library
    def __init__(self,filename=None,*args,**kws):
        self.filename = filename
        self.proline = []
        self.tlist = []

    def read(self,filename=None):
        """read ligand pdbqt file
        
        Return:
            proline : List[(x,y,z,charge,atomname,atomtype), ...]
            tlist   : List[[branch_begin_index, branch_end_index, i1, i2, ...], ...]
            m       : int   :   TORSDOF: number of freedom of torsions

        Note:
            when error happens, return values will be: `[], [], 0`

        Test:
            case 1: R, B, B, EB, EB, ER     (good)
            case 2: R, ER, B, EB, B, EB     (good)
            case 3: R, B, ER, EB, ...       (wrong, ER cannot be in-between)
            case 4: B, EB, ...              (wrong, no ROOT found)
            case 5: B, EB, R, ER, ...       (wrong, ROOT should be in first)
            case 6: EB, B, ...              (wrong, no ROOT found)
            case 7: R, B, EB, EB, ...       (wrong, BRANCH closed multiple times)
            case 8: R, EB, B, ...           (wrong, BRANCH closed earlier)
            case 9: R, ER, ER, ...          (wrong, ROOT closed multiple times)
        """
        if not filename: filename = self.filename
        if isinstance(filename, list):
            flines = filename
        elif os.path.isfile(filename):
            with open(filename,'rt') as f:
                flines = f.readlines()
        else:
            print(f'Error: file is not set: either can be a file or 2D-list')
            return [],[],0
        proline = []
        r = -1          # index for `ROOT` and `ENDROOT`
        k = -1          # index for `BRANCH` and `ENDBRANCH`
        b = e = None    # anchors for `BRANCH` and `ENDBRANCH`
        n = -1          # index for `ATOM` and `HETATM`
        tdict = {}      # temporary helper list
        tlist = []      # results
        m = None        # number for `TORSDOF`
        boroot = False
        ndxlist = []    # for check of ROOT/ENDROOT in-between
        for l in flines:
            if len(l) >= 4:
                s = l[:4].lower()
                if s == 'root':
                    r += 1
                    boroot = True
                    ndxlist.append('R')
                elif s == 'endr':
                    r -= 1
                    ndxlist.append('ER')
                elif s in ['atom','heta']:
                    n += 1
                    try:
                        atomname = l[12:16].strip()
                        x = float(l[30:38])
                        y = float(l[38:46])
                        z = float(l[46:54])
                        c = float(l[68:76])
                        atomtype = l[77:79].strip().upper()
                    except (ValueError,IndexError):
                        print(f'Error: line input: {l}')
                        return [],[],0
                    else:
                        #TODO when atomtype not found, get it from atomname
                        proline.append((x,y,z,c,atomname,atomtype))
                elif s == 'bran':
                    ndxlist.append('B')
                    k += 1
                    ls = l.split()
                    try:
                        if len(ls) != 3:
                            raise ValueError('Error: BRANCH: wrong settings')
                        b = int(ls[1])
                        e = int(ls[2])
                        if b > n+1 or e < n+1:
                            raise ValueError('Error: BRANCH: values are not correspondant')
                    except (ValueError,IndexError) as e:
                        print(f'Error: line input: BRANCH should have anchors: {l}\n  -> {e}')
                        return [],[],0
                elif s == 'endb':
                    ndxlist.append('EB')
                    k -= 1
                elif s == 'tors':
                    if m:
                        print(f'Error: duplicate line: {l}')
                        return [],[],0
                    ls = l.split()
                    try:
                        if len(ls) != 2:
                            raise ValueError('Error: TORSDOF: wrong setting')
                        m = int(ls[1])
                    except (ValueError,IndexError) as e:
                        print(f'Error: line input: {l}\n  -> {e}')
                        return [],[],0
                if k < -1:
                    print(f'Error: double ends: {l}')
                    return [],[],0
                elif k == -1:
                    if tdict:
                        v = 0
                        while v in tdict:
                            q = v + 1
                            ls = tdict[v]
                            while q in tdict:
                                ls.extend(tdict[q][2:])     # in-place edition
                                q += 1
                            v += 1
                        for v in sorted(tdict.keys()):
                            tlist.append(
                                tdict[v][:2] + sorted(list(set(tdict[v][2:])))
                            )
                        tdict = {}      # reset
                elif k in tdict:
                    tdict[k].append(n)
                else:
                    tdict[k] = [b-1,e-1]      # k controlled `b` and `e` values
                if r > 0:
                    print(f'Error: ROOT cannot be nested: {l}')
                    return [],[],0
                elif r < -1:
                    print(f'Error: ENDROOT: not-included: {l}')
                    return [],[],0
        if not tlist:
            print('Warning: no entry is found')
            return [],[],0
        if not boroot:
            print(f'Error: no ROOT entry is found')
            return [],[],0
        if r != -1:     # sequence matters, should be after `boroot`
            print(f'Error: ROOT is not closed')
            return [],[],0
        if m is None:
            print('Warning: TORSDOF: not set: number of freedom of torsions')
        elif m != len(tlist):
            print('Warning: number of freedom of torsions is not equal to defined entries')
            print(f'    ->: TORSDOF:{m}   !=   Defines:{len(tlist)}')
        for e in tlist:
            if e[1] not in e[2:]:
                print(f'Error: BRANCH: anchor is not correct: {e[0]}  {e[1]}')
                return [],[],0
        while True:
            if 'R' not in ndxlist:
                break
            i = ndxlist.index('R')
            j = ndxlist.index('ER')
            g = ndxlist[i:j+1].count('B')
            if g and (j-i-1) // g != 2:
                print('Error: ROOT/ENDROOT: cannot be in-between BRANCH')
                return [],[],0
            ndxlist.pop(i)
            ndxlist.pop(j-1)
        self.proline = proline
        self.tlist = tlist
        return proline,tlist,len(tlist)

    def info(self):
        if not self.proline: return
        ct = sum([i[3] for i in self.proline])
        if abs(ct-round(ct,0)) > 0.01:
            print(f'Warning: total charge {ct} is not an integer number')
        tdict = {}
        for l in self.proline:
            if l[5] in tdict:
                tdict[l[5]] += 1
            else:
                tdict[l[5]] = 1
        print()
        for k in sorted(tdict.keys()):
            print('Note: number of entries for atomtype {:}: {:}'.format(k,tdict[k]))
        print('\nTorsion Tree:\n #  Atom1--Atom2  Moved  List of Atoms Moved')
        print('___ ____________ ______ ______________________________________')
        for i,t in enumerate(self.tlist):
            print(
                '{:^3} {:^5}--{:^5} {:^6} {:}'.format(i+1,t[0],t[1],len(t[2:]),
                ', '.join([self.proline[j][4] for j in t[2:]]))
            )
        print()






