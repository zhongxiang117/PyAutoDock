import os

class ReadMap:
    def __init__(
        self,mapfile=None,npts=None,*args,**kws
    ):
        self.mapfile = mapfile
        self.npts = npts

    def read_map_from_autodock(self,mapfile=None,npts=None):
        if not mapfile: mapfile = self.mapfile
        if not os.path.isfile(mapfile):
            print(f'Error: not a file: {mapfile}')
            return []
        # try to find npts if not defined
        if not npts:
            with open(mapfile,'rt') as f:
                for line in f:
                    ltmp = line.split()
                    if len(ltmp) == 4 and ltmp[0].lower() == 'nelements':
                        try:
                            x = int(ltmp[1])
                            y = int(ltmp[2])
                            z = int(ltmp[3])
                        except ValueError:
                            pass
                        else:
                            npts = [x,y,z]
                            break
            if not npts:
                print('Error: npts: NELEMENTS not defined')
                return []
        data = []
        with open(mapfile,'rt') as f:
            # dump first 6 lines
            for i in range(6):
                f.readline()
            try:
                data = [float(l.strip()) for l in f]
            except ValueError:
                data = []
        if not data:
            print('Error: not all values are the number, double check it')
            return []
        if len(data) != (npts[0]+1)*(npts[1]+1)*(npts[2]+1):
            print('Error: number of data inside map is not equal to defined points')
            return []
        print(f'Note: read total number of {len(data)} points')
        return data


    def _split(self,data,npts):
        









