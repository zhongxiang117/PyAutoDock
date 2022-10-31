from PyAutoDock.read_pdbqt import ReadPDBQTLigand
from utils import allclose

import os

CWD = os.path.split(os.path.abspath(__file__))[0]

Info = ['good']
Mark = [
    [0,1,1,2,3,4,5,6,7,8,9],
    [2,3,3,4,5,6,7,8,9],
    [3,5,4,5,6,7,8,9],
    [5,6,6,7,8,9],
    [6,9,7,8,9],
    [0,10,10],
    [0,11,11,12,13],
    [0,13,12,13],
]
File = ["""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
ENDROOT
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ATOM      3  H11 MID     1      30.593  13.997  22.470  1.00  0.00     0.164 HD
BRANCH   3   4
ATOM      4  C19 MID     1      29.999  13.608  24.785  1.00 20.01     0.196 C 
BRANCH   4   6
ATOM      5  N18 MID     1      30.002  13.424  26.198  1.00 23.25    -0.290 N 
ATOM      6  H19 MID     1      30.729  13.877  26.718  1.00  0.00     0.174 HD
BRANCH   6   7
ATOM      7  S22 MID     1      29.437  12.014  26.802  1.00 24.09     0.261 S 
BRANCH   7   10
ATOM      8  A26 MID     1      30.856  11.006  27.018  1.00 19.26     0.103 A 
ATOM      9  A25 MID     1      30.841   9.779  26.364  1.00 19.50     0.013 A 
ATOM     10  A27 MID     1      31.920  11.408  27.767  1.00 15.97     0.013 A 
ENDBRANCH
ENDBRANCH
ENDBRANCH
ENDBRANCH
ENDBRANCH
BRANCH   1  11
ATOM     11  C9  MID     1      32.876  13.242  21.334  1.00 16.44     0.254 C 
ENDBRANCH
BRANCH   1  12
ATOM     12  C11 MID     1      32.632  15.671  21.590  1.00  6.47     0.058 C 
BRANCH   1  14
ATOM     13  A12 MID     1      32.636  16.760  22.572  1.00  6.00    -0.020 A 
ATOM     14  A14 MID     1      31.443  17.235  23.050  1.00  7.63    -0.004 A 
ENDBRANCH
ENDBRANCH
TORSDOF 10
"""]

Info.append('good')
Mark.append([
    [0,1,1,2],
    [2,3,3,4,5],
    [3,5,4,5],
])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
ENDROOT
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ATOM      3  H11 MID     1      30.593  13.997  22.470  1.00  0.00     0.164 HD
ENDBRANCH
BRANCH   3   4
ATOM      4  C19 MID     1      29.999  13.608  24.785  1.00 20.01     0.196 C 
BRANCH   4   6
ATOM      5  N18 MID     1      30.002  13.424  26.198  1.00 23.25    -0.290 N 
ATOM      6  H19 MID     1      30.729  13.877  26.718  1.00  0.00     0.174 HD
ENDBRANCH
ENDBRANCH
TORSDOF 10
""")

Info.append('good')
Mark.append([
    [0,1,1,2],
    [2,3,3,4,5],
    [3,5,4,5],
])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ATOM      3  H11 MID     1      30.593  13.997  22.470  1.00  0.00     0.164 HD
ENDBRANCH
BRANCH   3   4
ATOM      4  C19 MID     1      29.999  13.608  24.785  1.00 20.01     0.196 C 
BRANCH   4   6
ATOM      5  N18 MID     1      30.002  13.424  26.198  1.00 23.25    -0.290 N 
ATOM      6  H19 MID     1      30.729  13.877  26.718  1.00  0.00     0.174 HD
ENDBRANCH
ENDBRANCH
ENDROOT
TORSDOF 10
""")

Info.append('bad: no ROOT/entries')
Mark.append([])
File.append("""
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
""")

Info.append('bad: BRANCH double ends')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
ENDROOT
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ATOM      3  H11 MID     1      30.593  13.997  22.470  1.00  0.00     0.164 HD
ENDBRANCH
ENDBRANCH
""")

Info.append('bad: ENDROOT in-between')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   2
ENDROOT
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ATOM      3  H11 MID     1      30.593  13.997  22.470  1.00  0.00     0.164 HD
ENDBRANCH
""")

Info.append('bad: ROOT should be in first place / BRANCH should have anchors')
Mark.append([])
File.append("""
BRANCH   1   2
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
ENDBRANCH
ROOT
ENDROOT
""")

Info.append('bad: BRANCH closed multiple times / double ends')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ENDBRANCH
ENDBRANCH
ENDROOT
""")

Info.append('bad: BRANCH closed earlier / double ends')
Mark.append([])
File.append("""
ROOT
ENDBRANCH
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ENDBRANCH
ENDROOT
""")

Info.append('bad: BRANCH anchor not correct')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   3
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ENDBRANCH
ENDROOT
""")

Info.append('bad: ROOT nested')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
ROOT
BRANCH   1   3
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ENDBRANCH
ENDROOT
ENDROOT
""")

Info.append('bad: ROOT closed multiple times / open ENDROOT')
Mark.append([])
File.append("""
ROOT
ATOM      1  C8  MID     1      32.710  14.342  22.327  1.00 12.67     0.164 C 
BRANCH   1   2
ATOM      2  N7  MID     1      31.432  14.076  23.007  1.00 17.28    -0.337 N 
ENDBRANCH
ENDROOT
ENDROOT
""")

Info.append('bad: TORSDOF multiple times')
Mark.append([])
File.append("""
ROOT
TORSDOF 1
TORSDOF 2
ENDROOT
""")

fn = ReadPDBQTLigand()
for i,f,cl in zip(Info,File,Mark):
    proline,tlist,m= fn.read(f.split('\n'))
    print(f'Check: {i}')
    if tlist and len(tlist) == len(cl):
        for t,c in zip(tlist,cl):
            assert t == c
    print()




