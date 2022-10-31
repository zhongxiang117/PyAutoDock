

# level for debug
#   a) for debug, please set this number less than 20, but not equals to 0
#   b) for normal info printout, please setup in range 10<level<=20
#   c) for warning, in range 20<level<=30
#   d) for critical, in range level>=50
#
# level 0 means not set, level 40(error) is not used
DEBUG_LEVEL = 20

# for interaction with grids, unit in Angstrom
NONBONDED_CUTOFF = 8

# highest energy number
ENERGY_CLAMP = 1000

# number of grid points to be used in one Angstrom
NUMBER_OF_GRIDS_IN_ONE_ANGSTROM = 100

# cutoff for internal energy interactions, index starts at 0
INTERNAL_ENERGY_CUTOFF = 2048 - 1





# avoid quotient being zero
QUOTIENT_APPROACHES_ZERO = 0.00000001




NTORDIVS = 256
MAX_RUNS = 2000
MAX_TORS = 32
MAX_NONBONDS = 524288



# for internal use only

# for identification of bonds
# values from Handbook of Chemistry and Physics
# 0=C, 1=N, 2=O, 3=H, 4=XX, 5=P, 6=S
# Caution: when using, they should be mutual, C-N === N-C
DISTANCE_MINMAX = [
    [0, 0, 1.20,   1.545],  # p.3510; p.3511
    [0, 1, 1.10,   1.479],  # p.3510; p.3511
    [0, 2, 1.15,   1.47 ],  # p.3510; p.3512
    [0, 3, 1.022,  1.12 ],  # p.3518; p.3517
    [0, 4, 0.90,   1.545],  # AutoDock 3 defaults; p.3511
    [0, 5, 1.85,   1.89 ],  # p.3510; p.3510
    [0, 6, 1.55,   1.835],  # p.3510; p.3512
    [1, 1, 1.0974, 1.128],  # p.3513; p.3515
    [1, 2, 1.0619, 1.25 ],  # p.3515; p.3515
    [1, 3, 1.004,  1.041],  # p.3516; p.3515
    [1, 4, 0.90,   1.041],  # AutoDock 3 defaults; p. 3515
    [1, 5, 1.4910, 1.491],  # p.3515; p.3515
    [1, 6, 1.580,  1.672],  # 1czm.pdb sulfonamide; J. Chem. SOC., Dalton Trans., 1996, Pages 4063-4069
    [2, 2, 1.208,  1.51 ],  # p.3513; p.3515
    [2, 3, 0.955, 1.0289],  # p.3515; p.3515
    [2, 4, 0.955,  2.1  ],  # AutoDock 3 defaults
    [2, 5, 1.36,   1.67 ],  # p.3516;  p. 3517
    [2, 6, 1.41,   1.47 ],  # p.3517; p. 3515
    [3, 3, 100.,  -100. ],  # impossible values to prevent such bonds from forming.
    [3, 4, 0.9,    1.5  ],  # AutoDock 4 defaults
    [3, 5, 1.40,   1.44 ],  # p.3515;  p. 3515
    [3, 6, 1.325, 1.3455],  # p.3518; p. 3516
    [4, 4, 0.9,    2.1  ],  # AutoDock 3 defaults
    [4, 5, 0.9,    2.1  ],  # AutoDock 3 defaults
    [4, 6, 1.325,  2.1  ],  # p.3518; AutoDock 3 defaults
    [5, 5, 2.18,   2.23 ],  # p.3513; p.3513
    [5, 6, 1.83,   1.88 ],  # p.3516; p.3515
    [6, 6, 2.03,   2.05 ],  # p.3515; p.3515
]




