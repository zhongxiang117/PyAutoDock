

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




