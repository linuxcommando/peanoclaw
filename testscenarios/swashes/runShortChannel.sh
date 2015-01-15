#!/bin/bash

#../../bin/peano-claw-fullswof2d channelPseudo2D 36 36 12 1000 20 short --usePeano
#../../bin/peano-claw-fullswof2d channelPseudo2D 36 36 12 1000 1 short --usePeano
mpirun -np 11 ../../bin/peano-claw-fullswof2d-parallel channelPseudo2D 36 36 12 1000 1 short --usePeano
