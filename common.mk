# Definitions common to all Makefiles
# This file is included from other Makefiles in the augustus project.
AUGVERSION = 3.3.3

# set ZIPINPUT to false if you do not require input gzipped input genome files,
# get compilation errors about the boost iostreams library or
# the required libraries libboost-iostreams-dev and lib1g-dev are not available
ZIPINPUT = true

# set COMPGENEPRED to false if you do not require the comparative gene prediction mode (CGP) or
# the required libraries libgsl-dev, libboost-all-dev, libsuitesparse-dev and liblpsolve55-dev are not available
COMPGENEPRED = true
