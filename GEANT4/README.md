<table align="center"><tr><td align="center" width="9999">

# Geant4 - an Object-Oriented Toolkit for Simulation in HEP
## ABALONE photosensor

</td></tr></table>

----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Project top-level directory layout

    GEANT4
    │  
    ├── externals                      # files needed for the simulation to work 
    │					└── especially the SiPM code using [G4SiPM](https://github.com/ntim/g4sipm)
    ├── include                        # header files
    ├── plots                          # files for plotting results
    ├── resources                      # data of several SiPM
    ├── run			           # python scripts
    ├── src			           # source files
    ├── ABsimulation.cc		   # main function file
    ├── pe.mac			   # macro simulating the emission of a electron in the ABALONE
    └── decay.mac                      # macro simulating the decay of 176Lu in LYSO

The thermal noise and conditions for the digitizer hits are currently disabled to allow test on the decay simulations.
