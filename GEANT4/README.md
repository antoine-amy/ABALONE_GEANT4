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
    ├── plots                          # G4SiPM dependency, files for plotting results
    ├── resources                      # G4SiPM dependency, data of several SiPM
    ├── run			           # G4SiPM dependency, python scripts
    ├── src			           # source files
    ├── ABsimulation.cc		   # main function file
    ├── pe.mac			   # macro simulating the emission of a electron in the ABALONE
    ├── decay.mac                      # macro simulating the decay of 176Lu in LYSO
    └── calibration.mac                # macro used to calibrate the ABALONE, generates blue light photons
