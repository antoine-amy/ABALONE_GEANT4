
<table align="center"><tr><td align="center" width="9999">

# Geant4 - an Object-Oriented Toolkit for Simulation in HEP
## ABALONE photosensor

</td></tr></table>

----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Building the simulation from source code:
```
mkdir build
cd build
sudo cmake -DGeant4_DIR=<Your Geant4 install Directory (like /usr/share/geant4/geant4.10.07-install/lib/Geant4-10.7.1/)> <Your Directory>ABALONE_GEANT4/GEANT4
sudo make install -jN && sudo make
```

The Geant4 version used for the development of the simulation is 10.07.1.

----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Folder structure :

    My Directory
    │  
    ├── build
    ├── GEANT4
    ├── E_field
    └── results
    	├── SiPM
		└── tracking

----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Execute the simulation:

- In 'interactive mode' with visualization: ``` ./ABsimulation```
- In 'batch' mode from macro files: ```./ABsimulation -m full_run.mac```
