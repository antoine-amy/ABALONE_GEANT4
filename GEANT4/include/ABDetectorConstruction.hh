//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file ABDetectorConstruction.hh
/// \brief Definition of the ABDetectorConstruction class

#ifndef ABDetectorConstruction_h
#define ABDetectorConstruction_h 1

//#include "G4GDMLParser.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4FieldManager.hh"
#include <vector>
#include <iostream>
#include <string.h>
#include "G4Types.hh"

#include "G4Sipm.hh"
#include "housing/G4SipmHousing.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
//class G4UniformMagfield;
class G4GenericMessenger;
class ABElectricField;

/// Detector construction class to define materials and geometry.

class ABDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ABDetectorConstruction(std::string modelName, std::string housingName);
    virtual ~ABDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  private:
	G4SipmHousing* housing;
	G4SipmModel* createSipmModel(std::string name) const;
    G4SipmHousing* createHousing(std::string name, G4Sipm* sipm) const;	
	
  private:
	  //G4UniformMagfield* magField;
	  
	  G4GenericMessenger* fMessenger;

	  static G4ThreadLocal ABElectricField* fElectricField;
	  static G4ThreadLocal G4FieldManager* fFieldMgr;
	  
	  G4LogicalVolume* fElectricLogical;
	  G4VPhysicalVolume* DefineVolumes();
	  
  public:
    G4SipmModel* getSipmModel() const;
    G4SipmHousing* getSipmHousing() const;
    
  protected:
    G4LogicalVolume*  fScoringVolume; // useless

  private:
	G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
	G4int   fNofLayers;     // number of layers

 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

