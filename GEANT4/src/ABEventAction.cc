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
/// \file ABEventAction.cc
/// \brief Implementation of the ABEventAction class

#include "ABEventAction.hh"
//#include "ABCalorimeterSD.hh"
//#include "ABCalorHit.hh"

#include "hit/G4SipmHit.hh"
#include "hit/G4SipmSensitiveDetector.hh"

#include "ABAnalysis.hh"

#include "ABRunAction.hh"

//#include "persistency/PersistencyHandler.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4DCtable.hh"
#include "Randomize.hh"
#include <iomanip>

#include "G4DigiManager.hh"
#include "digi/G4SipmDigi.hh"
#include "digi/G4SipmVoltageTraceDigitizer.hh"

#include "persistency/PersistencyHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABEventAction::ABEventAction()
	: G4UserEventAction(),
	fAbsHCID(-1),
	fCoatingHCID(-1),
	fbaseHCID(-1),
	fringHCID(-1),
	fallHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABEventAction::~ABEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//ABCalorHitsCollection*ABEventAction::GetHitsCollection(G4int hcID,const G4Event* event) const {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABEventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength) const //G4double gapEdep, G4double gapTrackLength
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABEventAction::EndOfEventAction(const G4Event* event)
{
    PersistencyHandler* persistency = ((ABRunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->getPersistencyHandler();
    G4DigiManager* digiManager = G4DigiManager::GetDMpointer(); 
    G4DCtable* dcTable = digiManager->GetDCtable();
    
    for (int i =0; i< dcTable->entries();i++){
        G4String dmName = dcTable->GetDMname(i);
        G4VDigitizerModule* dm = digiManager->FindDigitizerModule(dmName);
        if (dm) {
            dm->Digitize();
        }
    }
    
    G4DCofThisEvent* dCof = event->GetDCofThisEvent();
       
	if (dCof != NULL) {
		for (int i = 0; i < dCof->GetCapacity(); ++i) {
			G4VDigiCollection* dc = dCof->GetDC(i);
			if (dc != NULL) {
				if (dynamic_cast<G4SipmDigiCollection*>(dc)) {
                    //persistency->persist((G4SipmDigiCollection*) dc);
				}
				
				if (dynamic_cast<G4SipmVoltageTraceDigiCollection*>(dc)) {
                    //std::cerr << "testing the voltage trace digi collection entries " << dc->GetSize() << std::endl;
                    persistency->persist((G4SipmVoltageTraceDigiCollection*) dc);
				}
			}
			
			else { std::cerr << " I am working here but the dc is NULL " << std::endl;}
		}
	}
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
