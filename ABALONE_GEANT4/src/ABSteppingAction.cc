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
/// \file ABSteppingAction.cc
/// \brief Implementation of the ABSteppingAction class

#include "ABSteppingAction.hh"
#include "ABEventAction.hh"
#include "ABDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABSteppingAction::ABSteppingAction() //ABEventAction* eventAction
: G4UserSteppingAction(),
  //fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABSteppingAction::~ABSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABSteppingAction::UserSteppingAction(const G4Step* step)
{

	auto analysisManager = G4AnalysisManager::Instance();

	G4String all_lv[16] = { "World" ,"GlassShellLV","AU_sealingLV","In_sealingLV","Cr_sealingLV","PhotocathodeLV","Base_logical","ring_logical1","ring_logical2","ring_logical3","ring_logical4" ,"ring_logical5","ring_logical6","Scint_sealing_Au_LV","Scint_sealing_Cr_LV","ScintillatorLV" };
	G4int logical_volume_pos= 99;
	G4int part_def_index = 99;

	auto log_volume_name = step->GetTrack()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
	for (int i = 0; i <= 15; i++) {
		if ( log_volume_name == all_lv[i]) {
			logical_volume_pos = i;
			break;
		}
	}
	
	G4String part_def_array[2] = { "e-" ,"opticalphoton" };
	auto part_def = step->GetTrack()->GetParticleDefinition()->GetParticleName();
	for (int i = 0; i <= 1; i++) {
		if (part_def == part_def_array[i]) {
			part_def_index = i;
			break;
		}
	}

	
	analysisManager->FillNtupleDColumn(0, static_cast<double>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()));
	analysisManager->FillNtupleDColumn(1, step->GetTrack()->GetParentID());
	analysisManager->FillNtupleDColumn(2, step->GetTrack()->GetTrackID());
	analysisManager->FillNtupleDColumn(3, part_def_index);
	analysisManager->FillNtupleDColumn(4, step->GetPreStepPoint()->GetPosition().getX()/mm);
	analysisManager->FillNtupleDColumn(5, step->GetPreStepPoint()->GetPosition().getY() / mm);
	analysisManager->FillNtupleDColumn(6, step->GetPreStepPoint()->GetPosition().getZ() / mm);
	analysisManager->FillNtupleDColumn(7, step->GetPostStepPoint()->GetProperTime());
	analysisManager->FillNtupleDColumn(8, step->GetPreStepPoint()->GetKineticEnergy()/keV);
	analysisManager->FillNtupleDColumn(9, step->GetTotalEnergyDeposit()/keV);
	analysisManager->FillNtupleDColumn(10, logical_volume_pos);
	analysisManager->AddNtupleRow();
	
	//G4cout << step->GetPreStepPoint()->GetKineticEnergy() / keV << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

