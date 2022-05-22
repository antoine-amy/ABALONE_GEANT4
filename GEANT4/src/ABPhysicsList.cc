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
/// \file ABPhysicsList.cc
/// \brief Implementation of the ABPhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ABPhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"

#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4PhysListFactory.hh"
#include "QGSP_BIC.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABPhysicsList::ABPhysicsList(): G4VModularPhysicsList()
{

  //add new units for radioActive decays
  // 
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);    

  //default physics
  fParticleList = new G4DecayPhysics();

  //default physics
  fRaddecayList = new G4RadioactiveDecayPhysics();

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics_option4();

  // Optical physics
  fOpPhysicsList = new G4OpticalPhysics();
	auto opticalParams=G4OpticalParameters::Instance();
	opticalParams->SetWLSTimeProfile("delta");
	opticalParams->SetScintYieldFactor(1.0);
	opticalParams->SetScintExcitationRatio(0.0);
	opticalParams->SetCerenkovMaxPhotonsPerStep(2000);
	opticalParams->SetCerenkovMaxBetaChange(100.0);
	opticalParams->SetCerenkovTrackSecondariesFirst(true);
	opticalParams->SetScintTrackSecondariesFirst(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABPhysicsList::~ABPhysicsList()
{
  delete fRaddecayList;
  delete fEmPhysicsList;
  delete fOpPhysicsList;
  delete fHadPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABPhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABPhysicsList::ConstructProcess()
{
  AddTransportation();
  // em
  fEmPhysicsList->ConstructProcess();
  // optical
  fOpPhysicsList->ConstructProcess();
  // decays
  //G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
  //G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();  
  //ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

  fParticleList->ConstructProcess();
  fRaddecayList->ConstructProcess();
  // had
  fHadPhysicsList->ConstructProcess();
  G4cout << "### ABPhysicsList::ConstructProcess is done" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABPhysicsList::SelectPhysicsList(const G4String& name)
{
  // default  Had physics
    fHadPhysicsList = new QGSP_BIC;
    fHadPhysicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABPhysicsList::SetCuts()
{
  SetCutsWithDefault();
}