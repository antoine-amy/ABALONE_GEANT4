
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
/// \file ABRunAction.cc
/// \brief Implementation of the ABRunAction class

#include "ABRunAction.hh"
#include "ABStackingAction.hh"

#include "ABAnalysis.hh"

#include "ABPrimaryGeneratorAction.hh"
#include "ABDetectorConstruction.hh"

#include <iostream>
#include <fstream>

// #include "ABRun.hh"
using namespace std;

#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4SipmUiMessenger.hh"
#include "ParticleSourceMessenger.hh"
#include "ABDetectorConstruction.hh"
#include "persistency/PersistencyHandler.hh"
#include "persistency/PersistVisitorFactory.hh"


#include <boost/filesystem.hpp>
#include <boost/lambda/bind.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ofstream myFile;
using namespace boost::filesystem;
using namespace boost::lambda;


ABRunAction::ABRunAction(): G4UserRunAction()
{
    std::string sep="_";
	
	fMessenger = new G4GenericMessenger(this, "/ABALONE/file/", "file name numbering");
	G4GenericMessenger::Command& fcCmd = fMessenger->DeclareProperty("file_number", file_count, "file numbers");
	fcCmd.SetParameterName("num", true);
	fcCmd.SetRange("num>=0");
	fcCmd.SetDefaultValue("1");
	
	persistencyHandler = new PersistencyHandler(PersistVisitorFactory::getInstance()->create(filename+sep+std::to_string(file_count)+sep+std::to_string(theta_source)));
    
	// set printing event number per each event
	//G4RunManager::GetRunManager()->SetPrintProgress(1);
	
	// Create analysis manager
	// The choice of analysis technology is done via selectin of a namespace
	// in ABAnalysis.hh	
	auto analysisManager = G4AnalysisManager::Instance();
	
	//photon_count_array.push_back(1);
	//photon_count_array.push_back(2);

	// Create directories 
	//analysisManager->SetHistoDirectoryName("histograms");
	//analysisManager->SetNtupleDirectoryName("ntuple");
	//analysisManager->SetVerboseLevel(1);
	//analysisManager->SetNtupleMerging(true);
	// Creating ntuple
	//analysisManager->CreateNtuple("ABALONE", "Tracking");
	//analysisManager->CreateH1("Hits","Photon hits", 1000, 0., 2500.);
	//analysisManager->CreateNtupleDColumn("Event_ID");
	//analysisManager->CreateNtupleDColumn("Parent_ID");
	//analysisManager->CreateNtupleDColumn("Track_ID");
	//analysisManager->CreateNtupleDColumn("Particle");
	//analysisManager->CreateNtupleDColumn("X");
	//analysisManager->CreateNtupleDColumn("Y");
	//analysisManager->CreateNtupleDColumn("Z");
	//analysisManager->CreateNtupleDColumn("Time");
	//analysisManager->CreateNtupleDColumn("KE");
	//analysisManager->CreateNtupleDColumn("DE");
	//analysisManager->CreateNtupleDColumn("Volume");
	//analysisManager->FinishNtuple();

	//std::cerr<< "Start of Run Action is Ready" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABRunAction::~ABRunAction()
{
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	
	
	G4cout<< "staring run " <<G4endl;
    std::string sepp = "_";
	persistencyHandler->open(filename+ sepp +std::to_string(nPrimeries)+ sepp+ std::to_string((int)theta_source)+sepp+"run_"+std::to_string(file_count));
    
	// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();
	// Open an output file
	G4String fileName = "../results/tracking/"+std::to_string(nPrimeries)+sepp+std::to_string((int)theta_source)+sepp+"track_"+std::to_string(file_count)+".root";
	analysisManager->OpenFile(fileName);
	//std::cerr<< fileName << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABRunAction::EndOfRunAction(const G4Run* /*run*/)
{
	//print histogram statistics
	auto analysisManager = G4AnalysisManager::Instance();
	// save histograms & ntuple
	analysisManager->Write();
	analysisManager->CloseFile();
	
    //persistencyHandler->persist(G4SipmUiMessenger::getInstance());

    //const ABDetectorConstruction* detector = (const ABDetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
	//persistencyHandler->persist(detector->getSipmHousing()->getSipm()->getModel());
	//persistencyHandler->persist(detector->getSipmHousing()->getSipm()->getModel()->getVoltageTraceModel());
    
    // Close output.
	persistencyHandler->close();
    	file_count++;
	//std::cerr << " at the end: " << photon_count_array.size() << std::endl; 
	//myFile.open("photon_number.csv");
	//for (int i = 0; i < photon_count_array.size(); i++) {
	//	myFile << photon_count_array[i] << "," << photon_sipm_hit[i] << endl;
	//}
	//myFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PersistencyHandler* ABRunAction::getPersistencyHandler() const {
	return persistencyHandler;
}

