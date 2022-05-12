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
/// \file ABsimulation.cc
/// \brief Main program of the AB simulation

#include "ABDetectorConstruction.hh"
#include "ABActionInitialization.hh"
#include "ABSteppingAction.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "time.h"
#include <stdlib.h>
#include <iostream>

#include "QGSP_BIC.hh"
#include "FTFP_BERT.hh"
#include "G4UImanager.hh"
#include "LBE.hh"
#include "QGSP_BIC.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"
//#include "G4DecayPhysics.hh"
#include "ABStackingAction.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4PhysListFactory.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

std::vector<int> photon_sipm_hit;

namespace {
	void PrintUsage() {
		//G4cerr << " Usage: " << G4endl;
		//G4cerr << " ABsimulation [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
		//G4cerr << "   note: -t option is available only for multi-threaded mode."
			//<< G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv){

	if (argc > 7) {
		PrintUsage();
		return 1;
	}

	std::cerr << "     ========================================================== " << std::endl;
	std::cerr << "     ======= A Geant4 Monte Carlo Simulation on ABALONE ======= " << std::endl;
	std::cerr << "     ========================================================== \n" << std::endl;
	
	std::cerr << "                           %%%%%%%%%%%                   /         " << std::endl; //18
	std::cerr << "                       %%%%%%%%%%%%%%%%%%%%             /         " << std::endl; //17
	std::cerr << "                     %%%%%%%%%%%%%%%%%%%%%%%%          /          " << std::endl; //16
	std::cerr << "                    %%%%%%%%%%%%%%%%%%%%%%%%%%        /           " << std::endl; //15
	std::cerr << "                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%      *            " << std::endl; //14
	std::cerr << "                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  " << std::endl; //13
	std::cerr << "                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 " << std::endl; //12
	std::cerr << "                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                " << std::endl; //11
	std::cerr << "               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               " << std::endl; //10
	std::cerr << "               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               " << std::endl; //9
	std::cerr << "              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              " << std::endl; //8
	std::cerr << "              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              " << std::endl; //7
	std::cerr << "              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              " << std::endl; //6
	std::cerr << "              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              " << std::endl; //5
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //4
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //3
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //2
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //1
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //-1
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //-2
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //-3
	std::cerr << "             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             " << std::endl; //-4
	std::cerr << " \n\n\n ----------- WELCOME TO ABALONE Simulation (ver 2.0.1)! ----------- \n" << std::endl;
	std::cerr << " ** Program Developed by: Sze Chun Yiu @ Stockholm University " << std::endl;
	std::cerr << " ** Electric Field Map Calculated via COMSOL multiphysics \n\n" << std::endl;
	std::cerr << " Reference:                            \n\n [1] Ferenc, Daniel, et al. ABALONETM Photosensors for the IceCube experiment. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment (2018). \n" << std::endl;
	std::cerr << " [2] Ferenc, Daniel, Andrew Chang, and Marija Šegedin Ferenc. The Novel ABALONE Photosensor Technology : 4 - Year Long Tests of Vacuum Integrity, Internal Pumpingand Afterpulsing.arXiv preprint arXiv:1703.04546 (2017). \n\n" << std::endl;
	

	G4String macro;
	G4String session;

#ifdef G4MULTITHREADED
    
	G4int nThreads = 4;
#endif

    for (G4int i = 1; i < argc; i = i + 2) {
		if (G4String(argv[i]) == "-m") macro = argv[i+ 1];
		else if (G4String(argv[i]) == "-u") session = argv[i + 1];

#ifdef G4MULTITHREADED
		else if (G4String(argv[i]) == "-t") {
			nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
		}
#endif
		else {
			PrintUsage();
			return 1;
		}
	}

	// Detect interactive mode (if no macro provided) and define UI session
	//
	G4UIExecutive* ui=0;
	if (!macro.size()){
		ui=new G4UIExecutive(argc, argv, session);
	}

#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	if (nThreads > 0) {
		runManager->SetNumberOfThreads(nThreads);
	}
#else
	G4RunManager* runManager=new G4RunManager;
#endif

	// Set mandatory initialization classes
	//
	std::cerr << "Initializing the Detector Consturction..." << std::endl;
	auto detConstruction = new ABDetectorConstruction("resources/sensl-microfc-60035-sm.properties","default");
	runManager->SetUserInitialization(detConstruction);

	std::cerr << "Initializing the Physics List..." << std::endl;
	
	G4VModularPhysicsList* physicsList=new QGSP_BIC;
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

	G4OpticalPhysics* opticalPhysics=new G4OpticalPhysics();
	auto opticalParams=G4OpticalParameters::Instance();
	opticalParams->SetWLSTimeProfile("delta");
	opticalParams->SetScintYieldFactor(1.0);
	opticalParams->SetScintExcitationRatio(0.0);
	opticalParams->SetCerenkovMaxPhotonsPerStep(2000);
	opticalParams->SetCerenkovMaxBetaChange(100.0);
	opticalParams->SetCerenkovTrackSecondariesFirst(true);
	opticalParams->SetScintTrackSecondariesFirst(true);

	//const G4RadioactiveDecay* theRadioactiveDecay=new G4RadioactiveDecay();
	//G4ProcessManager* pmanager = pmanager ->AddProcess(theRadioactiveDecay);


	//physicsList->RegisterPhysics( new G4RadioactiveDecay());
	physicsList->RegisterPhysics(opticalPhysics);
	
	physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics());

	//G4RadioactiveDecay* RadioactiveDecay=new G4RadioactiveDecay();
	//auto RadioactiveDecayParams=G4OpticalParameters::Instance();
	//RadioactiveDecayParams->SetAnalogueMonteCarlo(false);

  	//physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics());

	runManager->SetUserInitialization(physicsList);

	std::cerr << "Initializing the Actions..." << std::endl;
	auto actionInitialization=new ABActionInitialization();
	runManager->SetUserInitialization(actionInitialization);


	runManager->Initialize();

	auto UImanager=G4UImanager::GetUIpointer();
	UImanager->ApplyCommand("/process/verbose 0");
	UImanager->ApplyCommand("/vis/verbose errors");
	UImanager->ApplyCommand("/gps/verbose 0");
	UImanager->ApplyCommand("/control/verbose 0");


	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
	time_t systime = time(NULL);
	long seed = (long)systime;
	CLHEP::HepRandom::setTheSeed(seed);

	CLHEP::HepRandomEngine* theEngine = CLHEP::HepRandom::getTheEngine();
	auto theSeed = theEngine->getSeed();
	G4cout << "Random seed check: " << theSeed << G4endl;
	

	// Initialize visualization
	std::cerr << "Initializing Visualization..." << std::endl;
	auto visManager = new G4VisExecutive;
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	std::cerr << "Initializing User Interface Manager..." << std::endl;
	
	

	// Process macro or start UI session
	//
	if (macro.size()) { //macro.size()
		// batch mode
		G4String command = "/control/execute ";
		G4String mac_file = macro;
		UImanager->ApplyCommand(command + mac_file);
	}
	else {
		// interactive mode : define UI session
		UImanager->ApplyCommand("/control/execute vis.mac");
		if (ui->IsGUI()) {
			UImanager->ApplyCommand("/control/execute gui.mac");
		}
		ui->SessionStart();
		delete ui;
	}

	delete visManager;
	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
