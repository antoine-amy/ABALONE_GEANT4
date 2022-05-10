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
/// \file ABRunAction.hh
/// \brief Definition of the ABRunAction class

#ifndef ABRunAction_h
#define ABRunAction_h 1

#include "G4UserRunAction.hh"
//#include "G4Accumulable.hh"
#include "globals.hh"
#include "ABStackingAction.hh"
#include <vector>
#include "ABPrimaryGeneratorAction.hh"
#include "G4GenericMessenger.hh"

class G4Run;

class PersistencyHandler;

extern std::vector<int> photon_sipm_hit;

class ABRunAction : public G4UserRunAction
{
private:
    PersistencyHandler* persistencyHandler;
    std::string filename = "/home/amy/ABALONE_GEANT4/results/SiPM/SiPM_readout";
    int file_count;
	G4GenericMessenger* fMessenger;
public:
	ABRunAction();
	virtual ~ABRunAction();

	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);
    PersistencyHandler* getPersistencyHandler() const;


};
#endif

