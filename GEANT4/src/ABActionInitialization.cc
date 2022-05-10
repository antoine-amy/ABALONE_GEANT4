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
/// \file ABActionInitialization.cc
/// \brief Implementation of the ABActionInitialization class

#include "ABActionInitialization.hh"
#include "ABPrimaryGeneratorAction.hh"
#include "ABRunAction.hh"
#include "ParticleSourceMessenger.hh"
#include "ABEventAction.hh"
#include "ABSteppingAction.hh"
#include "ABStackingAction.hh"
#include "G4SipmUiMessenger.hh"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include "persistency/PersistencyHandlerMessenger.hh"

ABActionInitialization::ABActionInitialization(): G4VUserActionInitialization(){}

ABActionInitialization::~ABActionInitialization(){}

void ABActionInitialization::BuildForMaster() const{
	SetUserAction(new ABRunAction);
}

void ABActionInitialization::Build() const{
    G4SipmUiMessenger::getInstance();
	ParticleSourceMessenger::getInstance();
	PersistencyHandlerMessenger::getInstance();
	SetUserAction(new ABStackingAction);
	SetUserAction(new ABPrimaryGeneratorAction);
	SetUserAction(new ABRunAction);
    SetUserAction(new ABEventAction);
	SetUserAction(new ABSteppingAction);
}