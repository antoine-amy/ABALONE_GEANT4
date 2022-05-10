
#include "ABStackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

#include "ABRunAction.hh"

#include <vector>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABStackingAction::ABStackingAction()
	: gammaCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABStackingAction::~ABStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
ABStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
	//if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	//{ // particle is optical photon
	//	if (aTrack->GetParentID() > 0)
//		{ // particle is secondary
//			gammaCounter++;
//		}
//	}
	return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABStackingAction::NewStage()
{
	//G4cout << "Number of optical photons produces in this event : " << gammaCounter << G4endl;
	//std::cerr<< " reading: " << photon_count_array.size() << std::endl;
	//photon_count_array.push_back(gammaCounter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABStackingAction::PrepareNewEvent()
{
	//gammaCounter = 0;
	//std::cerr<< " checking the number of photon " << runaction_variable.photon_count.size() << std::endl;
}
