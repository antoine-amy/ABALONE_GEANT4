#ifndef ABStackingAction_H
#define ABStackingAction_H 1

#include "globals.hh"
#include "ABRunAction.hh"
#include "G4UserStackingAction.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern std::vector<int> photon_count_array;

class ABStackingAction : public G4UserStackingAction
{
 
public:
	ABStackingAction();
	~ABStackingAction();

public:
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	void NewStage();
	void PrepareNewEvent();
	
private:
	
	G4int gammaCounter;	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
