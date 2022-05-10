/// \file ABEventAction.hh
/// \brief Definition of the ABEventAction class

#ifndef ABEventAction_h
#define ABEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "ABRunAction.hh"
//#include "ABStackingAction.hh"

#include <G4HCofThisEvent.hh>
#include <G4DCofThisEvent.hh>
#include <G4VHitsCollection.hh>
#include <G4VDigiCollection.hh>
#include <G4Event.hh>



class ABEventAction : public G4UserEventAction
{

ABStackingAction ABStackingvariables_event;

public:
	ABEventAction();
	virtual ~ABEventAction();

	virtual void  BeginOfEventAction(const G4Event* event);
	virtual void    EndOfEventAction(const G4Event* event);

private:
	// methods
	//ABCalorHitsCollection* GetHitsCollection(G4int hcID,const G4Event* event) const;
	void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength) const;

	// data members                   
	G4int  fAbsHCID;
	G4int  fCoatingHCID;
	G4int  fringHCID;
	G4int  fbaseHCID;
	G4int fallHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
