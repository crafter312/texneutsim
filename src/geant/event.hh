#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsCollection.hh"
#include "G4SDManager.hh"

#include "G4AnalysisManager.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction(MyRunAction*);
    ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void AddEdep(G4double edep);
		void ToggleNeutron() { fHasNeut = true; }

  private:
    G4double fEdep;
		bool fHasNeut;
		G4int fFirstFrontInd; // index in neutron hit vectors of the first hit in the front layer, in hit order
		G4int fMaxEdepInd;    // index in neutron hit vectors of hit with largest energy deposition
		G4int fMinTimeind;    // index in neutron hit vectors of hit with earliest time
		MyRunAction* runAction;
};

#endif
