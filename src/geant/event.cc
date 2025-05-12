#include "event.hh"
#include "scinthit.hh"

#include <vector>

MyEventAction::MyEventAction(MyRunAction* _runAction)
{
  fEdep = 0.;
	fHasNeut = false;
	runAction = _runAction;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
	runAction->Clear();
  fEdep = 0.;
	fHasNeut = false;
}

void MyEventAction::EndOfEventAction(const G4Event *event)
{
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) return;

  // get the hits collection
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitsCollection");
  auto hitsCollection = static_cast<G4THitsCollection<ScintillatorHit>*>(hce->GetHC(hcID));
  if (!hitsCollection) return;


  // get the analysis manager
  G4AnalysisManager *man = G4AnalysisManager::Instance();



  // loop over all hits and write them out to disk
	int Nneuts = 0;
  for (size_t i = 0; i < hitsCollection->entries(); ++i) {
    ScintillatorHit* hit = (*hitsCollection)[i];

    // print some useful information
    //G4cout << "Time: " << hit->GetTime()
    //        << ", Energy: " << hit->GetEnergy()
    //        << ", Particle: " << hit->GetParticleName()
    //        //<< ", HitPosition " << hit->GetHitPosition()
    //        //<< ", DetPosition " << hit->GetDetectorPosition()
    //        //<< ", CopyNo " << hit->GetDetCopyNumber()
    //        << G4endl;



    if(hit->GetParticleName() != ""){ 

      man->FillNtupleIColumn(2,0,hit->GetEvent());
      man->FillNtupleIColumn(2,1,hit->GetIsPrimary());
      man->FillNtupleDColumn(2,2,hit->GetTime());
      man->FillNtupleDColumn(2,3,hit->GetInitialEnergy());
      man->FillNtupleDColumn(2,4,hit->GetDetEnergy());
      man->FillNtupleDColumn(2,5,hit->GetHitPosition()[0]);
      man->FillNtupleDColumn(2,6,hit->GetHitPosition()[1]);
      man->FillNtupleDColumn(2,7,hit->GetHitPosition()[2]);
      man->FillNtupleDColumn(2,8,hit->GetDetectorPosition()[0]);
      man->FillNtupleDColumn(2,9,hit->GetDetectorPosition()[1]);
      man->FillNtupleDColumn(2,10,hit->GetDetectorPosition()[2]);
      man->FillNtupleSColumn(2,11,hit->GetParticleName());
      man->FillNtupleIColumn(2,12,hit->GetDetCopyNumber());


      
      man->AddNtupleRow(2);
    }

		if(hit->GetParticleName() == "neutron") {
			Nneuts++;

			// fill vectors with neutron information
			runAction->FillVectors(
				hit->GetIsPrimary(),
				hit->GetTime(),
				hit->GetInitialEnergy(),
				hit->GetDetEnergy(),
				hit->GetHitPosition()[0],
				hit->GetHitPosition()[1],
				hit->GetHitPosition()[2],
				hit->GetDetectorPosition()[0],
				hit->GetDetectorPosition()[1],
				hit->GetDetectorPosition()[2]
			);
		}

  }

	// write hit wise information
	man->FillNtupleIColumn(3,0,(int)(Nneuts>0));
	man->FillNtupleIColumn(3,1,Nneuts);
	man->FillNtupleIColumn(3,2,(int)fHasNeut);
	man->AddNtupleRow(3);

}

void MyEventAction::AddEdep(G4double edep)
{ 
  fEdep += edep;
}
