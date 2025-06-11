#include "event.hh"
#include "scinthit.hh"
#include "scinthit_psum.hh"

#include <limits>
#include <vector>

MyEventAction::MyEventAction(MyRunAction* _runAction)
{
  fEdep = 0.;
	fHasNeut = false;
	fFirstFrontInd = -1;
	fMaxEdepInd = -1;
	fMinTimeind = -1;
	runAction = _runAction;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
	fHasNeut = false;
	fFirstFrontInd = -1;
	fMaxEdepInd = -1;
	fMinTimeind = -1; 
	runAction->Clear();
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
  for (size_t i = 0; i < hitsCollection->entries(); ++i) {
    ScintillatorHit* hit = (*hitsCollection)[i];

		if(hit->GetParticleName() == "") continue;

		// save all hits
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

	// get the summed proton hits collection
  hcID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorProtonSummedHitsCollection");
  auto hitsCollectionPSum = static_cast<G4THitsCollection<ScintillatorHitPSum>*>(hce->GetHC(hcID));
  if (!hitsCollectionPSum) return;

	// loop over summed proton hits
	int Nneuts = 0;
	double detPosZMin = std::numeric_limits<double>::max();
	double EdepMax = 0.;
	double timeMin = std::numeric_limits<double>::max();
  for (size_t i = 0; i < hitsCollectionPSum->entries(); ++i) {
    ScintillatorHitPSum* hit = (*hitsCollectionPSum)[i];

		// save event-wise information
		if(hit->GetDetectorPosition()[2]<detPosZMin) {
			detPosZMin = hit->GetDetectorPosition()[2];
			fFirstFrontInd = Nneuts;
		}
		if(hit->GetDetEnergy()>EdepMax) {
			EdepMax = hit->GetDetEnergy();
			fMaxEdepInd = Nneuts;
		}
		if(hit->GetTime()<timeMin) {
			timeMin = hit->GetTime();
			fMinTimeind = Nneuts;
		}

		// save neutron hits in vectors
		runAction->FillVectors(
			hit->GetEvent(),
			hit->GetTime(),
			hit->GetDetEnergy(),
			hit->GetDetCopyNumber(),
			hit->GetDetectorPosition()[0],
			hit->GetDetectorPosition()[1],
			hit->GetDetectorPosition()[2]
		);

		Nneuts++;
	}

	// write hit wise information
	man->FillNtupleIColumn(3,0,(int)(Nneuts>0));
	man->FillNtupleIColumn(3,1,Nneuts);
	man->FillNtupleIColumn(3,2,(int)fHasNeut);
	man->FillNtupleIColumn(3,3,fFirstFrontInd);
	man->FillNtupleIColumn(3,4,fMaxEdepInd);
	man->FillNtupleIColumn(3,5,fMinTimeind);
	man->AddNtupleRow(3);
}

void MyEventAction::AddEdep(G4double edep)
{ 
  fEdep += edep;
}
