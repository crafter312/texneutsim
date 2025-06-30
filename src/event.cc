#include "event.hh"

#include "scinthit.hh"
#include "scinthit_psum.hh"

#include "G4RunManager.hh"

#include <limits>
#include <vector>

MyEventAction::MyEventAction(MyRunAction* _runAction, Li6sim_alphapn& sim, G4double dist, RootOutput& out) : fLi6Sim(sim), fOutput(out), fTexNeutDistance(dist)
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

void MyEventAction::BeginOfEventAction(const G4Event *event)
{
  fEdep = 0.;
	fHasNeut = false;
	fFirstFrontInd = -1;
	fMaxEdepInd = -1;
	fMinTimeind = -1; 
	runAction->Clear();

	// Print progress
	G4int index = event->GetEventID();
	G4int Nevents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
	if(index%1000 == 0) G4cout <<	index << " of " << Nevents << G4endl;
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
	int NNeut = 0;
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

		if(hit->GetParticleName() == "neutron") NNeut++;
  }

	// get the summed proton hits collection
  hcID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorProtonSummedHitsCollection");
  auto hitsCollectionPSum = static_cast<G4THitsCollection<ScintillatorHitPSum>*>(hce->GetHC(hcID));
  if (!hitsCollectionPSum) return;

	// loop over summed proton hits
	int NProtSum = 0;
	double detPosZMin = std::numeric_limits<double>::max();
	double EdepMax = 0.;
	double timeMin = std::numeric_limits<double>::max();
  for (size_t i = 0; i < hitsCollectionPSum->entries(); ++i) {
    ScintillatorHitPSum* hit = (*hitsCollectionPSum)[i];
		if(hit->GetDetEnergy() < 0.3) continue; // hit energy is in MeV

		// save event-wise information
		if(hit->GetDetectorPosition()[2]<detPosZMin) {
			detPosZMin = hit->GetDetectorPosition()[2];
			fFirstFrontInd = NProtSum;
		}
		if(hit->GetDetEnergy()>EdepMax) {
			EdepMax = hit->GetDetEnergy();
			fMaxEdepInd = NProtSum;
		}
		if(hit->GetTime()<timeMin) {
			timeMin = hit->GetTime();
			fMinTimeind = NProtSum;
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

		NProtSum++;

		// Diagnostic output
		if(NNeut == 0) {
		//	if(i == 0) G4cout << "No neutron hits, getting proton hit origin information: " << G4endl;
		//	G4cout << "\torigin volume name: " << hit->GetTrackOriginVolumeName() << G4endl;
		//	G4cout << "\tcreation process name: " << hit->GetCreatorProcessName() << G4endl;
		//	G4cout << "\tparent particle ID: " << hit->GetParentID() << G4endl;
		}
	}

	// More diagnostic output
	if((NNeut == 0) && (NProtSum>0)) {
		bool hasOutputHeader = false;
		for (size_t i = 0; i < hitsCollection->entries(); ++i) {
			ScintillatorHit* hit2 = (*hitsCollection)[i];
			if ((hit2->GetParticleName() == "gamma") && (hit2->GetParentID() != -1)) {
				/*if (!hasOutputHeader) {
					G4cout << "Gamma information for production of no-neutron protons:" << G4endl;
					hasOutputHeader = true;
				}
				G4cout << "\torigin volume name: " << hit2->GetTrackOriginVolumeName() << G4endl;
				G4cout << "\tcreation process name: " << hit2->GetCreatorProcessName() << G4endl;
				G4cout << "\tparent particle ID: " << hit2->GetParentID() << G4endl;*/
			}
		}
	}

	// Note that this value is for direct neutron hits, and not for
	// summed proton hits. All the other event-wise information here
	// is related to the summed proton hits, as this is closer to
	// what might be measured in an actual experiment. However, I've
	// saved this value for diagnostic purposes.
	man->FillNtupleIColumn(3,0,NNeut);

	// write hit wise information
	man->FillNtupleIColumn(3,1,(int)(NProtSum>0));
	man->FillNtupleIColumn(3,2,NProtSum);
	man->FillNtupleIColumn(3,3,(int)fHasNeut);
	man->FillNtupleIColumn(3,4,fFirstFrontInd);
	man->FillNtupleIColumn(3,5,fMaxEdepInd);
	man->FillNtupleIColumn(3,6,fMinTimeind);
	man->AddNtupleRow(3);

	// Pass neutron information into charged particle simulation
	if(fMinTimeind>=0)
	fLi6Sim.SetExternalNeutronValues(
		runAction->GetTime(fMinTimeind)/ns,
		runAction->GetX(fMinTimeind)/cm,
		runAction->GetY(fMinTimeind)/cm,
		(runAction->GetZ(fMinTimeind) + fTexNeutDistance)/cm
	);
	fLi6Sim.DoSingleEventPostNeutron(fOutput);
}

void MyEventAction::AddEdep(G4double edep)
{ 
  fEdep += edep;
}
