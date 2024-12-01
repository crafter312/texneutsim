#include "event.hh"
#include "scinthit.hh"

MyEventAction::MyEventAction(MyRunAction*)
{
  fEdep = 0.;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event *event)
{
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) return;

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitsCollection");
  auto hitsCollection = static_cast<G4THitsCollection<ScintillatorHit>*>(hce->GetHC(hcID));
  if (!hitsCollection) return;


  // loop over all hits and write them out to disk
  for (size_t i = 0; i < hitsCollection->entries(); ++i) {
    ScintillatorHit* hit = (*hitsCollection)[i];
    G4cout << "Time: " << hit->GetTime()
            << ", Energy: " << hit->GetEnergy()
            << ", Particle: " << hit->GetParticleName()
            << ", HitPosition " << hit->GetHitPosition()
            << ", DetPosition " << hit->GetDetectorPosition()
            << ", CopyNo " << hit->GetDetCopyNumber()
            << G4endl;
  }



  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->FillNtupleDColumn(2,1,fEdep);
  man->AddNtupleRow(2);
}

void MyEventAction::AddEdep(G4double edep)
{ 
  fEdep += edep;
}