#include "run.hh"

MyRunAction::MyRunAction()
{

  G4AnalysisManager *man = G4AnalysisManager::Instance();


  // setup the root  file output 
  // first root tree
  man->CreateNtuple("Photons","Photons");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->CreateNtupleDColumn("fWlen");
  man->CreateNtupleDColumn("fTime");
  man->FinishNtuple(0);

  // second root tree
  man->CreateNtuple("PhotonHits","PhotonHits");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->FinishNtuple(1);

  // third root tree
  man->CreateNtuple("Scoring","Scoring");
  man->CreateNtupleIColumn("fEvent");     // 0 
  man->CreateNtupleIColumn("fIsPrimary"); // 1
  man->CreateNtupleDColumn("fTime");      // 2
  man->CreateNtupleDColumn("fInitialEnergy");    // 3
  man->CreateNtupleDColumn("fDetEnergy");    // 3
  man->CreateNtupleDColumn("fHitPosX");   // 4
  man->CreateNtupleDColumn("fHitPosY");   // 5
  man->CreateNtupleDColumn("fHitPosZ");   // 6
  man->CreateNtupleDColumn("fDetPosX");   // 7
  man->CreateNtupleDColumn("fDetPosY");   // 8
  man->CreateNtupleDColumn("fDetPosZ");   // 9
  man->CreateNtupleSColumn("fParticleName");  // 10
  man->CreateNtupleIColumn("fCopyNumber"); // 11
  man->FinishNtuple(2);

	// fourth root tree (for event-wise parameters)
	man->CreateNtuple("Events","Events");
	man->CreateNtupleIColumn("isNeutDet");                     // 0
	man->CreateNtupleIColumn("neutHitMult");                   // 1
	man->CreateNtupleIColumn("isNeutInVolume");                // 2
	man->CreateNtupleIColumn("firstFrontInd");                 // 3
	man->CreateNtupleIColumn("maxEdepInd");                    // 4
	man->CreateNtupleIColumn("minTimeind");                    // 5
	man->CreateNtupleIColumn("fEvent", event);                 // 6
  man->CreateNtupleDColumn("fTime", time);                   // 7
  man->CreateNtupleDColumn("fDetEnergy", detEnergy);         // 8
	man->CreateNtupleIColumn("fCopyNumber", copyNumber);       // 9
  man->CreateNtupleDColumn("fDetPosX", detPosX);             // 10
  man->CreateNtupleDColumn("fDetPosY", detPosY);             // 11
  man->CreateNtupleDColumn("fDetPosZ", detPosZ);             // 12
	man->FinishNtuple(3);
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  G4int runID = run->GetRunID();

  std::stringstream strRunID;
  strRunID << runID;

  man->OpenFile("texneutsim-output_run" + strRunID.str() + ".root");

}

void MyRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->Write();
  man->CloseFile();

}

void MyRunAction::Clear()
{
	event.clear();
	time.clear();
	detEnergy.clear();
	copyNumber.clear();
	detPosX.clear();
	detPosY.clear();
	detPosZ.clear();
}

// Note the order of the input arguments, make sure it is done correctly if you use this function
void MyRunAction::FillVectors(G4int ev, G4double t, G4double detE, G4int cpy, G4double dPosX, G4double dPosY, G4double dPosZ)
{
	event.push_back(ev);
	time.push_back(t);
	detEnergy.push_back(detE);
	copyNumber.push_back(cpy);
	detPosX.push_back(dPosX);
	detPosY.push_back(dPosY);
	detPosZ.push_back(dPosZ);
}
