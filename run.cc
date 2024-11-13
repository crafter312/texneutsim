#include "run.hh"

MyRunAction::MyRunAction()
{

  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->CreateNtuple("Photons","Photons");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->CreateNtupleDColumn("fWlen");
  man->CreateNtupleDColumn("fTime");
  man->FinishNtuple(0);

  man->CreateNtuple("PhotonHits","PhotonHits");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->FinishNtuple(1);

  

  man->CreateNtuple("Scoring","Scoring");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fEdep");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->CreateNtupleDColumn("fXneutron");
  man->CreateNtupleDColumn("fYneutron");
  man->CreateNtupleDColumn("fZneutron");
  man->CreateNtupleDColumn("fXgamma");
  man->CreateNtupleDColumn("fYgamma");
  man->CreateNtupleDColumn("fZgamma");
  man->FinishNtuple(2);

}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  G4int runID = run->GetRunID();

  std::stringstream strRunID;
  strRunID << runID;

  man->OpenFile("examplesim-output_run" + strRunID.str() + ".root");

}

void MyRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->Write();
  man->CloseFile();

}