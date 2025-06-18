#include "run.hh"

#include "G4AblaInterface.hh"
#include "G4AnalysisManager.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterface.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4RunManager.hh"

MyRunAction::MyRunAction(HistoManager& histo) : fHistoManager(histo)
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
  man->CreateNtupleIColumn("fEvent");         // 0 
  man->CreateNtupleIColumn("fIsPrimary");     // 1
  man->CreateNtupleDColumn("fTime");          // 2
  man->CreateNtupleDColumn("fInitialEnergy"); // 3
  man->CreateNtupleDColumn("fDetEnergy");     // 4
  man->CreateNtupleDColumn("fHitPosX");       // 5
  man->CreateNtupleDColumn("fHitPosY");       // 6
  man->CreateNtupleDColumn("fHitPosZ");       // 7
  man->CreateNtupleDColumn("fDetPosX");       // 8
  man->CreateNtupleDColumn("fDetPosY");       // 9
  man->CreateNtupleDColumn("fDetPosZ");       // 10
  man->CreateNtupleSColumn("fParticleName");  // 11
  man->CreateNtupleIColumn("fCopyNumber");    // 12
  man->FinishNtuple(2);

	// fourth root tree (for event-wise parameters)
	man->CreateNtuple("Events","Events");
	man->CreateNtupleIColumn("nMult");                   // 0
	man->CreateNtupleIColumn("isNeutDet");               // 1
	man->CreateNtupleIColumn("pSumMult");                // 2
	man->CreateNtupleIColumn("isNeutInVolume");          // 3
	man->CreateNtupleIColumn("firstFrontInd");           // 4
	man->CreateNtupleIColumn("maxEdepInd");              // 5
	man->CreateNtupleIColumn("minTimeind");              // 6
	man->CreateNtupleIColumn("fEvent", event);           // 7
  man->CreateNtupleDColumn("fTime", time);             // 8
  man->CreateNtupleDColumn("fDetEnergy", detEnergy);   // 9
	man->CreateNtupleIColumn("fCopyNumber", copyNumber); // 10
  man->CreateNtupleDColumn("fDetPosX", detPosX);       // 11
  man->CreateNtupleDColumn("fDetPosY", detPosY);       // 12
  man->CreateNtupleDColumn("fDetPosZ", detPosZ);       // 13
	man->FinishNtuple(3);
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  // Get hold of pointers to the INCL++ model interfaces
	std::vector<G4HadronicInteraction*> interactions = G4HadronicInteractionRegistry::Instance()->FindAllModels(G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName());
	for(std::vector<G4HadronicInteraction*>::const_iterator iInter = interactions.begin(), e = interactions.end(); iInter!=e; ++iInter) {
		G4INCLXXInterface *theINCLInterface = static_cast<G4INCLXXInterface*>(*iInter);
		if (theINCLInterface) {

			// Instantiate the ABLA model
			G4HadronicInteraction* interaction = G4HadronicInteractionRegistry::Instance()->FindModel("ABLA");
			G4AblaInterface* theAblaInterface = static_cast<G4AblaInterface*>(interaction);
			if(!theAblaInterface)
				theAblaInterface = new G4AblaInterface;

			// Couple INCL++ to ABLA
			G4cout << "Coupling INCLXX to ABLA" << G4endl;
			theINCLInterface->SetDeExcitation(theAblaInterface);
		}
	}

	//inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  G4int runID = run->GetRunID();

  std::stringstream strRunID;
  strRunID << runID;

  man->OpenFile("texneutsim-output_run" + strRunID.str() + ".root");

	fHistoManager.Book();
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->Write();
  man->CloseFile();

	fHistoManager.PrintStatistic();
  fHistoManager.Save();
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
