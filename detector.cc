#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
  quEff = new G4PhysicsOrderedFreeVector();

  std::ifstream datafile;
  datafile.open("efficiency.dat");

  while(1)
  {
    G4double wlen, queff;

    datafile>>wlen>>queff;

    if(datafile.eof())
      break;

    //G4cout <<wlen << " " << queff << G4endl;

    quEff->InsertValues(wlen,queff);
  }

  datafile.close();

  quEff->SetSpline(false);

}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{ 
  G4Track *track = aStep->GetTrack(); // this allows us to track our particle in the sensitive detector

  track->SetTrackStatus(fStopAndKill);

  G4StepPoint *preStepPoint = aStep->GetPreStepPoint(); // when the photon enters the detector
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint(); // when the photon leaves the detector

  // get the position of the photon when it hits the detector. 
  G4ThreeVector posPhoton = preStepPoint->GetPosition();
  G4ThreeVector momPhoton = preStepPoint->GetMomentum();

  G4double wlen = (1.239841939*eV/momPhoton.mag())*1E+03; // calculate the wavelength from the momentum

  const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

  // get the detector ID for the detector that was hit
  G4int copyNo = touchable->GetCopyNumber();

  // get the position of the detector that was hit
  G4VPhysicalVolume *physVol = touchable->GetVolume();
  G4ThreeVector posDetector = physVol->GetTranslation();

  //#ifdef G4MULTITHREADED
  ////G4cout << "Photon position: " << posPhoton << G4endl;
  ////G4cout << "Copy number: " << copyNo << G4endl;
  ////G4cout << "Detector position: " << posDetector << G4endl;
  //#endif

  G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->FillNtupleIColumn(0, 0, evt);
  man->FillNtupleDColumn(0, 1, posPhoton[0]);
  man->FillNtupleDColumn(0, 2, posPhoton[1]);
  man->FillNtupleDColumn(0, 3, posPhoton[2]);
  man->FillNtupleDColumn(0, 4, wlen);
  man->AddNtupleRow(0);

  if(G4UniformRand() < quEff->Value(wlen)){
    man->FillNtupleIColumn(1, 0, evt);
    man->FillNtupleDColumn(1, 1, posDetector[0]);
    man->FillNtupleDColumn(1, 2, posDetector[1]);
    man->FillNtupleDColumn(1, 3, posDetector[2]);
    man->AddNtupleRow(1);
  }

  return true;
}