#include "scintillatorSD.hh"

ScintillatorSD::ScintillatorSD(const G4String name) : G4VSensitiveDetector(name), fHitsCollection(nullptr)
{
  collectionName.insert("ScintillatorHitsCollection");
	collectionName.insert("ScintillatorProtonSummedHitsCollection");
}


ScintillatorSD::~ScintillatorSD()
{}


void ScintillatorSD::Initialize(G4HCofThisEvent *hce)
{
	// default hits
  fHitsCollection = new G4THitsCollection<ScintillatorHit>(SensitiveDetectorName, collectionName[0]);
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);

	// summed (per crystal) proton hits
	fPSumHitsCollection = new G4THitsCollection<ScintillatorHitPSum>(SensitiveDetectorName, collectionName[1]);
	hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]);
  hce->AddHitsCollection(hcID, fPSumHitsCollection);

  // clear hitID vector
  hitIDs.clear();
}


G4bool ScintillatorSD::ProcessHits(G4Step *step, G4TouchableHistory *ROhist)
{ 
	// Get the event number
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	// Get the current track
	G4Track* track = step->GetTrack();

	// Get the particle name
	auto particleName = track->GetDefinition()->GetParticleName();

	// Get the vector of secondary particles
	const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();

	/******** DIAGNOSTIC INFORMATION ********/

	if (particleName == "neutron") {
		//std::cout << "Neutron interaction detected!" << std::endl;
		const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
		if (process) {
			G4String processName = process->GetProcessName();
			//std::cout << "Process: " << processName << std::endl;

			//if (processName == "neutronInelastic" || processName == "nCapture" || processName == "hadElastic") {
			//  std::cout << "Neutron interacted via:" << std::endl;
			//	process->DumpInfo();
			//}
		}

		// Check if any secondaries are produced in this step
		if (secondaries->size() > 0) {
			//std::cout << "Number of secondaries: " << secondaries->size() << std::endl;
			for (const auto& secondary : *secondaries) {
				auto secondaryName = secondary->GetDefinition()->GetParticleName();
				auto secondaryEnergy = secondary->GetKineticEnergy();
				//std::cout << "Secondary: " << secondaryName 
				//          << ", Energy: " << secondaryEnergy / MeV << " MeV" << std::endl;
			}
		} 
	}
	else if (particleName == "proton") {
		const G4VProcess* process = step->GetPreStepPoint()->GetProcessDefinedStep();
		if (process) {
			//G4String processName = process->GetProcessName();
			//std::cout << "Process: " << processName << std::endl;
		}
	}
	else if (particleName != "opticalphoton") {
		//G4cout << "Particle name: " << particleName << G4endl;
	}

  // If nothing happened during the hit return 0, unless the
	// particle is a gamma (for diagnostics) or if there are secondaries
  G4double energyDeposit = step->GetTotalEnergyDeposit();
  if ((energyDeposit == 0) && (particleName != "gamma") && (secondaries->size() == 0)) return false;

  // define a new hit in the scintillator
  ScintillatorHit* newHit = new ScintillatorHit();

  // Get particle and parent information
  G4int trackID = track->GetTrackID();   // ID of this particle
  G4int parentID = track->GetParentID(); // Parent ID
  //G4Track *decayTrack = track->GetParentTrack(parentID);
  //G4cout << "trackID " << trackID << ", parentID " << parentID  << ", particleName " << particleName << G4endl;
  
  // Get the copy number of the volume
  G4int copyNumber = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber(); // Copy number of this volume

  // Get the position of the current volume
  G4ThreeVector HitPosition = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector DetPosition = fCopyPositions[copyNumber]; 

  // get the position of origination of particle
  G4ThreeVector originationPosition = track->GetVertexPosition();

  // get the time of originiation of particle
  //G4double originationTime = track->GetVertexTime();

  // get the energy of origination of partilce
  G4double originationEnergy = track->GetVertexKineticEnergy(); 

  G4double psptime = step->GetPreStepPoint()->GetGlobalTime();
  G4double gtime = track->GetGlobalTime();

  //if(parentID==0) {
  /*
    G4cout << "ParticleName: " << particleName
              << ", PSP time: " << psptime << " ns"
              << ", track gtime: " << gtime << " ns"
              << ", trackID: " << trackID
              << ", Production Position: " << originationPosition
              //<< ", trackID: " << track->GetTrackID() 
              << ", parentID: " << parentID
            //<< ", Production Time: " << originationTime << " ns"
              << ", Production Energy: " << originationEnergy / CLHEP::MeV << " MeV" << G4endl;
  //}*/


  // write the data from the hit to the hit collection
  newHit->SetEvent(evt);
  newHit->SetIsPrimary(parentID);
  newHit->SetTime(track->GetGlobalTime());
  newHit->SetInitialEnergy(originationEnergy/MeV);
  newHit->SetDetEnergy(energyDeposit/MeV);
  newHit->SetParticleName(particleName);
  newHit->SetCopyNumber(copyNumber);
  newHit->SetHitPosition(HitPosition);
  newHit->SetDetectorPosition(DetPosition);
  //newHit->SetParticleOriginTime();
  //newHit->SetParticleOriginPos();

	// Save gamma diagnostic information
	bool hasSecondaryProton = false;
  for (const auto& secondary : *secondaries) {
		G4String secondaryName = secondary->GetDefinition()->GetParticleName();
		if (secondaryName == "proton") hasSecondaryProton = true;
  }
  if ((particleName == "gamma") && hasSecondaryProton) {
		newHit->SetTrackOriginVolumeName(step->GetTrack()->GetLogicalVolumeAtVertex()->GetName());
		newHit->SetCreatorProcessName(step->GetTrack()->GetCreatorProcess()->GetProcessName());
		newHit->SetParentID(step->GetTrack()->GetParentID());
	}

  // add hit to hit collection
  fHitsCollection->insert(newHit);

  /******** PROTON SUMMING HITS ********/

  // return early if not proton
  if (particleName != "proton") return true;

	// search for already existing hit, increment if found
	for (size_t i = 0; i < hitIDs.size(); i++) {
		if (copyNumber != hitIDs[i]) continue;
		((ScintillatorHitPSum*)fPSumHitsCollection->GetHit(i))->IncrementE(energyDeposit/MeV);
		return true;
	}

  // make new hit if one with the relevant copy number doesn't already exist
	ScintillatorHitPSum* newHitSum = new ScintillatorHitPSum();
	newHitSum->SetEvent(evt);
  newHitSum->SetTime(track->GetGlobalTime());
  newHitSum->SetDetEnergy(energyDeposit/MeV);
  newHitSum->SetCopyNumber(copyNumber);
  newHitSum->SetDetectorPosition(DetPosition);

	// save diagnostic information to hit (this is not transfered to the output file)
	newHitSum->SetTrackOriginVolumeName(step->GetTrack()->GetLogicalVolumeAtVertex()->GetName());
	newHitSum->SetCreatorProcessName(step->GetTrack()->GetCreatorProcess()->GetProcessName());
	newHitSum->SetParentID(step->GetTrack()->GetParentID());

	// Add new hit
	hitIDs.push_back(copyNumber);
	fPSumHitsCollection->insert(newHitSum);

  return true;

  /////////////////////////////////////*

  /*
  G4Track *track = Step->GetTrack(); // this allows us to track our particle in the sensitive detector

  //track->SetTrackStatus(fStopAndKill);

  G4StepPoint *preStepPoint = Step->GetPreStepPoint(); // when the photon enters the detector
  G4StepPoint *postStepPoint = Step->GetPostStepPoint(); // when the photon leaves the detector

  G4double time = preStepPoint->GetGlobalTime();

  const G4VTouchable *touchable = Step->GetPreStepPoint()->GetTouchable();

  // get the detector ID for the detector that was hit
  G4int copyNo = touchable->GetCopyNumber();

  // get the position of the detector that was hit
  G4VPhysicalVolume *physVol = touchable->GetVolume();
  G4ThreeVector posDetector = physVol->GetTranslation();

  //#ifdef G4MULTITHREADED
  ////G4cout << "Copy number: " << copyNo << G4endl;
  ////G4cout << "Detector position: " << posDetector << G4endl;
  //#endif

  G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4AnalysisManager *man = G4AnalysisManager::Instance();

  
  */
  

  // handle optical photons
  /*
  if(track->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
    // get the position of the photon when it hits the detector. 
    G4ThreeVector posPhoton = preStepPoint->GetPosition();
    G4ThreeVector momPhoton = preStepPoint->GetMomentum();

    G4double wlen = (1.239841939*eV/momPhoton.mag())*1E+03; // calculate the wavelength from the momentum

    //#ifdef G4MULTITHREADED
    ////G4cout << "Photon position: " << posPhoton << G4endl;
    //#endif

    man->FillNtupleIColumn(0, 0, evt);
    man->FillNtupleDColumn(0, 1, posPhoton[0]);
    man->FillNtupleDColumn(0, 2, posPhoton[1]);
    man->FillNtupleDColumn(0, 3, posPhoton[2]);
    man->FillNtupleDColumn(0, 4, wlen);
    man->FillNtupleDColumn(0, 5, time);
    man->AddNtupleRow(0);

    if(G4UniformRand() < quEff->Value(wlen)){
      man->FillNtupleIColumn(1, 0, evt);
      man->FillNtupleDColumn(1, 1, posDetector[0]);
      man->FillNtupleDColumn(1, 2, posDetector[1]);
      man->FillNtupleDColumn(1, 3, posDetector[2]);
      man->AddNtupleRow(1);
    }

  }
  */



/*
  // handle neutrons
  if(track->GetParticleDefinition()->GetParticleName() == "neutron"){
    // get the position of the neutron when it hits the detector. 
    G4ThreeVector posNeutron = preStepPoint->GetPosition();
    G4ThreeVector momNeutron = preStepPoint->GetMomentum();

    //#ifdef G4MULTITHREADED
    ////G4cout << "Neutron position: " << posNeutron << G4endl;
    //#endif
    man->FillNtupleIColumn(2, 0, evt);
    man->FillNtupleDColumn(2, 1,0);
    man->FillNtupleDColumn(2, 2, posNeutron[0]);
    man->FillNtupleDColumn(2, 3, posNeutron[1]);
    man->FillNtupleDColumn(2, 4, posNeutron[2]);
    man->FillNtupleDColumn(2, 5, posNeutron[0]);
    man->FillNtupleDColumn(2, 6, posNeutron[1]);
    man->FillNtupleDColumn(2, 7, posNeutron[2]);
    man->FillNtupleDColumn(2,8,0);
    man->FillNtupleDColumn(2,9,0);
    man->FillNtupleDColumn(2,10,0);
    man->AddNtupleRow(2);
    
  }







  // handle gamma rays
  if(track->GetParticleDefinition()->GetParticleName() == "gamma"){
    // get the position of the neutron when it hits the detector. 
    G4ThreeVector posGamma = preStepPoint->GetPosition();
    G4ThreeVector momGamma = preStepPoint->GetMomentum();
    //#ifdef G4MULTITHREADED
    ////G4cout << "Gamma position: " << posNeutron << G4endl;
    //#endif
    man->FillNtupleIColumn(2, 0, evt);
    man->FillNtupleDColumn(2, 1,0);
    man->FillNtupleDColumn(2, 2, posGamma[0]);
    man->FillNtupleDColumn(2, 3, posGamma[1]);
    man->FillNtupleDColumn(2, 4, posGamma[2]);
    man->FillNtupleDColumn(2,5,0);
    man->FillNtupleDColumn(2,6,0);
    man->FillNtupleDColumn(2,7,0);
    man->FillNtupleDColumn(2, 8, posGamma[0]);
    man->FillNtupleDColumn(2, 9, posGamma[1]);
    man->FillNtupleDColumn(2, 10, posGamma[2]);
    man->AddNtupleRow(2);

  }

  

  return true;
  */
}
