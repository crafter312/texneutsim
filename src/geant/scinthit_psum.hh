#ifndef SCINTHITPSUM_HH
#define SCINTHITPSUM_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

class ScintillatorHitPSum : public G4VHit 
{
  public:
    ScintillatorHitPSum();
    virtual ~ScintillatorHitPSum();

    // Setters
    void SetEvent(const G4int& event) { fEvent = event; }
    void SetTime(const G4double& time) { fTime = time; }
    void SetDetEnergy(const G4double& energy) { fDetEnergy = energy; }
    void SetCopyNumber(const G4int& number) { fCopyNumber = number; }
    void SetDetectorPosition(const G4ThreeVector& detpos) { fDetectorPosition = detpos; }

		void SetTrackOriginVolumeName(const G4String& name) { fOriginVolumeName = name; }
		void SetCreatorProcessName(const G4String& name) { fCreatorProcessName = name; }
		void SetParentID(const G4int& id) { fParentID = id; }

    // Getters
    G4int GetEvent() const { return fEvent; }
    G4double GetTime() const { return fTime; }
    G4double GetDetEnergy() const { return fDetEnergy; }
    G4int GetDetCopyNumber() const { return fCopyNumber; }
    G4ThreeVector GetDetectorPosition() const { return fDetectorPosition; }

		G4String GetTrackOriginVolumeName() const { return fOriginVolumeName; }
		G4String GetCreatorProcessName() const { return fCreatorProcessName; }
		G4int GetParentID() const { return fParentID; }

    // Increment energy
    void IncrementE(G4double de) { fDetEnergy += de; };

  private:
    G4int fEvent;
    G4double fTime;
    G4double fDetEnergy;
    G4int fCopyNumber;
    G4ThreeVector fDetectorPosition;

		// Diagnostic information
		G4String fOriginVolumeName;
		G4String fCreatorProcessName;
		G4int fParentID;
};

#endif
