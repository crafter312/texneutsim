#ifndef SCINTHIT_HH
#define SCINTHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

class ScintillatorHit : public G4VHit 
{
  public:
    ScintillatorHit();
    virtual ~ScintillatorHit();

    // Setters
    void SetTime(G4double time) { fTime = time; }
    void SetEnergy(G4double energy) { fEnergy = energy; }
    void SetParticleName(const G4String& name) { fParticleName = name; }
    void SetCopyNumber(G4int number) { fCopyNumber = number; }
    void SetDetectorPosition(const G4ThreeVector& detpos) { fDetectorPosition = detpos; }
    void SetHitPosition(const G4ThreeVector& hitpos) { fHitPosition = hitpos; }

    // Getters
    G4double GetTime() const { return fTime; }
    G4double GetEnergy() const { return fEnergy; }
    G4String GetParticleName() const { return fParticleName; }
    G4int GetDetCopyNumber() const { return fCopyNumber; }
    G4ThreeVector GetDetectorPosition() const { return fDetectorPosition; }
    G4ThreeVector GetHitPosition() const { return fHitPosition; }
    
    

  private:
    G4double fTime;
    G4double fEnergy;
    G4ThreeVector fPosition;
    G4String fParticleName;
    G4int fCopyNumber;
    G4ThreeVector fDetectorPosition;
    G4ThreeVector fHitPosition;
};

#endif