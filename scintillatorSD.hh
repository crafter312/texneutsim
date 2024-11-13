#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "g4root.hh"
#include "run.hh"

class ScintillatorSD : public G4VSensitiveDetector
{
  public:
    ScintillatorSD(G4String);
    ~ScintillatorSD();

  private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

    G4PhysicsOrderedFreeVector *quEff;
};

#endif