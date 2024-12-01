#include "scinthit.hh"


ScintillatorHit::ScintillatorHit()
{
  fTime = 0.;
  fEnergy = 0.;
  fHitPosition = G4ThreeVector();
  fParticleName = "";
  fCopyNumber = 0;
  fDetectorPosition = G4ThreeVector();
  
}

ScintillatorHit::~ScintillatorHit()
{}