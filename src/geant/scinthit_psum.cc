#include "scinthit_psum.hh"

ScintillatorHitPSum::ScintillatorHitPSum()
{
  fEvent = -1;
  fTime = 0.;
  fDetEnergy = 0.;
  fCopyNumber = 0;
  fDetectorPosition = G4ThreeVector();

	fOriginVolumeName = "";
}

ScintillatorHitPSum::~ScintillatorHitPSum()
{}
