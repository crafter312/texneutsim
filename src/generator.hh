#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "rootoutput.h"
#include "Li6sim_alphapn.h"

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

#include "Randomize.hh"

#include <cmath>


class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
  public:
    MyPrimaryGenerator(Li6sim_alphapn&, G4double, G4double, RootOutput&);
    ~MyPrimaryGenerator();

    virtual void GeneratePrimaries(G4Event*);

  private:
		Li6sim_alphapn& fLi6Sim;
		RootOutput& fOutput;

		G4double fTexNeutDistance; // distance between center of target and center of first layer of TexNeut
		G4double targetThickness;  // in cm

    G4ParticleGun *fParticleGun;
    G4ThreeVector fConeAxis;
    G4double fConeAngle;
    G4double fConeApexRadius;
    G4String fParticleName;
};

#endif
