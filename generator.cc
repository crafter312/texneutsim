#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);// define the number of primary particles created per event

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  G4ParticleDefinition *particle = particleTable->FindParticle("geantino");
  G4ThreeVector pos(0.,0.,0.); // the position of the particle gun
  G4ThreeVector mom(0.,0.,1.); // the particle momentum
  // set particle gun parameters
  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(0.*GeV);
  fParticleGun->SetParticleDefinition(particle);
  
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();

  if(particle==G4Geantino::Geantino())
  {
      G4int Z = 11;
      G4int A = 22;

      G4double charge = 0.*eplus; // give charge if needed
      G4double energy = 0.*keV;

      G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z,A,energy);

      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(charge);
  }

  fParticleGun->GeneratePrimaryVertex(anEvent); // tell geant4 to generate the primary vertex
}
