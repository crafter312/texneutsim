#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);// define the number of primary particles created per event

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  G4ParticleDefinition *particle = particleTable->FindParticle("neutron");
  //G4ParticleDefinition *particle = particleTable->FindParticle("geantino");
  
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
    G4int Z = 98;
    G4int A = 252;

    G4double charge = 0.*eplus; // give charge if needed
    G4double energy = 0.*keV;

    G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z,A,energy);
    
    G4ThreeVector pos(0.,0.,0.); // the position of the particle gun
    G4ThreeVector mom(0.,0.,1.); // the particle momentum

    // set particle gun parameters
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(0.*GeV);
    fParticleGun->SetParticleDefinition(particle);


    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(charge);
  }

  if(particle!=G4Geantino::Geantino())
  {
    fConeAxis = G4ThreeVector(0., 0., 1.);
    fConeAngle = 30.*deg;
    fConeApexRadius = 0.*cm;

    // Generate random direction within the cone
    G4double cosTheta = std::cos(fConeAngle);
    G4double z = G4UniformRand() * (1. - cosTheta) + cosTheta;  // Random z within cone
    G4double phi = G4UniformRand() * 2. * M_PI;                   // Random azimuthal angle
    G4double r = std::sqrt(1. - z * z);
    G4double x = r * std::cos(phi);
    G4double y = r * std::sin(phi);

    // Rotate the direction vector to align with fConeAxis
    G4ThreeVector randomDirection = G4ThreeVector(x,y,z).rotateUz(fConeAxis);

    // Generate random position within the apex radius (if non-zero)
    G4ThreeVector sourcePosition = G4ThreeVector(0.,0.,-0.5*m);

    if (fConeApexRadius > 0.) {
      G4double radius = fConeApexRadius * std::sqrt(G4UniformRand());
      G4double theta = G4UniformRand() * 2. * M_PI;
      sourcePosition = G4ThreeVector(radius * std::cos(theta), radius * std::sin(theta), 0.);
    }

    // Set particle position
    fParticleGun->SetParticlePosition(sourcePosition);
    // Set particle momentum direction
    fParticleGun->SetParticleMomentumDirection(randomDirection);
    
    fParticleGun->SetParticleEnergy(1.0*MeV);
    fParticleGun->SetParticleDefinition(particle);



  }

  fParticleGun->GeneratePrimaryVertex(anEvent); // tell geant4 to generate the primary vertex
}
