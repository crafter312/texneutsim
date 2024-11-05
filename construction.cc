#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
  fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");

  fMessenger->DeclareProperty("nCols", nCols, "Number of columns");
  fMessenger->DeclareProperty("nRows", nRows, "Number of rows");
  fMessenger->DeclareProperty("isCherenkov", isCherenkov, "Toggle Cherenkov Detector");
    fMessenger->DeclareProperty("isScintillator", isScintillator, "Toggle Scintillator Detector");


  nCols = 10;
  nRows = 10;

  DefineMaterials();

  isCherenkov=false;
  isScintillator=true;
}

MyDetectorConstruction::~MyDetectorConstruction()
{

}

void MyDetectorConstruction::DefineMaterials()
{
  G4NistManager *nist = G4NistManager::Instance(); // this will call a nist database which can be used to source materials

  // make our first material, fused silica
  SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
  SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
  SiO2->AddElement(nist->FindOrBuildElement("O"), 2);

  // another definition for water
  H2O = new G4Material("H2O", 1.000*g/cm3, 2);
  H2O->AddElement(nist->FindOrBuildElement("H"), 2);
  H2O->AddElement(nist->FindOrBuildElement("O"), 1);

  // definition for carbon
  C = nist->FindOrBuildElement("C");

  // create a compount material made of elements and compounds
  Aerogel = new G4Material("Aerogel", 0.200*g/cm3, 3);
  Aerogel->AddMaterial(SiO2, 62.5*perCent);
  Aerogel->AddMaterial(H2O, 37.4*perCent);
  Aerogel->AddElement(C, 0.1*perCent);

  worldMat = nist->FindOrBuildMaterial("G4_AIR");

  // setup the refractive properties of the aerogel
  // 1.239841939 is the conversion from nm to eV
  G4double energy[2] = {1.239841939*eV/0.2, 1.239841939*eV/0.9};
  G4double rindexAerogel[2] = {1.1,1.1};
  G4double rindexWorld[2] = {1.,1.};

  G4MaterialPropertiesTable *mptAerogel = new G4MaterialPropertiesTable(); // make a material property for the aerogel and add the above properties to the table
  mptAerogel->AddProperty("RINDEX",energy,rindexAerogel,2);
  Aerogel->SetMaterialPropertiesTable(mptAerogel); // now add the material properties table to the material
  
  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX",energy,rindexWorld,2);
  worldMat->SetMaterialPropertiesTable(mptWorld);

  Na = nist->FindOrBuildElement("Na");
  I = nist->FindOrBuildElement("I");
  NaI = new G4Material("NaI",3.67*g/cm3,2);
  NaI->AddElement(Na,1);
  NaI->AddElement(I,1);
}


void MyDetectorConstruction::ConstructCherenkov()
{
  // here we can construction our aerogel
  solidRadiator = new G4Box("solidRadiator",0.4*m, 0.4*m, 0.01*m);
  logicRadiator = new G4LogicalVolume(solidRadiator,Aerogel,"logicRadiator");
  physRadiator = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.25*m),logicRadiator,"physRadiator", logicWorld, false, 0, true);
  fScoringVolume = logicRadiator;

  // now we will make a sensitive detector
  solidDetector = new G4Box("solidDetector",xWorld/nRows,yWorld/nCols,0.01*m);
  logicDetector = new G4LogicalVolume(solidDetector,worldMat,"logicDetector"); // we make it from air but it can still detect photons

  // here we generate many detectors at the same time. The value of rows and
  // columns is initialized in the header file but can be changed at runtime
  // by using the commands of the detector class
  for(G4int i=0; i<nRows; i++)// for loop to produce many detectors
  {
    for(G4int j=0; j<nCols; j++)
    {
      physDetector = new G4PVPlacement(0,
                                            G4ThreeVector(-0.5*m+(i+0.5)*m/nRows,-0.5*m+(j+0.5)*m/nCols,0.49*m),
                                            logicDetector,
                                            "physDetector",
                                            logicWorld,
                                            false,
                                            j+i*nCols, // this part gives a unique identity to each detector and is important for repeated devices
                                            true);
    }
  }
}


void MyDetectorConstruction::ConstructScintillator()
{
  solidScintillator = new G4Tubs("solidScintillator",10*cm,20*cm,30*cm,0*deg,360*deg);

  logicScintillator = new G4LogicalVolume(solidScintillator, NaI, "logicalScintillator");

  physScintillator = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicScintillator,"physScintillator",logicWorld, false, 0, true);

  fScoringVolume = logicScintillator;
 
  VisAttributes();
}


G4VPhysicalVolume *MyDetectorConstruction::Construct()
{

  // size of the world volume
  xWorld = 0.5*m;
  yWorld = 0.5*m;
  zWorld = 0.5*m;


  // Every material in geant4 has 3 parts
  //  The "solid", defines the size
  //  The "logical volume" which includes the material
  //  The "physical volume" places the volume in the geant4 simulation with coordinates and rotations, etc, and interacts with particles
  // The physical volume inherits the shape and logical volume from "logicWorld"
  solidWorld = new G4Box("solidWorld",xWorld,yWorld,zWorld); // creates a volume with half length x=0.5,y=0.5,z=0.5, giving a volumes 1x1x1 m^3
  logicWorld = new G4LogicalVolume(solidWorld,worldMat,"locigWorld");
  physWorld = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicWorld,"physWorld",0,false,0,true);

  if(isCherenkov) ConstructCherenkov();
  if(isScintillator) ConstructScintillator();


  return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{


  // here we make a sensative detector and we tell the logical detector that it is this sensitive detector
  MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

  if(isCherenkov) logicDetector->SetSensitiveDetector(sensDet);
  
}


void MyDetectorConstruction::VisAttributes()
{
  G4VisAttributes *scint_va = new G4VisAttributes(G4Color(0.8,0.4,0.4,0.4));
  scint_va->SetForceSolid(true);
  logicScintillator->SetVisAttributes(scint_va);

}
