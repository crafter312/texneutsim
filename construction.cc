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


  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;


  //***Elements
  H = new G4Element("H", "H", z=1., a=1.01*g/mole);
  C = new G4Element("C", "C", z=6., a=12.01*g/mole);
  N = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  O = new G4Element("O", "O", z=8., a= 16.00*g/mole);

  //***Materials
  //p-Terphenyl
  pTerp = new G4Material("pTerp",density=1.23*g/cm3,2);// this is really pterphenyl, not LXe
  pTerp->AddElement(H,14);
  pTerp->AddElement(C,18);
  //EJ-560 Silicone Optical Pad
  OpticalPadSilicone = new G4Material("OpticalPad_Silicone",density=1.03*g/cm3,2);
  OpticalPadSilicone->AddElement(H,6);
  OpticalPadSilicone->AddElement(C,2);
  //Liquid Xenon
  //fLXe = new G4Material("LXe",z=54.,a=131.29*g/mole,density=3.020*g/cm3);
  //Aluminum
  Al = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);
  //Air
  Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  //Glass
  Glass = new G4Material("Glass", density=1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);
 
  //***Material properties tables


  G4double lxe_Energy[]    = { 7.0*eV , 7.07*eV, 7.14*eV };
  G4int lxenum = sizeof(lxe_Energy)/sizeof(G4double);
  G4double lxe_SCINT[] = { 0.1, 1.0, 0.1 };
  assert(sizeof(lxe_SCINT) == sizeof(lxe_Energy));
  G4double lxe_RIND[]  = { 1.65 , 1.65, 1.65 };
  assert(sizeof(lxe_RIND) == sizeof(lxe_Energy));
  G4double lxe_ABSL[]  = { 10.*cm, 10.*cm, 10.*cm};
  assert(sizeof(lxe_ABSL) == sizeof(lxe_Energy));
  G4MaterialPropertiesTable *pTerp_mt = new G4MaterialPropertiesTable();
  pTerp_mt->AddProperty("FASTCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  pTerp_mt->AddProperty("SLOWCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  pTerp_mt->AddProperty("RINDEX",        lxe_Energy, lxe_RIND,  lxenum);
  pTerp_mt->AddProperty("ABSLENGTH",     lxe_Energy, lxe_ABSL,  lxenum);
  pTerp_mt->AddConstProperty("SCINTILLATIONYIELD",32000./MeV);
  pTerp_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  pTerp_mt->AddConstProperty("FASTTIMECONSTANT",3.3*ns);
  pTerp_mt->AddConstProperty("SLOWTIMECONSTANT",31.*ns);
  pTerp_mt->AddConstProperty("YIELDRATIO",1.0);
  pTerp->SetMaterialPropertiesTable(pTerp_mt);
  // Set the Birks Constant for the LXe scintillator
  pTerp->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
 


  // glass settings
  G4double glass_RIND[]={1.49,1.49,1.49};
  assert(sizeof(glass_RIND) == sizeof(lxe_Energy));
  G4double glass_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(glass_AbsLength) == sizeof(lxe_Energy));
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",lxe_Energy,glass_AbsLength,lxenum);
  glass_mt->AddProperty("RINDEX",lxe_Energy,glass_RIND,lxenum);
  Glass->SetMaterialPropertiesTable(glass_mt);

  // silicone pad settings
  G4double opticalpadsilicone_RIND[]={1.49,1.49,1.49};
  assert(sizeof(opticalpadsilicone_RIND) == sizeof(lxe_Energy));
  G4double opticalpadsilicone_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(opticalpadsilicone_AbsLength) == sizeof(lxe_Energy));
  G4MaterialPropertiesTable *opticalpadsilicone_mt = new G4MaterialPropertiesTable();
  opticalpadsilicone_mt->AddProperty("ABSLENGTH",lxe_Energy,opticalpadsilicone_AbsLength,lxenum);
  opticalpadsilicone_mt->AddProperty("RINDEX",lxe_Energy,opticalpadsilicone_RIND,lxenum);
  OpticalPadSilicone->SetMaterialPropertiesTable(opticalpadsilicone_mt);

  G4double vacuum_Energy[]={2.0*eV,7.0*eV,7.14*eV};
  const G4int vacnum = sizeof(vacuum_Energy)/sizeof(G4double);
  G4double vacuum_RIND[]={1.,1.,1.};
  assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND,vacnum);
  Vacuum->SetMaterialPropertiesTable(vacuum_mt);
  Air->SetMaterialPropertiesTable(vacuum_mt);//Give air the same rindex

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
  G4int imax=16; //16 //Rows
  G4int jmax=8;//8 //Columns 
  G4double xspacing=2.86258*cm;
  G4double zspacing=4.13766*cm;

  G4double xoffset = (imax-1)*xspacing/2.;
  G4double yoffset = 0.;
  G4double zoffset = jmax*zspacing/2.;

  G4double fD_mtl = 0.0635*cm;

  G4double fCube_x = 2.0*cm;
  G4double fCube_y = 2.0*cm;
  G4double fCube_z = 2.0*cm;

  G4double fPad_x = 0.1*cm;
  G4double fPad_y = fCube_y;
  G4double fPad_z = fCube_z;

  G4double fCube_mult = 6;

  G4double fScint_x = fCube_x*fCube_mult + (fCube_mult+1)*fPad_x;
  G4double fScint_y = fCube_y;
  G4double fScint_z = fCube_z;


  G4double housing_x = fScint_x-fPad_x*2;
  G4double housing_y = fScint_y+2.*fD_mtl;
  G4double housing_z = fScint_z+2.*fD_mtl;
 

  ///////////////////////////// Housing
  G4Box* outerBox = new G4Box("Outer Box",housing_x/2.,housing_y/2.,housing_z/2.);
  G4Box* innerBox = new G4Box("Inner Box",housing_x/2.,fScint_y/2.,fScint_z/2.);
  fHousing_box = new G4SubtractionSolid("housing_box",outerBox,innerBox);
  fHousing_log = new G4LogicalVolume(fHousing_box,G4Material::GetMaterial("Al"),"housing_log",0,0,0);
  G4RotationMatrix* housing_rot = new G4RotationMatrix();
  housing_rot->rotateX(0*deg);
  housing_rot->rotateY(0*deg);
  housing_rot->rotateZ(90*deg);	
  G4ThreeVector housing_loc = G4ThreeVector(0.,0.,0.);
    
  //Housing_Multiplier
  //i=Rows, j=Columns
  for(G4int i=0; i<imax; i++){ 
    for(G4int j=0; j<jmax; j++){
      housing_loc.setX(i*xspacing-xoffset);
      housing_loc.setZ(j*zspacing-zoffset);
      new G4PVPlacement(housing_rot,housing_loc,fHousing_log,"phys_Housing", logicWorld,false,i+j*(jmax-1),true);	
    }
  }

  
  // define the housing skin
  G4OpticalSurface *housingSkin = new G4OpticalSurface("housingSkin");
  G4LogicalSkinSurface *housingSurface = new G4LogicalSkinSurface("housingSkin",fHousing_log,housingSkin);
  housingSkin->SetType(dielectric_metal);
  housingSkin->SetFinish(polishedvm2000air);
  housingSkin->SetModel(glisur);


  
  G4double pp[] = {1.0*eV, 10.*eV};
  G4int num = sizeof(pp)/sizeof(G4double);
  G4double reflectivity[] = {1., 1.};
  assert(sizeof(reflectivity) == sizeof(pp));
  G4double efficiency[] = {0.0, 0.0};
  assert(sizeof(efficiency) == sizeof(pp));
  
  G4MaterialPropertiesTable* housingSkinProperty = new G4MaterialPropertiesTable();
  housingSkinProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
  housingSkinProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
  housingSkin->SetMaterialPropertiesTable(housingSkinProperty);

 
  // create the optical pads
  fPad_box = new G4Box("pad_box",fPad_x/2,fPad_y/2,fPad_z/2);
 
  fPad_log = new G4LogicalVolume(fPad_box,G4Material::GetMaterial("OpticalPad_Silicone"),"pad_log",0,0,0);
  for(G4int i=0; i<imax; i++){ 
    for(G4int j=0; j<jmax; j++){
      for(int k=0; k<fCube_mult+1; k++) {
        G4RotationMatrix* Rotation = new G4RotationMatrix();
        Rotation->rotateX(0*deg);
        Rotation->rotateY(0*deg);
        Rotation->rotateZ(90*deg);	
   
        G4ThreeVector pad_loc = G4ThreeVector(i*xspacing-xoffset,(-fScint_x+fPad_x)/2+(fCube_x+fPad_x)*k,j*zspacing-zoffset);
        G4VPhysicalVolume* pad = new G4PVPlacement(Rotation,pad_loc,fPad_log,
                                "pad_phys", logicWorld,
                                false,k+j*(fCube_mult+1)+i*(jmax*(fCube_mult+1)),true);

        //Surface properties for the optical pad
        G4OpticalSurface* padSkin = new G4OpticalSurface("PadSkin");

        G4LogicalSkinSurface* padSurface = new G4LogicalSkinSurface("PadSkin",fPad_log,padSkin);

        padSkin->SetType(dielectric_dielectric);
        padSkin->SetFinish(polished);
        padSkin->SetModel(glisur);

        G4double pp[] = {2.0*eV, 3.5*eV};
        G4int num = sizeof(pp)/sizeof(G4double);
        G4double reflectivity[] = {1., 1.};
        assert(sizeof(reflectivity) == sizeof(pp));
        G4double efficiency[] = {0.0, 0.0};
        assert(sizeof(efficiency) == sizeof(pp));
        
        G4MaterialPropertiesTable *scintWrapProperty = new G4MaterialPropertiesTable();
        scintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
        scintWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
        padSkin->SetMaterialPropertiesTable(scintWrapProperty);
      }
    }
  }

  

  // Create the discrete scintillator cubes
  fScint_box = new G4Box("scint_box",fCube_x/2.,fScint_y/2.,fScint_z/2.); 
  fScint_log = new G4LogicalVolume(fScint_box,G4Material::GetMaterial("pTerp"),"scint_log",0,0,0);
  

  for(G4int i=0; i<imax; i++){ 
    for(G4int j=0; j<jmax; j++){
      for(int k=0; k<fCube_mult; k++){
        //Create a Rotation Matrix
        G4RotationMatrix* Rotation = new G4RotationMatrix();
        Rotation->rotateX(90*deg);
        Rotation->rotateY(0*deg);
        Rotation->rotateZ(0*deg);	
      
        G4ThreeVector cube_loc = G4ThreeVector(i*xspacing-xoffset,-fScint_x/2+fCube_x/2+fPad_x+(fCube_x+fPad_x)*k,j*zspacing-zoffset);
        G4VPhysicalVolume* cube_phys = new G4PVPlacement(0,cube_loc,fScint_log,"scintillator",logicWorld,false,k+j*(fCube_mult+1)+i*(jmax*(fCube_mult+1)),true);

        //Surface properties for the scintillator cubes
        G4OpticalSurface* cubeSkin = new G4OpticalSurface("cubeSkin");

        //G4LogicalSkinSurface* abc = new G4LogicalSkinSurface("cubeSkin",fScint_log,cubeSkin);
        G4LogicalBorderSurface *cubeSurface = new G4LogicalBorderSurface("cubeSkin",cube_phys,physWorld,cubeSkin);

        cubeSkin->SetType(dielectric_dielectric);
        cubeSkin->SetFinish(groundair);
        cubeSkin->SetModel(glisur);

        G4double pp[] = {2.0*eV, 3.5*eV};
        G4int num = sizeof(pp)/sizeof(G4double);
        G4double reflectivity[] = {1., 1.};
        assert(sizeof(reflectivity) == sizeof(pp));
        G4double efficiency[] = {0.0, 0.0};
        assert(sizeof(efficiency) == sizeof(pp));
        
        G4MaterialPropertiesTable* cubeSkinProperty = new G4MaterialPropertiesTable();

        cubeSkinProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
        cubeSkinProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
        cubeSkin->SetMaterialPropertiesTable(cubeSkinProperty);
      }
    }
  } 


  ////////////////////Build PMTs
  G4double innerRadius_pmt = 0.*cm;
  G4double fOuterRadius_pmt = 17.5/2*mm;
  G4double height_pmt = 88./2*mm; //30*fD_mtl/2.;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;
  G4double cathode_depth = 0.4*cm;
  G4double cathode_thickness = 0.02*cm;
 
  //the "photocathode" is a metal slab at the back of the glass that
  //is only a very rough approximation of the real thing since it only
  //absorbs or detects the photons based on the efficiency set below
  fPmt = new G4Tubs("pmt_tube",innerRadius_pmt,fOuterRadius_pmt,
                    height_pmt,startAngle_pmt,spanningAngle_pmt);

  fPhotocath = new G4Tubs("photocath_tube",innerRadius_pmt,fOuterRadius_pmt,
                          cathode_thickness,startAngle_pmt,spanningAngle_pmt);

  fPmt_log = new G4LogicalVolume(fPmt,G4Material::GetMaterial("Glass"),
                                "pmt_log");
  fPhotocath_log = new G4LogicalVolume(fPhotocath,
                                      G4Material::GetMaterial("Al"),
                                      "photocath_log");
  //G4RotationMatrix* Rotation = new G4RotationMatrix();
  //Rotation->rotateX(0*deg);
  //Rotation->rotateY(0*deg);
  //Rotation->rotateZ(90*deg);	
  
  G4ThreeVector photocath_loc = G4ThreeVector(0.,0.,height_pmt-cathode_depth);
  
  new G4PVPlacement(0,photocath_loc,fPhotocath_log,"photocath",fPmt_log,false,0,true);	
  
  // define the photocathode skin
  G4double ephoton[] = {7.0*eV, 7.14*eV};
  num = sizeof(ephoton)/sizeof(G4double);
  G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[]={1.92,1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[]={1.69,1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf= new G4OpticalSurface("photocath_opsurf",glisur,polished,dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);
  new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);
  

  // make the mu-metal shields (in the simulation it's the grey tube)
  G4double fmushield_z = 9.4*cm;
  G4double fmushield_rmax = 2./2*cm;
  G4double fmushield_rmin = 1.892/2*cm;

  ///////////////Mu_Metal Shield Multiplier ////////////
  for(G4int i=0; i<imax; i++){ 
    for(G4int j=0; j<jmax; j++){ 
      for(int k=0; k<2; k++){
        fMuShield_tub = new G4Tubs("MuShield",fmushield_rmin,fmushield_rmax,fmushield_z/2.,0.*deg,360.*deg);
        int sign=0;
        if(k==0) sign=1;
        else sign=-1;
      
        G4ThreeVector fMuMetal_loc = G4ThreeVector(i*xspacing-xoffset,
                                                    -sign*fScint_x+sign*fCube_x-sign*fPad_x*(fCube_mult/2+0.5),
                                                    j*zspacing-zoffset);

        G4RotationMatrix* mu_rotation = new G4RotationMatrix();
        mu_rotation->rotateX(90*deg);
        mu_rotation->rotateY(0*deg);
        mu_rotation->rotateZ(0*deg);
        fMuShield_log = new G4LogicalVolume(fMuShield_tub,G4Material::GetMaterial("Al"),"fMuShield_log",0,0,0);

        G4VPhysicalVolume* fmuMetal_phys = new G4PVPlacement(mu_rotation,fMuMetal_loc,fMuShield_log,"MuShield",
                                logicWorld,false,k+j*(fCube_mult+1)+i*(jmax*(fCube_mult+1)),true);


        // define the mumetal skin
        G4OpticalSurface* muSkin = new G4OpticalSurface("muSkin");

        G4LogicalSkinSurface *muSurface = new G4LogicalSkinSurface("muSkin",fMuShield_log,muSkin);
        muSkin->SetType(dielectric_metal);
        muSkin->SetFinish(polished);
        muSkin->SetModel(glisur);

        G4double pp[] = {1.0*eV, 10.*eV};
        num = sizeof(pp)/sizeof(G4double);
        G4double reflectivity[] = {1., 1.};
        assert(sizeof(reflectivity) == sizeof(pp));
        G4double efficiency[] = {0.0, 0.0};
        assert(sizeof(efficiency) == sizeof(pp));

        G4MaterialPropertiesTable* muSkinProperty 
        = new G4MaterialPropertiesTable();

        muSkinProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
        muSkinProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
        muSkin->SetMaterialPropertiesTable(muSkinProperty);

      }
    }
  }

  ///////////////Arrange pmts around the outside of housing/////////////

  
  /////////////Left PMT/////////////
  G4RotationMatrix* pmt_rm1 = new G4RotationMatrix();
  pmt_rm1->rotateX(90*deg);
  G4RotationMatrix* pmt_rm2 = new G4RotationMatrix();
  pmt_rm2->rotateX(-90*deg);
  G4ThreeVector fPmt_loc1;
  G4ThreeVector fPmt_loc2;
  
   
  ///////////////PMT illustration and Multiplier/////////////
  for(G4int i=0; i<imax; i++){ 
    for(G4int j=0; j<jmax; j++){ 
      
      fPmt_loc1 = G4ThreeVector(i*xspacing-xoffset,-fScint_x/2.-height_pmt,j*zspacing-zoffset);
      fPmt_loc2 = G4ThreeVector(i*xspacing-xoffset, fScint_x/2.+height_pmt,j*zspacing-zoffset); 
      
      new G4PVPlacement(pmt_rm1,fPmt_loc1,fPmt_log,"pmt",logicWorld,
                            false,0+j*(fCube_mult+1)+i*(jmax*(fCube_mult+1)),true);

      new G4PVPlacement(pmt_rm2,fPmt_loc2,fPmt_log,"pmt",logicWorld,
                            false,1+j*(fCube_mult+1)+i*(jmax*(fCube_mult+1)),true);
              
        //fPmtPositions.push_back(fPmt_loc);	
      
    }
  }

  //solidScintillator = new G4Tubs("solidScintillator",10*cm,20*cm,30*cm,0*deg,360*deg);
  //logicScintillator = new G4LogicalVolume(solidScintillator, NaI, "logicalScintillator");
  //physScintillator = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicScintillator,"physScintillator",logicWorld, false, 0, true);
  //fScoringVolume = logicScintillator;
  
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
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8,1.0));
  housing_va->SetForceSolid(false);
  fHousing_log->SetVisAttributes(housing_va);

  G4VisAttributes* scint_va = new G4VisAttributes(G4Colour(0.4,0.4,0.8,0.4));
  scint_va->SetForceSolid(true);
  fScint_log->SetVisAttributes(scint_va);

  G4VisAttributes* cathode_va = new G4VisAttributes(G4Colour(0.8,0.4,0.4));
  cathode_va->SetForceSolid(true);
  fPhotocath_log->SetVisAttributes(cathode_va);

  G4VisAttributes* pmt_va = new G4VisAttributes(G4Colour(0.6,1.0,0.6,0.3));
  pmt_va->SetForceSolid(true);
  fPmt_log->SetVisAttributes(pmt_va);

  //G4VisAttributes* mushield_va = new G4VisAttributes(G4Colour(0.6,0.6,.6,0.3));
  //mushield_va->SetForceSolid(true);
  //fMuShield_log->SetVisAttributes(mushield_va);

  G4VisAttributes* pad_va = new G4VisAttributes(G4Colour(0.4,1.,0.4,0.4));
  pad_va->SetForceSolid(true);
  fPad_log->SetVisAttributes(pad_va);
}
