#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH


#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"

#include "detector.hh"


// Here we define a new detector construction which is inherited from the
// G4UserDetectorConstruction class supplied by Geant4.
class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction(); // constructor
    ~MyDetectorConstruction(); // destructor

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct(); // this is the main function that constructs the detector geometry

  private:
    G4Box *solidWorld, *solidRadiator, *solidDetector;
    G4Tubs *solidScintillator;
    G4LogicalVolume *logicWorld, *logicRadiator, *logicDetector, 
                    *logicScintillator;
    G4VPhysicalVolume *physWorld, *physRadiator, *physDetector, 
                      *physScintillator;
    G4int nCols,nRows;

    G4Material *SiO2, *H2O, *Aerogel, *worldMat, *NaI;
    G4Element *C, *Na, *I;

    G4GenericMessenger *fMessenger;

    G4LogicalVolume *fScoringVolume;

    G4double xWorld,yWorld,zWorld;

    bool isCherenkov;
    bool isScintillator;
        
    void DefineMaterials();
    void ConstructCherenkov();
    void ConstructScintillator();
    virtual void ConstructSDandField();
    void VisAttributes();
};

#endif
