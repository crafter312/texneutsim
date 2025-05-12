#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"

#include "G4AnalysisManager.hh"

#include <vector>

class MyRunAction : public G4UserRunAction
{
  public:
    MyRunAction();
    ~MyRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

		// Vector management functions
		void Clear();
		void FillVectors(G4int, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double);

	private:
		std::vector<G4int> isPrimary;
		std::vector<G4double> time;
		std::vector<G4double> initialEnergy;
		std::vector<G4double> detEnergy;
		std::vector<G4double> hitPosX;
		std::vector<G4double> hitPosY;
		std::vector<G4double> hitPosZ;
		std::vector<G4double> detPosX;
		std::vector<G4double> detPosY;
		std::vector<G4double> detPosZ;
};

#endif
