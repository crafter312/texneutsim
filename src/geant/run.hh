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
		void FillVectors(G4int, G4double, G4double, G4int, G4double, G4double, G4double);
		size_t GetSize() { return event.size(); }

	private:
		std::vector<G4int> event;
		std::vector<G4double> time;
		std::vector<G4double> detEnergy;
		std::vector<G4int> copyNumber;
		std::vector<G4double> detPosX;
		std::vector<G4double> detPosY;
		std::vector<G4double> detPosZ;
};

#endif
