#ifndef RUN_HH
#define RUN_HH

#include "Li6sim_alphapn.h"

#include "G4UserRunAction.hh"
#include "G4Run.hh"

#include <vector>

class MyRunAction : public G4UserRunAction
{
  public:
    MyRunAction(Li6sim_alphapn&);
    ~MyRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

		G4double GetTime(int i) { return time[i]; }
		G4double GetX(int i) { return detPosX[i]; }
		G4double GetY(int i) { return detPosY[i]; }
		G4double GetZ(int i) { return detPosZ[i]; }

		// Vector management functions
		void Clear();
		void FillVectors(G4int, G4double, G4double, G4int, G4double, G4double, G4double);
		size_t GetSize() { return event.size(); }

	private:
		Li6sim_alphapn& fLi6Sim;

		std::vector<G4int> event;
		std::vector<G4double> time;
		std::vector<G4double> detEnergy;
		std::vector<G4int> copyNumber;
		std::vector<G4double> detPosX;
		std::vector<G4double> detPosY;
		std::vector<G4double> detPosZ;
};

#endif
