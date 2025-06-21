#include <iostream>
#include <string>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"

#include "FTFP_INCLXX.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "MENATE_R.hh"

#include "construction.hh"
#include "generator.hh"
#include "Li6sim_alphapn.h"
#include "rootoutput.h"
#include "run.hh"
#include "stepping.hh"

using namespace std;

int main(int argc, char** argv) {

	G4RunManager *runManager = new G4RunManager();
	runManager->SetUserInitialization(new MyDetectorConstruction());

	// Initialize physics list
	G4VModularPhysicsList* physicsList = new FTFP_INCLXX;
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
	physicsList->ReplacePhysics(new MENATE_R());
	runManager->SetUserInitialization(physicsList);

	/******** CHARGED PARTICLE SIMULATION SETUP ********/

	// Total incoming beam energy in MeV, also used for the Fresco simulation.
	// If this changes, make sure to redo the Fresco simulations!
	double Ebeam = 63;

	// Other default values
	double Ex                 = 5.366;     // excitation energy of parent fragment in MeV
	double gamma              = 0.541;     // width of excited state of parent fragment in MeV
	double distanceFromTarget = 150;       // distance of Gobbi from the target in mm
	string suffix             = "alphapn"; // output file suffix

	// Initialize main simulation class
	Li6sim_alphapn sim(Ebeam, distanceFromTarget, Ex, gamma, suffix);
	//sim.AddExtraSuffix("perfTarg_noResolution2");

	// See Li6sim.h for default experiment parameters, which can
	// be changed via "Set..." commands as desired here.

	// Complete initialization of simulation class
	sim.Init();

	// Doesn't print absolutely everything, but can be changed if you wish
	sim.PrintSettings();

	// Initialize output manager
	RootOutput output(sim.GetSuffix(), sim.GetNFrags()-1); // -1 because neutron output handled separately

	// Single threaded setting of user action classes
	MyPrimaryGenerator *generator = new MyPrimaryGenerator();
	runManager->SetUserAction(generator);

	MyRunAction *runAction = new MyRunAction(sim);
	runManager->SetUserAction(runAction);

	MyEventAction *eventAction = new MyEventAction(runAction, sim, output);
	runManager->SetUserAction(eventAction);

	MySteppingAction *steppingAction = new MySteppingAction(eventAction);
	runManager->SetUserAction(steppingAction);


	// UI and visualization
	G4UIExecutive* ui = 0;
	if(argc == 1) ui = new G4UIExecutive(argc, argv);

	G4VisManager* visManager = new G4VisExecutive();
	visManager->Initialize();

	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if(ui) {
		UImanager->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
	}
	else {
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command + fileName);
	}

	return 0;
}
