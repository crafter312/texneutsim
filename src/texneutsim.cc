#include <iostream>
#include <string>

#include "TSystem.h"
#include "TInterpreter.h"

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

	// Make sure dictionary is properly linked and stuff
	G4cout << "Loading simlib shared library from: " << string(SOFILE) << G4endl;
	G4cout << "Loading simlib ROOT dictionary from: " << string(PCMFILE) << G4endl;
	gSystem->Load(SOFILE);
	TInterpreter::Instance()->AddIncludePath(PCMFILE);

	// Some simulation parameters
	const G4double inch = 2.54*cm;     // custom inch conversion
	G4double flangeDist = 18.969*inch; // distance between downstream side of target frame and upstream side of flange cover

	// NOTE THAT flangeDist NEEDS TO BE CORRECTED FOR THE TINY DISTANCE BETWEEN THE DOWNSTREAM SIDE OF THE TARGET
	// FRAME AND THE CENTER OF THE TARGET, ONCE THAT DISTANCE IS KNOWN. This does not affect the geometry of the
	// simulation, I'm pretty sure, except for the placement of the target relative to everything else will be
	// ever so slightly off.

	float thickness = 0.;                     // 12C target thickness in mg/cm^2 (3.026 is copied from Nic's experiment)
	float thick_cm  = (thickness / 2260.)*cm; // 12C target thickness in cm (divide by graphite density of 2260 mg/cm^3)

	G4RunManager *runManager = new G4RunManager();
	MyDetectorConstruction *detectorConstruction = new MyDetectorConstruction(flangeDist, thick_cm);
	runManager->SetUserInitialization(detectorConstruction);

	// Initialize physics list
	G4VModularPhysicsList* physicsList = new FTFP_INCLXX;
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
	physicsList->ReplacePhysics(new MENATE_R());
	runManager->SetUserInitialization(physicsList);

	/******** CHARGED PARTICLE SIMULATION SETUP ********/

	// Total incoming beam energy in MeV, also used for the Fresco simulation.
	// If this changes, make sure to redo the Fresco simulations!
	double Ebeam = 56;

	// Other default values
	double Ex                 = 5.366;     // excitation energy of parent fragment in MeV
	double gamma              = 0.541;     // width of excited state of parent fragment in MeV (should be 0.541)
	double distanceFromTarget = 90;        // distance of Gobbi from the target in mm
	string suffix             = "alphapn"; // output file suffix
	const bool hasNeutron     = true;      // flag to tell charged particle simulation that the neutron is simulated elsewhere

	// Initialize main simulation class
	Li6sim_alphapn sim(Ebeam, distanceFromTarget, Ex, gamma, suffix);
	sim.AddExtraSuffix("neutsigma0-1ns");

	// See Li6sim.h for default experiment parameters, which can
	// be changed via "Set..." commands as desired here.
	sim.SetEnableExternalNeutron(hasNeutron);
	sim.SetTargetThickness(thickness);
	sim.SetNeutTRes(0.1); // 0.5 by default
	sim.SetUseRealP(true); // false by default, does perfect charged particle reconstruction

	// Complete initialization of simulation class
	sim.Init();

	// Doesn't print absolutely everything, but can be changed if you wish
	sim.PrintSettings();

	// Initialize output manager
	RootOutput output(sim.GetSuffix(), sim.GetNFrags(), hasNeutron);

	// Single threaded setting of user action classes
	MyPrimaryGenerator *generator = new MyPrimaryGenerator(sim, detectorConstruction->GetTexNeutDist(), thick_cm, output);
	runManager->SetUserAction(generator);

	MyRunAction *runAction = new MyRunAction(sim);
	runManager->SetUserAction(runAction);

	MyEventAction *eventAction = new MyEventAction(runAction, sim, detectorConstruction->GetTexNeutDist(), output);
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
