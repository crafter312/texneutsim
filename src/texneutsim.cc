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

	// Distance between downstream side of target frame and center of first layer of TexNeut. This is flangeDist + a 3 inch
	// margin for the beam pipe joint, the scintillator housing wall thickness, and half the scintillator dimension. NOTE 
	// THAT flangeDist MUST BE CHANGED TO REFLECT THE SMALL DISTANCE BETWEEN THE CENTER OF THE TARGET AND THE DOWNSTREAM
	// SIDE OF THE TARGET FRAME, ONCE THAT VALUE IS KNOWN. See construction.cc for detailed TexNeut dimensions.
	G4double texNeutDistance = flangeDist + 3*inch + 1.0635*cm;

	float thickness = 3.026;                  // 12C target thickness in mg/cm^2 (copied from Nic's experiment)
	float thick_cm  = (thickness / 2260.)*cm; // 12C target thickness in cm (divide by graphite density of 2260 mg/cm^3)

	G4RunManager *runManager = new G4RunManager();
	runManager->SetUserInitialization(new MyDetectorConstruction(flangeDist, texNeutDistance + (thick_cm/2)));

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
	double gamma              = 0.;        // width of excited state of parent fragment in MeV 0.541
	double distanceFromTarget = 90;        // distance of Gobbi from the target in mm
	string suffix             = "alphapn"; // output file suffix
	const bool hasNeutron     = true;      // flag to tell charged particle simulation that the neutron is simulated elsewhere

	// Initialize main simulation class
	Li6sim_alphapn sim(Ebeam, distanceFromTarget, Ex, gamma, suffix);
	//sim.AddExtraSuffix("perfTarg_noResolution2");

	// See Li6sim.h for default experiment parameters, which can
	// be changed via "Set..." commands as desired here.
	sim.SetEnableExternalNeutron(hasNeutron);
	sim.SetTargetThickness(thickness);

	// Complete initialization of simulation class
	sim.Init();

	// Doesn't print absolutely everything, but can be changed if you wish
	sim.PrintSettings();

	// Initialize output manager
	RootOutput output(sim.GetSuffix(), sim.GetNFrags(), hasNeutron);

	// Single threaded setting of user action classes
	MyPrimaryGenerator *generator = new MyPrimaryGenerator(sim, texNeutDistance, thick_cm, output);
	runManager->SetUserAction(generator);

	MyRunAction *runAction = new MyRunAction(sim);
	runManager->SetUserAction(runAction);

	MyEventAction *eventAction = new MyEventAction(runAction, sim, texNeutDistance, output);
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
