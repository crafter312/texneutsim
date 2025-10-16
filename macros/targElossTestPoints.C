#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void targElossTestPoints() {

	// Read file and TTree
	string path = "";
	string ifname = "sim_alphapn_56MeV_90mm_zeroWidth_neutsigma0-1ns_diamond_targrecon_real_realNeut";
	size_t numentries;
	TFile *file = TFile::Open((path + ifname + ".root").c_str());
	if (!file || file->IsZombie()) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	// Check if tree exists in the file
	{
		TTree *tree = (TTree*)file->Get("t");
		if (!tree) {
			std::cerr << "Tree 't' not found in the file!" << std::endl;
			file->Close();
			return;
		}
		numentries = tree->GetEntries();
		cout << "Total entries in tree: " << numentries << endl;
	}

	// Create TTreeReader
	TTreeReader myReader("t", file);
	TTreeReaderValue<double> dETarg(myReader, "targEloss");
	TTreeReaderValue<double> inthick(myReader, "inthick");
	TTreeReaderValue<double> inthickrecimp(myReader, "inthickrecimp");
	TTreeReaderValue<bool> isFragDet(myReader, "isFragDet");
	TTreeReaderValue<bool> isNeutHit(myReader, "isNeutHit");
	TTreeReaderValue<vector<double>> dETests(myReader, "dETests");

	TCanvas* canv = new TCanvas("c", "Canvas", 800, 600);

	// Skip to desired entry
	for (int i = 0;;) {
		if (!myReader.Next()) break;
		if (!(*isFragDet) || !(*isNeutHit)) continue; // || (((*inthick)-(*inthickrecimp)) > -0.5)
		if (i == 0) break;
		i++; 
	}

	// 12C target thickness in mg/cm^2 (3.026 is copied from Nic's experiment), currently diamond
	float thickness = 17.575;

	// Make TGraph (make sure that n, x, and y all match as expected for how many test points were used)
	size_t n = (*dETests).size();
	double h = thickness / (double)(n - 1);
	double x[n], y[n];
	for (int i = 0; i < n; i++) {
		x[i] = h * (double)i;
		y[i] = (*dETests)[i];
	}

	TGraph g(n, x, y);
	g.SetTitle("Target Eloss Test Points;Reaction position (mg/cm^{2});Energy loss (MeV)");
	g.Draw("AC*");

	// Quadratic fit to target total energy loss function
	TF1 elossFit("quadratic", "([0]*x*x)+([1]*x)+[2]");
	elossFit.SetRange(-0.1, thickness+0.1);
	elossFit.SetParameter(0, 1);
	elossFit.SetParLimits(0, 0, 1);

	g.Fit(&elossFit, "NRQ", "", -0.1, thickness+0.1);
	elossFit.Draw("same");

	// Invert target energy loss function and calculate target reaction position
	//double inthickreconimproved = ((*dETarg) - elossFit.GetParameter(1)) / elossFit.GetParameter(0); // for linear fit

	double a = elossFit.GetParameter(0);
	cout << "a: " << a << endl;
	double b = elossFit.GetParameter(1);
	double c = elossFit.GetParameter(2) - (*dETarg);
	double disc = (b*b) - (4. * a * c);
	if (disc < 0) {
		cout << "Discriminant < 0, something went very wrong!" << endl;
		return;
	}
	double inthickreconimproved = 0.5 * (-b + sqrt(disc)) / a;
	inthickreconimproved = min(max(inthickreconimproved, 0.), (double)thickness); // clamp value to within target dimensions
	cout << "Real target Eloss: " << (*dETarg) << endl;
	cout << "Sampled target reaction position: " << (*inthick) << endl;
	cout << "Simulation reconstructed position: " << (*inthickrecimp) << endl;
	cout << "Macro reconstructed position: " << inthickreconimproved << endl;

	canv->SaveAs("target_eloss_test_points.png");

}
