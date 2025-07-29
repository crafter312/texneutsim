{
	// Load dictionary files
	gSystem->Load("libLI6SIM.so");
	TInterpreter::Instance()->AddIncludePath("libLI6SIM_rdict.pcm");

	gStyle->SetPalette(kBird);
	//std::unique_ptr<TFile> myFile(TFile::Open("texneutsim-output_run0.root", "UPDATE"));
	TBrowser b;

	gPad->SetTickx();
	gPad->SetTicky();
}
