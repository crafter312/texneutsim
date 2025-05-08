{
	gStyle->SetPalette(kBird);
	std::unique_ptr<TFile> myFile(TFile::Open("texneutsim-output_run0.root", "UPDATE"));
	TBrowser b;

	gPad->SetTickx();
	gPad->SetTicky();
}
