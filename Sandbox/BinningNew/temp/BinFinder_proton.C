#define BinFinder_proton_cxx
#include "BinFinder_proton.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void BinFinder_proton::Loop()
{
//   In a ROOT session, you can do:
//      root> .L BinFinder_proton.C
//      root> BinFinder_proton t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	gStyle->SetPalette(kBird);
	gROOT->SetBatch(kTRUE);

	TH2D * h_true_reco_mom = new TH2D("h_true_reco_mom", ";Muon Momentum (MCS) [GeV];Muon Momentum (Truth) [GeV]", 120, 0, 2.5, 120, 0, 2.5);
	double bins_mumom[13] = {0.30, 0.36, 0.41, 0.44, 0.49, 0.53, 0.56, 0.59, 0.63, 0.73, 0.81,1.27, 1.50};

	TH2D * h_true_reco_mom_rightbin = new TH2D("h_true_reco_mom_rightbin", ";Muon Momentum (MCS) [GeV];Muon Momentum (Truth) [GeV]", 6, bins_mumom, 6, bins_mumom);

	TH1D * h_reco_mom_pre_truth = new TH1D("h_reco_mom_pre_truth", ";Muon Momentum (Truth) [GeV];", 60, 0, 4);

	TH1D * h_reco_mom_pre_truth_firstbin = new TH1D("h_reco_mom_pre_truth_firstbin", ";Muon Momentum (Truth) [GeV];", 60, 0, 4);


	Long64_t nentries = fChain->GetEntriesFast();

	double lower_bin = 0.3;
	double upper_bin = 1.5;

	std::vector<double> bin_v;
	bin_v.push_back(lower_bin);

	std::vector<double> points_v, error_x_v, mean_v, std_v;

	bool exit_flag = false;

        bool first = true;

	while (!exit_flag) {

		std::cout << ">>> Using lower_bin = " << lower_bin << ", and upper_bin = " << upper_bin << std::endl;

		Long64_t nbytes = 0, nb = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;

			double true_mom = pmom_true;
			double mcs_mom = pmom_reco;

			//h_true_reco_mom->Fill(true_mom, mcs_mom);
			//h_true_reco_mom_rightbin->Fill(true_mom, mcs_mom);
			if (first) {
			  h_true_reco_mom->Fill(mcs_mom, true_mom);
			  h_true_reco_mom_rightbin->Fill(mcs_mom, true_mom);
                        }

			//if (true_mom > lower_bin && true_mom < upper_bin) {
				//h_reco_mom_pre_truth->Fill(mcs_mom);
			//}
 
			if (mcs_mom > lower_bin && mcs_mom < upper_bin) {
                                h_reco_mom_pre_truth->Fill(true_mom);
			}

			if (first && mcs_mom > 0.3 && mcs_mom < 0.45) {
				h_reco_mom_pre_truth_firstbin->Fill(true_mom);
			}

		}  //end of loop over all the entries
                std::cout<<"Filled the histograms of the recovstrue and reco per truth"<<std::endl;
                double array[2];
                //array is the array of the fitting parameters
		if (first) {
			DoFit(h_reco_mom_pre_truth_firstbin, "test", "test", array, false);
			bin_v.push_back(lower_bin);
			points_v.push_back((lower_bin + 0)/2);
			error_x_v.push_back((lower_bin - 0)/2);
			mean_v.push_back(array[0]);
			std_v.push_back(array[1]);
		}
		first = false;
                //==================================================================================
                //fit the distribution of the whole bin
		DoFit(h_reco_mom_pre_truth, "test", "test", array, false);
		if (lower_bin > 0.29 && lower_bin < 0.31 && upper_bin > 0.99 && upper_bin < 1.1){ //why
			ostringstream oss;
			oss << "Reco. Momentum: " << lower_bin << "-" << upper_bin;
			std::string label = oss.str();
			ostringstream oss2;
			oss2 << "reco_momentum_" << lower_bin << "-" << upper_bin;
			std::string label2 = oss2.str();
			DoFit(h_reco_mom_pre_truth, label, label2, array, true);
		}
		std::cout << "mean " << array[0] << ", std " << array[1] << std::endl;
		std::cout << "reco bin: " << upper_bin - lower_bin << std::endl;
		std::cout << "true bin: " << 2*array[1] << std::endl;
		double truth_bin = upper_bin - lower_bin;
		double reco_bin = 2*array[1];
                //========================================================================================= 
		if (std::abs(truth_bin - reco_bin) < 0.005) {
			bin_v.push_back(upper_bin);
			points_v.push_back((upper_bin + lower_bin)/2);
			error_x_v.push_back((upper_bin - lower_bin)/2);
			mean_v.push_back(array[0]);
			std_v.push_back(array[1]);
			// Rerun fitting just to save the plot
			ostringstream oss;
			oss << "Reco Momentum: " << lower_bin << "-" << upper_bin;
			std::string label = oss.str();
			ostringstream oss2;
			oss2 << "reco_momentum_" << lower_bin << "_" << upper_bin;
			std::string label2 = oss2.str();
			DoFit(h_reco_mom_pre_truth, label, label2, array, true);

			lower_bin = upper_bin;
			upper_bin = lower_bin + 2;
		} else {
			upper_bin-=0.005 ;
		}
		if (lower_bin > 1.0) exit_flag = true;
		h_reco_mom_pre_truth->Reset();
		//exit_flag = true;
	} //end of while loop

	std::cout << "Bins: " << std::endl;
	for (auto b : bin_v)
		std::cout << "\t" << b << std::endl;



	TCanvas * c0 = new TCanvas;
	c0->SetRightMargin(0.13);
	TLine *line = new TLine(0,0,2.5,2.5);
	line->SetLineColor(kRed);
	gStyle->SetPalette(kBird);
	h_true_reco_mom->Draw("colz");
	line->Draw();

	/*((TPaletteAxis*)h_true_reco_mom->GetListOfFunctions()->FindObject("palette"))->SetX1NDC(0.875);
        ((TPaletteAxis*)h_true_reco_mom->GetListOfFunctions()->FindObject("palette"))->SetX1NDC(0.92);
        gPad->Modified();
        gPad->Update();*/
        //TPaletteAxis* palette = (TPaletteAxis*)h_true_reco_mom->GetListOfFunctions()->FindObject("palette");
        //palette->SetX1NDC(0.875);
        //palette->SetX2NDC(0.92);
        //palette->SetY1NDC(0.2);
        //palette->SetY2NDC(0.8);
        //gPad->Modified();
        //gPad->Update();

	TString n = "final_nopoints";

	c0->SaveAs("output_proton/" + n + ".pdf");
	c0->SaveAs("output_proton/" + n + ".C");


	TCanvas * c00 = new TCanvas;
	c00->SetRightMargin(0.13);
	c00->Modified();
        c00->Update();
	TLine *line2 = new TLine(0,0,2.5,2.5);
	line2->SetLineColor(kRed);
	h_true_reco_mom_rightbin->Draw("colz text");
	line2->Draw();

        //TPaletteAxis* palette00 = (TPaletteAxis*)h_true_reco_mom_rightbin->GetListOfFunctions()->FindObject("palette");
        //palette00->SetX1NDC(0.875);
        // palette00->SetX2NDC(0.92);
        //gPad->Modified();
        //gPad->Update();

	n = "final_nopoints_rightbin";

	c00->SaveAs("output_proton/" + n + ".pdf");
	c00->SaveAs("output_proton/" + n + ".C");



	double* points = &points_v[0];
	double* mean_array = &mean_v[0];
	double* error_x = &error_x_v[0];
	double* sigma_array = &std_v[0];
	int pts = points_v.size();

	TCanvas * c = new TCanvas;
	c->SetRightMargin(0.13);
        c->Modified();
        c->Update();

	TGraphErrors* gr = new TGraphErrors(pts,points,mean_array,error_x,sigma_array);
	gr->SetMarkerStyle(kFullCircle);
	gr->SetMarkerSize(1.0);
	gr->SetLineWidth(2.0);

	gStyle->SetPalette(kBird);
	h_true_reco_mom->Draw("colz");
	gr->Draw("P");

	line->Draw();



        TPaletteAxis* palette = (TPaletteAxis*)h_true_reco_mom->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.875);
        palette->SetX2NDC(0.92);
        //palette->SetY1NDC(0.2);
        //palette->SetY2NDC(0.8);
        gPad->Modified();
        gPad->Update();


        //gr->Fit("pol1", "", "", 0.1, 0.9);
        //TF1 *myfunc = gr->GetFunction("pol1");

	n = "final";

	c->SaveAs("output_proton/" + n + ".pdf");
	c->SaveAs("output_proton/" + n + ".C");
}


void BinFinder_proton::DoFit(TH1D * the_histo, std::string name, std::string save_name, double *array, bool save_plot) {

	TCanvas * c = new TCanvas;

	the_histo->Fit("gaus");
	TF1 *myfunc = the_histo->GetFunction("gaus");
	the_histo->Draw("histo");
	myfunc->Draw("same");

	std::cout << "mean: " << myfunc->GetParameter(1) << std::endl;
	std::cout << "sigma: " << myfunc->GetParameter(2) << std::endl;
	array[0] = myfunc->GetParameter(1);
	array[1] = myfunc->GetParameter(2);

	std::string j = name + " GeV";

	TLatex* range = new TLatex(0.88,0.86, j.c_str());
	range->SetTextColor(kGray+2);
	range->SetNDC();
	range->SetNDC();
        //range->SetTextSize(2/30.);
        range->SetTextSize(0.038);
	range->SetTextAlign(32);
	range->Draw();

	TString n = save_name;
	n = n + "_gaus";

	if (save_plot) {
		c->SaveAs("output_proton/" + n + ".pdf");
		c->SaveAs("output_proton/" + n + ".C");
	}

  //_txtfile << name << " GeV |" << array[0] << "|" << array[1] << std::endl;


}










