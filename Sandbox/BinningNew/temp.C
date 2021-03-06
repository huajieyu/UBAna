void temp()
{
//=========Macro generated from canvas: c1_n1094/c1_n1094
//=========  (Mon Feb 19 16:18:24 2018) by ROOT version6.06/06
   TCanvas *c1_n1094 = new TCanvas("c1_n1094", "c1_n1094",452,241,700,502);
   gStyle->SetOptStat(0);
   c1_n1094->Range(-0.5,-46.33125,4.5,416.9813);
   c1_n1094->SetFillColor(10);
   c1_n1094->SetBorderMode(0);
   c1_n1094->SetBorderSize(2);
   c1_n1094->SetFrameLineWidth(2);
   c1_n1094->SetFrameBorderMode(0);
   c1_n1094->SetFrameLineWidth(2);
   c1_n1094->SetFrameBorderMode(0);
   
   TH1D *h_reco_mom_pre_truth__3__1 = new TH1D("h_reco_mom_pre_truth__3__1","",60,0,4);
   h_reco_mom_pre_truth__3__1->SetBinContent(1,3);
   h_reco_mom_pre_truth__3__1->SetBinContent(2,12);
   h_reco_mom_pre_truth__3__1->SetBinContent(3,31);
   h_reco_mom_pre_truth__3__1->SetBinContent(4,55);
   h_reco_mom_pre_truth__3__1->SetBinContent(5,159);
   h_reco_mom_pre_truth__3__1->SetBinContent(6,353);
   h_reco_mom_pre_truth__3__1->SetBinContent(7,323);
   h_reco_mom_pre_truth__3__1->SetBinContent(8,291);
   h_reco_mom_pre_truth__3__1->SetBinContent(9,212);
   h_reco_mom_pre_truth__3__1->SetBinContent(10,91);
   h_reco_mom_pre_truth__3__1->SetBinContent(11,57);
   h_reco_mom_pre_truth__3__1->SetBinContent(12,41);
   h_reco_mom_pre_truth__3__1->SetBinContent(13,32);
   h_reco_mom_pre_truth__3__1->SetBinContent(14,16);
   h_reco_mom_pre_truth__3__1->SetBinContent(15,18);
   h_reco_mom_pre_truth__3__1->SetBinContent(16,5);
   h_reco_mom_pre_truth__3__1->SetBinContent(17,9);
   h_reco_mom_pre_truth__3__1->SetBinContent(18,9);
   h_reco_mom_pre_truth__3__1->SetBinContent(19,5);
   h_reco_mom_pre_truth__3__1->SetBinContent(20,3);
   h_reco_mom_pre_truth__3__1->SetBinContent(21,3);
   h_reco_mom_pre_truth__3__1->SetBinContent(22,2);
   h_reco_mom_pre_truth__3__1->SetBinContent(23,2);
   h_reco_mom_pre_truth__3__1->SetBinContent(25,2);
   h_reco_mom_pre_truth__3__1->SetBinContent(27,1);
   h_reco_mom_pre_truth__3__1->SetBinContent(28,2);
   h_reco_mom_pre_truth__3__1->SetBinContent(30,1);
   h_reco_mom_pre_truth__3__1->SetBinContent(34,1);
   h_reco_mom_pre_truth__3__1->SetBinContent(61,1);
   h_reco_mom_pre_truth__3__1->SetEntries(1740);
   h_reco_mom_pre_truth__3__1->SetStats(0);
   
   TF1 *gaus1 = new TF1("gaus","gaus",0,4);
   gaus1->SetFillColor(10);
   gaus1->SetFillStyle(0);
   gaus1->SetLineColor(2);
   gaus1->SetLineWidth(2);
   gaus1->SetChisquare(162.2047);
   gaus1->SetNDF(25);
   gaus1->GetXaxis()->SetNdivisions(506);
   gaus1->GetXaxis()->SetLabelFont(42);
   gaus1->GetXaxis()->SetTitleSize(0.055);
   gaus1->GetXaxis()->SetTitleOffset(0.8);
   gaus1->GetXaxis()->SetTitleFont(42);
   gaus1->GetYaxis()->SetNdivisions(506);
   gaus1->GetYaxis()->SetLabelFont(42);
   gaus1->GetYaxis()->SetTitleSize(0.055);
   gaus1->GetYaxis()->SetTitleOffset(0.9);
   gaus1->GetYaxis()->SetTitleFont(42);
   gaus1->SetParameter(0,324.7977);
   gaus1->SetParError(0,11.59541);
   gaus1->SetParLimits(0,0,0);
   gaus1->SetParameter(1,0.4501394);
   gaus1->SetParError(1,0.003441674);
   gaus1->SetParLimits(1,0,0);
   gaus1->SetParameter(2,0.1291438);
   gaus1->SetParError(2,0.003274946);
   gaus1->SetParLimits(2,0,2.058316);
   h_reco_mom_pre_truth__3__1->GetListOfFunctions()->Add(gaus1);
   h_reco_mom_pre_truth__3__1->SetLineWidth(2);
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetTitle("Muon Momentum (Truth) [GeV]");
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetNdivisions(506);
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetLabelFont(42);
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetTitleSize(0.055);
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetTitleOffset(0.8);
   h_reco_mom_pre_truth__3__1->GetXaxis()->SetTitleFont(42);
   h_reco_mom_pre_truth__3__1->GetYaxis()->SetNdivisions(506);
   h_reco_mom_pre_truth__3__1->GetYaxis()->SetLabelFont(42);
   h_reco_mom_pre_truth__3__1->GetYaxis()->SetTitleSize(0.055);
   h_reco_mom_pre_truth__3__1->GetYaxis()->SetTitleOffset(0.9);
   h_reco_mom_pre_truth__3__1->GetYaxis()->SetTitleFont(42);
   h_reco_mom_pre_truth__3__1->GetZaxis()->SetNdivisions(506);
   h_reco_mom_pre_truth__3__1->GetZaxis()->SetLabelFont(42);
   h_reco_mom_pre_truth__3__1->GetZaxis()->SetTitleSize(0.055);
   h_reco_mom_pre_truth__3__1->GetZaxis()->SetTitleOffset(0.8);
   h_reco_mom_pre_truth__3__1->GetZaxis()->SetTitleFont(42);
   h_reco_mom_pre_truth__3__1->Draw("histo");
   
   TF1 *gaus2 = new TF1("gaus","gaus",0,4);
   gaus2->SetFillColor(10);
   gaus2->SetFillStyle(0);
   gaus2->SetLineColor(2);
   gaus2->SetLineWidth(2);
   gaus2->SetChisquare(162.2047);
   gaus2->SetNDF(25);
   gaus2->GetXaxis()->SetNdivisions(506);
   gaus2->GetXaxis()->SetLabelFont(42);
   gaus2->GetXaxis()->SetTitleSize(0.055);
   gaus2->GetXaxis()->SetTitleOffset(0.8);
   gaus2->GetXaxis()->SetTitleFont(42);
   gaus2->GetYaxis()->SetNdivisions(506);
   gaus2->GetYaxis()->SetLabelFont(42);
   gaus2->GetYaxis()->SetTitleSize(0.055);
   gaus2->GetYaxis()->SetTitleOffset(0.9);
   gaus2->GetYaxis()->SetTitleFont(42);
   gaus2->SetParameter(0,324.7977);
   gaus2->SetParError(0,11.59541);
   gaus2->SetParLimits(0,0,0);
   gaus2->SetParameter(1,0.4501394);
   gaus2->SetParError(1,0.003441674);
   gaus2->SetParLimits(1,0,0);
   gaus2->SetParameter(2,0.1291438);
   gaus2->SetParError(2,0.003274946);
   gaus2->SetParLimits(2,0,2.058316);
   gaus2->Draw("same");
   TLatex *   tex = new TLatex(0.88,0.86,"Reco Momentum: 0.3-0.56 GeV");
tex->SetNDC();
   tex->SetTextAlign(32);
   tex->SetTextColor(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03789474);
   tex->SetLineWidth(2);
   tex->Draw();
   c1_n1094->Modified();
   c1_n1094->cd();
   c1_n1094->SetSelected(c1_n1094);
}
