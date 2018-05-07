#ifndef __BASE_PLOTTINGTOOLS_CXX__
#define __BASE_PLOTTINGTOOLS_CXX__

#include "PlottingTools.h"

namespace Base {

  void PlottingTools::DrawPreliminary() {
    TLatex* prelim = new TLatex(0.94,0.93, "MicroBooNE Preliminary");
    prelim->SetTextFont(62);
    prelim->SetTextColor(kGray+2);
    prelim->SetNDC();
    prelim->SetTextSize(1/30.);
    prelim->SetTextAlign(32);
    //prelim->SetTextSize(0.04631579);
    prelim->Draw();
  }

  void PlottingTools::DrawPreliminaryXSec() {
    TLatex* prelim = new TLatex(0.9,0.93, "MicroBooNE Preliminary");
    prelim->SetTextFont(62);
    prelim->SetTextColor(kGray+2);
    prelim->SetNDC();
    prelim->SetTextSize(1/30.);
    prelim->SetTextAlign(32);
    //prelim->SetTextSize(0.04631579);
    prelim->Draw();
  }

  void PlottingTools::DrawSimulation() {
    TLatex* prelim = new TLatex(0.9,0.93, "MicroBooNE Simulation");
    prelim->SetTextColor(kGray+1);
    prelim->SetNDC();
    prelim->SetTextSize(2/30.);
    prelim->SetTextAlign(32);
    prelim->SetTextSize(0.04631579);
    prelim->Draw();
  }

  void PlottingTools::DrawOverlay() {
    TLatex* prelim = new TLatex(0.94,0.93, "BNB#nu MC + Cosmic Data Overlay");
    prelim->SetTextFont(62);
    prelim->SetTextColor(kRed+1);
    prelim->SetNDC();
    prelim->SetTextSize(1/30.);
    prelim->SetTextAlign(32);
    //prelim->SetTextSize(0.04631579);
    prelim->Draw();
  }

	void PlottingTools::DrawPOT(double pot) {
  
  std::stringstream sstm2;
  sstm2 << "Accumulated POT: " << pot;
  std::string str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str()); 
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}





void PlottingTools::DrawSimPOT(double pot, double target) {
  
 std::stringstream sstm;
  sstm << "Simulated POT: " << pot;
  std::string str = sstm.str();
  
  TLatex* pot_latex = new TLatex(.10, .96, str.c_str());
  pot_latex->SetTextColor(kGray+2);
  pot_latex->SetNDC();
  pot_latex->SetTextSize(1/30.);
  pot_latex->SetTextAlign(10); //left adjusted
  pot_latex->Draw();
  
  
  std::stringstream sstm2;
  sstm2 << "Scaled to POT: " << target;
  str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}

void PlottingTools::DrawPOTRatio(double pot) {
  
  std::stringstream sstm2;
  sstm2 << "Accumulated POT: " << pot;
  std::string str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.13, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}




TLegend* PlottingTools::DrawTHStack(THStack *hs_trklen,
                   double pot_scaling,
                   bool _breakdownPlots,
                   std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    if (iter.second == NULL || iter.first == "intimecosmic" || iter.first == "beam-off") continue;
    iter.second->Scale(pot_scaling);
  }
  
  if (themap["beam-off"] != NULL) {
    themap["beam-off"]->SetLineColor(kBlue+2);
    themap["beam-off"]->SetFillColor(kBlue+2);
    themap["beam-off"]->SetFillStyle(3004);
    hs_trklen->Add(themap["beam-off"]);
    themap["total"]->Add(themap["beam-off"]);
  }

  if (themap["intimecosmic"] != NULL) {
    themap["intimecosmic"]->SetLineColor(kBlue+2);
    themap["intimecosmic"]->SetFillColor(kBlue+2);
    themap["intimecosmic"]->SetFillStyle(3004);
    hs_trklen->Add(themap["intimecosmic"]);
  }

  if (_breakdownPlots) {
    themap["cosmic_nostopmu"]->SetLineColor(kBlue+2);
    themap["cosmic_nostopmu"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic_nostopmu"]);
    themap["cosmic_stopmu"]->SetLineColor(kBlue);
    themap["cosmic_stopmu"]->SetFillColor(kBlue);
    hs_trklen->Add(themap["cosmic_stopmu"]);
    themap["outfv_nostopmu"]->SetLineColor(kGreen+3);
    themap["outfv_nostopmu"]->SetFillColor(kGreen+3);
    hs_trklen->Add(themap["outfv_nostopmu"]);
    themap["outfv_stopmu"]->SetLineColor(kGreen+2);
    themap["outfv_stopmu"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv_stopmu"]);
    themap["nc_proton"]->SetLineColor(kGray+2);
    themap["nc_proton"]->SetFillColor(kGray+2);
    hs_trklen->Add(themap["nc_proton"]);
    themap["nc_pion"]->SetLineColor(kGray+1);
    themap["nc_pion"]->SetFillColor(kGray+1);
    hs_trklen->Add(themap["nc_pion"]);
    themap["nc_other"]->SetLineColor(kGray);
    themap["nc_other"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc_other"]);
  }
  else {
    themap["cosmic"]->SetLineColor(kBlue+2);
    themap["cosmic"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic"]);
    themap["outfv"]->SetLineColor(kGreen+2);
    themap["outfv"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv"]);
    themap["nc"]->SetLineColor(kGray);
    themap["nc"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc"]);
  }
  
  themap["anumu"]->SetLineColor(kOrange-3);
  themap["anumu"]->SetFillColor(kOrange-3);
  hs_trklen->Add(themap["anumu"]);
  themap["nue"]->SetLineColor(kMagenta+1);
  themap["nue"]->SetFillColor(kMagenta+1);
  hs_trklen->Add(themap["nue"]);
  
  if (_breakdownPlots) {
    themap["signal_nostopmu"]->SetLineColor(kRed+2);
    themap["signal_nostopmu"]->SetFillColor(kRed+2);
    hs_trklen->Add(themap["signal_nostopmu"]);
    themap["signal_stopmu"]->SetLineColor(kRed);
    themap["signal_stopmu"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal_stopmu"]);
  }
  else {
    themap["signal"]->SetLineColor(kRed);
    themap["signal"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal"]);
  }
  hs_trklen->Draw("hist");
  
  //h_trklen_total->DrawCopy("hist");
  
  //gStyle->SetHatchesLineWidth(1);
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");
  
  
  
  
  
  
  
  TLegend* leg2;
  if (_breakdownPlots){
    // leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
    leg2 = new TLegend(0.6015038,0.3101075,0.9235589,0.8468817,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breakdownPlots) {
    sstm << "#nu_{#mu} CC (stopping #mu), " << std::setprecision(2)  << themap["signal_stopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["signal_stopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    sstm << "#nu_{#mu} CC (other), " << std::setprecision(2)  << themap["signal_nostopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["signal_nostopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["signal_stopmu"],"#nu_{#mu} CC (stopping #mu)","f");
    // leg2->AddEntry(themap["signal_nostopmu"],"#nu_{#mu} CC (other)","f");
  } else {
    sstm << "#nu_{#mu} CC (signal), " << std::setprecision(2)  << themap["signal"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["signal"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  // nue
  sstm << "#nu_{e}, #bar{#nu}_{e} CC, " << std::setprecision(2)  << themap["nue"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["nue"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // anumu
  sstm << "#bar{#nu}_{#mu} CC, " << std::setprecision(2)  << themap["anumu"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["anumu"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // nc, outfv, cosmic
  if (_breakdownPlots) {
    sstm << "NC (other), " << std::setprecision(2)  << themap["nc_other"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc_other"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["nc_other"],"NC (other)","f");

    sstm << "NC (pion), " << std::setprecision(2)  << themap["nc_pion"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc_pion"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["nc_pion"],"NC (pion)","f");

    sstm << "NC (proton), " << std::setprecision(2)  << themap["nc_proton"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc_proton"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["nc_proton"],"NC (proton)","f");

    sstm << "OUTFV (stopping #mu), " << std::setprecision(2)  << themap["outfv_stopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["outfv_stopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["outfv_stopmu"],"OUTFV (stopping #mu)","f");

    sstm << "OUTFV (other), " << std::setprecision(2)  << themap["outfv_nostopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["outfv_nostopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["outfv_nostopmu"],"OUTFV (other)","f");

    sstm << "Cosmic (stopping #mu), " << std::setprecision(2)  << themap["cosmic_stopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["cosmic_stopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["cosmic_stopmu"],"Cosmic (stopping #mu)","f");

    sstm << "Cosmic (other), " << std::setprecision(2)  << themap["cosmic_nostopmu"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["cosmic_nostopmu"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["cosmic_nostopmu"],"Cosmic (other)","f");

    if (themap["intimecosmic"] != NULL) {
      // leg2->AddEntry(themap["intimecosmic"],"In-time cosmics","f");
      sstm << "n-time cosmics, " << std::setprecision(2)  << themap["intimecosmic"]->Integral() / themap["total"]->Integral()*100. << "%";
      leg2->AddEntry(themap["intimecosmic"],sstm.str().c_str(),"f");
      sstm.str("");
    }
  } else {
    sstm << "NC, " << std::setprecision(2)  << themap["nc"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "OUTFV, " << std::setprecision(2)  << themap["outfv"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["outfv"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "Cosmic, " << std::setprecision(2)  << themap["cosmic"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["cosmic"],sstm.str().c_str(),"f");
    sstm.str("");
  }

  if (themap["beam-off"] != NULL) {
    sstm << "Data (Beam-off), " << std::setprecision(2)  << themap["beam-off"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["beam-off"],sstm.str().c_str(),"f");
    sstm.str("");
    // leg2->AddEntry(themap["beam-off"],"Data (Beam-off)","f");
  }
  
  leg2->Draw();
  
  return leg2;

}



//**********************************************************

TLegend* PlottingTools::DrawTHStack2D(THStack *hs_trklen,
                   double pot_scaling,
                   bool _breakdownPlots,
                   std::map<std::string,TH2D*> themap){
  
  
  for (auto iter : themap) {
    if (iter.second == NULL || iter.first == "intimecosmic" || iter.first == "beam-off") continue;
    iter.second->Scale(pot_scaling);
    iter.second->SetLineWidth(1);
  }
  
  if (themap["beam-off"] != NULL) {
    themap["beam-off"]->SetLineColor(kBlack);
    themap["beam-off"]->SetFillColor(kOrange);
    //themap["beam-off"]->SetFillStyle(3004);
    hs_trklen->Add(themap["beam-off"]);
  }

  if (themap["intimecosmic"] != NULL) {
    themap["intimecosmic"]->SetLineColor(kBlack);
    themap["intimecosmic"]->SetFillColor(kBlue+2);
    themap["intimecosmic"]->SetFillStyle(3004);
    hs_trklen->Add(themap["intimecosmic"]);
  }

  if (false/*_breakdownPlots*/) {
    themap["cosmic_nostopmu"]->SetLineColor(kBlack);
    themap["cosmic_nostopmu"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic_nostopmu"]);
    themap["cosmic_stopmu"]->SetLineColor(kBlack);
    themap["cosmic_stopmu"]->SetFillColor(kBlue);
    hs_trklen->Add(themap["cosmic_stopmu"]);
    themap["outfv_nostopmu"]->SetLineColor(kBlack);
    themap["outfv_nostopmu"]->SetFillColor(kGreen+3);
    hs_trklen->Add(themap["outfv_nostopmu"]);
    themap["outfv_stopmu"]->SetLineColor(kBlack);
    themap["outfv_stopmu"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv_stopmu"]);
    themap["nc_proton"]->SetLineColor(kBlack);
    themap["nc_proton"]->SetFillColor(kGray+2);
    hs_trklen->Add(themap["nc_proton"]);
    themap["nc_pion"]->SetLineColor(kBlack);
    themap["nc_pion"]->SetFillColor(kGray+1);
    hs_trklen->Add(themap["nc_pion"]);
    themap["nc_other"]->SetLineColor(kBlack);
    themap["nc_other"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc_other"]);
  }
  else {
    themap["cosmic"]->SetLineColor(kBlack);
    themap["cosmic"]->SetFillColor(kBlue+2);
    hs_trklen->Add(themap["cosmic"]);
    themap["outfv"]->SetLineColor(kBlack);
    themap["outfv"]->SetFillColor(kGreen+2);
    hs_trklen->Add(themap["outfv"]);
    themap["nc"]->SetLineColor(kBlack);
    themap["nc"]->SetFillColor(kGray);
    hs_trklen->Add(themap["nc"]);
  }
  
  themap["anumu"]->SetLineColor(kBlack);
  themap["anumu"]->SetFillColor(kOrange-3);
  hs_trklen->Add(themap["anumu"]);
  themap["nue"]->SetLineColor(kBlack);
  themap["nue"]->SetFillColor(kMagenta+1);
  hs_trklen->Add(themap["nue"]);
  
  if (false/*_breakdownPlots*/) {
    themap["signal_nostopmu"]->SetLineColor(kBlack);
    themap["signal_nostopmu"]->SetFillColor(kRed+2);
    hs_trklen->Add(themap["signal_nostopmu"]);
    themap["signal_stopmu"]->SetLineColor(kBlack);
    themap["signal_stopmu"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal_stopmu"]);
  }
  else {
    themap["signal"]->SetLineColor(kBlack);
    themap["signal"]->SetFillColor(kRed);
    hs_trklen->Add(themap["signal"]);
  }



  hs_trklen->Draw("lego1");
  
  //h_trklen_total->DrawCopy("hist");
  
  //gStyle->SetHatchesLineWidth(1);
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  //themap["total"]->Draw("E2 same");
  
  
  
  
  
  
  
  TLegend* leg2;
  if (_breakdownPlots){
    leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breakdownPlots) {
    leg2->AddEntry(themap["signal_stopmu"],"#nu_{#mu} CC (stopping #mu)","f");
    leg2->AddEntry(themap["signal_nostopmu"],"#nu_{#mu} CC (other)","f");
  } else {
    sstm << "#nu_{#mu} CC (signal), " << std::setprecision(2)  << themap["signal"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["signal"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  // nue
  sstm << "#nu_{e}, #bar{#nu}_{e} CC, " << std::setprecision(2)  << themap["nue"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["nue"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // anumu
  sstm << "#bar{#nu}_{#mu} CC, " << std::setprecision(2)  << themap["anumu"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["anumu"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // nc, outfv, cosmic
  if (_breakdownPlots) {
    leg2->AddEntry(themap["nc_other"],"NC (other)","f");
    leg2->AddEntry(themap["nc_pion"],"NC (pion)","f");
    leg2->AddEntry(themap["nc_proton"],"NC (proton)","f");
    leg2->AddEntry(themap["outfv_stopmu"],"OUTFV (stopping #mu)","f");
    leg2->AddEntry(themap["outfv_nostopmu"],"OUTFV (other)","f");
    leg2->AddEntry(themap["cosmic_stopmu"],"Cosmic (stopping #mu)","f");
    leg2->AddEntry(themap["cosmic_nostopmu"],"Cosmic (other)","f");
    if (themap["intimecosmic"] != NULL) {
      leg2->AddEntry(themap["intimecosmic"],"In-time cosmics","f");
    }
  } else {
    sstm << "NC, " << std::setprecision(2)  << themap["nc"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["nc"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "OUTFV, " << std::setprecision(2)  << themap["outfv"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["outfv"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "Cosmic, " << std::setprecision(2)  << themap["cosmic"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["cosmic"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  //leg2->AddEntry(themap["total"],"MC Stat Unc.","f");
  //leg2->Draw();
  
  return leg2;

}



//**********************************************************
TLegend* PlottingTools::DrawTHStack2(THStack *hs_trklen,
                  double pot_scaling,
                  bool _breakdownPlots,
                  std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    if (iter.second == NULL || iter.first == "intimecosmic" || iter.first == "beam-off") continue;
    iter.second->Scale(pot_scaling);
  }
  
  if (themap["beam-off"] != NULL) {
    themap["beam-off"]->SetLineColor(kBlue+2);
    themap["beam-off"]->SetFillColor(kBlue+2);
    themap["beam-off"]->SetFillStyle(3004);
    themap["total"]->Add(themap["beam-off"]);
    hs_trklen->Add(themap["beam-off"]);
  }
  
  themap["background"]->SetLineColor(kBlue+2);
  themap["background"]->SetFillColor(kBlue+2);
  hs_trklen->Add(themap["background"]);
  
  themap["signal"]->SetLineColor(kRed+2);
  themap["signal"]->SetFillColor(kRed+2);
  hs_trklen->Add(themap["signal"]);
  
  hs_trklen->Draw("hist");
  
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");
  
  
  
  
  TLegend* leg2;
  // leg2 = new TLegend(0.13,0.69,0.45,0.87,NULL,"brNDC");
  leg2 = new TLegend(0.1582915,0.6798623,0.5,0.8760757,NULL,"brNDC");

  std::stringstream sstm;
  
  sstm << "Signal, " << std::setprecision(2)  << themap["signal"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["signal"],sstm.str().c_str(),"f");
  sstm.str("");
  
  if (themap["beam-off"] != NULL) {
    sstm << "Background (MC), " << std::setprecision(2)  << themap["background"]->Integral() / themap["total"]->Integral()*100. << "%";
  } else {
    sstm << "Background, " << std::setprecision(2)  << themap["background"]->Integral() / themap["total"]->Integral()*100. << "%";
  }
  leg2->AddEntry(themap["background"],sstm.str().c_str(),"f");
  sstm.str("");

  if (themap["beam-off"] != NULL) {
    sstm << "Background (Off-Beam), " << std::setprecision(2)  << themap["beam-off"]->Integral() / themap["total"]->Integral()*100. << "%";
    leg2->AddEntry(themap["beam-off"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  //leg2->AddEntry(themap["total"],"MC Stat Unc.","f");
  //leg2->Draw();
  
  return leg2;

}

//**********************************************************
TLegend* PlottingTools::DrawTHStack3(THStack *hs_trklen,
                      double pot_scaling,
                      bool _breakdownPlots,
                      std::map<std::string,TH1D*> themap){
  
  
  for (auto iter : themap) {
    iter.second->Scale(pot_scaling);
  }
  
  themap["proton"]->SetLineColor(kBlue+2);
  themap["proton"]->SetFillColor(kBlue+2);
  hs_trklen->Add(themap["proton"]);
  
  themap["pion"]->SetLineColor(kGreen+2);
  themap["pion"]->SetFillColor(kGreen+2);
  hs_trklen->Add(themap["pion"]);
  
  themap["photon"]->SetLineColor(kOrange-3);
  themap["photon"]->SetFillColor(kOrange-3);
  hs_trklen->Add(themap["photon"]);
  
  themap["electron"]->SetLineColor(kMagenta+1);
  themap["electron"]->SetFillColor(kMagenta+1);
  hs_trklen->Add(themap["electron"]);
  
  themap["else"]->SetLineColor(kGray+2);
  themap["else"]->SetFillColor(kGray+2);
  hs_trklen->Add(themap["else"]);
  
  themap["muon"]->SetLineColor(kRed+2);
  themap["muon"]->SetFillColor(kRed+2);
  hs_trklen->Add(themap["muon"]);
  
  hs_trklen->Draw("histo");
  
  themap["total"]->SetFillColor(kBlack);
  themap["total"]->SetFillStyle(3005);
  themap["total"]->Draw("E2 same");
  
  
  
  
  TLegend* leg2;
  
  leg2 = new TLegend(0.6475645,0.5136842,0.8767908,0.8336842,NULL,"brNDC");
  
  std::stringstream sstm;
  
  sstm << "Muon, " << std::setprecision(2)  << themap["muon"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["muon"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Proton, " << std::setprecision(2)  << themap["proton"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["proton"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Pion, " << std::setprecision(2)  << themap["pion"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["pion"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Photon, " << std::setprecision(2)  << themap["photon"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["photon"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Electron, " << std::setprecision(2)  << themap["electron"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["electron"],sstm.str().c_str(),"f");
  sstm.str("");
  
  sstm << "Other, " << std::setprecision(2)  << themap["else"]->Integral() / themap["total"]->Integral()*100. << "%";
  leg2->AddEntry(themap["else"],sstm.str().c_str(),"f");
  sstm.str("");
  
  leg2->AddEntry(themap["total"],"MC Stat Unc.","f");
  leg2->Draw();
  
  return leg2;
  
}


void PlottingTools::DrawDataHisto(TH1D* histo) {

  histo->SetMarkerStyle(kFullCircle);
  histo->SetMarkerSize(0.9);

  histo->Draw("E1 same");
  
}

void PlottingTools::DrawDataHisto2D(TH2D* histo) {

  histo->SetMarkerStyle(kFullCircle);
  histo->SetMarkerSize(0.9);

  histo->Draw("E1 same");
  
}
}

#endif
