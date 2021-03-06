#ifndef __BASE_REWEIGHTINGPLOTTER_CXX__
#define __BASE_REWEIGHTINGPLOTTER_CXX__

#include "ReweightingPlotter.h"

namespace Base {

	void ReweightingPlotter::SetEventBootstrapMap(std::map<std::string, BootstrapTH1D> map_bs)
  {
    _map_bs = map_bs;

    _configured = true;
  }

  void ReweightingPlotter::SetEfficiencyBootstraps(BootstrapTH1D eff_num, BootstrapTH1D eff_den)
  {
    _bs_eff_num = eff_num;
    _bs_eff_den = eff_den;
  }

  void ReweightingPlotter::SetXSecBootstrap(BootstrapTH1D bs_xsec) 
  {

    _bs_xsec = bs_xsec;

  }




  void ReweightingPlotter::MakePlots(int variable, bool normalised, bool makeLaTeX) 
  {

    if (!_configured) {
      std::cout << _name << "Not configured." << std::endl;
      throw std::exception();
    }

    TH1D *histo; // the nominal histogram
    std::map<std::string, std::vector<TH1D>> histo_map; // a map from function name to 2 TH1D (+ and - 1 sigma)
    std::map<std::string, std::vector<TH1D>> histo_map_den; // (only for eff plot) a map from function name to 2 TH1D (+ and - 1 sigma)
    TH1D *histo_p1; // for each function, this will reoresent the p1 histo
    TH1D *histo_m1; // for each function, this will reoresent the m1 histo

    double efficiency_nominal;
    double efficiency_p1;
    double efficiency_m1;

    if (variable == 2 || variable == 3 || variable == 12 || variable == 13 || variable == 102) {

      BootstrapTH1D bs = _map_bs["total"];

      histo = (TH1D*) (bs.GetNominal().Clone("histo_hateroot"));

      histo_map = bs.UnpackPMHisto();

    }

    if (variable == 0 || variable == 1 || variable == 10 || variable == 11 || variable == 100) {

      // Construct nominal histogram
      //histo->Reset();
      histo = (TH1D*) _bs_eff_num.GetNominal().Clone("nominal_eff");
      TH1D * temp_eff_den = (TH1D*)_bs_eff_den.GetNominal().Clone("temp_eff_den");
      histo->Divide(temp_eff_den);
      efficiency_nominal = _bs_eff_num.GetNominal().Integral() / _bs_eff_den.GetNominal().Integral();

      // Numerator for efficiency calculation 
      histo_map = _bs_eff_num.UnpackPMHisto();
      histo_map_den = _bs_eff_den.UnpackPMHisto();

    }



    // Make a directory to store the plots
    if (variable == 0) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 1) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 2) system("mkdir ./EvtWgtEventPlots");
    if (variable == 3) system("mkdir ./EvtWgtEventPlots");
  
    if (variable == 10) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 11) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 12) system("mkdir ./EvtWgtEventPlots");
    if (variable == 13) system("mkdir ./EvtWgtEventPlots");

    if (variable == 100) system("mkdir ./EvtWgtEfficiencyPlots");
    if (variable == 102) system("mkdir ./EvtWgtEventPlots");

    // Avoid root to dislay the canvases
    //gROOT->SetBatch(kTRUE);
  
    // Opening a text file to write the integrals
    std::ofstream outfile;
    if (variable == 0) outfile.open("./EvtWgtEfficiencyPlots/IntegralsPmu.txt");
    if (variable == 1) outfile.open("./EvtWgtEfficiencyPlots/IntegralsCosThetaMu.txt");
    if (variable == 2) outfile.open("./EvtWgtEventPlots/IntegralsPmu.txt");
    if (variable == 3) outfile.open("./EvtWgtEventPlots/IntegralsCosThetaMu.txt");
  
    if (variable == 10) outfile.open("./EvtWgtEfficiencyPlots/IntegralsPproton.txt");
    if (variable == 11) outfile.open("./EvtWgtEfficiencyPlots/IntegralsCosThetaproton.txt");
    if (variable == 12) outfile.open("./EvtWgtEventPlots/IntegralsPproton.txt");
    if (variable == 13) outfile.open("./EvtWgtEventPlots/IntegralsCosThetaproton.txt");
  
    if (variable == 100) outfile.open("./EvtWgtEfficiencyPlots/IntegralsThetaMuP.txt");
    if (variable == 102) outfile.open("./EvtWgtEventPlots/IntegralsThetaMuP.txt");
  
    // Open LaTeX file to write the table
    std::ofstream latexFile;
    std::ofstream latexFile3;
    if(makeLaTeX) {
      if (variable == 0) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPmu.tex");
      if (variable == 1) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyCosThetaMu.tex");
      if (variable == 2) latexFile.open("./EvtWgtEventPlots/evtwgtEventPmu.tex");
      if (variable == 3) latexFile.open("./EvtWgtEventPlots/evtwgtEventCosThetaMu.tex");

      if (variable == 10) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPproton.tex");
      if (variable == 11) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyCosThetaproton.tex");
      if (variable == 12) latexFile.open("./EvtWgtEventPlots/evtwgtEventPproton.tex");
      if (variable == 13) latexFile.open("./EvtWgtEventPlots/evtwgtEventCosThetaproton.tex");

      if (variable == 100) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyThetaMuP.tex");
      if (variable == 102) latexFile.open("./EvtWgtEventPlots/evtwgtEfficiencyThetaMuP.tex");

      latexFile << "\\begin{table}[]" << std::endl;
      latexFile << "\\caption{}" << std::endl;
      latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << std::endl;
      latexFile << "\\label{tab:}" << std::endl;
      latexFile << "\\centering" << std::endl;
      latexFile << "\\begin{tabular}{ccc}" << std::endl;
      latexFile << "\\toprule" << std::endl;
      if (variable == 0 || variable == 1 || variable == 10 || variable == 11) latexFile << "  &  Efficiency  &  Difference (\\%) \\\\" << std::endl;
      else latexFile << "  &  Integral  &  Difference (\\%) \\\\" << std::endl;
      latexFile << "\\midrule" << std::endl;

      if (variable == 0) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPmu_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 1) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyCosThetaMu_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 2) latexFile3.open("./EvtWgtEventPlots/evtwgtEventPmu_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 3) latexFile3.open("./EvtWgtEventPlots/evtwgtEventCosThetaMu_figures.tex", std::ofstream::out | std::ofstream::trunc);

      if (variable == 10) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPproton_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 11) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyCosThetaproton_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 12) latexFile3.open("./EvtWgtEventPlots/evtwgtEventPproton_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 13) latexFile3.open("./EvtWgtEventPlots/evtwgtEventCosThetaproton_figures.tex", std::ofstream::out | std::ofstream::trunc);

      if (variable == 100) latexFile3.open("./EvtWgtEfficiencyPlots/evtwgtEventThetaMuP_figures.tex", std::ofstream::out | std::ofstream::trunc);
      if (variable == 102) latexFile3.open("./EvtWgtEventPlots/evtwgtEventThetaMuP_figures.tex", std::ofstream::out | std::ofstream::trunc);

      latexFile3 << "\\begin{figure}[t]" << std::endl;
      latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
      latexFile3 << "\\centering" << std::endl;
    }
    

    bool is_first = true;

    int function_counter = 0;

    for(auto iter : histo_map) {

      function_counter++;

      std::string function_name = iter.first;

      //std::cout << "This is function name: " << function_name << std::endl;


      if (variable == 2 || variable == 3 || variable == 12 || variable == 13 || variable == 102) {

        histo_p1 = (TH1D*) iter.second.at(1).Clone("histo_p1_hateroot");
        histo_m1 = (TH1D*) iter.second.at(0).Clone("histo_m1_hateroot");

      }

      if (variable == 0 || variable == 1 || variable == 10 || variable == 11 || variable == 100) {

        //histo_p1->Reset();
        histo_p1 = (TH1D*) iter.second.at(1).Clone("p1");
        TH1D * temp_eff_den_p1 = (TH1D*)histo_map_den[function_name].at(1).Clone("temp_eff_den_p1");
        histo_p1->Divide(temp_eff_den_p1);

        efficiency_p1 = iter.second.at(1).Integral() / histo_map_den[function_name].at(1).Integral();

        //histo_m1->Reset();
        histo_m1 = (TH1D*) iter.second.at(0).Clone("m1");
        TH1D * temp_eff_den_m1 = (TH1D*)histo_map_den[function_name].at(0).Clone("temp_eff_den_m1");
        histo_m1->Divide(temp_eff_den_m1);

        efficiency_m1 = iter.second.at(0).Integral() / histo_map_den[function_name].at(0).Integral();

      }



      TString SaveName;
      if(variable == 0 || variable == 2) SaveName = "Pmu_"+function_name;
      if(variable == 1 || variable == 3) SaveName = "CosThetaMu_"+function_name;

      if(variable == 10 || variable == 12) SaveName = "Pproton_"+function_name;
      if(variable == 11 || variable == 13) SaveName = "CosThetaproton_"+function_name;

      if(variable == 100 || variable == 102) SaveName = "ThetaMuP_"+function_name;

      TString LegName  = GetLegendName(function_name);

      if(normalised) SaveName += "_normalised";


      double histo_Int = histo->Integral();
      double histo_p1_Int = histo_p1->Integral();
      double histo_m1_Int = histo_m1->Integral();


      if (normalised) {
        histo->Scale(1./histo_Int);
        histo_p1->Scale(1./histo_p1_Int);
        histo_m1->Scale(1./histo_m1_Int);
      }

      // Calculate integrals
      outfile << function_name << std::endl;
      outfile << "Integral Nominal:    " << histo->Integral() << std::endl;
      outfile << "Integral nominal_p1: " << histo_p1->Integral() << std::endl;
      outfile << "Integral nominal_m1: " << histo_m1->Integral() << std::endl;

      outfile << "Difference w.r.t. Nominal (%):" << std::endl;
      outfile << "nominal_p1: " << (histo_p1->Integral()-histo->Integral())/(histo->Integral())*100. << std::endl;
      outfile << "nominal_m1: " << (histo_m1->Integral()-histo->Integral())/(histo->Integral())*100. << std::endl;
      outfile << "--------------------------------------" << std::endl << std::endl;


      if(makeLaTeX) {
        if ( (variable == 2) || (variable == 3 || variable == 12 || variable ==13 || variable == 102) ) {
          if (is_first) latexFile << "Nominal" << " & " << histo->Integral() << " & 0" << "\\\\" << std::endl;
          latexFile << "\\midrule" << std::endl;
          latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ & " << histo_p1->Integral() << " & " << (histo_p1->Integral()-histo->Integral())/(histo->Integral())*100. << "\\\\" << std::endl;
          latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ & " << histo_m1->Integral() << " & " << (histo_m1->Integral()-histo->Integral())/(histo->Integral())*100. << "\\\\" << std::endl;
        }

        if (variable == 0 || variable == 1 || variable == 10 || variable == 11 || variable == 100) {  // save the value of the efficiency to the LaTeX file

          double eff_nom = efficiency_nominal * 100.;
          double eff_p1  = efficiency_p1      * 100.;
          double eff_m1  = efficiency_m1      * 100.;

          if (is_first) latexFile << "Nominal" << " & " << std::setprecision(4) << eff_nom << " & 0" << "\\\\" << std::endl;
          latexFile << "\\midrule" << std::endl;
        
          latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ & "
          << std::setprecision(4) << eff_p1 << " & "
          << std::setprecision(4) << (eff_p1-eff_nom)/eff_nom*100. << "\\\\" << std::endl;
        
          latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ & "
          << std::setprecision(4) << eff_m1 << " & "
          << std::setprecision(4) << (eff_m1-eff_nom)/eff_nom*100. << "\\\\" << std::endl;
        }

      }


      // Define the Canvas
      TCanvas *c = new TCanvas("c", "canvas", 800, 1200);



      // Upper plot will be in pad1
      TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1.0, 1.0);
      pad1->SetBottomMargin(0); // Upper and lower plot are joined
      pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();               // pad1 becomes the current pad
      if (variable == 0 || variable == 1 || variable == 10 || variable == 11) histo_p1->SetMaximum(1.);
      histo_p1->SetStats(0);          // No statistics on upper plot
      histo_m1->SetStats(0);          // No statistics on upper plot
      histo->SetStats(0);          // No statistics on upper plot
      histo_p1->Draw("histo");               // Draw h1
      histo->Draw("histo same");         // Draw h2 on top of h1
      histo_m1->Draw("histo same");
    
      //uBooNESimulation();
    
      if (normalised) {
        // TLatex
        double x;
        x = 0.87;
        if (variable == 3 || variable  == 13) x = 0.46;
        double y = 0.52;
        double size = 28;
        int color = 1;
        int font = 43;
        int align = 32;
        TLatex *latex = new TLatex( x, y, "Area Normalised" );
        latex->SetNDC();
        latex->SetTextSize(size);
        latex->SetTextColor(color);
        latex->SetTextFont(font);
        latex->SetTextAlign(align);
      
        latex->Draw();
      
      }

      // Do not draw the Y axis label on the upper plot and redraw a small
      // axis instead, in order to avoid the first label (0) to be clipped.
      histo->GetYaxis()->SetLabelSize(0.);
      TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
      axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      axis->SetLabelSize(15);
      axis->Draw();
      
      // Legend for the upper plot
      TLegend* leg;
      if (variable == 0 || variable == 2 || variable == 10 || variable == 12 || variable == 100 || variable == 102) leg = new TLegend(0.65,0.6,.85,0.87);
      if (variable == 1 || variable == 3 || variable == 11 || variable == 13) leg = new TLegend(0.216792,0.5723502,0.4172932,0.843318,NULL,"brNDC");
      leg->SetTextFont(42);
      leg->SetBorderSize(0);
      //leg->SetHeader("");
      //leg->SetTextFont(42);
      //leg->AddEntry(histo, "BBA2005 (Nominal)");
      //leg->AddEntry(histo_p1, "Dipole");
      leg->AddEntry(histo, "Nominal");
      leg->AddEntry(histo_p1, LegName + " + 1#sigma");
      leg->AddEntry(histo_m1, LegName + " - 1#sigma");
      leg->Draw();
    
      // lower plot will be in pad
      c->cd();          // Go back to the main canvas before defining pad2
      TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.5);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.5);
      pad2->SetGridx(); // vertical grid
      //pad2->SetGridy(); // orizontal grid
      pad2->Draw();
      pad2->cd();       // pad2 becomes the current pad
    
    
      // Define the first ratio plot
      TH1D *ratio_p1 = (TH1D*)histo_p1->Clone("ratio_p1");
      //ratio_p1->SetMinimum(0.92);  // Define Y ..
      //ratio_p1->SetMaximum(1.08); // .. range
      //ratio_p1->Sumw2();
      ratio_p1->SetStats(0);      // No statistics on lower plot
      ratio_p1->Add(histo, -1);
      ratio_p1->Divide(histo);
      ratio_p1->SetLineWidth(2);
      ratio_p1->SetLineColor(kRed+1);
      //ratio_p1->Draw("hist");       // Draw the ratio plot
    
      // Define the second ratio plot
      TH1D *ratio_m1 = (TH1D*)histo_m1->Clone("ratio_m1");
      //ratio_m1->SetMinimum(0.9);  // Define Y ..
      //ratio_m1->SetMaximum(1.1); // .. range
      //ratio_m1->Sumw2();
      ratio_m1->SetStats(0);      // No statistics on lower plot
      ratio_m1->Add(histo, -1);
      ratio_m1->Divide(histo);
      ratio_m1->SetLineWidth(2);
      ratio_m1->SetLineColor(kGreen+2);
      //ratio_m1->Draw("hist same");       // Draw the ratio plot
    
    
    
      // Try to set the Y range for the ratio plots
      double max = ratio_p1->GetMaximum();
      double min = ratio_p1->GetMinimum();
      if (ratio_m1->GetMaximum() > max) max = ratio_m1->GetMaximum();
      if (ratio_m1->GetMinimum() < min) min = ratio_m1->GetMinimum();
      //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
      //ratio_p1->SetMinimum(max);  // Define Y ..
      //ratio_p1->SetMaximum(min); // .. range
      //ratio_m1->SetMinimum(max);  // Define Y ..
      //ratio_m1->SetMaximum(min); // .. range

    
      // Draw the ratio plot
      //ratio_p1->Draw("hist");
      //ratio_m1->Draw("hist same");
    
      THStack *hs = new THStack("hs","");
      hs->Add(ratio_p1);
      hs->Add(ratio_m1);
      hs->SetMaximum(hs->GetMaximum("nostack")+0.01);
      hs->SetMinimum(hs->GetMinimum("nostack")-0.01);
      //std::cout << "hs->GetMinimum(): " << hs->GetMinimum("nostack") << std::endl;
      hs->Draw("NOSTACK histo");
    
    
      //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
    
      histo_p1->SetMinimum(0.0001); // Otherwise 0 label overlaps (need to do it after the THStack, otherwise sets the minimum)


      //**********************
      //
      // Settings
      //
      //**********************

      // h1 settings
      histo->SetLineColor(kBlack);
      histo->SetLineWidth(2);

      // Y axis h1 plot settings
      histo_p1->GetYaxis()->SetTitle("Selected Events");
      if (variable == 0 || variable == 1 || variable == 10 || variable ==11 || variable == 100) histo_p1->GetYaxis()->SetTitle("Efficiency");
      histo_p1->GetYaxis()->CenterTitle();
      histo_p1->GetYaxis()->SetTitleSize(25);
      histo_p1->GetYaxis()->SetTitleFont(43);
      histo_p1->GetYaxis()->SetTitleOffset(1.55);

      // h2 settings
      histo_p1->SetLineColor(kRed+1);
      histo_p1->SetLineWidth(2);

      // h3 settings
      histo_m1->SetLineColor(kGreen+2);
      histo_m1->SetLineWidth(2);

      // Ratio plot (ratio_p1) settings
      ratio_p1->SetTitle(""); // Remove the ratio title
    
    
      hs->GetYaxis()->SetTitle("Ratio");
      hs->GetYaxis()->CenterTitle();
      hs->GetYaxis()->SetNdivisions(505);
      hs->GetYaxis()->SetTitleSize(25);
      hs->GetYaxis()->SetTitleFont(43);
      hs->GetYaxis()->SetTitleOffset(1.2);
      hs->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      hs->GetYaxis()->SetLabelSize(15);
    
    
      if(variable == 0 || variable == 2) hs->GetXaxis()->SetTitle("p_{#mu} [GeV]");
      if(variable == 1 || variable == 3) hs->GetXaxis()->SetTitle("cos#theta_{#mu}");
      if(variable == 10 || variable == 12) hs->GetXaxis()->SetTitle("p_{p} [GeV]");
      if(variable == 11 || variable == 13) hs->GetXaxis()->SetTitle("cos#theta_{p}");
      if(variable == 100 || variable ==102) hs->GetXaxis()->SetTitle("#theta_{#mu p}");
      hs->GetXaxis()->CenterTitle();
      hs->GetXaxis()->SetTitleSize(25);
      hs->GetXaxis()->SetTitleFont(43);
      hs->GetXaxis()->SetTitleOffset(3.5);
      hs->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      hs->GetXaxis()->SetLabelSize(20);
    
      // Draw linea at 1 in ratio plot
      TLine *line;
      if (variable == 0 || variable == 2) line = new TLine(0,0,2.5,0);
      if (variable == 10 || variable == 12) line = new TLine(0,0,1.5,0);
      if (variable ==100 || variable ==102) line = new TLine(0, 0, 3.14, 0);
      if (variable == 1 || variable == 3 || variable == 11 || variable == 13) line = new TLine(-1,1,1,1);
      line->SetLineColor(kBlack);
      line->SetLineStyle(9); // dashed
      line->Draw();
    
      if (variable == 0 || variable == 1 || variable == 10 || variable == 11 || variable == 100) {
        c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".C");
         c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".pdf");
      }
      if (variable == 2 || variable == 3 || variable == 12 || variable == 13 || variable == 102) {
        c->Print("./EvtWgtEventPlots/" + SaveName + ".C");
        c->Print("./EvtWgtEventPlots/" + SaveName + ".pdf");
      }

      if (makeLaTeX) {
        latexFile3 << "\\subfloat[][]" << std::endl;
        latexFile3 << "   {\\includegraphics[width=.35\\textwidth]{images/EvtWgtEfficiencyPlots/" << SaveName << "}" << std::endl;
        latexFile3 << "   \\label{fig:" << "EvtWgtEfficiencyPlots_" << SaveName << "}} \\quad" << std::endl;
      }

      if (makeLaTeX && function_counter % 12 == 0) {
        latexFile3 << "\\caption{Efficiency Plots}" << std::endl;
        latexFile3 << "\\label{fig:EvtWgtEfficiencyPlots}" << std::endl;
        latexFile3 << "\\end{adjustwidth}" << std::endl;
        latexFile3 << "\\end{figure}" << std::endl;

        latexFile3 << "\\begin{figure}[t]" << std::endl;
        latexFile3 << "\\ContinuedFloat" << std::endl;
        latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
        latexFile3 << "\\centering" << std::endl;

      }
      //--------------------------------------------------------------------
      //c_all->cd();   
      //if(function_counter==1) ratio_p1->Draw("hist");
      //if(function_counter>1) ratio_p1->Draw("hist, same");
      //c_all->Close();
      //--------------------------------------------------------------------

      is_first = false;

    } // end loop functions
    
    if(makeLaTeX) {

      latexFile << "\\bottomrule" << std::endl;
      latexFile << "\\end{tabular}" << std::endl;
      latexFile << "\\end{table}" << std::endl;


      latexFile3 << "\\caption{Efficiency Plots}" << std::endl;
      latexFile3 << "\\label{fig:EvtWgtXsecDiffPlots}" << std::endl;
      latexFile3 << "\\end{adjustwidth}" << std::endl;
      latexFile3 << "\\end{figure}" << std::endl;
    }



  }

  void ReweightingPlotter::MakeBackgroundPlots(int variable, bool normalised, bool makeLaTeX) {
  

      // Pmu: variable == 0
  // CosThetaMu: variable == 1


  

  
  // Create output directory
  system("mkdir ./EvtWgtBackgroundPlots");
  system("mkdir ./EvtWgtBackgroundPlots_reducedLegend"); 
  system("mkdir ./EvtWgtBackground_RelativeUncertainty"); 
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  
  double scaleFactor = 1.;//5.3e19/1.22e20;
  
  std::cout << "numu:   " << _map_bs["signal"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "anumu:  " << _map_bs["anumu"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "nue:    " << _map_bs["nue"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "nc:     " << _map_bs["nc"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "outfv:  " << _map_bs["outfv"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "cosmic: " << _map_bs["cosmic"].GetNominal().Integral() * scaleFactor << std::endl;
  std::cout << "cc_other: "<<_map_bs["cc_other"].GetNominal().Integral() * scaleFactor << std::endl;
  std::ofstream latexFile;
  std::ofstream latexFile3;

  double numu_nominal = _map_bs["signal"].GetNominal().Integral();
  double anumu_nominal = _map_bs["anumu"].GetNominal().Integral();
  double nue_nominal = _map_bs["nue"].GetNominal().Integral();
  double nc_nominal = _map_bs["nc"].GetNominal().Integral();
  double outfv_nominal = _map_bs["outfv"].GetNominal().Integral();
  double cosmic_nominal = _map_bs["cosmic"].GetNominal().Integral();
  double ccother_nominal = _map_bs["cc_other"].GetNominal().Integral();

  



  if(makeLaTeX) {
    if (variable == 0) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPmu.tex");
    if (variable == 1) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundCosThetaMu.tex");
    if (variable == 2) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPproton.tex");
    if (variable == 3) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundCosThetaproton.tex");
    if (variable == 4) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundThetamup.tex");

    latexFile << "\\begin{table}[]" << std::endl;
    latexFile << "\\begin{adjustwidth}{-2.1cm}{-1cm}" << std::endl;
    latexFile << "\\caption{}" << std::endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << std::endl;
    latexFile << "\\label{tab:}" << std::endl;
    latexFile << "\\centering" << std::endl;
    latexFile << "\\tiny" << std::endl;
    latexFile << "\\begin{tabular}{c   c c  c c  c c  c c  c c  c c }" << std::endl;
    latexFile << "\\toprule" << std::endl;
    latexFile << "  &  \\multicolumn{2}{ c }{$\\nu_\\mu$ CC} &  \\multicolumn{2}{ c }{$\\bar{\\nu}_\\mu$ CC}   & \\multicolumn{2}{ c }{$\\nu_e$, $\\bar{\\nu}_e$ CC}   &  \\multicolumn{2}{ c }{NC} &  \\multicolumn{2}{ c }{OUTFV} &  \\multicolumn{2}{ c }{Cosmic} \\\\" << std::endl;
    latexFile << "  &  Events   &  Diff. (\\%) &  Events  &  Diff. (\\%) & Events  &  Diff. (\\%) &  Events &  Diff. (\\%)  &  Events &  Diff. (\\%)  &  Events &  Diff. (\\%) \\\\" << std::endl;
    latexFile << "\\midrule" << std::endl;
    latexFile << "Nominal & " << _map_bs["signal"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["anumu"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["nue"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["nc"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["outfv"].GetNominal().Integral() << " & 0 "
    << " & " << _map_bs["cosmic"].GetNominal().Integral() << " & 0 " 
    << " & " << _map_bs["cc_other"].GetNominal().Integral() << " & 0 "<< "\\\\" << std::endl;

    if (variable == 0) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPmu_figures.tex", std::ofstream::out | std::ofstream::trunc);
    if (variable == 1) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundCosThetaMu_figures.tex", std::ofstream::out | std::ofstream::trunc);
    if (variable == 0) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPproton_figures.tex", std::ofstream::out | std::ofstream::trunc);
    if (variable == 1) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundCosThetaproton_figures.tex", std::ofstream::out | std::ofstream::trunc);
    if (variable == 0) latexFile3.open("./EvtWgtBackgroundPlots/evtwgtBackgroundThetamup_figures.tex", std::ofstream::out | std::ofstream::trunc);




    latexFile3 << "\\begin{figure}[t]" << std::endl;
    latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
    latexFile3 << "\\centering" << std::endl;
  }
  
  
  
  
  
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  
  
  std::map<std::string, std::vector<TH1D>> map_bs_signal = _map_bs["signal"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_anumu = _map_bs["anumu"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_nue = _map_bs["nue"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_nc = _map_bs["nc"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_outfv = _map_bs["outfv"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_cosmic = _map_bs["cosmic"].UnpackPMHisto();
  std::map<std::string, std::vector<TH1D>> map_bs_ccother = _map_bs["cc_other"].UnpackPMHisto();

  
  
  
  /*std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<number of Universes is "<<map_bs_ccother.GetNUniverses()<<std::endl; 
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<integral of the m1(ccother) is "<<map_bs_ccother.at(75).Integral()<<std::endl; 
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<integral of the p1(ccother) is "<<map_bs_ccother.at(75).Integral()<<std::endl; 

  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<integral of the nominal(cosmics) is "<<map_bs_cosmic.at(75).Integral()<<std::endl; 
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<integral of the m1(cosmics) is "<<map_bs_cosmic.at(75).Integral()<<std::endl; 
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<integral of the p1(cosmics) is "<<map_bs_cosmic.at(75).Integral()<<std::endl; 
  */
 

  int function_counter = 0;

  for (auto iter : map_bs_signal) {

    function_counter++;

    std::string function_name = iter.first;

    auto ii = map_bs_anumu.find(function_name);

    std::vector<TH1D> anumu_v = ii->second;
    const TH1D* anumu_nom = &_map_bs["anumu"].GetNominal();
    ii = map_bs_nue.find(function_name);
    std::vector<TH1D> nue_v = ii->second;
    const TH1D* nue_nom = &_map_bs["nue"].GetNominal();
    ii = map_bs_nc.find(function_name);
    std::vector<TH1D> nc_v = ii->second;
    const TH1D* nc_nom = &_map_bs["nc"].GetNominal();
    ii = map_bs_outfv.find(function_name);
    std::vector<TH1D> outfv_v = ii->second;
    const TH1D* outfv_nom = &_map_bs["outfv"].GetNominal();
    ii = map_bs_cosmic.find(function_name);
    std::vector<TH1D> cosmic_v = ii->second;
    const TH1D* cosmic_nom = &_map_bs["cosmic"].GetNominal();
    ii = map_bs_ccother.find(function_name);
    std::vector<TH1D> ccother_v = ii->second;
    const TH1D* ccother_nom = &_map_bs["cc_other"].GetNominal();


    
    TString SaveName;
    if (variable == 0) SaveName = "Pmu_" + function_name;
    if (variable == 1) SaveName = "CosThetaMu_" + function_name;
    if (variable == 2) SaveName = "Pproton_" + function_name;
    if (variable == 3) SaveName = "CosThetaproton_" + function_name;
    if (variable == 4) SaveName = "Thetamup_" + function_name;
      
    c->SetLogy();                                 // Log scale
    //c->SetGridy();                                // Horizontal grid
    
    
    TH1D* histo_numu;
    TH1D* histo_numu_p1;
    TH1D* histo_numu_m1;
    TH1D* histo_anumu;
    TH1D* histo_anumu_p1;
    TH1D* histo_anumu_m1;
    TH1D* histo_nue;
    TH1D* histo_nue_p1;
    TH1D* histo_nue_m1;
    TH1D* histo_nc;
    TH1D* histo_nc_p1;
    TH1D* histo_nc_m1;
    TH1D* histo_outfv;
    TH1D* histo_outfv_p1;
    TH1D* histo_outfv_m1;
    TH1D* histo_cosmic;
    TH1D* histo_cosmic_p1;
    TH1D* histo_cosmic_m1;
    TH1D* histo_ccother;
    TH1D* histo_ccother_p1;
    TH1D* histo_ccother_m1;

    
    if (variable == 0 || variable == 1 || variable == 2 || variable ==3 || variable == 4) {
      histo_numu        = (TH1D*)_map_bs["signal"].GetNominal().Clone("histo_numu");
      histo_numu_p1     = (TH1D*)iter.second.at(1).Clone("histo_numu_p1");
      histo_numu_m1     = (TH1D*)iter.second.at(0).Clone("histo_numu_m1");
      
      histo_anumu       = (TH1D*) (anumu_nom->Clone("histo_anumu"));
      histo_anumu_p1    = (TH1D*) (anumu_v.at(1).Clone("histo_anumu_p1"));
      histo_anumu_m1    = (TH1D*) (anumu_v.at(0).Clone("histo_anumu_m1"));
      
      histo_nue         = (TH1D*) (nue_nom->Clone("histo_nue"));
      histo_nue_p1      = (TH1D*) (nue_v.at(1).Clone("histo_nue_p1"));
      histo_nue_m1      = (TH1D*) (nue_v.at(0).Clone("histo_nue_m1"));
      
      histo_nc          = (TH1D*) (nc_nom->Clone("histo_nc"));
      histo_nc_p1       = (TH1D*) (nc_v.at(1).Clone("histo_nc_p1"));
      histo_nc_m1       = (TH1D*) (nc_v.at(0).Clone("histo_nc_m1"));

      histo_outfv       = (TH1D*) (outfv_nom->Clone("histo_outfv"));
      histo_outfv_p1    = (TH1D*) (outfv_v.at(1).Clone("histo_outfv_p1"));
      histo_outfv_m1    = (TH1D*) (outfv_v.at(0).Clone("histo_outfv_m1"));

      histo_cosmic      = (TH1D*) (cosmic_nom->Clone("histo_cosmic"));
      histo_cosmic_p1   = (TH1D*) (cosmic_v.at(1).Clone("histo_cosmic_p1"));
      histo_cosmic_m1   = (TH1D*) (cosmic_v.at(0).Clone("histo_cosmic_m1"));

      histo_ccother      = (TH1D*) (ccother_nom->Clone("histo_ccother"));
      histo_ccother_p1   = (TH1D*) (ccother_v.at(1).Clone("histo_ccother_p1"));
      histo_ccother_m1   = (TH1D*) (ccother_v.at(0).Clone("histo_ccother_m1"));
      //std::cout<<"libo test ccother<<<<<<<<<<<<<<<"<<histo_ccother_p1->Integral()<<"  "<<histo_ccother_m1->Integral()<<std::endl;
      //std::cout<<"libo test cosmic<<<<<<<<<<<<<<<"<<histo_cosmic_p1->Integral()<<"  "<<histo_cosmic_m1->Integral()<<std::endl;


    }
    else {
      std::cout << "Invalid option. Exit." << std::endl;
      exit(0);
    }
    
    
    // Settings
    histo_numu->SetStats(0);          // No statistics on upper plot
    if (variable == 1) histo_numu->SetMinimum(1);
    if (variable == 1) histo_numu->SetMaximum(1e5);

    if (variable == 0) histo_numu->GetXaxis()->SetTitle("Reconstructed p_{#mu} [GeV]");
    if (variable == 1) histo_numu->GetXaxis()->SetTitle("Reconstructed cos#theta_{#mu}");

    if (variable == 2) histo_numu->GetXaxis()->SetTitle("Reconstructed p_{proton} [GeV]");
    if (variable == 3) histo_numu->GetXaxis()->SetTitle("Reconstructed cos#theta_{p}");
    if (variable == 4) histo_numu->GetXaxis()->SetTitle("Reconstructed #theta_{#mu, p}");

    histo_numu->GetXaxis()->CenterTitle();
    histo_numu->GetXaxis()->SetTitleSize(25);
    histo_numu->GetXaxis()->SetTitleFont(43);
    histo_numu->GetXaxis()->SetTitleOffset(1.45);
    
    histo_numu->GetYaxis()->SetTitle("Events");
    histo_numu->GetYaxis()->CenterTitle();
    histo_numu->GetYaxis()->SetTitleSize(25);
    histo_numu->GetYaxis()->SetTitleFont(43);
    histo_numu->GetYaxis()->SetTitleOffset(1.55);
    
    double lineWidth = 1;
    
    // Nominal
    histo_numu   ->Draw("histo");
    histo_anumu  ->Draw("histo same");
    histo_nue    ->Draw("histo same");
    histo_nc     ->Draw("histo same");
    histo_outfv  ->Draw("histo same");
    histo_cosmic ->Draw("histo same");
    histo_ccother ->Draw("histo same");

    histo_numu   ->SetLineColor(kRed+1);
    histo_anumu  ->SetLineColor(kOrange+1);
    histo_nue    ->SetLineColor(kViolet+1);
    histo_nc     ->SetLineColor(kGray+1);
    histo_outfv  ->SetLineColor(kGreen+1);
    histo_cosmic ->SetLineColor(kBlue+1);
    histo_ccother ->SetLineColor(kCyan+1);
    
    histo_numu   ->SetLineWidth(lineWidth+1);
    histo_anumu  ->SetLineWidth(lineWidth+1);
    histo_nue    ->SetLineWidth(lineWidth+1);
    histo_nc     ->SetLineWidth(lineWidth+1);
    histo_outfv  ->SetLineWidth(lineWidth+1);
    histo_cosmic ->SetLineWidth(lineWidth+1);
    histo_ccother ->SetLineWidth(lineWidth+1);

    // +- 1 sigma
    histo_numu_p1   ->Draw("histo same");
    histo_anumu_p1  ->Draw("histo same");
    histo_nue_p1    ->Draw("histo same");
    histo_nc_p1     ->Draw("histo same");
    histo_outfv_p1  ->Draw("histo same");
    histo_cosmic_p1 ->Draw("histo same");
    histo_ccother_p1 ->Draw("histo same");

    histo_numu_m1   ->Draw("histo same");
    histo_anumu_m1  ->Draw("histo same");
    histo_nue_m1    ->Draw("histo same");
    histo_nc_m1     ->Draw("histo same");
    histo_outfv_m1  ->Draw("histo same");
    histo_cosmic_m1 ->Draw("histo same");
    histo_ccother_m1 ->Draw("histo same");

    histo_numu_p1  ->SetLineColor(kRed+1);
    histo_anumu_p1 ->SetLineColor(kOrange+1);
    histo_nue_p1   ->SetLineColor(kViolet+1);
    histo_nc_p1    ->SetLineColor(kGray+1);
    histo_outfv_p1 ->SetLineColor(kGreen+1);
    histo_cosmic_p1->SetLineColor(kBlue+1);
    histo_ccother_p1 ->SetLineColor(kCyan+1);

    histo_numu_m1  ->SetLineColor(kRed+1);
    histo_anumu_m1 ->SetLineColor(kOrange+1);
    histo_nue_m1   ->SetLineColor(kViolet+1);
    histo_nc_m1    ->SetLineColor(kGray+1);
    histo_outfv_m1 ->SetLineColor(kGreen+1);
    histo_cosmic_m1->SetLineColor(kBlue+1);
    histo_ccother_m1 ->SetLineColor(kCyan+1);

    histo_numu_p1  ->SetLineStyle(7);
    histo_anumu_p1 ->SetLineStyle(7);
    histo_nue_p1   ->SetLineStyle(7);
    histo_nc_p1    ->SetLineStyle(7);
    histo_outfv_p1 ->SetLineStyle(7);
    histo_cosmic_p1->SetLineStyle(7);
    histo_ccother_p1->SetLineStyle(7);

    histo_numu_m1  ->SetLineStyle(3);
    histo_anumu_m1 ->SetLineStyle(3);
    histo_nue_m1   ->SetLineStyle(3);
    histo_nc_m1    ->SetLineStyle(3);
    histo_outfv_m1 ->SetLineStyle(3);
    histo_cosmic_m1->SetLineStyle(3);
    histo_ccother_m1 ->SetLineStyle(3);

    histo_numu_p1  ->SetLineWidth(lineWidth);
    histo_anumu_p1 ->SetLineWidth(lineWidth);
    histo_nue_p1   ->SetLineWidth(lineWidth);
    histo_nc_p1    ->SetLineWidth(lineWidth);
    histo_outfv_p1 ->SetLineWidth(lineWidth);
    histo_cosmic_p1->SetLineWidth(lineWidth);
    histo_ccother_p1->SetLineWidth(lineWidth);

    histo_numu_m1  ->SetLineWidth(lineWidth);
    histo_anumu_m1 ->SetLineWidth(lineWidth);
    histo_nue_m1   ->SetLineWidth(lineWidth);
    histo_nc_m1    ->SetLineWidth(lineWidth);
    histo_outfv_m1 ->SetLineWidth(lineWidth);
    histo_cosmic_m1->SetLineWidth(lineWidth);
    histo_ccother_m1 ->SetLineWidth(lineWidth);

    
    // Legend
    TLegend* leg;
    if (variable == 0 || variable == 2) leg = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1 || variable == 3 || variable ==4) leg = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    
    leg->AddEntry(histo_numu,  "Signal (nominal)");
    leg->AddEntry(histo_anumu, "#bar{#nu}_{#mu} CC (nominal)");
    leg->AddEntry(histo_nue,   "#nu_{e} CC (nominal)");
    leg->AddEntry(histo_nc,    "NC all flavours (nominal)");
    leg->AddEntry(histo_outfv,   "OUTFV (nominal)");
    leg->AddEntry(histo_cosmic,    "Cosmic (nominal)");
    leg->AddEntry(histo_ccother, "CC-Other (nominal");

    leg->AddEntry(histo_numu_p1,   "Signal (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_anumu_p1,  "#bar{#nu}_{#mu} CC (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_nue_p1,    "#nu_{e} CC (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_nc_p1,     "NC all flavours (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_outfv_p1,  "OUTFV (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_cosmic_p1, "Cosmic (" + GetLegendName(function_name) + " + 1 #sigma)");
    leg->AddEntry(histo_ccother_p1, "CC-Other (" +GetLegendName(function_name) + " +1 #sigma)");

    leg->AddEntry(histo_numu_m1,   "Signal (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_anumu_m1,  "#bar{#nu}_{#mu} CC (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_nue_m1,    "#nu_{e} CC (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_nc_m1,     "NC all flavours (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_outfv_m1,  "OUTFV (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_cosmic_m1, "Cosmic (" + GetLegendName(function_name) + " - 1 #sigma)");
    leg->AddEntry(histo_ccother_m1, "CC-Other (" + GetLegendName(function_name) + "- 1 #sigma");

    leg->Draw();
    
    
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".pdf");

    

























    

    // Reduced Legend Plots
    histo_numu->GetXaxis()->SetTitleOffset(1.10);
    histo_numu->GetYaxis()->SetTitleOffset(1.30);
    histo_numu->GetXaxis()->SetTitleSize(30);
    histo_numu->GetYaxis()->SetTitleSize(30);

    histo_numu  ->Draw("histo");
    histo_anumu ->Draw("histo same");
    histo_nue   ->Draw("histo same");
    histo_nc    ->Draw("histo same");
    histo_outfv ->Draw("histo same");
    histo_cosmic->Draw("histo same");
    histo_ccother -> Draw("histo same");

    histo_numu_p1  ->Draw("histo same");
    histo_anumu_p1 ->Draw("histo same");
    histo_nue_p1   ->Draw("histo same");
    histo_nc_p1    ->Draw("histo same");
    histo_outfv_p1 ->Draw("histo same");
    histo_cosmic_p1->Draw("histo same");
    histo_ccother_p1 ->Draw("histo same");

    histo_numu_m1  ->Draw("histo same");
    histo_anumu_m1 ->Draw("histo same");
    histo_nue_m1   ->Draw("histo same");
    histo_nc_m1    ->Draw("histo same");
    histo_outfv_m1 ->Draw("histo same");
    histo_cosmic_m1->Draw("histo same");
    histo_ccother_m1 ->Draw("histo same");

    TLegend* leg2;
    if (variable == 0 || variable == 2 ) leg2 = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1 || variable == 3 || variable == 4) leg2 = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);

    leg2->AddEntry(histo_numu,  "Signal");
    leg2->AddEntry(histo_anumu, "#bar{#nu}_{#mu} CC");
    leg2->AddEntry(histo_nue,   "#nu_{e}, #bar{#nu}_{e} CC");
    leg2->AddEntry(histo_nc,    "NC all flavours");
    leg2->AddEntry(histo_outfv, "OUTFV");
    leg2->AddEntry(histo_cosmic,"Cosmic");
    leg2->AddEntry(histo_ccother, "CC-Other");

    leg2->Draw();

    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".pdf");

    //try to plot the relative uncertainty here

    c->SetLogy(false);

     
    TH1D * histo_totalbkg_nominal;
    TH1D * histo_totalbkg_p1;
    TH1D * histo_totalbkg_m1;

    histo_totalbkg_nominal=(TH1D*)histo_anumu->Clone();
    histo_totalbkg_nominal->Add(histo_nue);
    histo_totalbkg_nominal->Add(histo_nc);
    histo_totalbkg_nominal->Add(histo_outfv);
    histo_totalbkg_nominal->Add(histo_cosmic);
    histo_totalbkg_nominal->Add(histo_ccother);

    histo_totalbkg_p1=(TH1D*)histo_anumu_p1->Clone();
    histo_totalbkg_p1->Add(histo_nue_p1);
    histo_totalbkg_p1->Add(histo_nc_p1);
    histo_totalbkg_p1->Add(histo_outfv_p1);
    histo_totalbkg_p1->Add(histo_cosmic_p1);
    histo_totalbkg_p1->Add(histo_ccother_p1);
    histo_totalbkg_p1->Add(histo_totalbkg_nominal, -1.0);
    histo_totalbkg_p1->Divide(histo_totalbkg_nominal);

    histo_totalbkg_m1=(TH1D*)histo_anumu_m1->Clone();
    histo_totalbkg_m1->Add(histo_nue_m1);
    histo_totalbkg_m1->Add(histo_nc_m1);
    histo_totalbkg_m1->Add(histo_outfv_m1);
    histo_totalbkg_m1->Add(histo_cosmic_m1);
    histo_totalbkg_m1->Add(histo_ccother_m1);
    histo_totalbkg_m1->Add(histo_totalbkg_nominal, -1.0);
    histo_totalbkg_m1->Divide(histo_totalbkg_nominal);


    histo_totalbkg_p1->SetLineColor(kRed);
    histo_totalbkg_m1->SetLineColor(kBlue);
    //histo_totalbkg_p1->SetMaximum();
    //histo_totalbkg_m1->SetMinimum(-0.5);
    histo_totalbkg_p1->Draw("hist");
    histo_totalbkg_m1->Draw("hist,same");


    c->Print("./EvtWgtBackground_RelativeUncertainty/" + SaveName + ".C");
    c->Print("./EvtWgtBackground_RelativeUncertainty/" + SaveName + ".pdf");










    if (makeLaTeX) {
        latexFile3 << "\\subfloat[][]" << std::endl;
        latexFile3 << "   {\\includegraphics[width=.35\\textwidth]{images/EvtWgtBackgroundPlots/" << SaveName << "}" << std::endl;
        latexFile3 << "   \\label{fig:" << "EvtWgtEfficiencyPlots_" << SaveName << "}} \\quad" << std::endl;
    }

    if (makeLaTeX && function_counter % 12 == 0) {
      latexFile3 << "\\caption{Signal and Background Plots}" << std::endl;
      latexFile3 << "\\label{fig:EvtWgtBackgroundPlotss}" << std::endl;
      latexFile3 << "\\end{adjustwidth}" << std::endl;
      latexFile3 << "\\end{figure}" << std::endl;

      latexFile3 << "\\begin{figure}[t]" << std::endl;
      latexFile3 << "\\ContinuedFloat" << std::endl;
      latexFile3 << "\\begin{adjustwidth}{-1cm}{-1cm}" << std::endl;
      latexFile3 << "\\centering" << std::endl;

    }




    double numu_p1;
    double numu_m1;
    double anumu_p1;
    double anumu_m1;
    double nue_p1;
    double nue_m1;
    double nc_p1;
    double nc_m1;
    double outfv_p1;
    double outfv_m1;    
    double cosmic_p1;
    double cosmic_m1;
    double ccother_p1;
    double ccother_m1;

    if (variable == 0 || variable == 1 || variable == 2 || variable == 3 || variable == 4) {
      numu_p1   = histo_numu_p1   ->Integral();
      numu_m1   = histo_numu_m1   ->Integral();
      anumu_p1  = histo_anumu_p1  ->Integral();
      anumu_m1  = histo_anumu_m1  ->Integral();
      nue_p1    = histo_nue_p1    ->Integral();
      nue_m1    = histo_nue_p1    ->Integral();
      nc_p1     = histo_nc_p1     ->Integral();
      nc_m1     = histo_nc_p1     ->Integral();
      outfv_p1  = histo_outfv_p1  ->Integral();
      outfv_m1  = histo_outfv_p1  ->Integral();
      cosmic_p1 = histo_cosmic_p1 ->Integral();
      cosmic_m1 = histo_cosmic_p1 ->Integral();
      ccother_p1 = histo_ccother_p1 ->Integral();
      ccother_m1 = histo_ccother_p1 ->Integral();
    } 
    
    if(makeLaTeX) {
      double numu_diff_p1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      double numu_diff_m1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      
      double anumu_diff_p1 = (anumu_p1 -anumu_nominal)/anumu_nominal * 100.;
      double anumu_diff_m1 = (anumu_m1 -anumu_nominal)/anumu_nominal * 100.;

      double nue_diff_p1 = (nue_p1 -nue_nominal)/nue_nominal * 100.;
      double nue_diff_m1 = (nue_m1 -nue_nominal)/nue_nominal * 100.;
      
      double nc_diff_p1 = (nc_p1 -nc_nominal)/nc_nominal * 100.;
      double nc_diff_m1 = (nc_p1 -nc_nominal)/nc_nominal * 100.;

      double outfv_diff_p1 = (outfv_p1 -outfv_nominal)/outfv_nominal * 100.;
      double outfv_diff_m1 = (outfv_p1 -outfv_nominal)/outfv_nominal * 100.;

      double cosmic_diff_p1 = (cosmic_p1 -cosmic_nominal)/cosmic_nominal * 100.;
      double cosmic_diff_m1 = (cosmic_p1 -cosmic_nominal)/cosmic_nominal * 100.;

      double ccother_diff_p1 = (ccother_p1 -ccother_nominal)/ccother_nominal * 100.;
      double ccother_diff_m1 = (ccother_p1 -ccother_nominal)/ccother_nominal * 100.;

      latexFile << "\\midrule" << std::endl;
      latexFile << "$" << GetLegendName(function_name) << " + 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_p1  << " & " << std::setprecision(2) << numu_diff_p1
      << " & " << std::setprecision(5) << anumu_p1 << " & " << std::setprecision(2) << anumu_diff_p1
      << " & " << std::setprecision(5) << nue_p1   << " & " << std::setprecision(2) << nue_diff_p1
      << " & " << std::setprecision(5) << nc_p1    << " & " << std::setprecision(2) << nc_diff_p1
      << " & " << std::setprecision(5) << outfv_p1 << " & " << std::setprecision(2) << outfv_diff_p1
      << " & " << std::setprecision(5) << cosmic_p1<< " & " << std::setprecision(2) << cosmic_diff_p1 
      << " & " << std::setprecision(5) << ccother_p1<< " & " << std::setprecision(2) << ccother_diff_p1 
      << "\\\\" << std::endl;
      
      latexFile << "$" << GetLegendName(function_name) << " - 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_m1  << " & " << std::setprecision(2) << numu_diff_m1
      << " & " << std::setprecision(5) << anumu_m1 << " & " << std::setprecision(2) << anumu_diff_m1
      << " & " << std::setprecision(5) << nue_m1   << " & " << std::setprecision(2) << nue_diff_m1
      << " & " << std::setprecision(5) << nc_m1    << " & " << std::setprecision(2) << nc_diff_m1
      << " & " << std::setprecision(5) << outfv_m1 << " & " << std::setprecision(2) << outfv_diff_m1
      << " & " << std::setprecision(5) << cosmic_m1<< " & " << std::setprecision(2) << cosmic_diff_m1 
      << " & " << std::setprecision(5) << ccother_m1<< " & " << std::setprecision(2) << ccother_diff_m1 
      << "\\\\" << std::endl;
    }
  }
  
  
  
  
  if(makeLaTeX) {
    latexFile << "\\bottomrule" << std::endl;
    latexFile << "\\end{tabular}" << std::endl;
    latexFile << "\\end{adjustwidth}" << std::endl;
    latexFile << "\\end{table}" << std::endl;

    latexFile3 << "\\caption{Signal and Background Plots}" << std::endl;
    latexFile3 << "\\label{fig:EvtWgtBackgroundPlots}" << std::endl;
    latexFile3 << "\\end{adjustwidth}" << std::endl;
    latexFile3 << "\\end{figure}" << std::endl;
  }
 

  }


  void ReweightingPlotter::MakeXsecDiffPlots(int variable, bool makeLaTeX) {

    //Get the background histograms and efficiency for nominal
    //Creat output directory
    /*system("mkdir, /EvtWgtXSecUnisimPlots");
    //Avoid root to display the cavases
    gROOT->SetBatch(kTRUE);
    
    double scaleFactor = 1.0;
    std::ofstream latexFile;
    std::ofstream latexFile3;
    
    double numu_nominal = _map_bs["signal"].GetNominal().Integral();
    double anumu_nominal = _map_bs["anumu"].GetNominal().Integral();
    double nue_nominal = _map_bs["nue"].GetNominal().Integral();
    double nc_nominal = _map_bs["nc"].GetNominal().Integral();
    double outfv_nominal = _map_bs["outfv"].GetNominal().Integral();
    double cosmic_nominal = _map_bs["cosmic"].GetNominal().Integral();
    double ccother_nominal = _map_bs["cc_other"].GetNominal().Integral();
   
    std::map<std::string, std::vector<TH1D>> map_bs_signal = _map_bs["signal"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_anumu = _map_bs["anumu"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_nue = _map_bs["nue"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_nc = _map_bs["nc"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_outfv = _map_bs["outfv"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_cosmic = _map_bs["cosmic"].UnpackPMHisto();
    std::map<std::string, std::vector<TH1D>> map_bs_ccother = _map_bs["cc_other"].UnpackPMHisto();

    int function_counter  = 0;
    for (auto iter : map_bs_signal){
        function_counter++;
        std::string function_name = iter.first;
        
        auto ii = map_bs_anumu.find(function_name);
        std::vector<TH1D> anumu_v = ii->second;
        const TH1D* anumu_nom = &_map_bs["anumu"].GetNominal();

        ii = map_bs_nue.find(function_name);
        std::vector<TH1D> nue_v = ii->second;
        const TH1D* nue_nom = &_map_bs["nue"].GetNominal();

        ii = map_bs_nc.find(function_name);
        std::vector<TH1D> nc_v = ii->second;
        const TH1D* nc_nom = &_map_bs["nc"].GetNominal();

        ii = map_bs_outfv.find(function_name);
        std::vector<TH1D> outfv_v = ii->second;
        const TH1D* outfv_nom = &_map_bs["outfv"].GetNominal();

        ii = map_bs_cosmic.find(function_name);
        std::vector<TH1D> cosmic_v = ii->second;
        const TH1D* cosmic_nom = &_map_bs["cosmic"].GetNominal();
        
        ii = map_bs_ccother.find(function_name);
        std::vector<TH1D> ccother_v = ii->second;
        const TH1D* ccother_nom = &_map_bs["cc_other"].GetNominal();

        TString SaveName; 
        if (variable == 0 ) SaveName = "Pmu_" + function_name;
        if (variable == 1 ) SaveName = "CosThetaMu_" + function_name;
        if (variable == 2 ) SaveName = "Proton_" + function_name;
        if (variable == 3 ) SaveName = "CosThetaP_" + function_name;
        if (variable == 4 ) SaveName = "Thetamup_" + function_name;
        //c->SetLogy();
        //c->SetGridy();

        TH1D* histo_numu;
        TH1D* histo_numu_p1;
        TH1D* histo_numu_m1;
        TH1D* histo_anumu;
        TH1D* histo_anumu_p1;
        TH1D* histo_anumu_m1;
        TH1D* histo_nue;
        TH1D* histo_nue_p1;
        TH1D* histo_nue_m1;
        TH1D* histo_nc;
        TH1D* histo_nc_p1;
        TH1D* histo_nc_m1;
        TH1D* histo_outfv;
        TH1D* histo_outfv_p1;
        TH1D* histo_outfv_m1;
        TH1D* histo_cosmic;
        TH1D* histo_cosmic_p1;
        TH1D* histo_cosmic_m1;
        TH1D* histo_ccother;
        TH1D* histo_ccother_p1;
        TH1D* histo_ccother_m1;





        if (variable == 0 || variable == 1 || variable == 2 || variable == 3 || variable == 4) {
        histo_numu        = (TH1D*)_map_bs["signal"].GetNominal().Clone("histo_numu");
        histo_numu_p1     = (TH1D*)iter.second.at(1).Clone("histo_numu_p1");
        histo_numu_m1     = (TH1D*)iter.second.at(0).Clone("histo_numu_m1");
      
        histo_anumu       = (TH1D*) (anumu_nom->Clone("histo_anumu"));
        histo_anumu_p1    = (TH1D*) (anumu_v.at(1).Clone("histo_anumu_p1"));
        histo_anumu_m1    = (TH1D*) (anumu_v.at(0).Clone("histo_anumu_m1"));
      
        histo_nue         = (TH1D*) (nue_nom->Clone("histo_nue"));
        histo_nue_p1      = (TH1D*) (nue_v.at(1).Clone("histo_nue_p1"));
        histo_nue_m1      = (TH1D*) (nue_v.at(0).Clone("histo_nue_m1"));
      
        histo_nc          = (TH1D*) (nc_nom->Clone("histo_nc"));
        histo_nc_p1       = (TH1D*) (nc_v.at(1).Clone("histo_nc_p1"));
        histo_nc_m1       = (TH1D*) (nc_v.at(0).Clone("histo_nc_m1"));

        histo_outfv       = (TH1D*) (outfv_nom->Clone("histo_outfv"));
        histo_outfv_p1    = (TH1D*) (outfv_v.at(1).Clone("histo_outfv_p1"));
        histo_outfv_m1    = (TH1D*) (outfv_v.at(0).Clone("histo_outfv_m1"));

        histo_cosmic      = (TH1D*) (cosmic_nom->Clone("histo_cosmic"));
        histo_cosmic_p1   = (TH1D*) (cosmic_v.at(1).Clone("histo_cosmic_p1"));
        histo_cosmic_m1   = (TH1D*) (cosmic_v.at(0).Clone("histo_cosmic_m1"));

        histo_ccother      = (TH1D*) (ccother_nom->Clone("histo_ccother"));
        histo_ccother_p1   = (TH1D*) (ccother_v.at(1).Clone("histo_ccother_p1"));
        histo_ccother_m1   = (TH1D*) (ccother_v.at(0).Clone("histo_ccother_m1"));
        }
        else {
        std::cout << "Invalid option. Exit." << std::endl;
        exit(0);
        }

        //Get the histogram of the denominator
        std::map<std::string, std::vector<TH1D>> histo_map_den; //a map from function name to 2 TH1D
        std::map<std::string, std::vector<TH1D>> histo_map; // only for eff plot, a map from function name
        //Get the histograms of data file 
        

        //Calculate efficiency and do the background subtraction
        TH1D *histo; TH1D *histo_p1; TH1D *histo_m1;

        histo = (TH1D*)_bs_eff_num.GetNominal().Clone("nominal_eff");
        TH1D * temp_eff_den = (TH1D*)_bs_eff_den.GetNominal().Clone("temp_eff_den");  
        histo->Divide(temp_eff_den); //get the efficiency of the nominal value
       

        histo_map = _bs_eff_num.UnpackPMHisto();
        histo_map_den = _bs_eff_den.UnpackPMHisto();

 
        //std::ofstream outfile_genieunisim;
        //outfile_genieunisim.open("./EvtWgtXSecUnisimPlots");
        
        //open latex file to write the table
        std::ofstream latexFile;
        if(makeLaTeX){
              if(variable == 0) latexFile.open("./EvtWgtXSecUnisimPlots/evtwgtXSecPmu.tex");
              if(variable == 1) latexFile.open("./EvtWgtXSecUnisimPlots/evtwgtXSecCosThetaMu.tex");
              if(variable == 2) latexFile.open("./EvtWgtXSecUnisimPlots/evtwgtXSecPproton.tex");
              if(variable == 3) latexFile.open("./EvtWgtXSecUnisimPlots/evtwgtXSecCosThetaP.tex");
              if(variable == 4) latexFile.open("./EvtWgtXSecUnisimPlots/evtwgtXSecThetaMuP.tex");
              latexFile <<"\\begin{table}[]"<<std::endl;
              latexFile <<"\\caption{}"<<std::endl;
              latexFile <<"\\captionsetup{format=hang,labelfont={sf,bf}}"<<std::endl;
              latexFile <<"\\label{tab:}"<<std::endl;
              latexFile <<"\\centering"<<std::endl;
              latexFile <<"\\begin{tabular}{ccc}"<<std::endl;
              latexFile <<"\\toprule"<<std::endl;
              latexFile <<"\\ & xsec & Difference (\\%) \\\\"<<std::endl;
              latexFile <<"\\midrule"<<std::endl;
                            
        }    
        histo_p1 = (TH1D*) iter.second.at(1).Clone("p1");
        histo_m1 = (TH1D*) iter.second.at(0).Clone("m1");
        TH1D * temp_eff_den_p1 = (TH1D*)histo_map_den[function_name].at(1).Clone("temp_eff_den_p1");
        histo_p1->Divide(temp_eff_den_p1);
        
        TH1D * temp_eff_den_m1 = (TH1D*)histo_map_den[function_name].at(0).Clone("temp_eff_den_m1");
        histo_m1->Divide(temp_eff_den_m1);
         

        TCanvas *c = new TCanvas("c", "canvas", 1600, 800);
        TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.52, 0.95, 0.97);
        pad1->SetBottomMargin(0); //upper and lower plot are joint"
        pad1->SetGridx();
        pad1->Draw();
        pad1->cd();
        histo_p1->SetStats(0);
        histo_m1->SetStats(0);
        histo_p1->Draw("histo");
        histo->Draw("histo same");
        histo_m1->Draw("histo same");
        */ 
        //Get the smearing matrix
        /*
        //calculate cross section with the efficiency, signal and background
        CrossSectionBootstrapCalculator1D _xsec_bs_calc_temp;
        _xsec_bs_calc_temp.SetFluxCorrectionWeight(flux_correction_weight);
        _xsec_bs_calc_temp.Reset();
        _xsec_bs_calc_temp.SetScaleFactors(scale_Factor_mc_bnbcosmic, scale_Factor_bnbon, scale_factor_extbnb, scale_factor_mc_dirt);
        _xsec_bs_calc_temp.SetPOT(bnbon_pot_meas); 
        _xsec_bs_calc_temp.MigrationMatrixDimensions(n_bins_temp, n_bins_temp);
        _xsec_bs_calc_temp.SetBkgToSubtrac(bkg_name);
        
        CrossSectionBootstrapCalculator1D _xsec_bs_calc_temp_p1;
        _xsec_bs_calc_temp_p1.SetFluxCorrectionWeight(flux_correction_weight);
        _xsec_bs_calc_temp_p1.Reset();
        _xsec_bs_calc_temp_p1.SetScaleFactors(scale_Factor_mc_bnbcosmic, scale_Factor_bnbon, scale_factor_extbnb, scale_factor_mc_dirt);
        _xsec_bs_calc_temp_p1.SetPOT(bnbon_pot_meas); 
        _xsec_bs_calc_temp_p1.MigrationMatrixDimensions(n_bins_temp, n_bins_temp);
        _xsec_bs_calc_temp_p1.SetBkgToSubtrac(bkg_name);

        CrossSectionBootstrapCalculator1D _xsec_bs_calc_temp_m1;
        _xsec_bs_calc_temp_m1.SetFluxCorrectionWeight(flux_correction_weight);
        _xsec_bs_calc_temp_m1.Reset();
        _xsec_bs_calc_temp_m1.SetScaleFactors(scale_Factor_mc_bnbcosmic, scale_Factor_bnbon, scale_factor_extbnb, scale_factor_mc_dirt);
        _xsec_bs_calc_temp_m1.SetPOT(bnbon_pot_meas); 
        _xsec_bs_calc_temp_m1.MigrationMatrixDimensions(n_bins_temp, n_bins_temp);
        _xsec_bs_calc_temp_m1.SetBkgToSubtrac(bkg_name);
        */
        /*double bins_mumom[7] = {0.00, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
        double bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
        int n_bins_mumom = 6;
        int n_bins_mucostheta = 12;


        double bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.50};
        double bins_pcostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};

        int n_bins_pmom = 10; 
        int n_bins_pcostheta = 9;

        double bins_muptheta[7] = {0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14};
        int n_bins_muptheta = 6;
  

        int n_bins_temp = -999;
        TString name_temp = "";
        TString var_temp = "";
        TString varx_temp = "";
        TString vary_temp = "";
        TH1D *h_truth_xsec_temp;
        */
        /*if(variable == 0) {
            n_bins_temp = n_bins_mumom;
            name_temp = "trkmom";
            var_temp = ";p_{#mu}^{reco}[GeV];SelectedEvents";
            varx_temp = "p_{#mu}^{reco} [GeV]";
            vary_temp = "d#sigma/dp_{#mu}^{reco} [10^{-38} cm^{2}/GeV]";
            h_truth_xsec_temp = (TH1D*)h_truth_xsec_mumom.Clone();
        }
        if(variable == 1) {
            n_bins_temp = n_bins_mucostheta;
            name_temp = "trkcostheta";
            var_temp = ";cos#theta_{#mu}^{reco};SelectedEvents";
            varx_temp = "cos#theta_{#mu}^{reco}";
            vary_temp = "d#sigma/dcos#theta_{#mu}^{reco} [10^{-38} cm^{2}]";
            h_truth_xsec_temp = (TH1D*)h_truth_xsec_mucostheta.Clone();
        }
        if(variable == 10) {
            n_bins_temp = n_bins_pmom;
            name_temp = "trkpmom";
            var_temp = ";p_{proton}^{reco}[GeV];SelectedEvents";
            varx_temp = "p_{proton}^{reco} [GeV]";
            vary_temp = "d#sigma/dp_{proton}^{reco} [10^{-38} cm^{2}/GeV]";
            h_truth_xsec_temp = (TH1D*)h_truth_xsec_pmom.Clone();
        }
        if(variable == 11) {
            n_bins_temp = n_bins_pcostheta;
            name_temp = "trkpcostheta";
            var_temp = ";cos#theta_{proton}^{reco}[GeV];SelectedEvents";
            varx_temp = "cos#theta_{proton}^{reco} [GeV]";
            vary_temp = "d#sigma/dcos#tehta_{proton}^{reco} [10^{-38} cm^{2}]";
            h_truth_xsec_temp = (TH1D*)h_truth_xsec_pangle.Clone();
        }
        if(variable == 100) {
            n_bins_temp = n_bins_muptheta;
            name_temp = "thetamup";
            var_temp = ";#theta_{#mu, p}^{reco};SelectedEvents";
            varx_temp = "#theta_{#mu, p}^{reco}";
            vary_temp = "d#sigma/d#theta_{#mu, p}^{reco} [10^{-38} cm^{2}/Rad]";
            h_truth_xsec_temp = (TH1D*)h_truth_xsec_thetamup.Clone();
        }






        TMatrix S_2d; S_2d.Clear(); S_2d.ResizeTo(n_bins_temp+1, n_bins_temp+1);
        MigrationMatrix2D migrationmatrix2d;
        migrationmatrix2d.SetNbins(n_bins_temp, n_bins_temp);
        */
        /*
        xsec_calc_temp.Reset();
        xsec_calc_temp.set_verbosity(Base::msg::kINFO);
        
        xsec_calc_temp.SetTruthHistograms(h_truth_xsec_temp);
        xsec_calc_temp.SetNameAndLabel(name_temp, var_temp);

        TH1D * xsec_varname_temp = _xsec_calc_temp.ExtractCrossSection(bkg_names, varx_temp;vary_temp);
        TH1D * xsec_varname_mc_temp.xsec_calc.GetMCCrossSection();
        xsec_varname_temp->Write(save_name.c_str());



        TH1D *temp_xsec_p1; TH1D *temp_xsec; TH1D *temp_xsec_m1;
          
        

        TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.02, 0.95, 0.47);
        pad2->SetBottomMargin(0);
        pad2->SetGridX();
        pad2->Draw();
        pad2->cd();
        temp_xsec_p1->Draw("histo");
        temp_xsec->Draw("histo same")
        temp_xsec_m1->Draw("histo same");
        */
     

    //} //end of loop over all the functions


  }

  TString ReweightingPlotter::GetLegendName(std::string fName) {
    
    TString legName = "null";
    if (fName.find("qema") != std::string::npos) legName = "M_{A}^{CCQE}";
    if (fName.find("ncelEta") != std::string::npos) legName = "#eta^{NCEL}";
    if (fName.find("ncelAxial") != std::string::npos) legName = "M_{A}^{NCEL}";
    if (fName.find("qevec") != std::string::npos) legName = "VecFF";
    if (fName.find("ccresAxial") != std::string::npos) legName = "M_{A}^{CCRES}";
    if (fName.find("ccresVector") != std::string::npos) legName = "M_{V}^{CCRES}";
    if (fName.find("ncresAxial") != std::string::npos) legName = "M_{A}^{NCRES}";
    if (fName.find("ncresVector") != std::string::npos) legName = "M_{V}^{NCRES}";
    if (fName.find("cohMA") != std::string::npos) legName = "M_{A}^{COH#pi}";
    if (fName.find("cohR0") != std::string::npos) legName = "R_{0}^{COH#pi}";
    if (fName.find("NonResRvp1pi") != std::string::npos) legName = "NonResRvp1pi";
    if (fName.find("NonResRvbarp1pi") != std::string::npos) legName = "NonResRvbarp1pi";
    if (fName.find("NonResRvp2pi") != std::string::npos) legName = "NonResRvp2pi";
    if (fName.find("NonResRvbarp2pi") != std::string::npos) legName = "NonResRvbarp2pi";
    if (fName.find("ResDecayGamma") != std::string::npos) legName = "BR(#gamma)";
    if (fName.find("ResDecayEta") != std::string::npos) legName = "BR(#eta)";
    if (fName.find("ResDecayTheta") != std::string::npos) legName = "Ang. distr.";
    if (fName.find("DISAth") != std::string::npos) legName = "A_{HT}^{BY}";
    if (fName.find("DISBth") != std::string::npos) legName = "B_{HT}^{BY}";
    if (fName.find("DISCv1u") != std::string::npos) legName = "C_{V1u}^{BY}";
    if (fName.find("DISCv2u") != std::string::npos) legName = "C_{V2u}^{BY}";
    if (fName.find("AGKYxF") != std::string::npos) legName = "x_{F}";
    if (fName.find("AGKYpT") != std::string::npos) legName = "p_{T}";
    if (fName.find("FormZone") != std::string::npos) legName = "FZ";
    if (fName.find("FermiGasModelKf") != std::string::npos) legName = "CCQE-PauliSup";
    if (fName.find("FermiGasModelSf") != std::string::npos) legName = "Fermi Gas to SF";
    if (fName.find("IntraNukeNmfp") != std::string::npos) legName = "x_{mfp}^{N}";
    if (fName.find("IntraNukeNcex") != std::string::npos) legName = "x_{cex}^{N}";
    if (fName.find("ntraNukeNel") != std::string::npos) legName = "x_{el}^{N}";
    if (fName.find("ntraNukeNinel") != std::string::npos) legName = "x_{inel}^{N}";
    if (fName.find("ntraNukeNabs") != std::string::npos) legName = "x_{abs}^{N}";
    if (fName.find("ntraNukeNpi") != std::string::npos) legName = "x_{pi}^{N}";
    if (fName.find("IntraNukePImfp") != std::string::npos) legName = "x_{mfp}^{PI}";
    if (fName.find("IntraNukePIcex") != std::string::npos) legName = "x_{cex}^{PI}";
    if (fName.find("ntraNukePIel") != std::string::npos) legName = "x_{el}^{PI}";
    if (fName.find("ntraNukePIinel") != std::string::npos) legName = "x_{inel}^{PI}";
    if (fName.find("ntraNukePIabs") != std::string::npos) legName = "x_{abs}^{PI}";
    if (fName.find("ntraNukePIpi") != std::string::npos) legName = "x_{pi}^{PI}";
  
    return legName;
  }
  
}

#endif
