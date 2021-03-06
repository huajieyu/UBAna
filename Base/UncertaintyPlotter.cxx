#ifndef __BASE_UNCERTAINTYPLOTTER_CXX__
#define __BASE_UNCERTAINTYPLOTTER_CXX__

#include "UncertaintyPlotter.h"

namespace Base {

  void UncertaintyPlotter::Reset()
  {
    _cov_names.clear();
    _frac_cov_v.clear();
  }

  void UncertaintyPlotter::AddFracCovarianceMatrix(std::string name, TH2D cov)
  {
    _cov_names.push_back(name);
    _frac_cov_v.push_back(cov);
  }

  void UncertaintyPlotter::MakePlot(std::string file_name, bool no_legend)
  {
    if (_2d_case) {
      MakePlot2D(file_name, no_legend);
    } else {
      MakePlot1D(file_name, no_legend);
    }
  }

  void UncertaintyPlotter::MakePlot2D(std::string file_name, bool no_legend)
  {

    std::vector<std::string> costhetamu_ranges = {"-1.00 #leq cos(#theta_{#mu}^{reco}) < -0.50",
                                                  "-0.50 #leq cos(#theta_{#mu}^{reco}) < 0.00",
                                                  "0.00 #leq cos(#theta_{#mu}^{reco}) < 0.27",
                                                  "0.27 #leq cos(#theta_{#mu}^{reco}) < 0.45",
                                                  "0.45 #leq cos(#theta_{#mu}^{reco}) < 0.62",
                                                  "0.62 #leq cos(#theta_{#mu}^{reco}) < 0.76",
                                                  "0.76 #leq cos(#theta_{#mu}^{reco}) < 0.86",
                                                  "0.86 #leq cos(#theta_{#mu}^{reco}) < 0.94",
                                                  "0.94 #leq cos(#theta_{#mu}^{reco}) < 1.00",
                                                  "nan #leq cos(#theta_{#mu}^{reco}) < nan",
                                                  "nan #leq cos(#theta_{#mu}^{reco}) < nan",};


    int x_bins = _xsec_2d->GetNBinsX();
    LOG_INFO() << "n bins x " << x_bins << std::endl;

    int horizontal_division = 2;
    int vertical_division = floor(x_bins / 2.);

    if (x_bins / 2. != floor(x_bins / 2.)) vertical_division++;

    LOG_INFO() << "Horizontal divisions " << horizontal_division << std::endl;
    LOG_INFO() << "Vertical divisions " << vertical_division << std::endl;

    TCanvas *c_xsec_split = new TCanvas("c_xsec_split", "multipads", 0, 0, 1144,746);
    c_xsec_split->Divide(3, 3, 0.01, 0.01);

    // TCanvas *c_xsec_split = new TCanvas("c_xsec_split", "multipads",0,45,1006,1150);
    // c_xsec_split->Divide(horizontal_division, vertical_division, 0.02, 0.01);

    std::vector<std::vector<int>> bin_numbers;


    std::vector<TH1D> _xsec_data_histos;
    std::vector<std::vector<TH1D>> _xsec_unc_histos; // Per systematic, per costheta bin
    _xsec_unc_histos.resize(_cov_names.size());


    for (int i = 0; i < x_bins; i++) {
      std::vector<int> b_numbers;
      _xsec_data_histos.emplace_back(*_xsec_2d->ProjectionY(("fuck_"+std::to_string(i)).c_str(), i+1, b_numbers));
      bin_numbers.emplace_back(b_numbers);
    }

    for (size_t c = 0; c < _cov_names.size(); c++) {
      for (int i = 0; i < x_bins; i++) {
        TString s;
        s.Form ("h_unc_%d", i); // = "h_unc_" + i;
        TH1D h = *(TH1D*)_xsec_data_histos.at(i).Clone(s);
        int counter = 1;
        for (int j = bin_numbers.at(i).front(); j <= bin_numbers.at(i).back(); j++) {
          // if (_cov_names.at(c) == "DETECTOR") std::cout << "bin number = " << j << ", " <<  << std::endl;
          h.SetBinContent(counter, std::sqrt(_frac_cov_v.at(c).GetBinContent(j, j)) * 100. );
          counter++;
        }
        _xsec_unc_histos.at(c).emplace_back(h);
      }
    }


    //
    // Draw
    //

    // TLegend *l = new TLegend(0.1448987,0.1039989,0.7901438,0.9024423,NULL,"brNDC");
    TLegend *l = new TLegend(0.3248686,0.3708945,0.9927898,0.8356745,NULL,"brNDC");
    l->SetFillColor(0);
    l->SetFillStyle(0);
    // l->SetTextSize(0.1140633);
    l->SetTextSize(0.06639715);

    for (size_t i = 0; i < _xsec_unc_histos.at(0).size(); i++) {

      c_xsec_split->cd(i+1);

      gPad->SetBottomMargin(0.16);
      gPad->SetLeftMargin(0.18);
      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.06);

      // gPad->SetRightMargin(0.14);
      // gPad->SetTopMargin(0.13);

      gStyle->SetTitleFontSize(0.07);
      gStyle->SetTitleStyle(0);

      double max = _xsec_unc_histos.at(_cov_names.size()-1).at(i).GetMaximum();

      for (size_t c = 0; c < _cov_names.size(); c++) {

        if (i==0) l->AddEntry(&_xsec_unc_histos.at(c).at(i), _cov_names.at(c).c_str(), "l");

        _xsec_unc_histos.at(c).at(i).SetTitle("");
        // _xsec_unc_histos.at(c).at(i).SetTitle(costhetamu_ranges.at(i).c_str());
        _xsec_unc_histos.at(c).at(i).GetXaxis()->SetTitle(_xaxis_title.c_str());
        _xsec_unc_histos.at(c).at(i).GetYaxis()->SetTitle("Relative Uncertainty [%]");
        _xsec_unc_histos.at(c).at(i).GetXaxis()->CenterTitle();
        _xsec_unc_histos.at(c).at(i).GetYaxis()->CenterTitle();
        _xsec_unc_histos.at(c).at(i).SetLineColor(_line_colors.at(c));
        _xsec_unc_histos.at(c).at(i).SetLineStyle(_line_styles.at(c));
        _xsec_unc_histos.at(c).at(i).SetFillColor(0); // white

        _xsec_unc_histos.at(c).at(i).GetXaxis()->SetTitleOffset(0.92);
        _xsec_unc_histos.at(c).at(i).GetXaxis()->SetTitleSize(0.07);
        _xsec_unc_histos.at(c).at(i).GetXaxis()->SetLabelSize(0.06);

        _xsec_unc_histos.at(c).at(i).GetYaxis()->SetTitleOffset(0.90);
        _xsec_unc_histos.at(c).at(i).GetYaxis()->SetTitleSize(0.06);
        _xsec_unc_histos.at(c).at(i).GetYaxis()->SetLabelSize(0.06);

        _xsec_unc_histos.at(c).at(i).SetMinimum(0);
        _xsec_unc_histos.at(c).at(i).SetMaximum(max * 1.2);
        // if (i < 5) {
        //   _xsec_unc_histos.at(c).at(i).SetMinimum(0);
        //   _xsec_unc_histos.at(c).at(i).SetMaximum(200);
        // }
        _xsec_unc_histos.at(c).at(i).Draw("histo same");
      }


      TLatex * costheta_label = new TLatex(0.9500818,0.8657913, costhetamu_ranges.at(i).c_str());
      costheta_label->SetNDC();
      costheta_label->SetTextSize(0.07);
      costheta_label->SetTextAlign(32);
      costheta_label->SetLineWidth(2);
      costheta_label->SetTextFont(42);
      costheta_label->Draw();

    }

    // Legend
    c_xsec_split->cd(_xsec_unc_histos.at(0).size() + 1);
    if (!no_legend) l->Draw();

    c_xsec_split->SaveAs((file_name + ".pdf").c_str());
    c_xsec_split->SaveAs((file_name + ".C").c_str());


  }

  

  void UncertaintyPlotter::MakePlot1D(std::string file_name, bool no_legend)
  {

    

    std::vector<TH1D> histos;
    histos.resize(_cov_names.size());

    for (size_t i = 0; i < _cov_names.size(); i++)
    {
      histos.at(i) = * (TH1D*)_xsec_1d.Clone((_cov_names.at(i) + "_h").c_str());
      histos.at(i).Reset();

      for (int j = 1; j <= _xsec_1d.GetNbinsX(); j++)
      {
        histos.at(i).SetBinContent(j, std::sqrt(_frac_cov_v.at(i).GetBinContent(j, j)) * 100.);
      }
    }


    TCanvas * c = new TCanvas();

    c->SetBottomMargin(0.12);

    TLegend *l = new TLegend(0.3825215,0.65,0.7363897,0.8547368,NULL,"brNDC");
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.03578947);


    for (size_t i = 0; i < _cov_names.size(); i++) 
    {
      if (i == 0) {
        histos.at(i).SetMinimum(0);
        histos.at(i).SetMaximum(80);
      }
      histos.at(i).SetLineColor(_line_colors.at(i));
      histos.at(i).GetXaxis()->SetTitle(_xaxis_title.c_str());
      histos.at(i).GetYaxis()->SetTitle("Relative Uncertainty [%]");
      histos.at(i).Draw("same");

      l->AddEntry(&histos.at(i), _cov_names.at(i).c_str(), "l");
    }

    if (!no_legend) l->Draw();

    PlottingTools::DrawSimulationXSec();

    c->SaveAs((file_name + ".pdf").c_str());
    c->SaveAs((file_name + ".C").c_str());

  }


}

#endif
