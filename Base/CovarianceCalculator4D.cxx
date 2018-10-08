#ifndef __BASE_COVARIANCECALCULATOR4D_CXX__
#define __BASE_COVARIANCECALCULATOR4D_CXX__

#include "CovarianceCalculator4D.h"
#include "PlottingTools.h"

namespace Base {
  

  void CovarianceCalculator4D::SetBootstrap(BootstrapTH2D bs)
  {
    LOG_DEBUG() << "Seting up via BootstrapTH2D." << std::endl;
    _bs = bs;
    _polybin_mode = false;
  }

  void CovarianceCalculator4D::SetBootstrap(BootstrapTH2DPoly bs)
  {
    LOG_DEBUG() << "Seting up via BootstrapTH2DPoly (poly bins)." << std::endl;
    _bs_poly = bs;
    _polybin_mode = true;
  }

  void CovarianceCalculator4D::SetPrefix(std::string prefix) 
  {
    _prefix = prefix;
  }

  void CovarianceCalculator4D::GetCovarianceMatrix(TH2D &h)
  {
    h = _M_h;
  }

  void CovarianceCalculator4D::GetFractionalCovarianceMatrix(TH2D &h)
  {
    h = _M_frac_h;
  }

  void CovarianceCalculator4D::AddExtraDiagonalUncertainty(double value) 
  {
    _extra_relative_uncertainty = value;
  }

  void CovarianceCalculator4D::CalculateCovarianceMatrix() 
  {
    if (_polybin_mode) {
      this->CalculateCovarianceMatrixPoly();
    }
    else {
      this->CalculateCovarianceMatrixNormal();
    }

  }


  void CovarianceCalculator4D::CalculateCovarianceMatrixNormal() 
  {

    int n_bins_x = _bs.GetNbinsX();
    int n_bins_y = _bs.GetNbinsY();


    // Resize the covariance matrix
    _M.resize(n_bins_x, 
              std::vector<std::vector<std::vector<double>>> (n_bins_y,
                                                             std::vector<std::vector<double>>(n_bins_x,
                                                                                             std::vector<double> (n_bins_y, 0.
                                                                                                                 )
                                                                                            )
                                                            )
              );

    // Resize the fractional covariance matrix
    _M_frac.resize(n_bins_x, 
              std::vector<std::vector<std::vector<double>>> (n_bins_y,
                                                             std::vector<std::vector<double>>(n_bins_x,
                                                                                             std::vector<double> (n_bins_y, 0.
                                                                                                                 )
                                                                                            )
                                                            )
              );

    // Resize the correlation matrix
    _RHO.resize(n_bins_x, 
              std::vector<std::vector<std::vector<double>>> (n_bins_y,
                                                             std::vector<std::vector<double>>(n_bins_x,
                                                                                             std::vector<double> (n_bins_y, 0.
                                                                                                                 )
                                                                                            )
                                                            )
              );



    int counter = 0;
    int loop_length = _bs.GetNbinsX() * _bs.GetNbinsX() * _bs.GetNbinsY() * _bs.GetNbinsY();

    _bs.ResetIterator();
    LOG_NORMAL() << " Calculating Cov Matrix with " << _bs.GetNUniverses() << " universes. Nominal histogram is excluded." << std::endl;

    double number_of_universes = (double)_bs.GetNUniverses() - 1;

    for (int i = 0; i < _bs.GetNbinsX(); i++) {

      for (int j = 0; j < _bs.GetNbinsY(); j++) {

        for (int m = 0; m < _bs.GetNbinsX(); m++) {

          for (int n = 0; n < _bs.GetNbinsY(); n++) {

            counter++;
            PlottingTools::DrawProgressBar((double)counter/(double)loop_length, 70);

            // std::cout << "i = " << i << ", j = " << j << ", m = " << m << ", n = " << n << std::endl;

            // Reset the matrix element
            _M[i][j][m][n] = 0.;
            _M_frac[i][j][m][n] = 0.;
            _RHO[i][j][m][n] = 0.;

            // Nominal cross section for bin ij and mn 
            double N_ij_cv = _bs.GetNominal().GetBinContent(i+1, j+1);
            double N_mn_cv = _bs.GetNominal().GetBinContent(m+1, n+1);

            _bs.ResetIterator();

            for (int s = 0; s < number_of_universes; s++) {

              double N_ij_s = _bs.NextUniverse().GetBinContent(i+1, j+1);
              double N_mn_s = _bs.SameUniverse().GetBinContent(m+1, n+1);

              _M[i][j][m][n] += (N_ij_s - N_ij_cv) * (N_mn_s - N_mn_cv) / number_of_universes;

            } // universe loop

            _M_frac[i][j][m][n] = _M[i][j][m][n] / (N_ij_cv * N_mn_cv);

            if ( (i == m) && (j == n)) { // Diagonal
              _M_frac[i][j][m][n] += _extra_relative_uncertainty * _extra_relative_uncertainty;
              _M[i][j][m][n] += (N_ij_cv * _extra_relative_uncertainty) * (N_mn_cv * _extra_relative_uncertainty);
            }

          } // bin n loop
        } // bin m loop
      } // bin j loop
    } // bin i loop




    // if (_verbose) {
    //   std::cout << _name << "Printing Covariance Matrix M = " << std::endl;
    //   for (int i = 0; i < _bs.GetNbinsX(); i++) {
    //     for (int j = 0; j < _bs.GetNbinsY(); j++) {
    //       for (int m = 0; m < _bs.GetNbinsX(); m++) {
    //         for (int n = 0; n < _bs.GetNbinsY(); n++) { 
    //           std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _M[i][j][m][n] << std::endl;
    //         }
    //       }
    //     }
    //   }
    //   std::cout << _name << "Printing Fractional Covariance Matrix M = " << std::endl;
    //   for (int i = 0; i < _bs.GetNbinsX(); i++) {
    //     for (int j = 0; j < _bs.GetNbinsY(); j++) {
    //       for (int m = 0; m < _bs.GetNbinsX(); m++) {
    //         for (int n = 0; n < _bs.GetNbinsY(); n++) { 
    //           std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _M_frac[i][j][m][n] << std::endl;
    //         }
    //       }
    //     }
    //   }
    // }





    for (int i = 0; i < _bs.GetNbinsX(); i++) {

      for (int j = 0; j < _bs.GetNbinsY(); j++) {

        for (int m = 0; m < _bs.GetNbinsX(); m++) {

          for (int n = 0; n < _bs.GetNbinsY(); n++) {

            LOG_DEBUG() << "This is " << i << ", " << j << ", " << m << ", " << n << ", cov matrix is " << _M[j][j][m][n] << std::endl;

            _RHO[i][j][m][n] += _M[i][j][m][n] / (std::sqrt(_M[i][j][i][j]) * std::sqrt(_M[m][n][m][n]));

            if (_RHO[i][j][m][n] < -1. || _RHO[i][j][m][n] > 1.) {
              std::cout << "WARNING!!! Corraltion Matrix rho is smaller than -1 or greater than +1, value: _RHO[" << i << "][" << j << "][" << m << "][" << n << "]" << _RHO[i][j][m][n] << std::endl;
            }
          } // bin n loop
        } // bin m loop
      } // bin i loop
    } // bin i loop


    

  }





  void CovarianceCalculator4D::CalculateCovarianceMatrixPoly() 
  {

    int n_bins = _bs_poly.GetNumberOfBins();

    LOG_DEBUG() << "Start covariance matrix calculation with " << n_bins << " bins." << std::endl;

    // Resize the covariance matrix
    _M_p.ResizeTo(n_bins, n_bins);

    // Resize the fractional covariance matrix
    _M_frac_p.ResizeTo(n_bins, n_bins);

    // Resize the correlation matrix
    _RHO_p.ResizeTo(n_bins, n_bins);


    LOG_DEBUG() << "Matrices have been resized to " << n_bins << " bins." << std::endl;
    
    int counter = 0;
    int loop_length = n_bins * n_bins;

    _bs_poly.ResetIterator();

    LOG_NORMAL() << "Calculating Cov Matrix with " << _bs_poly.GetNUniverses() -1 << " universes." << std::endl;

    double number_of_universes = (double)_bs_poly.GetNUniverses() - 1;

    for (int i = 0; i < n_bins; i++) {

      for (int j = 0; j < n_bins; j++) {

        counter++;
        PlottingTools::DrawProgressBar((double)counter/(double)loop_length, 70);

        // Reset the matrix element
        _M_p[i][j] = 0.;
        _M_frac_p[i][j] = 0.;
        _RHO_p[i][j] = 0.;

        // Nominal cross section for bin ij and mn 
        double N_i_cv = _bs_poly.GetNominal()->GetBinContent(i+1);
        double N_j_cv = _bs_poly.GetNominal()->GetBinContent(j+1);

    
        _bs_poly.ResetIterator();

        for (int s = 0; s < number_of_universes; s++) {

          double N_i_s = _bs_poly.NextUniverse().GetBinContent(i+1);
          double N_j_s = _bs_poly.SameUniverse().GetBinContent(j+1);


          _M_p[i][j] += (N_i_s - N_i_cv) * (N_j_s - N_j_cv) / number_of_universes;
          

        } // universe loop

        _M_frac_p[i][j] = _M_p[i][j] / (N_i_cv * N_j_cv);

        if (i == j) { // Diagonal
          _M_frac_p[i][j] += _extra_relative_uncertainty * _extra_relative_uncertainty;
          _M_p[i][j] += (N_i_cv * _extra_relative_uncertainty) * (N_j_cv * _extra_relative_uncertainty);
        }

      } // bin j loop
    } // bin i loop
    std::cout << std::endl;


    

    LOG_DEBUG() << _name << "Printing Covariance Matrix M = " << std::endl;
    for (int i = 0; i < n_bins; i++) {
      for (int j = 0; j < n_bins; j++) {
        LOG_DEBUG() << "(" << i << ", " << j << ") => " << _M_p[i][j] << std::endl;
      }
    }

    LOG_DEBUG() << _name << "Printing Fractional Covariance Matrix M = " << std::endl;
    for (int i = 0; i < n_bins; i++) {
      for (int j = 0; j < n_bins; j++) {
        LOG_DEBUG() << "(" << i << ", " << j << ") => " << _M_frac_p[i][j] << std::endl;
      }
    }





    for (int i = 0; i < n_bins; i++) {

      for (int j = 0; j < n_bins; j++) {

        _RHO_p[i][j] += _M_p[i][j] / (std::sqrt(_M_p[i][i]) * std::sqrt(_M_p[j][j]));

        // if (_RHO_p[i][j] < -1 || _RHO_p[i][j] > 1) {
        //   std::cout << "WARNING!!! Corraltion Matrix rho is smaller than -1 or greater than +1, value: _RHO[" << i << "][" << j << "]" << _RHO_p[i][j] << std::endl;
        // }

      } // bin i loop
    } // bin i loop


  }



  void CovarianceCalculator4D::PlotMatrices()
  {
    if (_polybin_mode) {
      this->PlotMatricesPoly();
    }
    else {
      this->PlotMatricesNormal();
    }
  }



  void CovarianceCalculator4D::PlotMatricesNormal()
  {

    int n_bins = _bs.GetNbinsX() * _bs.GetNbinsY();

    TH2D * cov_matrix_histo = new TH2D("cov_matrix_histo", "",           n_bins, 0, n_bins, n_bins, 0, n_bins);
    TH2D * frac_cov_matrix_histo = new TH2D("frac_cov_matrix_histo", "", n_bins, 0, n_bins, n_bins, 0, n_bins);
    TH2D * corr_matrix_histo = new TH2D("corr_matrix_histo", "",         n_bins, 0, n_bins, n_bins, 0, n_bins);


    for (int i = 0; i < _bs.GetNbinsX(); i++) {
      for (int j = 0; j < _bs.GetNbinsY(); j++) {
        for (int m = 0; m < _bs.GetNbinsX(); m++) {
          for (int n = 0; n < _bs.GetNbinsY(); n++) { 

            int a = j + i * _bs.GetNbinsY() + 1;
            int b = n + m * _bs.GetNbinsY() + 1;

            // std::cout << "a: " << a << ", b: " << b << "(" << i << ", " << j << ", " << m << ", " << n << ")" << std::endl;

            cov_matrix_histo->SetBinContent(a, b, _M[i][j][m][n]);
            frac_cov_matrix_histo->SetBinContent(a, b, _M_frac[i][j][m][n]);
            corr_matrix_histo->SetBinContent(a, b, _RHO[i][j][m][n]);
            
          }
        }
      } 
    }

    this->PlotMatricesBase(cov_matrix_histo, frac_cov_matrix_histo, corr_matrix_histo);

  }


  void CovarianceCalculator4D::PlotMatricesPoly()
  {

    int n_bins = _bs_poly.GetNumberOfBins();

    TH2D * cov_matrix_histo = new TH2D("cov_matrix_histo", "",           n_bins, 0, n_bins, n_bins, 0, n_bins);
    TH2D * frac_cov_matrix_histo = new TH2D("frac_cov_matrix_histo", "", n_bins, 0, n_bins, n_bins, 0, n_bins);
    TH2D * corr_matrix_histo = new TH2D("corr_matrix_histo", "",         n_bins, 0, n_bins, n_bins, 0, n_bins);


    for (int i = 0; i < n_bins; i++) {
      for (int j = 0; j < n_bins; j++) {

        cov_matrix_histo->SetBinContent(i+1, j+1, _M_p[i][j]);
        frac_cov_matrix_histo->SetBinContent(i+1, j+1, _M_frac_p[i][j]);
        corr_matrix_histo->SetBinContent(i+1, j+1, _RHO_p[i][j]);

      } 
    }

    this->PlotMatricesBase(cov_matrix_histo, frac_cov_matrix_histo, corr_matrix_histo);

  }


  void CovarianceCalculator4D::PlotMatricesBase(TH2D * cov_matrix_histo, TH2D * frac_cov_matrix_histo, TH2D * corr_matrix_histo)
  {

    int n_bins = cov_matrix_histo->GetNbinsX();


    TString name;


    const Int_t NCont = 100;
    const Int_t NRGBs = 5;
    Double_t mainColour[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 1.00 };
    Double_t otherColour[NRGBs]   = { 0.99,0.80, 0.60, 0.40, 0.20 };
      //Double_t otherOtherColour[NRGBs]   = { 0.9,0.80, 0.80, 0.80, 0.80 };
    Double_t stops[NRGBs] = { 0.00, 0.05, 0.1, 0.4, 1.00 };

    TColor::CreateGradientColorTable(NRGBs, stops, mainColour, otherColour, otherColour, NCont);
    gStyle->SetNumberContours(NCont);


    //
    // Create axis lables
    //

    TH2F *h = new TH2F("h", "", cov_matrix_histo->GetNbinsX(), 0, cov_matrix_histo->GetNbinsX(),
      cov_matrix_histo->GetNbinsY(), 0, cov_matrix_histo->GetNbinsY());

    // h->SetMaximum(1);

    if (!_polybin_mode) {
      int i_label_number = 0;
      int j_label_number = 0;
      for (int i = 0; i <  cov_matrix_histo->GetNbinsX()+1; i++) {
        std::ostringstream oss;
        oss << i_label_number << "," << j_label_number;
        if (j_label_number % _bs.GetNbinsY() == 0) {
          i_label_number ++;
          j_label_number = 0;
        }
        j_label_number++;
        std::string label = oss.str();
        if (i == 0) continue;
        h->GetXaxis()->SetBinLabel(i,label.c_str());
        h->GetYaxis()->SetBinLabel(i,label.c_str());
      }
    } else {
      for (int i = 0; i < cov_matrix_histo->GetNbinsX(); i++) {
        std::ostringstream oss;
        oss << i + 1;
        std::string label = oss.str();
        h->GetXaxis()->SetBinLabel(i+1,label.c_str());
        h->GetYaxis()->SetBinLabel(i+1,label.c_str());
      }
    }

    h->GetXaxis()->SetLabelOffset(0.004);
    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelOffset(0.004);
    h->GetYaxis()->SetLabelSize(0.03);
    h->GetXaxis()->SetTitle("Bin Number");
    h->GetYaxis()->SetTitle("Bin Number");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();


    //
    // Create lines to divide primary bins
    //

    std::vector<TLine*> lines;

    if (_polybin_mode) {

      Int_t * separators = _bs_poly.GetNominal()->GetSeparators();
      Int_t separators_length = _bs_poly.GetNominal()->GetSeparatorsLength();

      if (separators_length > 50) separators_length = 0;

      // Vertical lines
      Int_t sum = 0;
      for (Int_t s = 0; s < separators_length - 1; s++) {
        TLine *line = new TLine(separators[s] + sum, 0, separators[s] + sum, n_bins);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        lines.emplace_back(line);
        sum += separators[s];
      }

      // Horizontal lines
      sum = 0;
      for (Int_t s = 0; s < separators_length - 1; s++) {
        TLine *line = new TLine(0, separators[s] + sum, n_bins, separators[s] + sum);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        lines.emplace_back(line);
        sum += separators[s];
      }
    }
    else {

      for (int i = 1; i < _bs.GetNbinsX(); i++) {
        TLine *line = new TLine(_bs.GetNbinsY()  * i, 0, _bs.GetNbinsY() * i, n_bins);
        line->SetLineColor(kGreen+2);
        line->SetLineWidth(2);
        lines.emplace_back(line);
      }        
      for (int i = 1; i < _bs.GetNbinsX(); i++) {
        TLine *line = new TLine(0, _bs.GetNbinsY() * i, n_bins, _bs.GetNbinsY() * i);
        line->SetLineColor(kGreen+2);
        line->SetLineWidth(2);
        lines.emplace_back(line);
      }
    }



    //
    // Create TLatex lables
    //

    TLatex* prelim = new TLatex(0.10,0.97, "i, m = cos(#theta_{#mu}) bins");
    prelim->SetTextColor(kBlack);
    prelim->SetTextFont(42);
    prelim->SetNDC();
    prelim->SetTextSize(1/30.);
    prelim->SetTextAlign(12);

    TLatex* prelim2 = new TLatex(0.10,0.93, "j, n = p_{#mu} bins");
    prelim2->SetTextColor(kBlack);
    prelim2->SetTextFont(42);
    prelim2->SetNDC();
    prelim2->SetTextSize(1/30.);
    prelim2->SetTextAlign(12);



    // 
    // Draw the proper matrices
    //

    TCanvas * cov_c = new TCanvas();
    cov_c->SetRightMargin(0.13);
    cov_c->SetFixedAspectRatio();
    cov_matrix_histo->SetMarkerColor(kBlack);
    cov_matrix_histo->SetMarkerSize(1.8);
    cov_matrix_histo->GetXaxis()->CenterTitle();
    cov_matrix_histo->GetYaxis()->CenterTitle();
    cov_matrix_histo->GetXaxis()->SetTitle("Bin Number");
    cov_matrix_histo->GetYaxis()->SetTitle("Bin Number");
    cov_matrix_histo->GetXaxis()->SetTickLength(0);
    cov_matrix_histo->GetYaxis()->SetTickLength(0);
    h->Draw();
      // cov_matrix_histo->Draw("colz text same");
    cov_matrix_histo->Draw("colz same");

    for (auto l : lines)
      l->Draw();

    if (!_polybin_mode) {
      prelim->Draw();
      prelim2->Draw();
    }

    PlottingTools::DrawSimulationXSec();
    name = _prefix + "_cov_matrix_2d";
    cov_c->SaveAs(name + ".pdf");
    cov_c->SaveAs(name + ".C","C");


    TCanvas * cov_frac_c = new TCanvas();
    cov_frac_c->SetRightMargin(0.13);
    cov_frac_c->SetFixedAspectRatio();
    frac_cov_matrix_histo->SetMarkerColor(kBlack);
    frac_cov_matrix_histo->SetMarkerSize(1.8);
    frac_cov_matrix_histo->GetXaxis()->CenterTitle();
    frac_cov_matrix_histo->GetYaxis()->CenterTitle();
    frac_cov_matrix_histo->GetXaxis()->SetTitle("Bin Number");
    frac_cov_matrix_histo->GetYaxis()->SetTitle("Bin Number");
    frac_cov_matrix_histo->GetXaxis()->SetTickLength(0);
    frac_cov_matrix_histo->GetYaxis()->SetTickLength(0);
    h->Draw();
      // frac_cov_matrix_histo->Draw("colz text same");
    frac_cov_matrix_histo->Draw("colz same");

    for (auto l : lines)
      l->Draw();

    if (!_polybin_mode) {
      prelim->Draw();
      prelim2->Draw();
    }

    PlottingTools::DrawSimulationXSec();
    name = _prefix + "_cov_frac_matrix_2d";
    cov_frac_c->SaveAs(name + ".pdf");
    cov_frac_c->SaveAs(name + ".C","C");


    TCanvas * corr_c = new TCanvas();
    corr_c->SetRightMargin(0.13);
    corr_c->SetFixedAspectRatio();
    corr_matrix_histo->SetMarkerColor(kBlack);
    corr_matrix_histo->SetMarkerSize(1.8);
    corr_matrix_histo->GetXaxis()->CenterTitle();
    corr_matrix_histo->GetYaxis()->CenterTitle();
    corr_matrix_histo->GetXaxis()->SetTitle("Bin Number");
    corr_matrix_histo->GetYaxis()->SetTitle("Bin Number");
    corr_matrix_histo->GetXaxis()->SetTickLength(0);
    corr_matrix_histo->GetYaxis()->SetTickLength(0);
    h->Draw();
      // corr_matrix_histo->Draw("colz text same");
    corr_matrix_histo->Draw("colz same");

    for (auto l : lines)
      l->Draw();

    if (!_polybin_mode) {
      prelim->Draw();
      prelim2->Draw();
    }

    PlottingTools::DrawSimulationXSec();
    name = _prefix + "_corr_matrix_2d";
    corr_c->SaveAs(name + ".pdf");
    corr_c->SaveAs(name + ".C","C");

    gStyle->SetPalette(kRainBow);




    _M_h = *cov_matrix_histo;
    _M_frac_h = *frac_cov_matrix_histo;
    _RHO_h = *corr_matrix_histo;

  }
}

#endif
