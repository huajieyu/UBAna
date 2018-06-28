#ifndef __BASE_MIGRATIONMATRIX4D_CXX__
#define __BASE_MIGRATIONMATRIX4D_CXX__

#include "MigrationMatrix4D.h"

namespace Base {

	void MigrationMatrix4D::SetScaleFactors(double bnbcosmic, double bnbon, double extbnb, double intimecosmic)
  {
  
    _configured = true;
  }

  void MigrationMatrix4D::SetOutputFileName(std::string name) 
  {
    _f_out.open(name, std::ios::out | std::ios::trunc);
  }

  void MigrationMatrix4D::SetTTree(TTree *t)
  {
    _tree = t;
  }

  void MigrationMatrix4D::UseWeights(std::string weight_name)
  {
    _use_weights = true;
    _weight_name = weight_name;
  }

  void MigrationMatrix4D::SetBins(double *var1_b, int n_var1_bins, double *var2_b, int n_var2_bins)
  {

    for (int i = 0; i < n_var1_bins; i++) 
    {
      _var1_bins.push_back(std::make_pair(var1_b[i], var1_b[i+1]));
    }

    for (int i = 0; i < n_var2_bins; i++) 
    {
      _var2_bins.push_back(std::make_pair(var2_b[i], var2_b[i+1]));
    }

    std::cout << _prefix << "Number of var1 bins: " << _var1_bins.size() << std::endl;
    std::cout << _prefix << "Number of var2 bins: " << _var2_bins.size() << std::endl;

    _reco_per_true = new TH2D("reco_per_true", "reco_per_true", n_var1_bins, var1_b, n_var2_bins, var2_b);
    
  }


  void MigrationMatrix4D::CalculateMigrationMatrix() 
  {
    Double_t        mom_true;
    Double_t        mom_mcs;

    Double_t        angle_true;
    Double_t        angle_reco;

    Double_t        event_weight;

    std::vector<std::string> wgtsnames_genie_multisim;
    std::vector<double> wgts_genie_multisim;

    std::vector<std::string> wgtsnames_flux_multisim;
    std::vector<double> wgts_flux_multisim;

    TBranch        *b_mom_true; 
    TBranch        *b_mom_mcs;

    TBranch        *b_angle_true;
    TBranch        *b_angle_reco;

    TBranch        *b_event_weight;

    TBranch        *b_wgtsnames_genie_multisim;
    TBranch        *b_wgts_genie_multisim;

    TBranch        *b_wgtsnames_flux_multisim;
    TBranch        *b_wgts_flux_multisim;

    _tree->SetMakeClass(1);

    _tree->SetBranchAddress("mom_true", &mom_true, &b_mom_true);
    _tree->SetBranchAddress("mom_mcs", &mom_mcs, &b_mom_mcs);

    _tree->SetBranchAddress("angle_true", &angle_true, &b_angle_true);
    _tree->SetBranchAddress("angle_reco", &angle_reco, &b_angle_reco);

    _tree->SetBranchAddress("event_weight", &event_weight, &b_event_weight);

    _tree->SetBranchAddress("wgtsnames_genie_multisim", &wgtsnames_genie_multisim, &b_wgtsnames_genie_multisim);
    _tree->SetBranchAddress("wgts_genie_multisim", &wgts_genie_multisim, &b_wgts_genie_multisim);

    _tree->SetBranchAddress("wgtsnames_flux_multisim", &wgtsnames_flux_multisim, &b_wgtsnames_flux_multisim);
    _tree->SetBranchAddress("wgts_flux_multisim", &wgts_flux_multisim, &b_wgts_flux_multisim);

    Long64_t nentries = _tree->GetEntriesFast();

    

    // Resize the smearing matrix
    _S.resize( _var1_bins.size(), 
              std::vector<std::vector<std::vector<double>>> (_var2_bins.size(),
                                                             std::vector<std::vector<double>>(_var1_bins.size(),
                                                                                             std::vector<double> (_var2_bins.size(), 0.
                                                                                                                 )
                                                                                            )
                                                            )
              );

    int counter = 0;

    for (int i = 0; i < _var1_bins.size(); i++) {
      for (int j = 0; j < _var2_bins.size(); j++) {
        for (int m = 0; m < _var1_bins.size(); m++) {
          for (int n = 0; n < _var2_bins.size(); n++) { 
            if(_verbose) std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _S[i][j][m][n] << std::endl;
            counter++;
          }
        }
      }
    }
    std::cout << _prefix << "Total entries: " << counter << std::endl;


    // True bin m, n
    //int m = 0, n = 0;

    for (int m = 0; m < _var1_bins.size(); m++) {
      for (int n = 0; n < _var2_bins.size(); n++) {

        std::cout << _prefix << "m = " << m << ", n = " << n << std::endl;

        auto v1_bin = _var1_bins.at(m);
        auto v2_bin = _var2_bins.at(n);

        std::cout << _prefix << "Here 1" << std::endl;

        if(_verbose) std::cout << _prefix << "b1: " << v1_bin.first << " - " << v1_bin.second << std::endl;
        if(_verbose) std::cout << _prefix << "b2: " << v2_bin.first << " - " << v2_bin.second << std::endl;

        _reco_per_true->Reset();

        std::cout << _prefix << "Here 2" << std::endl;

        for (Long64_t jentry=0; jentry < nentries;jentry++) {
          _tree->GetEntry(jentry);


          // Weights starts
          // double evt_weight = event_weight;

          // if (_use_weights) {

          //   std::cout << _prefix << "Using weight with name: " << _weight_name << std::endl;

          //   for (size_t i = 0; i < wgtsnames_genie_multisim.size(); i++) {
          //     if (wgtsnames_genie_multisim.at(i) == _weight_name){
          //       evt_weight *= wgts_genie_multisim.at(i);
          //       std::cout << _prefix << "Weight with name: " << _weight_name << " found." << std::endl;
          //       break;
          //     }
          //   }

          //   for (size_t i = 0; i < wgtsnames_flux_multisim.size(); i++) {
          //     if (wgtsnames_flux_multisim.at(i) == _weight_name){
          //       evt_weight *= wgts_flux_multisim.at(i);
          //       std::cout << _prefix << "Weight with name: " << _weight_name << " found." << std::endl;
          //       break;
          //     }
          //   }

          // }
          // Weights ends

          // std::cout << _prefix << "mom_true: " << mom_true << ", mom_mcs: " << mom_mcs << ", evt_weight: " << evt_weight << std::endl;

         
          if (  angle_true > v1_bin.first && angle_true < v1_bin.second
             && mom_true > v2_bin.first   && mom_true < v2_bin.second) {

            // Filling reco bin i, j
            _reco_per_true->Fill(angle_reco, mom_mcs/*, evt_weight*/);
          }
        }

        std::cout << _prefix << "Here 3" << std::endl;

        // Normalize to get a probability
        _reco_per_true->Scale(1./_reco_per_true->Integral());

        // Set values to matrix
        TCanvas *c = new TCanvas();
        _reco_per_true->Draw("colz text");
        for (int i = 0; i < _var1_bins.size(); i++) {
          for (int j = 0; j < _var2_bins.size(); j++) {
            if(_verbose) std::cout << "(" << i << ", " << j << ")" << _reco_per_true->GetBinContent(i+1, j+1) << std::endl;

            double value = _reco_per_true->GetBinContent(i+1, j+1);
            if (std::isnan(value))
              value = 0.;
  
            _S[i][j][m][n] = value;
          }
        }

        std::cout << _prefix << "Here 4" << std::endl;

        // Saving the plot
        std::stringstream sstm;
        sstm << "True Bin (" << m << ", " << n << ")";
        std::string str = sstm.str();
        _reco_per_true->SetTitle(str.c_str());
        _reco_per_true->GetXaxis()->SetTitle("cos(#theta_{#mu}) [Reco Bin i]");
        _reco_per_true->GetYaxis()->SetTitle("p_{#mu} (GeV) [Reco Bin j]");
        _reco_per_true->GetYaxis()->SetTitleOffset(0.8);

        sstm.str("");
        sstm << "smearing_matrix_true_" << m << "_" << n;

        std::cout << _prefix << "Here 5" << std::endl;


        TString name = _folder + sstm.str();
        c->SaveAs(name + ".pdf");
        c->SaveAs(name + ".C","C");

        std::cout << _prefix << "Here 6" << std::endl;
      }
    }

    std::cout << _prefix << "Here 7" << std::endl;


    // if(_verbose) {
      for (int i = 0; i < _var1_bins.size(); i++) {
        for (int j = 0; j < _var2_bins.size(); j++) {
          for (int m = 0; m < _var1_bins.size(); m++) {
            for (int n = 0; n < _var2_bins.size(); n++) { 
              std::cout << "(" << i << ", " << j << ", " << m << ", " << n << ") => " << _S[i][j][m][n] << std::endl;
            }
          }
        }
      }
    // }

    std::cout << _prefix << "Here 8" << std::endl;

    std::cout << _prefix << "_S.size() " << _S.size() << std::endl;

    // return _S;

  }

  void MigrationMatrix4D::PlotMatrix()
  {

    int n_bins = _var1_bins.size() * _var2_bins.size();
    TH2D *h_sm = new TH2D("h_sm", "", n_bins, 0, n_bins, n_bins, 0, n_bins);


    for (int n = 0; n < _var2_bins.size(); n++) {   // pmu true
      for (int m = 0; m < _var1_bins.size(); m++) {  // theta true
        for (int j = 0; j < _var2_bins.size(); j++) {  // pmu reco
          for (int i = 0; i < _var1_bins.size(); i++) {  // theta reco
        
            int reco_bin = i + j * _var1_bins.size() + 1;
            int true_bin = m + n * _var1_bins.size() + 1;
            h_sm->SetBinContent(reco_bin, true_bin, _S[i][j][m][n]);
            if(_verbose) std::cout << "(i, j, m, n) = (" << i << ", " << j << ", " << m << ", " << n << "   reco_bin: " << reco_bin << ", true_bin: " << true_bin << ", S: " << _S[i][j][m][n] << std::endl;

          }
        }
      }
    }

    std::vector<std::string> bin_labels;

    for (int j = 0; j < _var2_bins.size(); j++) {  
      for (int i = 0; i < _var1_bins.size(); i++) {  

        //int bin = i + j * _var1_bins.size() + 1;

        std::stringstream sstm;
        sstm << i;
        std::string str = sstm.str();

        bin_labels.emplace_back(str);

      }
    }

    for (int i = 0; i < n_bins; i++) {
      h_sm->GetXaxis()->SetBinLabel(i+1, bin_labels.at(i).c_str());
      h_sm->GetYaxis()->SetBinLabel(i+1, bin_labels.at(i).c_str());
    }

    std::vector<TLine*> lines;

    for (int i = 1; i < _var2_bins.size(); i++) {
      TLine *line = new TLine(_var1_bins.size() * i, 0, _var1_bins.size() * i, n_bins);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      lines.emplace_back(line);
    }

    for (int i = 1; i < _var2_bins.size(); i++) {
      TLine *line = new TLine(0, _var1_bins.size() * i, n_bins, _var1_bins.size() * i);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      lines.emplace_back(line);
    }


    TCanvas *c_sm = new TCanvas();
    h_sm->Draw("colz");
 
    for (auto l : lines)
     l->Draw();


/*
    // Other axis
    gPad->Update();
    Double_t xmin = gPad->GetUxmin();
    Double_t ymin = 0;
    Double_t xmax = gPad->GetUxmax();
    Double_t ymax = 5;

    TGaxis *axis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin()-1,gPad->GetUxmax(),gPad->GetUymin()-1,ymin,ymax,506,"+L");
    axis->SetLineColor(kRed+1);
    axis->SetLabelColor(kRed+1);
    axis->SetTextFont(42);
    axis->SetLabelSize(.04);
    axis->SetLabelOffset(.005);


    axis->Draw();
    */


    TString name = _folder + "full_migration_matrix_4d";
    c_sm->SaveAs(name + ".pdf");
    c_sm->SaveAs(name + ".C","C");

  }


  void MigrationMatrix4D::PrintSmearingMatrixLatex()
  {

    for (int m = 0; m < _var1_bins.size(); m++) {
      for (int n = 0; n < _var2_bins.size(); n++) {
        this->PrintSmearingMatrixLatex(m, n);
      }
    }

  }



  void MigrationMatrix4D::PrintSmearingMatrixLatex(int true_m, int true_n)
  {


    if (!_f_out.is_open()) {
      std::cout << "File not opened." << std::endl;
      return;
    }

    _f_out << "\\begin{equation}" << std::endl;
    _f_out << "S_{ij" << true_m << true_n << "} =" << std::endl;
    _f_out << "\\begin{bmatrix}" << std::endl;

    for (int i = 0; i < _var1_bins.size(); i++) {
      for (int j = 0; j < _var2_bins.size(); j++) {

        _f_out << _S[i][j][true_m][true_n] << "  &  ";

      }

      _f_out << " \\\\" << std::endl;
    }

    _f_out << "\\end{bmatrix}" << std::endl;
    _f_out << "\\end{equation}" << std::endl << std::endl;

  }
}

#endif
