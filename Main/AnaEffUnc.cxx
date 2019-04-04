#ifndef __MAIN_ANAEFFUNC_CXX__
#define __MAIN_ANAEFFUNC_CXX__

#include "AnaEffUnc.h"

namespace Main{
  void AnaEffUnc::SetBNBCosmicFile(std::string f) {
   
    mc_bnbcosmic_file_name = f;

  }
  void AnaEffUnc::SetInputBNBCosmicFile(std::string f) {
   
    mc_bnbcosmic_input_file_name = f;

  }


  void AnaEffUnc::SetPrefix(std::string p) {
    _prefix = p;
  }
  void AnaEffUnc::SetOutputFile(std::string f)
  {
    mc_bnbcosmic_output_file_name = f;
  }


  void AnaEffUnc::GetRMS(std::map<std::string, TH1D*>* bnbcosmic){
    
       std::cout<<""<<bnbcosmic->size()<<std::endl;
  }




  void AnaEffUnc::DoAnaEffUnc(){


  std::cout<<"Start open the input files "<<std::endl;
  TFile* mc_bnbcosmic_file = TFile::Open(mc_bnbcosmic_file_name.c_str(), "READ");


  TFile* mc_bnbcosmic_input_file = TFile::Open(mc_bnbcosmic_input_file_name.c_str(), "READ");

 
  std::cout<<"<<<<<<<<<<<<<beging of DoAnaEffAna<<<<<<<<<<<<<<<<<"<<std::endl;


  TFile *file_out = new TFile(mc_bnbcosmic_output_file_name.c_str(),"RECREATE");
  if ( file_out->IsOpen() ) {
    LOG_NORMAL() << "File opened successfully." << std::endl;
  } else {
    LOG_CRITICAL() << "Failed to open file." << std::endl;
    throw std::exception();
  }
 


  std::cout<<"<<<<<<<<<<<<<Get the GENIE reweight factors<<<<<<<<<<"<<std::endl;
  TChain *chain_ubxsec;
  chain_ubxsec = new TChain("UBXSec/tree");
  chain_ubxsec ->Add(mc_bnbcosmic_input_file_name.c_str());
  UBXSecEvent * t = new UBXSecEvent(chain_ubxsec);
  int evts = chain_ubxsec -> GetEntries();
  std::cout<<"total number of entries is "<<evts<<std::endl;     
  chain_ubxsec -> GetEntry(0);
   
  std::vector<double> wgts_genie_multisim;



  std::vector<std::string> fname_flux_multisim;
  fname_flux_multisim.clear();
  fname_flux_multisim.resize(t->evtwgt_flux_multisim_nweight.at(1));

  std::vector<std::string> fname_extra_syst;
  fname_extra_syst.clear();
  fname_extra_syst.resize(100);
  std::cout<<t->evtwgt_extra_syst_multisim_nweight.at(1)<<std::endl;


 
  std::ostringstream oss;

  ofstream outfile;
  outfile.open("rewght_genie_rwfac");
      for (size_t i = 0; i < t->evtwgt_genie_multisim_weight.size(); i++) {
         wgts_genie_multisim.push_back(t->evtwgt_genie_multisim_weight.at(i).at(0));
         wgts_genie_multisim.push_back(t->evtwgt_genie_multisim_weight.at(i).at(1));
         //outfile<<std::setprecision(60)<<t->evtwgt_genie_multisim_weight.at(i).at(0)<<std::setprecision(60)<<t->evtwgt_genie_multisim_weight.at(i).at(1)<<std::endl;
       }
      std::cout<<t->evtwgt_flux_multisim_weight.size()<<std::endl;
      std::cout<<"size of the flux syst is "<<t->evtwgt_flux_multisim_funcname.size()<<std::endl;


      for(size_t i_wgt = 0; i_wgt < t->evtwgt_flux_multisim_funcname.size(); i_wgt++){
        oss.str("");
        oss<< "universe"<<i_wgt;  
        //std::string func_name = t->evtwgt_flux_multisim_funcname.at(i_wgt);
        //std::cout<<i_wgt<< " th function name is"<<func_name<<std::endl; 
        
        fname_flux_multisim.at(i_wgt) = oss.str();         
        std::cout<<fname_flux_multisim.at(i_wgt)<<std::endl; 
      }
      std::vector<double> wgts_flux_multisim;
      wgts_flux_multisim.clear();
      wgts_flux_multisim.resize(fname_flux_multisim.size(), 1.);

      for(size_t i_func = 0; i_func < t->evtwgt_flux_multisim_funcname.size(); i_func++){
         std::string func_name = t->evtwgt_flux_multisim_funcname.at(i_func);
         //size_t found = std::string::npos;
        
         if(func_name =="bnbcorrection_FluxHist") continue;
         for (size_t i_wgt = 0; i_wgt < fname_flux_multisim.size(); i_wgt++) {
           wgts_flux_multisim.at(i_wgt) *= t->evtwgt_flux_multisim_weight.at(i_func).at(i_wgt);
         }
      }

      std::cout<<"size of the extra syst is "<<t->evtwgt_extra_syst_multisim_funcname.size()<<std::endl;

      std::vector<double> wgts_extra_syst;

      wgts_extra_syst.clear();
      wgts_extra_syst.resize(fname_extra_syst.size(), 1.);

      for(size_t i_func=0; i_func<t->evtwgt_extra_syst_multisim_funcname.size(); i_func++){
        std::string func_name=t->evtwgt_extra_syst_multisim_funcname.at(i_func);
        std::cout<<i_func<<" th function name of the extra syst is "<<func_name<<std::endl;
        if(func_name=="bnbcorrection_FluxHist") continue;
        for(size_t i_wgt=0; i_wgt<fname_extra_syst.size(); i_wgt++){
            wgts_extra_syst.at(i_wgt) *=t->evtwgt_extra_syst_multisim_weight.at(i_func).at(i_wgt);
            //std::cout<<t->evtwgt_extra_syst_multisim_weight.at(i_func).at(i_wgt)<<std::endl;
        }
      }



   outfile.close();
   for(unsigned int i_flux=0; i_flux<wgts_flux_multisim.size(); i_flux++){
       std::cout<<wgts_flux_multisim[i_flux]<<std::endl;

   }




  for(unsigned int i_extra=0; i_extra<wgts_extra_syst.size(); i_extra++){
       std::cout<<wgts_extra_syst[i_extra]<<std::endl;

  }





  

  UBXSecEventHisto1D * _event_histo_1d_mc = 0;
  mc_bnbcosmic_file->GetObject("UBXSecEventHisto1D", _event_histo_1d_mc);

  /*for(auto it : _event_histo_1d_mc->bs_genie_multisim_eff_mumom_num){
     std::cout<<it.first()<<std::endl;
     std::map<std::string, TH1D*> bs_map_num=it.second();
  }*/




















  //========================================================================== 
  int n_genie=wgts_genie_multisim.size();
  TH1D* histo_num[n_genie+1];
  TH1D* histo_den[n_genie+1];
  
  
  histo_num[0]=_event_histo_1d_mc->h_eff_mumom_num;  
  histo_den[0]=_event_histo_1d_mc->h_eff_mumom_den;



  histo_num[0]->Write();
  histo_den[0]->Write();

  for(unsigned int i=0; i<wgts_genie_multisim.size(); i++){
  
 
  histo_num[i+1]=(TH1D*)histo_num[0]->Clone();
  histo_num[i+1]->Scale(wgts_genie_multisim[i]);
  
  histo_den[i+1]=(TH1D*)histo_den[0]->Clone();
  histo_den[i+1]->Scale(wgts_genie_multisim[i]);

  std::cout<<i<<" th weight factor is "<<wgts_genie_multisim[i]<<std::endl;

  histo_num[i+1]->Write();
  histo_den[i+1]->Write();

  }


  TCanvas *c1=new TCanvas("c","c");
  TEfficiency * pEff=new TEfficiency(*histo_num[1], *histo_den[1]);
  pEff->Draw("AP");
  TEfficiency * pEff2=new TEfficiency(*histo_num[2], *histo_den[2]);
  pEff2->Draw("same");
  c1->SaveAs("test.png");
  LOG_NORMAL() << "Checkpoint 1" << std::endl;

  TH1D* histo_muangle_num[n_genie+1];
  TH1D* histo_muangle_den[n_genie+1];
  
  histo_muangle_num[0]=_event_histo_1d_mc->h_eff_muangle_num;  
  histo_muangle_den[0]=_event_histo_1d_mc->h_eff_muangle_den;

  histo_muangle_num[0]->Write();
  histo_muangle_den[0]->Write();

  for(unsigned int i=0; i<wgts_genie_multisim.size(); i++){
  
 
  histo_muangle_num[i+1]=(TH1D*)histo_muangle_num[0]->Clone();
  histo_muangle_num[i+1]->Scale(wgts_genie_multisim[i]);
  
  histo_muangle_den[i+1]=(TH1D*)histo_muangle_den[0]->Clone();
  histo_muangle_den[i+1]->Scale(wgts_genie_multisim[i]);

  std::cout<<i<<" th weight factor is "<<wgts_genie_multisim[i]<<std::endl;

  histo_muangle_num[i+1]->Write();
  histo_muangle_den[i+1]->Write();

  }












  LOG_NORMAL() << "Checkpoint 2" << std::endl;
  TH1D* histo_pangle_num[n_genie+1];
  TH1D* histo_pangle_den[n_genie+1];
  
  histo_pangle_num[0]=_event_histo_1d_mc->h_eff_pangle_num;  
  histo_pangle_den[0]=_event_histo_1d_mc->h_eff_pangle_den;

  histo_pangle_num[0]->Write();
  histo_pangle_den[0]->Write();

  for(unsigned int i=0; i<wgts_genie_multisim.size(); i++){
  
 
  histo_pangle_num[i+1]=(TH1D*)histo_pangle_num[0]->Clone();
  histo_pangle_num[i+1]->Scale(wgts_genie_multisim[i]);
  
  histo_pangle_den[i+1]=(TH1D*)histo_pangle_den[0]->Clone();
  histo_pangle_den[i+1]->Scale(wgts_genie_multisim[i]);

  std::cout<<i<<" th weight factor is "<<wgts_genie_multisim[i]<<std::endl;

  histo_pangle_num[i+1]->Write();
  histo_pangle_den[i+1]->Write();

  }



  LOG_NORMAL() << "Checkpoint 3" << std::endl;
  TH1D* histo_pmom_num[n_genie+1];
  TH1D* histo_pmom_den[n_genie+1];
  
  histo_pmom_num[0]=_event_histo_1d_mc->h_eff_pmom_num;  
  histo_pmom_den[0]=_event_histo_1d_mc->h_eff_pmom_den;

  histo_pmom_num[0]->Write();
  histo_pmom_den[0]->Write();

  for(unsigned int i=0; i<wgts_genie_multisim.size(); i++){
  
 
  histo_pmom_num[i+1]=(TH1D*)histo_pmom_num[0]->Clone();
  histo_pmom_num[i+1]->Scale(wgts_genie_multisim[i]);
  
  histo_pmom_den[i+1]=(TH1D*)histo_pmom_den[0]->Clone();
  histo_pmom_den[i+1]->Scale(wgts_genie_multisim[i]);

  std::cout<<i<<" th weight factor is "<<wgts_genie_multisim[i]<<std::endl;

  histo_pmom_num[i+1]->Write();
  histo_pmom_den[i+1]->Write();

  }




  LOG_NORMAL() << "Checkpoint 4" << std::endl;
  //loop over all the histograms and plot the efficiency



  LOG_NORMAL() << "Checkpoint 5" << std::endl;
  
 
  file_out->cd();

  LOG_NORMAL() << "Checkpoint 6" << std::endl;







  /*LOG_NORMAL() << "Checkpoint 8" << std::endl;
  std::map<std::string,std::map<std::string,TH1D*>>* hmap_temp_genie_multisim_bs;
  mc_bnbcosmic_file->GetObject("hmap_trkmom_genie_multisim_bs", hmap_temp_genie_multisim_bs);


  TH1D *histo_num;
  TH1D *histo_den;
  histo_num = (TH1D*) bs_genie_multisim_eff_num.GetNominal().Clone();
  histo_den = (TH1D*) bs_genie_multisim_eff_den.GetNominal().Clone();
  histo_num->SetLineColor(kRed);
  histo_den->SetLineColor(kBlue);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetCanvasSize(1500, 1500);
  c1->SetWindowSize(500, 500);
  c1->cd();
  histo_den->Draw();
  histo_num->Draw("same");
  c1->SaveAs("test.png");
  */

  //std::map<std::string, std::vector<TH1D>> histo_map;
  //std::map<std::string, std::vector<TH1D>> histo_map_den;



  //histo_map=bs_genie_multisim_eff_num.UnpackPMHisto();
  //histo_map_den=bs_genie_multisim_eff_den.UnpackPMHisto();
  //std::cout<<"GetN Universes= "<<bs_genie_multisim_eff_num.GetNWeights()<<std::endl;



  /*TCanvas *temp_c = new TCanvas("temp_c", "temp_c");
  TH1D* histo_p1;
  TH1D* histo_m1; 
  temp_c->cd();

  int ind=0;
  for(auto iter: histo_map) {

     std::string function_name = iter.first;
     //std::cout<<function_name<<std::endl;i
     std::cout<<"ind= "<<ind<<std::endl;
     histo_p1 = (TH1D*) iter.second.at(1).Clone();
     histo_m1 = (TH1D*) iter.second.at(0).Clone();
     
     histo_p1->Draw();
     histo_m1->Draw("same");
     temp_c->SaveAs("test1.png");
     ind++;
    
  }
  */




  


  //std::cout<<histo_den->GetNbinsX()<<" total number of entries is "<<histo_den->Integral()<<std::endl;
    
 
  file_out->Close();





  std::cout<<"<<<<<<<<<<<<<end of DoAnaEffAna<<<<<<<<<<<<<<<<<"<<std::endl;


  }//End of DoAnaEffUnc

}
#endif
