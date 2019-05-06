void makeClass() {
  
  TFile *f = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root");
  std::cout<<"opened the input file "<<std::endl;
  TTree *v = (TTree*)f->Get("true_reco_tree");
  std::cout<<"get a tree content called true_reco_tree"<<std::endl;
  v->MakeClass("RewgtFac");
}
