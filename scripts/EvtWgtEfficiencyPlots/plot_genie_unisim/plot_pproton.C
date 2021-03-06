void plot_pproton(){
   TCanvas *c = new TCanvas("c", "canvas",0,0,800,800);
   gStyle->SetOptStat(0);
   c->SetHighLightColor(2);
   c->Range(0,0,1,1);
   c->SetFillColor(10);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameLineWidth(2);
   c->SetFrameBorderMode(0);
  
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("");
   hs->SetMinimum(-0.0);
   hs->SetMaximum(0.04);
   Double_t xAxis2332[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1F *hs_stack_116 = new TH1F("hs_stack_116","",10, xAxis2332);
   hs_stack_116->SetMinimum(-0.02804466);
   hs_stack_116->SetMaximum(0.025164);
   hs_stack_116->SetDirectory(0);
   hs_stack_116->SetStats(0);
   hs_stack_116->SetLineWidth(2);
   hs_stack_116->GetXaxis()->SetTitle("p_{p} [GeV]");
   hs_stack_116->GetXaxis()->CenterTitle(true);
   hs_stack_116->GetXaxis()->SetNdivisions(506);
   hs_stack_116->GetXaxis()->SetLabelFont(43);
   hs_stack_116->GetXaxis()->SetLabelSize(20);
   hs_stack_116->GetXaxis()->SetTitleSize(25);
   hs_stack_116->GetXaxis()->SetTitleOffset(1.2);
   hs_stack_116->GetXaxis()->SetTitleFont(43);
   hs_stack_116->GetYaxis()->SetTitle("Relative Uncertainty");
   hs_stack_116->GetYaxis()->CenterTitle(true);
   hs_stack_116->GetYaxis()->SetNdivisions(505);
   hs_stack_116->GetYaxis()->SetLabelFont(43);
   hs_stack_116->GetYaxis()->SetLabelSize(15);
   hs_stack_116->GetYaxis()->SetTitleSize(25);
   hs_stack_116->GetYaxis()->SetTitleOffset(1.2);
   hs_stack_116->GetYaxis()->SetTitleFont(43);
   hs_stack_116->GetZaxis()->SetNdivisions(506);
   hs_stack_116->GetZaxis()->SetLabelFont(42);
   hs_stack_116->GetZaxis()->SetTitleSize(0.055);
   hs_stack_116->GetZaxis()->SetTitleOffset(0.8);
   hs_stack_116->GetZaxis()->SetTitleFont(42);
   hs->SetHistogram(hs_stack_116);
   Double_t xAxis2166[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_m1__2078 = new TH1D("ratio_m1__2078","bs_genie_pm1_eff_pmom_num_title",10, xAxis2166);
   ratio_m1__2078->SetBinContent(1,0.01516669);
   ratio_m1__2078->SetBinContent(2,0.007177945);
   ratio_m1__2078->SetBinContent(3,0.006472644);
   ratio_m1__2078->SetBinContent(4,0.002958813);
   ratio_m1__2078->SetBinContent(5,0.003242375);
   ratio_m1__2078->SetBinContent(6,0.0005566654);
   ratio_m1__2078->SetBinContent(7,-0.001046542);
   ratio_m1__2078->SetBinContent(8,0.002298257);
   ratio_m1__2078->SetBinContent(9,0.0004287917);
   ratio_m1__2078->SetBinContent(10,0.001606918);
   ratio_m1__2078->SetBinContent(11,0.008681631);
   ratio_m1__2078->SetBinError(1,0.01617107);
   ratio_m1__2078->SetBinError(2,0.01450416);
   ratio_m1__2078->SetBinError(3,0.01583794);
   ratio_m1__2078->SetBinError(4,0.01654145);
   ratio_m1__2078->SetBinError(5,0.01764004);
   ratio_m1__2078->SetBinError(6,0.01940555);
   ratio_m1__2078->SetBinError(7,0.02187321);
   ratio_m1__2078->SetBinError(8,0.02439185);
   ratio_m1__2078->SetBinError(9,0.03227455);
   ratio_m1__2078->SetBinError(10,0.02527879);
   ratio_m1__2078->SetBinError(11,0.1893304);
   ratio_m1__2078->SetEntries(0.3402996);
   ratio_m1__2078->SetStats(0);

   ci = TColor::GetColor("#cc0000");
   ratio_m1__2078->SetLineColor(ci);
   ratio_m1__2078->SetLineWidth(2);
   ratio_m1__2078->GetXaxis()->SetNdivisions(506);
   ratio_m1__2078->GetXaxis()->SetLabelFont(42);
   ratio_m1__2078->GetXaxis()->SetTitleSize(0.055);
   ratio_m1__2078->GetXaxis()->SetTitleOffset(0.8);
   ratio_m1__2078->GetXaxis()->SetTitleFont(42);
   ratio_m1__2078->GetYaxis()->SetNdivisions(506);
   ratio_m1__2078->GetYaxis()->SetLabelFont(42);
   ratio_m1__2078->GetYaxis()->SetTitleSize(0.055);
   ratio_m1__2078->GetYaxis()->SetTitleOffset(0.9);
   ratio_m1__2078->GetYaxis()->SetTitleFont(42);
   ratio_m1__2078->GetZaxis()->SetNdivisions(506);
   ratio_m1__2078->GetZaxis()->SetLabelFont(42);
   ratio_m1__2078->GetZaxis()->SetTitleSize(0.055);
   ratio_m1__2078->GetZaxis()->SetTitleOffset(0.8);
   ratio_m1__2078->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_m1__2078,"");

   Double_t xAxis2177[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_p1__2087 = new TH1D("ratio_p1__2087","",10, xAxis2177);
   ratio_p1__2087->SetBinContent(1,0.02699475);
   ratio_p1__2087->SetBinContent(2,0.0163978);
   ratio_p1__2087->SetBinContent(3,0.01555183);
   ratio_p1__2087->SetBinContent(4,0.006107463);
   ratio_p1__2087->SetBinContent(5,0.008083202);
   ratio_p1__2087->SetBinContent(6,0.002726967);
   ratio_p1__2087->SetBinContent(7,-0.00220253);
   ratio_p1__2087->SetBinContent(8,0.003577695);
   ratio_p1__2087->SetBinContent(9,0.001642419);
   ratio_p1__2087->SetBinContent(10,-0.002014959);
   ratio_p1__2087->SetBinContent(11,-0.007455829);
   ratio_p1__2087->SetBinError(1,0.01647435);
   ratio_p1__2087->SetBinError(2,0.01476904);
   ratio_p1__2087->SetBinError(3,0.01613113);
   ratio_p1__2087->SetBinError(4,0.01679587);
   ratio_p1__2087->SetBinError(5,0.01792351);
   ratio_p1__2087->SetBinError(6,0.01968552);
   ratio_p1__2087->SetBinError(7,0.02213311);
   ratio_p1__2087->SetBinError(8,0.02470413);
   ratio_p1__2087->SetBinError(9,0.0326334);
   ratio_p1__2087->SetBinError(10,0.02547213);
   ratio_p1__2087->SetBinError(11,0.1886323);
   ratio_p1__2087->SetEntries(1.297202);
   ratio_p1__2087->SetStats(0);

   ci = TColor::GetColor("#009900");
   ratio_p1__2087->SetLineColor(ci);
   ratio_p1__2087->SetLineWidth(2);
   ratio_p1__2087->GetXaxis()->SetNdivisions(506);
   ratio_p1__2087->GetXaxis()->SetLabelFont(42);
   ratio_p1__2087->GetXaxis()->SetTitleSize(0.055);
   ratio_p1__2087->GetXaxis()->SetTitleOffset(0.8);
   ratio_p1__2087->GetXaxis()->SetTitleFont(42);
   ratio_p1__2087->GetYaxis()->SetNdivisions(506);
   ratio_p1__2087->GetYaxis()->SetLabelFont(42);
   ratio_p1__2087->GetYaxis()->SetTitleSize(0.055);
   ratio_p1__2087->GetYaxis()->SetTitleOffset(0.9);
   ratio_p1__2087->GetYaxis()->SetTitleFont(42);
   ratio_p1__2087->GetZaxis()->SetNdivisions(506);
   ratio_p1__2087->GetZaxis()->SetLabelFont(42);
   ratio_p1__2087->GetZaxis()->SetTitleSize(0.055);
   ratio_p1__2087->GetZaxis()->SetTitleOffset(0.8);
   ratio_p1__2087->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_p1__2087,"");
   Double_t xAxis2201[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_p1__2107 = new TH1D("ratio_p1__2107","",10, xAxis2201);
   ratio_p1__2107->SetBinContent(1,-0.002237403);
   ratio_p1__2107->SetBinContent(2,-0.007541357);
   ratio_p1__2107->SetBinContent(3,-0.009694624);
   ratio_p1__2107->SetBinContent(4,-0.01009345);
   ratio_p1__2107->SetBinContent(5,-0.01056139);
   ratio_p1__2107->SetBinContent(6,-0.009020988);
   ratio_p1__2107->SetBinContent(7,-0.007884449);
   ratio_p1__2107->SetBinContent(8,-0.01012272);
   ratio_p1__2107->SetBinContent(9,-0.01097028);
   ratio_p1__2107->SetBinContent(10,-0.01010987);
   ratio_p1__2107->SetBinContent(11,-0.02263437);
   ratio_p1__2107->SetBinError(1,0.01597383);
   ratio_p1__2107->SetBinError(2,0.01437178);
   ratio_p1__2107->SetBinError(3,0.01570128);
   ratio_p1__2107->SetBinError(4,0.01642967);
   ratio_p1__2107->SetBinError(5,0.01751793);
   ratio_p1__2107->SetBinError(6,0.01931174);
   ratio_p1__2107->SetBinError(7,0.02179849);
   ratio_p1__2107->SetBinError(8,0.02424513);
   ratio_p1__2107->SetBinError(9,0.03210198);
   ratio_p1__2107->SetBinError(10,0.02515964);
   ratio_p1__2107->SetBinError(11,0.1868222);
   ratio_p1__2107->SetEntries(0.08823653);
   ratio_p1__2107->SetStats(0);

   ci = TColor::GetColor("#990099");
   ratio_p1__2107->SetLineColor(ci);
   ratio_p1__2107->SetLineWidth(2);
   ratio_p1__2107->GetXaxis()->SetNdivisions(506);
   ratio_p1__2107->GetXaxis()->SetLabelFont(42);
   ratio_p1__2107->GetXaxis()->SetTitleSize(0.055);
   ratio_p1__2107->GetXaxis()->SetTitleOffset(0.8);
   ratio_p1__2107->GetXaxis()->SetTitleFont(42);
   ratio_p1__2107->GetYaxis()->SetNdivisions(506);
   ratio_p1__2107->GetYaxis()->SetLabelFont(42);
   ratio_p1__2107->GetYaxis()->SetTitleSize(0.055);
   ratio_p1__2107->GetYaxis()->SetTitleOffset(0.9);
   ratio_p1__2107->GetYaxis()->SetTitleFont(42);
   ratio_p1__2107->GetZaxis()->SetNdivisions(506);
   ratio_p1__2107->GetZaxis()->SetLabelFont(42);
   ratio_p1__2107->GetZaxis()->SetTitleSize(0.055);
   ratio_p1__2107->GetZaxis()->SetTitleOffset(0.8);
   ratio_p1__2107->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_p1__2107,"");

   Double_t xAxis2219[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_p1__2122 = new TH1D("ratio_p1__2122","",10, xAxis2219);
   ratio_p1__2122->SetBinContent(1,0.002079068);
   ratio_p1__2122->SetBinContent(2,0.006830282);
   ratio_p1__2122->SetBinContent(3,0.009220248);
   ratio_p1__2122->SetBinContent(4,0.008309796);
   ratio_p1__2122->SetBinContent(5,0.009258665);
   ratio_p1__2122->SetBinContent(6,0.008148527);
   ratio_p1__2122->SetBinContent(7,0.007545939);
   ratio_p1__2122->SetBinContent(8,0.008703272);
   ratio_p1__2122->SetBinContent(9,0.01025636);
   ratio_p1__2122->SetBinContent(10,0.009773401);
   ratio_p1__2122->SetBinContent(11,0.01874762);
   ratio_p1__2122->SetBinError(1,0.01600495);
   ratio_p1__2122->SetBinError(2,0.01447278);
   ratio_p1__2122->SetBinError(3,0.01584718);
   ratio_p1__2122->SetBinError(4,0.01657644);
   ratio_p1__2122->SetBinError(5,0.01768247);
   ratio_p1__2122->SetBinError(6,0.01946914);
   ratio_p1__2122->SetBinError(7,0.0219614);
   ratio_p1__2122->SetBinError(8,0.02447047);
   ratio_p1__2122->SetBinError(9,0.03243975);
   ratio_p1__2122->SetBinError(10,0.02541262);
   ratio_p1__2122->SetBinError(11,0.1909207);
   ratio_p1__2122->SetEntries(1.438681);
   ratio_p1__2122->SetStats(0);

   ci = TColor::GetColor("#000099");
   ratio_p1__2122->SetLineColor(ci);
   ratio_p1__2122->SetLineWidth(2);
   ratio_p1__2122->GetXaxis()->SetNdivisions(506);
   ratio_p1__2122->GetXaxis()->SetLabelFont(42);
   ratio_p1__2122->GetXaxis()->SetTitleSize(0.055);
   ratio_p1__2122->GetXaxis()->SetTitleOffset(0.8);
   ratio_p1__2122->GetXaxis()->SetTitleFont(42);
   ratio_p1__2122->GetYaxis()->SetNdivisions(506);
   ratio_p1__2122->GetYaxis()->SetLabelFont(42);
   ratio_p1__2122->GetYaxis()->SetTitleSize(0.055);
   ratio_p1__2122->GetYaxis()->SetTitleOffset(0.9);
   ratio_p1__2122->GetYaxis()->SetTitleFont(42);
   ratio_p1__2122->GetZaxis()->SetNdivisions(506);
   ratio_p1__2122->GetZaxis()->SetLabelFont(42);
   ratio_p1__2122->GetZaxis()->SetTitleSize(0.055);
   ratio_p1__2122->GetZaxis()->SetTitleOffset(0.8);
   ratio_p1__2122->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_p1__2122,"");
   Double_t xAxis2334[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_m1__2218 = new TH1D("ratio_m1__2218","bs_genie_pm1_eff_pmom_num_title",10, xAxis2334);
   ratio_m1__2218->SetBinContent(1,-0.007625935);
   ratio_m1__2218->SetBinContent(2,-0.008313792);
   ratio_m1__2218->SetBinContent(3,-0.009755521);
   ratio_m1__2218->SetBinContent(4,-0.01203256);
   ratio_m1__2218->SetBinContent(5,-0.01222449);
   ratio_m1__2218->SetBinContent(6,-0.009601826);
   ratio_m1__2218->SetBinContent(7,-0.009612976);
   ratio_m1__2218->SetBinContent(8,-0.01121044);
   ratio_m1__2218->SetBinContent(9,-0.01033002);
   ratio_m1__2218->SetBinContent(10,0.009653015);
   ratio_m1__2218->SetBinContent(11,-0.02493564);
   ratio_m1__2218->SetBinError(1,0.0159353);
   ratio_m1__2218->SetBinError(2,0.014361);
   ratio_m1__2218->SetBinError(3,0.01568894);
   ratio_m1__2218->SetBinError(4,0.01640442);
   ratio_m1__2218->SetBinError(5,0.01749635);
   ratio_m1__2218->SetBinError(6,0.01930815);
   ratio_m1__2218->SetBinError(7,0.02179082);
   ratio_m1__2218->SetBinError(8,0.02425426);
   ratio_m1__2218->SetBinError(9,0.03215221);
   ratio_m1__2218->SetBinError(10,0.02544394);
   ratio_m1__2218->SetBinError(11,0.1867886);
   ratio_m1__2218->SetEntries(0.08105455);
   ratio_m1__2218->SetStats(0);

   ci = TColor::GetColor("#00ccc7");
   ratio_m1__2218->SetLineColor(ci);
   ratio_m1__2218->SetLineWidth(2);
   ratio_m1__2218->GetXaxis()->SetNdivisions(506);
   ratio_m1__2218->GetXaxis()->SetLabelFont(42);
   ratio_m1__2218->GetXaxis()->SetTitleSize(0.055);
   ratio_m1__2218->GetXaxis()->SetTitleOffset(0.8);
   ratio_m1__2218->GetXaxis()->SetTitleFont(42);
   ratio_m1__2218->GetYaxis()->SetNdivisions(506);
   ratio_m1__2218->GetYaxis()->SetLabelFont(42);
   ratio_m1__2218->GetYaxis()->SetTitleSize(0.055);
   ratio_m1__2218->GetYaxis()->SetTitleOffset(0.9);
   ratio_m1__2218->GetYaxis()->SetTitleFont(42);
   ratio_m1__2218->GetZaxis()->SetNdivisions(506);
   ratio_m1__2218->GetZaxis()->SetLabelFont(42);
   ratio_m1__2218->GetZaxis()->SetTitleSize(0.055);
   ratio_m1__2218->GetZaxis()->SetTitleOffset(0.8);
   ratio_m1__2218->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_m1__2218,"");
   

   TLegend *leg=new TLegend (0.6, 0.6, 0.9, 0.9);
   leg->AddEntry(ratio_m1__2078, "IntraNukeNabs");
   leg->AddEntry(ratio_p1__2087, "IntraNukeNel");
   leg->AddEntry(ratio_p1__2107, "IntraNukePIabs");
   leg->AddEntry(ratio_p1__2122, "IntraNukePIinel");
   leg->AddEntry(ratio_m1__2218, "QE Ma");
   hs->Draw("nostack histo");
   leg->Draw("same");
   TLine *line = new TLine(-1,1,1,1);
   line->SetLineStyle(9);
   line->Draw();
   c->cd();
   c->Modified();
   //c->cd();
   //c->SetSelected(c);
   c->Print("GENIE_unisim_pproton.png");

}
