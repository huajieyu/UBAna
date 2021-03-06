void Pproton_genie_IntraNukeNel_Genie()
{
//=========Macro generated from canvas: c/canvas
//=========  (Wed May 29 09:33:53 2019) by ROOT version6.06/08
   TCanvas *c = new TCanvas("c", "canvas",0,0,800,1200);
   gStyle->SetOptStat(0);
   c->SetHighLightColor(2);
   c->Range(0,0,1,1);
   c->SetFillColor(10);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameLineWidth(2);
   c->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "pad1",0,0.5,1,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(0.15,0.0001,1.65,1.1111);
   pad1->SetFillColor(10);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetGridx();
   pad1->SetBottomMargin(0);
   pad1->SetFrameLineWidth(2);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameLineWidth(2);
   pad1->SetFrameBorderMode(0);
   Double_t xAxis2173[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *p1__2084 = new TH1D("p1__2084","bs_genie_pm1_eff_pmom_num_title",10, xAxis2173);
   p1__2084->SetBinContent(1,0.1848015);
   p1__2084->SetBinContent(2,0.3155587);
   p1__2084->SetBinContent(3,0.3496736);
   p1__2084->SetBinContent(4,0.3541776);
   p1__2084->SetBinContent(5,0.3466555);
   p1__2084->SetBinContent(6,0.336982);
   p1__2084->SetBinContent(7,0.3181852);
   p1__2084->SetBinContent(8,0.2845372);
   p1__2084->SetBinContent(9,0.249451);
   p1__2084->SetBinContent(10,0.130088);
   p1__2084->SetBinContent(11,0.02917074);
   p1__2084->SetBinError(1,0.002157756);
   p1__2084->SetBinError(2,0.00332194);
   p1__2084->SetBinError(3,0.004020974);
   p1__2084->SetBinError(4,0.004260808);
   p1__2084->SetBinError(5,0.004444055);
   p1__2084->SetBinError(6,0.00475594);
   p1__2084->SetBinError(7,0.005057221);
   p1__2084->SetBinError(8,0.005029569);
   p1__2084->SetBinError(9,0.005820416);
   p1__2084->SetBinError(10,0.002370461);
   p1__2084->SetBinError(11,0.00392351);
   p1__2084->SetMinimum(0.0001);
   p1__2084->SetMaximum(1);
   p1__2084->SetEntries(45075.2);
   p1__2084->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc0000");
   p1__2084->SetLineColor(ci);
   p1__2084->SetLineWidth(2);
   p1__2084->GetXaxis()->SetNdivisions(506);
   p1__2084->GetXaxis()->SetLabelFont(42);
   p1__2084->GetXaxis()->SetTitleSize(0.055);
   p1__2084->GetXaxis()->SetTitleOffset(0.8);
   p1__2084->GetXaxis()->SetTitleFont(42);
   p1__2084->GetYaxis()->SetTitle("Efficiency");
   p1__2084->GetYaxis()->CenterTitle(true);
   p1__2084->GetYaxis()->SetNdivisions(506);
   p1__2084->GetYaxis()->SetLabelFont(42);
   p1__2084->GetYaxis()->SetTitleSize(25);
   p1__2084->GetYaxis()->SetTitleOffset(1.55);
   p1__2084->GetYaxis()->SetTitleFont(43);
   p1__2084->GetZaxis()->SetNdivisions(506);
   p1__2084->GetZaxis()->SetLabelFont(42);
   p1__2084->GetZaxis()->SetTitleSize(0.055);
   p1__2084->GetZaxis()->SetTitleOffset(0.8);
   p1__2084->GetZaxis()->SetTitleFont(42);
   p1__2084->Draw("histo");
   Double_t xAxis2174[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *nominal_eff__2085 = new TH1D("nominal_eff__2085","bs_genie_pm1_eff_pmom_num_title",10, xAxis2174);
   nominal_eff__2085->SetBinContent(1,0.1799439);
   nominal_eff__2085->SetBinContent(2,0.3104677);
   nominal_eff__2085->SetBinContent(3,0.3443188);
   nominal_eff__2085->SetBinContent(4,0.3520276);
   nominal_eff__2085->SetBinContent(5,0.3438759);
   nominal_eff__2085->SetBinContent(6,0.3360656);
   nominal_eff__2085->SetBinContent(7,0.3188875);
   nominal_eff__2085->SetBinContent(8,0.2835229);
   nominal_eff__2085->SetBinContent(9,0.2490419);
   nominal_eff__2085->SetBinContent(10,0.1303507);
   nominal_eff__2085->SetBinContent(11,0.02938986);
   nominal_eff__2085->SetBinError(1,0.002032018);
   nominal_eff__2085->SetBinError(2,0.003160235);
   nominal_eff__2085->SetBinError(3,0.003831178);
   nominal_eff__2085->SetBinError(4,0.004099247);
   nominal_eff__2085->SetBinError(5,0.004270532);
   nominal_eff__2085->SetBinError(6,0.004598628);
   nominal_eff__2085->SetBinError(7,0.004923352);
   nominal_eff__2085->SetBinError(8,0.004874604);
   nominal_eff__2085->SetBinError(9,0.005672053);
   nominal_eff__2085->SetBinError(10,0.002324941);
   nominal_eff__2085->SetBinError(11,0.003916605);
   nominal_eff__2085->SetEntries(47572.1);
   nominal_eff__2085->SetStats(0);
   nominal_eff__2085->SetLineWidth(2);
   nominal_eff__2085->GetXaxis()->SetNdivisions(506);
   nominal_eff__2085->GetXaxis()->SetLabelFont(42);
   nominal_eff__2085->GetXaxis()->SetTitleSize(0.055);
   nominal_eff__2085->GetXaxis()->SetTitleOffset(0.8);
   nominal_eff__2085->GetXaxis()->SetTitleFont(42);
   nominal_eff__2085->GetYaxis()->SetNdivisions(506);
   nominal_eff__2085->GetYaxis()->SetLabelFont(42);
   nominal_eff__2085->GetYaxis()->SetLabelSize(0);
   nominal_eff__2085->GetYaxis()->SetTitleSize(0.055);
   nominal_eff__2085->GetYaxis()->SetTitleOffset(0.9);
   nominal_eff__2085->GetYaxis()->SetTitleFont(42);
   nominal_eff__2085->GetZaxis()->SetNdivisions(506);
   nominal_eff__2085->GetZaxis()->SetLabelFont(42);
   nominal_eff__2085->GetZaxis()->SetTitleSize(0.055);
   nominal_eff__2085->GetZaxis()->SetTitleOffset(0.8);
   nominal_eff__2085->GetZaxis()->SetTitleFont(42);
   nominal_eff__2085->Draw("histo same");
   Double_t xAxis2175[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *m1__2086 = new TH1D("m1__2086","bs_genie_pm1_eff_pmom_num_title",10, xAxis2175);
   m1__2086->SetBinContent(1,0.1751664);
   m1__2086->SetBinContent(2,0.3050947);
   m1__2086->SetBinContent(3,0.3385413);
   m1__2086->SetBinContent(4,0.3492137);
   m1__2086->SetBinContent(5,0.3412064);
   m1__2086->SetBinContent(6,0.3352381);
   m1__2086->SetBinContent(7,0.3198701);
   m1__2086->SetBinContent(8,0.2829704);
   m1__2086->SetBinContent(9,0.2487809);
   m1__2086->SetBinContent(10,0.1305518);
   m1__2086->SetBinContent(11,0.02961796);
   m1__2086->SetBinError(1,0.002049621);
   m1__2086->SetBinError(2,0.003225124);
   m1__2086->SetBinError(3,0.003914458);
   m1__2086->SetBinError(4,0.004222417);
   m1__2086->SetBinError(5,0.004400654);
   m1__2086->SetBinError(6,0.004758596);
   m1__2086->SetBinError(7,0.005111244);
   m1__2086->SetBinError(8,0.005027495);
   m1__2086->SetBinError(9,0.005834775);
   m1__2086->SetBinError(10,0.002389818);
   m1__2086->SetBinError(11,0.003990369);
   m1__2086->SetEntries(44162.26);
   m1__2086->SetStats(0);

   ci = TColor::GetColor("#009900");
   m1__2086->SetLineColor(ci);
   m1__2086->SetLineWidth(2);
   m1__2086->GetXaxis()->SetNdivisions(506);
   m1__2086->GetXaxis()->SetLabelFont(42);
   m1__2086->GetXaxis()->SetTitleSize(0.055);
   m1__2086->GetXaxis()->SetTitleOffset(0.8);
   m1__2086->GetXaxis()->SetTitleFont(42);
   m1__2086->GetYaxis()->SetNdivisions(506);
   m1__2086->GetYaxis()->SetLabelFont(42);
   m1__2086->GetYaxis()->SetTitleSize(0.055);
   m1__2086->GetYaxis()->SetTitleOffset(0.9);
   m1__2086->GetYaxis()->SetTitleFont(42);
   m1__2086->GetZaxis()->SetNdivisions(506);
   m1__2086->GetZaxis()->SetLabelFont(42);
   m1__2086->GetZaxis()->SetTitleSize(0.055);
   m1__2086->GetZaxis()->SetTitleOffset(0.8);
   m1__2086->GetZaxis()->SetTitleFont(42);
   m1__2086->Draw("histo same");
   TGaxis *gaxis = new TGaxis(-5,20,-5,220,20,220,510,"");
   gaxis->SetLabelOffset(0.005);
   gaxis->SetLabelSize(15);
   gaxis->SetTickSize(0.03);
   gaxis->SetGridLength(0);
   gaxis->SetTitleOffset(1);
   gaxis->SetTitleSize(0.04);
   gaxis->SetTitleColor(1);
   gaxis->SetTitleFont(62);
   gaxis->SetLabelFont(43);
   gaxis->Draw();
   
   TLegend *leg = new TLegend(0.65,0.6,0.85,0.87,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("nominal_eff","Nominal","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("p1","x_{el}^{N} + 1#sigma","lpf");
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#cc0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("m1","x_{el}^{N} - 1#sigma","lpf");
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#009900");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.155201,0.9162116,0.844799,0.9837884,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("bs_genie_pm1_eff_pmom_num_title");
   pt->Draw();
   pad1->Modified();
   c->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "pad2",0,0.05,1,0.5);
   pad2->Draw();
   pad2->cd();
   pad2->Range(0.15,-0.1100948,1.65,0.03699475);
   pad2->SetFillColor(10);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetGridx();
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.5);
   pad2->SetFrameLineWidth(2);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameLineWidth(2);
   pad2->SetFrameBorderMode(0);
   
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("");
   hs->SetMinimum(-0.03655002);
   hs->SetMaximum(0.03699475);
   Double_t xAxis2176[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1F *hs_stack_90 = new TH1F("hs_stack_90","",10, xAxis2176);
   hs_stack_90->SetMinimum(-0.03655002);
   hs_stack_90->SetMaximum(0.03699475);
   hs_stack_90->SetDirectory(0);
   hs_stack_90->SetStats(0);
   hs_stack_90->SetLineWidth(2);
   hs_stack_90->GetXaxis()->SetTitle("p_{p} [GeV]");
   hs_stack_90->GetXaxis()->CenterTitle(true);
   hs_stack_90->GetXaxis()->SetNdivisions(506);
   hs_stack_90->GetXaxis()->SetLabelFont(43);
   hs_stack_90->GetXaxis()->SetLabelSize(20);
   hs_stack_90->GetXaxis()->SetTitleSize(25);
   hs_stack_90->GetXaxis()->SetTitleOffset(3.5);
   hs_stack_90->GetXaxis()->SetTitleFont(43);
   hs_stack_90->GetYaxis()->SetTitle("Ratio");
   hs_stack_90->GetYaxis()->CenterTitle(true);
   hs_stack_90->GetYaxis()->SetNdivisions(505);
   hs_stack_90->GetYaxis()->SetLabelFont(43);
   hs_stack_90->GetYaxis()->SetLabelSize(15);
   hs_stack_90->GetYaxis()->SetTitleSize(25);
   hs_stack_90->GetYaxis()->SetTitleOffset(1.2);
   hs_stack_90->GetYaxis()->SetTitleFont(43);
   hs_stack_90->GetZaxis()->SetNdivisions(506);
   hs_stack_90->GetZaxis()->SetLabelFont(42);
   hs_stack_90->GetZaxis()->SetTitleSize(0.055);
   hs_stack_90->GetZaxis()->SetTitleOffset(0.8);
   hs_stack_90->GetZaxis()->SetTitleFont(42);
   hs->SetHistogram(hs_stack_90);
   
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

   ci = TColor::GetColor("#cc0000");
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
   hs->Add(ratio_p1,"");
   Double_t xAxis2178[11] = {0.3, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.5}; 
   
   TH1D *ratio_m1__2088 = new TH1D("ratio_m1__2088","bs_genie_pm1_eff_pmom_num_title",10, xAxis2178);
   ratio_m1__2088->SetBinContent(1,-0.02655002);
   ratio_m1__2088->SetBinContent(2,-0.01730611);
   ratio_m1__2088->SetBinContent(3,-0.01677942);
   ratio_m1__2088->SetBinContent(4,-0.007993182);
   ratio_m1__2088->SetBinContent(5,-0.0077628);
   ratio_m1__2088->SetBinContent(6,-0.002462378);
   ratio_m1__2088->SetBinContent(7,0.003081127);
   ratio_m1__2088->SetBinContent(8,-0.0019485);
   ratio_m1__2088->SetBinContent(9,-0.001048317);
   ratio_m1__2088->SetBinContent(10,0.001543036);
   ratio_m1__2088->SetBinContent(11,0.007761161);
   ratio_m1__2088->SetBinError(1,0.01604214);
   ratio_m1__2088->SetBinError(2,0.01454481);
   ratio_m1__2088->SetBinError(3,0.01590876);
   ratio_m1__2088->SetBinError(4,0.01671756);
   ratio_m1__2088->SetBinError(5,0.0178327);
   ratio_m1__2088->SetBinError(6,0.01969119);
   ratio_m1__2088->SetBinError(7,0.02225484);
   ratio_m1__2088->SetBinError(8,0.02469883);
   ratio_m1__2088->SetBinError(9,0.03267471);
   ratio_m1__2088->SetBinError(10,0.02557836);
   ratio_m1__2088->SetBinError(11,0.1902493);
   ratio_m1__2088->SetEntries(0.07722656);
   ratio_m1__2088->SetStats(0);

   ci = TColor::GetColor("#009900");
   ratio_m1__2088->SetLineColor(ci);
   ratio_m1__2088->SetLineWidth(2);
   ratio_m1__2088->GetXaxis()->SetNdivisions(506);
   ratio_m1__2088->GetXaxis()->SetLabelFont(42);
   ratio_m1__2088->GetXaxis()->SetTitleSize(0.055);
   ratio_m1__2088->GetXaxis()->SetTitleOffset(0.8);
   ratio_m1__2088->GetXaxis()->SetTitleFont(42);
   ratio_m1__2088->GetYaxis()->SetNdivisions(506);
   ratio_m1__2088->GetYaxis()->SetLabelFont(42);
   ratio_m1__2088->GetYaxis()->SetTitleSize(0.055);
   ratio_m1__2088->GetYaxis()->SetTitleOffset(0.9);
   ratio_m1__2088->GetYaxis()->SetTitleFont(42);
   ratio_m1__2088->GetZaxis()->SetNdivisions(506);
   ratio_m1__2088->GetZaxis()->SetLabelFont(42);
   ratio_m1__2088->GetZaxis()->SetTitleSize(0.055);
   ratio_m1__2088->GetZaxis()->SetTitleOffset(0.8);
   ratio_m1__2088->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_m1,"");
   hs->Draw("nostack histo");
   TLine *line = new TLine(0,0,1.5,0);
   line->SetLineStyle(9);
   line->Draw();
   pad2->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
