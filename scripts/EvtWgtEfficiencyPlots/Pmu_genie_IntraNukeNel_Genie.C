void Pmu_genie_IntraNukeNel_Genie()
{
//=========Macro generated from canvas: c/canvas
//=========  (Wed May 29 09:31:57 2019) by ROOT version6.06/08
   TCanvas *c = new TCanvas("c", "canvas",51,67,800,1200);
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
   pad1->Range(-0.3125,0.0001,2.8125,1.1111);
   pad1->SetFillColor(10);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetGridx();
   pad1->SetBottomMargin(0);
   pad1->SetFrameLineWidth(2);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameLineWidth(2);
   pad1->SetFrameBorderMode(0);
   Double_t xAxis67[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1D *p1__56 = new TH1D("p1__56","bs_genie_pm1_eff_mumom_num_title",6, xAxis67);
   p1__56->SetBinContent(1,0.07295402);
   p1__56->SetBinContent(2,0.2443323);
   p1__56->SetBinContent(3,0.2948866);
   p1__56->SetBinContent(4,0.3156348);
   p1__56->SetBinContent(5,0.3004829);
   p1__56->SetBinContent(6,0.2745997);
   p1__56->SetBinContent(7,0.1813466);
   p1__56->SetBinError(1,0.002364327);
   p1__56->SetBinError(2,0.003099929);
   p1__56->SetBinError(3,0.002676683);
   p1__56->SetBinError(4,0.002529491);
   p1__56->SetBinError(5,0.002632647);
   p1__56->SetBinError(6,0.003509123);
   p1__56->SetBinError(7,0.009175752);
   p1__56->SetMinimum(0.0001);
   p1__56->SetMaximum(1);
   p1__56->SetEntries(47048.63);
   p1__56->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc0000");
   p1__56->SetLineColor(ci);
   p1__56->SetLineWidth(2);
   p1__56->GetXaxis()->SetNdivisions(506);
   p1__56->GetXaxis()->SetLabelFont(42);
   p1__56->GetXaxis()->SetTitleSize(0.055);
   p1__56->GetXaxis()->SetTitleOffset(0.8);
   p1__56->GetXaxis()->SetTitleFont(42);
   p1__56->GetYaxis()->SetTitle("Efficiency");
   p1__56->GetYaxis()->CenterTitle(true);
   p1__56->GetYaxis()->SetNdivisions(506);
   p1__56->GetYaxis()->SetLabelFont(42);
   p1__56->GetYaxis()->SetTitleSize(25);
   p1__56->GetYaxis()->SetTitleOffset(1.55);
   p1__56->GetYaxis()->SetTitleFont(43);
   p1__56->GetZaxis()->SetNdivisions(506);
   p1__56->GetZaxis()->SetLabelFont(42);
   p1__56->GetZaxis()->SetTitleSize(0.055);
   p1__56->GetZaxis()->SetTitleOffset(0.8);
   p1__56->GetZaxis()->SetTitleFont(42);
   p1__56->Draw("histo");
   Double_t xAxis68[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1D *nominal_eff__57 = new TH1D("nominal_eff__57","bs_genie_pm1_eff_mumom_num_title",6, xAxis68);
   nominal_eff__57->SetBinContent(1,0.0733842);
   nominal_eff__57->SetBinContent(2,0.2413326);
   nominal_eff__57->SetBinContent(3,0.2907289);
   nominal_eff__57->SetBinContent(4,0.3118239);
   nominal_eff__57->SetBinContent(5,0.2958878);
   nominal_eff__57->SetBinContent(6,0.2707195);
   nominal_eff__57->SetBinContent(7,0.1787536);
   nominal_eff__57->SetBinError(1,0.002292434);
   nominal_eff__57->SetBinError(2,0.002959182);
   nominal_eff__57->SetBinError(3,0.002559558);
   nominal_eff__57->SetBinError(4,0.002423418);
   nominal_eff__57->SetBinError(5,0.002513036);
   nominal_eff__57->SetBinError(6,0.003351046);
   nominal_eff__57->SetBinError(7,0.008765173);
   nominal_eff__57->SetEntries(50064.41);
   nominal_eff__57->SetStats(0);
   nominal_eff__57->SetLineWidth(2);
   nominal_eff__57->GetXaxis()->SetNdivisions(506);
   nominal_eff__57->GetXaxis()->SetLabelFont(42);
   nominal_eff__57->GetXaxis()->SetTitleSize(0.055);
   nominal_eff__57->GetXaxis()->SetTitleOffset(0.8);
   nominal_eff__57->GetXaxis()->SetTitleFont(42);
   nominal_eff__57->GetYaxis()->SetNdivisions(506);
   nominal_eff__57->GetYaxis()->SetLabelFont(42);
   nominal_eff__57->GetYaxis()->SetLabelSize(0);
   nominal_eff__57->GetYaxis()->SetTitleSize(0.055);
   nominal_eff__57->GetYaxis()->SetTitleOffset(0.9);
   nominal_eff__57->GetYaxis()->SetTitleFont(42);
   nominal_eff__57->GetZaxis()->SetNdivisions(506);
   nominal_eff__57->GetZaxis()->SetLabelFont(42);
   nominal_eff__57->GetZaxis()->SetTitleSize(0.055);
   nominal_eff__57->GetZaxis()->SetTitleOffset(0.8);
   nominal_eff__57->GetZaxis()->SetTitleFont(42);
   nominal_eff__57->Draw("histo same");
   Double_t xAxis69[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1D *m1__58 = new TH1D("m1__58","bs_genie_pm1_eff_mumom_num_title",6, xAxis69);
   m1__58->SetBinContent(1,0.07337579);
   m1__58->SetBinContent(2,0.238099);
   m1__58->SetBinContent(3,0.2863986);
   m1__58->SetBinContent(4,0.3075716);
   m1__58->SetBinContent(5,0.2906537);
   m1__58->SetBinContent(6,0.2663608);
   m1__58->SetBinContent(7,0.1766526);
   m1__58->SetBinError(1,0.002389897);
   m1__58->SetBinError(2,0.003034695);
   m1__58->SetBinError(3,0.002611396);
   m1__58->SetBinError(4,0.002476538);
   m1__58->SetBinError(5,0.00256064);
   m1__58->SetBinError(6,0.003424622);
   m1__58->SetBinError(7,0.008961065);
   m1__58->SetEntries(46335.7);
   m1__58->SetStats(0);

   ci = TColor::GetColor("#009900");
   m1__58->SetLineColor(ci);
   m1__58->SetLineWidth(2);
   m1__58->GetXaxis()->SetNdivisions(506);
   m1__58->GetXaxis()->SetLabelFont(42);
   m1__58->GetXaxis()->SetTitleSize(0.055);
   m1__58->GetXaxis()->SetTitleOffset(0.8);
   m1__58->GetXaxis()->SetTitleFont(42);
   m1__58->GetYaxis()->SetNdivisions(506);
   m1__58->GetYaxis()->SetLabelFont(42);
   m1__58->GetYaxis()->SetTitleSize(0.055);
   m1__58->GetYaxis()->SetTitleOffset(0.9);
   m1__58->GetYaxis()->SetTitleFont(42);
   m1__58->GetZaxis()->SetNdivisions(506);
   m1__58->GetZaxis()->SetLabelFont(42);
   m1__58->GetZaxis()->SetTitleSize(0.055);
   m1__58->GetZaxis()->SetTitleOffset(0.8);
   m1__58->GetZaxis()->SetTitleFont(42);
   m1__58->Draw("histo same");
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
   
   TPaveText *pt = new TPaveText(0.15,0.9163265,0.85,0.9836735,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("bs_genie_pm1_eff_mumom_num_title");
   pt->Draw();
   pad1->Modified();
   c->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "pad2",0,0.05,1,0.5);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-0.3125,-0.08090928,2.8125,0.02552971);
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
   hs->SetMinimum(-0.02768979);
   hs->SetMaximum(0.02552971);
   Double_t xAxis70[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1F *hs_stack_12 = new TH1F("hs_stack_12","",6, xAxis70);
   hs_stack_12->SetMinimum(-0.02768979);
   hs_stack_12->SetMaximum(0.02552971);
   hs_stack_12->SetDirectory(0);
   hs_stack_12->SetStats(0);
   hs_stack_12->SetLineWidth(2);
   hs_stack_12->GetXaxis()->SetTitle("p_{#mu} [GeV]");
   hs_stack_12->GetXaxis()->CenterTitle(true);
   hs_stack_12->GetXaxis()->SetNdivisions(506);
   hs_stack_12->GetXaxis()->SetLabelFont(43);
   hs_stack_12->GetXaxis()->SetLabelSize(20);
   hs_stack_12->GetXaxis()->SetTitleSize(25);
   hs_stack_12->GetXaxis()->SetTitleOffset(3.5);
   hs_stack_12->GetXaxis()->SetTitleFont(43);
   hs_stack_12->GetYaxis()->SetTitle("Ratio");
   hs_stack_12->GetYaxis()->CenterTitle(true);
   hs_stack_12->GetYaxis()->SetNdivisions(505);
   hs_stack_12->GetYaxis()->SetLabelFont(43);
   hs_stack_12->GetYaxis()->SetLabelSize(15);
   hs_stack_12->GetYaxis()->SetTitleSize(25);
   hs_stack_12->GetYaxis()->SetTitleOffset(1.2);
   hs_stack_12->GetYaxis()->SetTitleFont(43);
   hs_stack_12->GetZaxis()->SetNdivisions(506);
   hs_stack_12->GetZaxis()->SetLabelFont(42);
   hs_stack_12->GetZaxis()->SetTitleSize(0.055);
   hs_stack_12->GetZaxis()->SetTitleOffset(0.8);
   hs_stack_12->GetZaxis()->SetTitleFont(42);
   hs->SetHistogram(hs_stack_12);
   
   Double_t xAxis71[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1D *ratio_p1__59 = new TH1D("ratio_p1__59","",6, xAxis71);
   ratio_p1__59->SetBinContent(1,-0.005862051);
   ratio_p1__59->SetBinContent(2,0.01242966);
   ratio_p1__59->SetBinContent(3,0.01430091);
   ratio_p1__59->SetBinContent(4,0.0122211);
   ratio_p1__59->SetBinContent(5,0.01552971);
   ratio_p1__59->SetBinContent(6,0.01433297);
   ratio_p1__59->SetBinContent(7,0.01450587);
   ratio_p1__59->SetBinError(1,0.04487678);
   ratio_p1__59->SetBinError(2,0.0177587);
   ratio_p1__59->SetBinError(3,0.01273932);
   ratio_p1__59->SetBinError(4,0.01123443);
   ratio_p1__59->SetBinError(5,0.01230108);
   ratio_p1__59->SetBinError(6,0.01792408);
   ratio_p1__59->SetBinError(7,0.07099219);
   ratio_p1__59->SetEntries(1.28236);
   ratio_p1__59->SetStats(0);

   ci = TColor::GetColor("#cc0000");
   ratio_p1__59->SetLineColor(ci);
   ratio_p1__59->SetLineWidth(2);
   ratio_p1__59->GetXaxis()->SetNdivisions(506);
   ratio_p1__59->GetXaxis()->SetLabelFont(42);
   ratio_p1__59->GetXaxis()->SetTitleSize(0.055);
   ratio_p1__59->GetXaxis()->SetTitleOffset(0.8);
   ratio_p1__59->GetXaxis()->SetTitleFont(42);
   ratio_p1__59->GetYaxis()->SetNdivisions(506);
   ratio_p1__59->GetYaxis()->SetLabelFont(42);
   ratio_p1__59->GetYaxis()->SetTitleSize(0.055);
   ratio_p1__59->GetYaxis()->SetTitleOffset(0.9);
   ratio_p1__59->GetYaxis()->SetTitleFont(42);
   ratio_p1__59->GetZaxis()->SetNdivisions(506);
   ratio_p1__59->GetZaxis()->SetLabelFont(42);
   ratio_p1__59->GetZaxis()->SetTitleSize(0.055);
   ratio_p1__59->GetZaxis()->SetTitleOffset(0.8);
   ratio_p1__59->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_p1,"");
   Double_t xAxis72[7] = {0, 0.18, 0.3, 0.48, 0.75, 1.14, 2.5}; 
   
   TH1D *ratio_m1__60 = new TH1D("ratio_m1__60","bs_genie_pm1_eff_mumom_num_title",6, xAxis72);
   ratio_m1__60->SetBinContent(1,-0.0001147081);
   ratio_m1__60->SetBinContent(2,-0.01339899);
   ratio_m1__60->SetBinContent(3,-0.01489469);
   ratio_m1__60->SetBinContent(4,-0.01363707);
   ratio_m1__60->SetBinContent(5,-0.01768979);
   ratio_m1__60->SetBinContent(6,-0.01610043);
   ratio_m1__60->SetBinContent(7,-0.01175354);
   ratio_m1__60->SetBinError(1,0.04512722);
   ratio_m1__60->SetBinError(2,0.01756428);
   ratio_m1__60->SetBinError(3,0.01257804);
   ratio_m1__60->SetBinError(4,0.01111253);
   ratio_m1__60->SetBinError(5,0.01212643);
   ratio_m1__60->SetBinError(6,0.0176999);
   ratio_m1__60->SetBinError(7,0.07012744);
   ratio_m1__60->SetEntries(0.07583567);
   ratio_m1__60->SetStats(0);

   ci = TColor::GetColor("#009900");
   ratio_m1__60->SetLineColor(ci);
   ratio_m1__60->SetLineWidth(2);
   ratio_m1__60->GetXaxis()->SetNdivisions(506);
   ratio_m1__60->GetXaxis()->SetLabelFont(42);
   ratio_m1__60->GetXaxis()->SetTitleSize(0.055);
   ratio_m1__60->GetXaxis()->SetTitleOffset(0.8);
   ratio_m1__60->GetXaxis()->SetTitleFont(42);
   ratio_m1__60->GetYaxis()->SetNdivisions(506);
   ratio_m1__60->GetYaxis()->SetLabelFont(42);
   ratio_m1__60->GetYaxis()->SetTitleSize(0.055);
   ratio_m1__60->GetYaxis()->SetTitleOffset(0.9);
   ratio_m1__60->GetYaxis()->SetTitleFont(42);
   ratio_m1__60->GetZaxis()->SetNdivisions(506);
   ratio_m1__60->GetZaxis()->SetLabelFont(42);
   ratio_m1__60->GetZaxis()->SetTitleSize(0.055);
   ratio_m1__60->GetZaxis()->SetTitleOffset(0.8);
   ratio_m1__60->GetZaxis()->SetTitleFont(42);
   hs->Add(ratio_m1,"");
   hs->Draw("nostack histo");
   TLine *line = new TLine(0,0,2.5,0);
   line->SetLineStyle(9);
   line->Draw();
   pad2->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
