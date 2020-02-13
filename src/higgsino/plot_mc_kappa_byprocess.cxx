///// plot_ratios: plots rMJ and rmT, and combinations of these

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace Functions;

namespace{
  bool debug = false;
  // options "zll", "qcd", "ttbar", "search"
  string sample_name = "ttbar";
  TString tag = "";
  bool boosted = false;
  float lumi=137;

  struct aplot{
    TString name;
    vector<string> abcd;
    NamedFunc baseline;
    vector<TString> bincuts;
  };
}


void plotRatio(vector<vector<vector<GammaParams> > > &allyields, aplot &plotdef,
               vector<vector<vector<int> > > &indices, vector<TString> &leglabels,
               vector<shared_ptr<Process> > &procs, torch::OrderedDict<string, string> abcdcuts);
// void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
//                 TString baseline, vector<shared_ptr<Process> > &procs);
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  set<int> years;
  years = {2016, 2017, 2018};

  string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  // string mc_skim_folder = "mc/merged_higmc_higtight/";
  string mc_skim_folder = "mc/merged_higmc_higtight/";
  if (sample_name=="ttbar") mc_skim_folder = "mc/merged_higmc_higlep1/";
  // if (sample_name=="ttbar") foldermc = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higlep1/";
  // if (sample_name=="zll") foldermc = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higlep2/";

  set<string> alltags = {"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                         "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", 
                         "*_ST_*.root", "*_WJetsToLNu*.root", "*_ZJet*.root", 
                         "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                         "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                         "*_QCD_HT2000toInf_*", "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                         "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"};

  set<string> ttxtags = {"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"};
  
  set<string> qcdtags = {"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                         "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                         "*_QCD_HT2000toInf_*"};

  NamedFunc base_filters = HigUtilities::pass_2016 && "stitch";
  vector<shared_ptr<Process> > procs;
  // All bkg. kappa
  procs.push_back(Process::MakeShared<Baby_pico>("All bkg.", Process::Type::background, 
    kBlack, attach_folder(mc_base_folder, years, mc_skim_folder,alltags), base_filters));
  // Sample specific kappa
  if (sample_name=="zll" || sample_name=="search") 
    procs.push_back(Process::MakeShared<Baby_pico>("Z#rightarrow#nu#nu", Process::Type::background, 
      kOrange+1, attach_folder(mc_base_folder, years, mc_skim_folder,{"*_ZJet*.root"}), base_filters));
  if (sample_name=="qcd") 
    procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, 
      colors("other"),attach_folder(mc_base_folder, years, mc_skim_folder,qcdtags), base_filters)); 
  if (sample_name=="ttbar" || sample_name=="search") 
    procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,
      colors("tt_1l"),attach_folder(mc_base_folder, years, mc_skim_folder,ttxtags), base_filters));

 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////// Defining cuts ///////////////////////////////////////////////
  // Baseline definitions
  // NamedFunc wgt = "w_lumi*w_isr";
  NamedFunc wgt = "w_lumi*w_isr";
  string baseline = "";
  if (boosted) {
    baseline += "ht>600 && nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && nbm>=2";
    baseline += "&& fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  } else {
    string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
    baseline += "njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  }

  torch::OrderedDict<string, NamedFunc> baselines; // all desired cut combinations
  // zll skim:  ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100))
  // nlep==2 && Max$(lep_pt)>40 
  if (boosted) {
    baselines.insert("vanilla", baseline);
  } else {
    if (sample_name=="zll") {
      baselines.insert("zll_CS",baseline+"&& nlep==2 && met<50");
    } else if (sample_name=="qcd") {
      baselines.insert("qcd_CS",baseline+"&& nvlep==0 && ntk==0 && lowDphiFix");
    } else if (sample_name=="ttbar") {
      baseline += "&& nlep==1 && mt<100";
      baselines.insert("ttX_CS_loose",baseline);
      baselines.insert("ttX_CS",baseline+"&& ntk==0 && !lowDphiFix");
      baselines.insert("ttX_CS_low_drbb", baseline+"&& hig_cand_drmax[0]<=1.1");
      baselines.insert("ttX_CS_high_drbb", baseline+"&& hig_cand_drmax[0]>1.1");
    } else if (sample_name=="search") {
        baseline += "&& nvlep==0 && ntk==0 && !lowDphiFix";
        baselines.insert("SR_vanilla", baseline);
        baselines.insert("SR_low_drbb", baseline+"&& hig_cand_drmax[0]<=1.1");
        baselines.insert("SR_high_drbb", baseline+"&& hig_cand_drmax[0]>1.1");
    }
  }

  vector<TString> metcuts;
  string metdef = "met";
  if (sample_name=="zll") metdef = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  if (sample_name=="zll" || sample_name=="ttbar") {
    metcuts.push_back(metdef+"<=50");
    metcuts.push_back(metdef+">50&&"+metdef+"<=100");
  }
  // if (sample_name!="qcd") metcuts.push_back(metdef+">100&&"+metdef+"<=150");
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  if (boosted) {
    metcuts.push_back(metdef+">300&&"+metdef+"<=500");
    metcuts.push_back(metdef+">500&&"+metdef+"<=700");
    metcuts.push_back(metdef+">700");
  } else {
    metcuts.push_back(metdef+">300&&"+metdef+"<=400");
    metcuts.push_back(metdef+">400");
  }

  string c_ab="nbm==2";
  string c_cd="nbm==3&&nbl==3";
  string c_ef="nbm>=3&&nbl>=4";
  if (sample_name=="qcd" || sample_name=="zll"){
    c_ab="nbm==0";
    c_cd="nbm==1";
    c_ef="nbt==1";  
  }

  string c_hig="hig_cand_am[0]>100 && hig_cand_am[0]<=140";
  string c_sbd="!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";

  //list of all cuts of interest
  torch::OrderedDict<string, string> abcdcuts; 
  if (boosted) {
    //Using 2016 method
    abcdcuts.insert("D", "boostedRegionIdx==0");
    abcdcuts.insert("C", "boostedRegionIdx==3");
    abcdcuts.insert("B1", "boostedRegionIdx==1");
    abcdcuts.insert("A1", "boostedRegionIdx==4");
    abcdcuts.insert("B2", "boostedRegionIdx==2");
    abcdcuts.insert("A2", "boostedRegionIdx==5");

    // Using average mass
    abcdcuts.insert("am_D", "amBoostedRegionIdx==0");
    abcdcuts.insert("am_C", "amBoostedRegionIdx==3");
    abcdcuts.insert("am_B1", "amBoostedRegionIdx==1");
    abcdcuts.insert("am_A1", "amBoostedRegionIdx==4");
    abcdcuts.insert("am_B2", "amBoostedRegionIdx==2");
    abcdcuts.insert("am_A2", "amBoostedRegionIdx==5");

  } else {
    abcdcuts.insert("sbd2b", c_ab+"&&"+c_sbd);
    abcdcuts.insert("hig2b", c_ab+"&&"+c_hig);
    abcdcuts.insert("sbd3b", c_cd+"&&"+c_sbd);
    abcdcuts.insert("hig3b", c_cd+"&&"+c_hig);
    abcdcuts.insert("sbd4b", c_ef+"&&"+c_sbd);
    abcdcuts.insert("hig4b", c_ef+"&&"+c_hig);
  }

  // Makes a plot for each vector in plots
  vector<aplot> plots;
  if (boosted) {
    for (auto &ibase: baselines) {
      plots.push_back({"1H_2016_"+ibase.key(),vector<string>({"D","C","B1","A1"}), ibase.value(), metcuts});
      plots.push_back({"2H_2016_"+ibase.key(),vector<string>({"D","C","B2","A2"}), ibase.value(), metcuts});
      plots.push_back({"1H_AM_"+ibase.key(),vector<string>({"am_D","am_C","am_B1","am_A1"}), ibase.value(), metcuts});
      plots.push_back({"2H_AM_"+ibase.key(),vector<string>({"am_D","am_C","am_B2","am_A2"}), ibase.value(), metcuts});
    }
  } else {
    for (auto &ibase: baselines) {
      plots.push_back({"3b_"+ibase.key(),vector<string>({"sbd2b","hig2b","sbd3b","hig3b"}), ibase.value(), metcuts});
      plots.push_back({"4b_"+ibase.key(),vector<string>({"sbd2b","hig2b","sbd4b","hig4b"}), ibase.value(), metcuts});
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////
  PlotMaker pm;
  vector<vector<vector<NamedFunc> > > allcuts(plots.size(), vector<vector<NamedFunc> > (4));
  for(size_t iplot=0; iplot<plots.size(); iplot++){
    for(size_t iabcd=0; iabcd<plots[iplot].abcd.size(); iabcd++){
      vector<TableRow> table_cuts;
      for(size_t ibin=0; ibin<plots[iplot].bincuts.size(); ibin++){
        NamedFunc totcut=plots[iplot].baseline && plots[iplot].bincuts[ibin] && abcdcuts[plots[iplot].abcd[iabcd]];
        table_cuts.push_back(TableRow("", totcut,0,0,wgt));
        allcuts[iplot][iabcd].push_back(totcut);
      } // Loop over bins
      TString tname = "rmj"; tname += iplot; tname += iabcd;
      pm.Push<Table>(tname.Data(),  table_cuts, procs, false, false);
    } // Loop over abcdcuts
  } // Loop over plots

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////

  for(size_t iplot=0; iplot<plots.size(); iplot++){
    vector<vector<vector<GammaParams> > > allyields(procs.size(), vector<vector<GammaParams> >(4));
    for(size_t iabcd=0; iabcd<plots[iplot].abcd.size(); iabcd++){
      Table * yield_table = static_cast<Table*>(pm.Figures()[iplot*4+iabcd].get());
      for(size_t ibkg=0; ibkg<procs.size(); ibkg++)
        allyields[ibkg][iabcd] = yield_table->Yield(procs[ibkg].get(), lumi);
    } // Loop over ABCD cuts

    //// Print MC/Data yields, cuts applied, kappas, preds
    // if(debug) printDebug(allcuts[iplot], allyields, baseline, procs);

    //// RMJ/RMJ
    vector<vector<vector<int> > > indices;
    vector<TString> leglabels;

    //// 3b kappa
    indices.clear(); leglabels.clear();
    for(int ibkg=0; ibkg<static_cast<int>(procs.size()); ibkg++) { 
      // defining kappa as (3/2) /(1/0)
      indices.push_back(vector<vector<int> >({{ibkg, 3, 1},  {ibkg, 2, -1}, {ibkg, 1, -1}, {ibkg, 0, 1}}));
      leglabels.push_back(procs[ibkg]->name_);
    }
    plotRatio(allyields, plots[iplot], indices, leglabels, procs, abcdcuts);
  } // Loop over plots


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<plots.size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plotRatio(vector<vector<vector<GammaParams> > > &allyields, aplot &plotdef,
               vector<vector<vector<int> > > &indices, vector<TString> &leglabels, 
               vector<shared_ptr<Process> > &procs, torch::OrderedDict<string, string> abcdcuts){

  size_t ngraphs = indices.size();
  size_t nbins = allyields[0][0].size();

  //// Finding all ratios for all graphs
  float val(1.), valup(1.), valdown(1.);
  vector<vector<vector<float> > > ratios(ngraphs);
  float maxr=-1., minr=1e6;
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    // Finding powers to calculate ratio
    vector<float> powers;
    for(size_t ipow=0; ipow<indices[igraph].size(); ipow++) powers.push_back(indices[igraph][ipow][2]);

    // Finding ratios for each bin
    for(size_t ibin=0; ibin<nbins; ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      for(size_t ind=0; ind<indices[igraph].size(); ind++) {
        size_t ibkg = indices[igraph][ind][0];
        size_t iabcd = indices[igraph][ind][1];
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[ibkg][iabcd][ibin].NEffective());
        weights.back().push_back(allyields[ibkg][iabcd][ibin].Weight());
      } // Loop over indices

      // Throwing toys to find ratios and uncertainties
      val = calcKappa(entries, weights, powers, valdown, valup);
      if(valdown<0) valdown = 0;
      ratios[igraph].push_back(vector<float>({val, valdown, valup}));
      if(maxr < val+valup) maxr = val+valup;
      if(minr > val-valdown) minr = val-valdown;
    } // Loop over bins
  } // Loop over graphs

  //// Finding ytitle
  TString ytitle="#kappa";
  if(indices[0].size()==2){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    if((int(ind0)==abcdcuts.idx("hig2b") && int(ind1)==abcdcuts.idx("sbd2b"))) ytitle = "Hig/Sbd";
  }
  if(indices[0].size()==4){
    // size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    // size_t ind2=indices[0][2][1], ind3=indices[0][3][1];
  }
  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Ratio");
  setPlotStyle(opts);

  //// Plotting kappas
  TCanvas can("can","");
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);

  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = maxr*1.2;
  if(maxy>3) maxy = 3;
  if(maxy<2) maxy = 2;
  TH1D histo("histo", plotdef.name, nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.GetYaxis()->SetTitleSize(0.06);
  histo.GetYaxis()->SetTitleOffset(0.8);
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  vector<vector<double> > vx(ngraphs), vexh(ngraphs), vexl(ngraphs);
  vector<vector<double> > vy(ngraphs), veyh(ngraphs), veyl(ngraphs);
  for(size_t ibin=0; ibin<nbins; ibin++){
    histo.GetXaxis()->SetBinLabel(ibin+1, CodeToRootTex(plotdef.bincuts[ibin].Data()).c_str());
    // xval is the x position of the first marker in the group
    double xval = ibin+1, minxb = 0.15, binw = 0;
    // If there is more than one point in the group, it starts minxb to the left of the center of the bin
    // binw is the distance between points in the njets group
    if(ngraphs>1) {
      xval -= minxb;
      binw = 2*minxb/(ngraphs-1);
    }
    for(size_t igraph=0; igraph<ngraphs; igraph++){
      vx[igraph].push_back(xval);
      xval += binw;
      vexl[igraph].push_back(0);
      vexh[igraph].push_back(0);
      vy[igraph]  .push_back(ratios[igraph][ibin][0]);
      veyl[igraph].push_back(ratios[igraph][ibin][1]);
      veyh[igraph].push_back(ratios[igraph][ibin][2]);
    } // Loop over TGraphs
  } // Loop over bin cuts

  //// Drawing legend and TGraphs
  double legX(opts.LeftMargin()+0.023), legY(1-opts.TopMargin()-0.03), legSingle = 0.05;
  double legW = 0.19*ngraphs, legH = legSingle;
  int Ncol = ngraphs;
  if(ngraphs>3) {
    legH *= 2;
    Ncol = (ngraphs+1)/2;
    legW = 0.25*Ncol;
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.2); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(Ncol);

  Palette colors("txt/colors.txt", "default");
  vector<int> mcolors, styles;
  for(size_t ibkg=0; ibkg<procs.size(); ibkg++) {
    mcolors.push_back(procs[ibkg]->color_);
    styles.push_back(19+ibkg);
  }
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    graph[igraph] = TGraphAsymmErrors(vx[igraph].size(), &(vx[igraph][0]), &(vy[igraph][0]),
                                    &(vexl[igraph][0]), &(vexh[igraph][0]), &(veyl[igraph][0]), &(veyh[igraph][0]));
    graph[igraph].SetMarkerStyle(styles[igraph]); 
    if(leglabels[igraph].Contains("All")) graph[igraph].SetMarkerStyle(21);
    graph[igraph].SetMarkerSize(1.4);
    graph[igraph].SetMarkerColor(mcolors[igraph]);
    graph[igraph].SetLineColor(mcolors[igraph]); graph[igraph].SetLineWidth(2);
    graph[igraph].Draw("p0 same");
    leg.AddEntry(&graph[igraph], leglabels[igraph], "p");
  } // Loop over TGraphs
  leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  //cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  //cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname = "plots/ratio_"+sample_name+"_"+plotdef.name+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

} // plotRatio



// allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
// void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
//                 TString baseline, vector<shared_ptr<Process> > &procs){
//   int digits = 3;
//   cout<<endl<<endl<<"============================ Printing cuts  ============================"<<endl;
//   cout<<"-- Baseline cuts: "<<baseline<<endl<<endl;
//   for(size_t ibin=0; ibin<allcuts[0].size(); ibin++){
//     for(size_t iabcd=0; iabcd<allcuts.size(); iabcd++){
//       for(size_t ibkg=0; ibkg<procs.size(); ibkg++)
//         cout<<procs[ibkg]->name_<<": "    <<setw(9)<<RoundNumber(allyields[ibkg][iabcd][ibin].Yield(), digits)<<", ";
//       cout<<"  - "<< allcuts[iabcd][ibin]<<endl;
//     } // Loop over ABCD cuts
//     cout<<endl;
//   } // Loop over bin cuts

// } // printDebug


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},      // Luminosity to normalize MC with (no data)
      {"sample", required_argument, 0, 's'},    // Which sample_name to use: standard, met150, 2015 data
      {"boo", no_argument, 0, 0},       // resolved or boosted?
      {"tag", required_argument, 0, 't'},    //
      {"debug", no_argument, 0, 'd'},           // Debug: prints yields and cuts used
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:s:t:o:d", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      sample_name = optarg;
      break;
    case 't':
      tag = optarg;
      break;
    case 'd':
      debug = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "boo"){
        boosted = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
