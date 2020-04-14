#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 1;
  string cr_sample = "search";
  bool do_data = false;
  bool do_metint = true; //Integrate in MET?
  bool do_nbint = false;  // Integrate all Nb?
  bool do_ge3b = false;  // Integrate 3b and 4b?
  // signal points to include and their colors
  vector<string> sigm = {"400","900"}; 
  vector<int> sig_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
}
  
int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  if (do_data) {
    log_lumi_info = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
    lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
    log_lumi = log_lumi().Bottom(BottomType::ratio);
    lin_lumi = lin_lumi().Bottom(BottomType::ratio);
  }
  vector<PlotOpt> linplot = {lin_lumi_info};
  vector<PlotOpt> logplot = {log_lumi_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years; 
  years = {2016, 2017, 2018};

  string base_dir(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string mc_skim_dir("mc/merged_higmc_preselect/"), data_skim_dir("mc/merged_higdata_higloose/");
  if (cr_sample=="ttbar")    {mc_skim_dir = "mc/merged_higmc_higlep1T/"; data_skim_dir = "merged_higdata_higlep1T/";} 
  else if (cr_sample=="zll") {mc_skim_dir = "mc/merged_higmc_higlep2T/"; data_skim_dir = "merged_higdata_higlep2T/";} 
  else if (cr_sample=="qcd") {mc_skim_dir = "mc/merged_higmc_higqcd/";  data_skim_dir = "merged_higdata_higqcd/";} 
  string sig_skim_dir("SMS-TChiHH_2D/mergednn_higmc_preselect/");


  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", 
                                    "*_TTZ*.root", "*_TTW*.root",
                                    "*_TTGJets*.root", "*_ttHTobb*.root","*_TTTT*.root"
                                  });
  mctags["vjets"]   = set<string>({ 
                                   "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"
                                 });
  mctags["qcd"]   = set<string>({
                                   "*QCD_HT200to300_Tune*", // these have too low weights
                                   "*QCD_HT300to500_Tune*", 
                                   "*QCD_HT500to700_Tune*",
                                   "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                   "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*",
                                 });
  mctags["other"]   = set<string>({
                                   "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*_ST_*.root"
                                 });
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baseline cuts ///////////////////////////////////////
  string baseline_s = "pass && stitch && met/mht<2 && met/met_calo<2 && njet>=4 && njet<=5";
  baseline_s += " && nvlep==0 && ntk==0 && !low_dphi_met";
  if (cr_sample=="ttbar") baseline_s += " && nlep==1 && mt<100";
  else if (cr_sample=="zll") baseline_s += " && nlep==2 && met<50";
  else if (cr_sample=="qcd") baseline_s += " && nvlep==0 && ntk==0 && low_dphi_met";
  NamedFunc baseline = baseline_s;

  // Baseline definitions
  NamedFunc wgt = "w_lumi*w_isr"*HigUtilities::w_years*Higfuncs::eff_higtrig;//Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  string cr_sample_label = "";
  if (cr_sample=="zll") cr_sample_label = "2L CS";
  if (cr_sample=="qcd") cr_sample_label = "Low #Delta#phi CS";
  if (cr_sample=="ttbar") cr_sample_label = "1L CS";

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(base_dir, years, mc_skim_dir,mctags["qcd"]),  baseline)); // kill the huge error band in the ratio plot 
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", 
      Process::Type::background, colors("tt_1l"),    attach_folder(base_dir, years, mc_skim_dir,mctags["ttx"]),  baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(base_dir, years, mc_skim_dir,mctags["vjets"]),baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("Other",      
    Process::Type::background, kGreen-2,           attach_folder(base_dir, years, mc_skim_dir,mctags["other"]),baseline));      

  if (do_data) {
    procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
      attach_folder(base_dir, years, data_skim_dir, {"*root"}),  baseline)); 
  }
  if (cr_sample == "search") {
    for (unsigned isig(0); isig<sigm.size(); isig++)
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        sig_colors[isig], attach_folder(base_dir, years, mc_skim_dir, {"*TChiHH_mGluino-"+sigm[isig]+"*.root"}), baseline));
  }

  PlotMaker pm;

  vector<string> metcuts;
  string metdef = "met";
  if (cr_sample=="zll") metdef = "ll_pt[0]";
  if (do_metint) {
    if (cr_sample=="zll") metcuts.push_back(metdef+">0");
    else if (cr_sample=="ttbar") metcuts.push_back(metdef+">0");
    else metcuts.push_back(metdef+">150");
  } else {
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300&&"+metdef+"<=400");
    metcuts.push_back(metdef+">400");
  }
  
  vector<string> nbcuts;
  if (cr_sample=="zll" || cr_sample=="qcd") {
    if (do_nbint) {
      nbcuts.push_back("(nbm==0 || nbm==1)");
    } else {
      nbcuts.push_back("nbm==0");
      nbcuts.push_back("nbm==1");
    }
  } 
  if (cr_sample=="ttbar" || cr_sample=="search") {
    if (do_nbint) {
      nbcuts.push_back("nbt>=2");
    } else {
      nbcuts.push_back("nbt==2&&nbm==2");
      if (do_ge3b) {
        nbcuts.push_back("nbt>=2&&nbm>=3");
      } else {
        nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
        nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
      }
    }
  }


  vector<TString> vc_drmax;
  vc_drmax.push_back("1");
  // vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  // vc_drmax.push_back("hig_cand_drmax[0]<=1.1");

  string hig_trim = "hig_cand_dm[0]<=40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<=2.2";

  pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), "1", procs, linplot).Weight(wgt).Tag(cr_sample);
  pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), "1", procs, linplot).Weight(wgt).Tag(cr_sample);
  pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), "1", procs, linplot).Weight(wgt).Tag(cr_sample);

  pm.Push<Hist1D>(Axis(20,0,2000,"ht", "H_{T} [GeV]"), "1", procs, logplot).Weight(wgt).Tag(cr_sample); 

  if (cr_sample=="zll") {
    pm.Push<Hist1D>(Axis(24,0,600,metdef, "p_{T}^{Z} [GeV]",{200,300,400}), 
      hig_trim, procs, logplot).Weight(wgt).Tag(cr_sample);
  } else {
    pm.Push<Hist1D>(Axis(9,150,600,"met", "E_{T}^{miss} [GeV]",{200,300,400}), 
      hig_trim, procs, linplot).Weight(wgt).Tag(cr_sample);
  } 
  pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_bcat, "b-tag category (TTML)"), 
    hig_bcat>0., procs, linplot).Weight(wgt).Tag(cr_sample);

  for(auto &imet: metcuts) { 
    for(auto &inb: nbcuts) {
      for(auto &idrmax: vc_drmax) {
        pm.Push<Hist1D>(Axis(24,0,240,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
          imet+"&&"+inb+"&&hig_cand_dm[0]<40", 
          procs, linplot).Weight(wgt).Tag(cr_sample).RightLabel({cr_sample_label});
        pm.Push<Hist1D>(Axis(24,0,240,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
          imet+"&&"+inb+"&&"+idrmax+"&&"+hig_trim, 
          procs, linplot).Weight(wgt).Tag(cr_sample).RightLabel({cr_sample_label});
        pm.Push<Hist1D>(Axis(15,0,150,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}), 
          imet+"&&"+inb+"&& hig_cand_drmax[0]<=2.2", procs, linplot).Weight(wgt).Tag(cr_sample)
          .RightLabel({cr_sample_label});
        pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
          imet+"&&"+inb+"&& hig_cand_dm[0]<=40", procs, linplot).Weight(wgt).Tag(cr_sample)
          .RightLabel({cr_sample_label});
      }
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    // Which cr_sample to use: standard, met150, 2015 data         
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      cr_sample = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      printf("Bad option! Found option name %s\n", optname.c_str());
      exit(1);
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
