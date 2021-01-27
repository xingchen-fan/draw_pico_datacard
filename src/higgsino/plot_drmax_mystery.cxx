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
#include "higgsino/hig_utilities.hpp"
#include "higgsino/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  string cr_sample = "search";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes)
    // .RatioMaximum(3.5)
    .LegendColumns(3);
  // PlotOpt log_shapes_info = lin_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years;
  years = {2016, 2017, 2018};

  string base_dir(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string mc_skim_dir("mc/merged_higmc_preselect/"), data_skim_dir("mc/merged_higdata_higloose/");
  if (cr_sample=="ttbar")    {mc_skim_dir = "mc/merged_higmc_higlep1T/"; data_skim_dir = "merged_higdata_higlep1/";} 
  else if (cr_sample=="zll") {mc_skim_dir = "mc/merged_higmc_higlep2T/"; data_skim_dir = "merged_higdata_higlep2/";} 
  else if (cr_sample=="qcd") {mc_skim_dir = "mc/merged_higmc_higqcd/";  data_skim_dir = "merged_higdata_higqcd/";} 

  set<string> alltags; 
  if (cr_sample=="ttbar" || cr_sample=="search") alltags = {"*TTJets_*Lept*"};//,
                                      // "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", 
                                      // "*ttHTobb*.root","*_TTTT*.root"};
  if (cr_sample=="zll") alltags = {"*DYJetsToLL*.root"};
  if (cr_sample=="qcd") alltags = {//"*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                //"*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"};

  set<string> allfiles = attach_folder(base_dir, years, mc_skim_dir,alltags);

  // Baseline definitions
  string baseline = "stitch && pass && njet>=4 && njet<=5"; // "&& ntk==0" - should not affect kinematics, so don't kill stats
  if (cr_sample=="search") baseline += "&& met>150 && nbt>=2 && nvlep==0 && ntk==0 && !low_dphi_met";
  if (cr_sample=="ttbar")  baseline += "&& nbt>=2 && nlep==1  && mt<=100";
  if (cr_sample=="zll")    baseline += "&& nlep==2  && met<50";
  if (cr_sample=="qcd")    baseline += "&& met>150 && nvlep==0 && low_dphi_met";
  
  baseline += "&& hig_cand_dm[0]<=40";

  ////// Nb cuts
  vector<string> nbcuts;
  if (cr_sample=="zll" || cr_sample=="qcd") {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
    nbcuts.push_back("nbm>=3");
  } else {
    nbcuts.push_back("nbm==2");
    nbcuts.push_back("nbm==3 && nbl==3");
    nbcuts.push_back("nbm>=3 && nbl>=4");    
  }

  string sname = "t#bar{t}+X";
  if (cr_sample=="qcd") sname = "QCD";
  if (cr_sample=="zll") sname = "Z#rightarrow ll";

  vector<int> colors = {kOrange+3, kRed+1, kOrange-3, kSpring+10};
  if (cr_sample=="ttbar" || cr_sample=="search") 
    colors = {kOrange, kAzure+1, kBlue+1};
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(0); inb<nbcuts.size(); inb++){
    procs.push_back(Process::MakeShared<Baby_pico>(sname+" "+CodeToRootTex(nbcuts[inb]), 
      Process::Type::background, colors[inb], allfiles, baseline+"&&"+nbcuts[inb]));
  }

  vector<string> nisr_cuts = {"nisr==0", "nisr==1", "nisr==2"};
  vector<int> colors_isr = {kBlack, kPink+5, kPink+1};
  vector<shared_ptr<Process> > procs_isr = vector<shared_ptr<Process> >();
  for (unsigned iisr(0); iisr<nisr_cuts.size(); iisr++){
    procs_isr.push_back(Process::MakeShared<Baby_pico>(nisr_cuts[iisr], 
      Process::Type::background, colors_isr[iisr], allfiles, baseline+"&&"+nisr_cuts[iisr]));
  }

  vector<TString> metcuts;
  string metdef = "met";
  if (cr_sample=="zll") metdef = "ll_pt[0]";
  metcuts.push_back("1");
  if (cr_sample=="ttbar" || cr_sample=="zll") {
    metcuts.push_back(metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
  }
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  metcuts.push_back(metdef+">300&&"+metdef+"<=400");
  metcuts.push_back(metdef+">400");

  vector<string> dr_cuts;
  dr_cuts.push_back("1");
  dr_cuts.push_back("hig_cand_drmax[0]<=1.1");
  dr_cuts.push_back("hig_cand_drmax[0]>1.1 && hig_cand_drmax[0]<=2.2");
  dr_cuts.push_back("hig_cand_drmax[0]>2.2");

  PlotMaker pm;
  NamedFunc wgt = "w_lumi*w_isr" * Higfuncs::w_years;//Higfuncs::weight_higd;

  for (auto idr: dr_cuts) {
    pm.Push<Hist1D>(Axis(13,0,260,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      idr+"&&"+baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");
    pm.Push<Hist1D>(Axis(13,0,260,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      idr+"&&"+baseline, procs_isr, plt_types).Weight(wgt).Tag(cr_sample+"_isrcats").RatioTitle("N_{ISR}","N_{ISR}=0");
  }

  for (auto imet: metcuts) 
    pm.Push<Hist1D>(Axis(13,0,260,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      imet+"&&"+baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");

  pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
    baseline, procs_isr, plt_types).Weight(wgt).Tag(cr_sample+"_isrcats").RatioTitle("N_{ISR}","N_{ISR}=0");
  for (auto iisr: nisr_cuts) 
    pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      iisr+"&&"+baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");

  pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
    baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");
  for (auto imet: metcuts) 
    pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
      imet+"&&"+baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");

  for (auto idr: dr_cuts)
    pm.Push<Hist1D>(Axis(20,0,400,"isr_tru_pt", "p_{T, ISR}", {}),
      idr+"&&"+baseline, procs, plt_types).Weight(wgt).Tag(cr_sample+"_bcats").RatioTitle("Bkg. nb","Bkg. 2b");

  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      cr_sample = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
