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
  bool do_allbkg = false;
  string cr_sample = "search";
  float lumi = 1;
  bool unblind = false;
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
  if (do_allbkg) { // don't include QCD MC unless in QCD control cr_sample
    if(cr_sample=="qcd") alltags = {"*TTJets_*Lept*", 
                                 "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
                                 "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
                                 "*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                 "*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
    else  alltags = {"*TTJets_*Lept*",
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"
          };
  }
  set<string> allfiles = attach_folder(base_dir, years, mc_skim_dir,alltags);

  // Baseline definitions
  string baseline = "stitch && pass && njet>=4 && njet<=5 && met>150 && met<=200"; // "&& ntk==0" - should not affect kinematics, so don't kill stats
  if (cr_sample=="search") baseline += "&& met>150 && nbt>=2 && nvlep==0 && ntk==0 && !low_dphi_met";
  if (cr_sample=="ttbar")  baseline += "&& met>150&& nbt>=2 && nlep==1  && mt<=100";
  if (cr_sample=="zll")    baseline += "&& nlep==2  && met<50";
  if (cr_sample=="qcd")    baseline += "&& met>150 && nvlep==0 && low_dphi_met";
  
  baseline += "&& hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200";

  ////// Nb cuts
  vector<string> nbcuts;
  unsigned firstnb = 0;
  nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm==1");
  nbcuts.push_back("nbm==2");
  if (cr_sample=="ttbar" || cr_sample=="search") {
    firstnb = 2;
    nbcuts.push_back("nbm==3 && nbl==3");
    nbcuts.push_back("nbm>=3 && nbl>=4");
  } 

  string sname = "t#bar{t}+X";
  if (cr_sample=="qcd") sname = "QCD";
  if (cr_sample=="zll") sname = "Z#rightarrow ll";
  if (do_allbkg) sname = "Bkg.";

  vector<int> colors = {kGreen+3, kGreen+1, kOrange, kAzure+1, kBlue+1};

  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    // if (cr_sample=="qcd" && inb==nbcuts.size()-1) continue;
    procs.push_back(Process::MakeShared<Baby_pico>(sname+" "+RoundNumber(inb,0).Data()+"b", 
      Process::Type::background, colors[inb], allfiles, baseline+"&&"+nbcuts[inb]));
  }

  vector<int> colors_trub = {kAzure-4, kTeal-8, kOrange-4, kPink+2, kMagenta-1};
  vector<shared_ptr<Process> > procs_trub = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    procs_trub.push_back(Process::MakeShared<Baby_pico>(sname+" "+RoundNumber(inb,0).Data()+" B-hadrons", 
        Process::Type::background, colors_trub[inb], allfiles, Functions::ntrub==inb && baseline));
  }

  string lumi_s=RoundNumber(lumi,1).Data();

  vector<vector<string>> combos, combos_labels;
  int color_data = kBlue-7;
  if (cr_sample=="zll") { // do 0b vs 1b
    color_data = kOrange+1;
    combos.push_back({"nbm==0","nbm==1"});
    combos.push_back({"nbm==1","nbm==2"});
  } else if (cr_sample=="qcd") { // do 0b vs 1b and 2b vs 3+b
    color_data = kOrange;
    combos.push_back({"nbm==0","nbm==1"});
    combos.push_back({"nbm==1","nbm==2"});
    combos.push_back({"nbm==2","nbm==3"});
    combos.push_back({"nbm==2", "nbm>=3 && nbl>=4"});
  } else if (cr_sample=="search" || cr_sample=="ttbar") {
    combos.push_back({"nbm==2", "nbm==3 && nbl==3"});
    combos.push_back({"nbm==2", "nbm>=3 && nbl>=4"});
  }
  for (unsigned ind(0); ind<combos.size(); ind++){
    vector<string> label;
    for(unsigned lab=0; lab<2; lab++){
      if (Contains(combos[ind][lab],"nbm==0")) label.push_back("0b");
      if (Contains(combos[ind][lab],"nbm==1")) label.push_back("1b");
      if (Contains(combos[ind][lab],"nbm==2")) label.push_back("2b");
      if (Contains(combos[ind][lab],"nbm==3")) label.push_back("3b");
      if (Contains(combos[ind][lab],"nbm>=3") && !Contains(combos[ind][lab],"nbl>=4")) label.push_back("#geq 3b");
      if (Contains(combos[ind][lab],"nbl>=4")) label.push_back("4b");
    }
    combos_labels.push_back(label);
  }

  vector<vector<shared_ptr<Process> >> procs_data;
  if (unblind) {
    for (unsigned ind(0); ind<combos.size(); ind++){
      vector<string> icomb = combos[ind], ilab = combos_labels[ind];
      procs_data.push_back(vector<shared_ptr<Process> >());
      procs_data.back().push_back(Process::MakeShared<Baby_pico>(ilab[0]+" Data "+lumi_s+" fb^{-1}", 
                       Process::Type::background, kBlack, 
                       attach_folder(base_dir, years, data_skim_dir,{"*root"}),
                       /*Higfuncs::trig_hig &&*/ icomb[0]));
      procs_data.back().back()->SetFillColor(color_data);
      procs_data.back().back()->SetLineColor(color_data);
      procs_data.back().back()->SetLineWidth(2);

      procs_data.back().push_back(Process::MakeShared<Baby_pico>(ilab[1]+" Data "+lumi_s+" fb^{-1}", 
                       Process::Type::data, kBlack, 
                       attach_folder(base_dir, years, data_skim_dir,{"*root"}),
                        /*Higfuncs::trig_hig &&*/ icomb[1]));
    }
  }


  string metcut = "met>300";
  if (cr_sample=="zll") metcut = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>0";
  else if (cr_sample=="ttbar") metcut = "1";

  vector<string> xcuts;
  xcuts.push_back("1");
  // xcuts.push_back("hig_cand_drmax[0]<=1.1");
  // xcuts.push_back("hig_cand_drmax[0]>1.1");

  PlotMaker pm;
  NamedFunc wgt = "w_lumi*w_isr" * Higfuncs::w_years*Higfuncs::eff_higtrig;//Higfuncs::weight_higd;

  for (unsigned ic(0); ic<xcuts.size(); ic++){
    if (unblind && cr_sample!="search") {
      for (unsigned i(0); i<combos.size(); i++)
        pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
          baseline+"&&"+xcuts[ic], procs_data[i], plt_types)
    .Tag(cr_sample+"_datavdata"+to_string(i)).RatioTitle(combos_labels[i][1],combos_labels[i][0]);
    }

    pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      baseline+"&&"+xcuts[ic], procs, plt_types).Weight(wgt).Tag(cr_sample+"_shape_bcats")
      .RatioTitle("Bkg. nb","Bkg. 2b");
    pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      baseline+"&&"+xcuts[ic], procs_trub, plt_types).Weight(wgt).Tag(cr_sample+"_shape_trub");

    pm.Push<Hist1D>(Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}),
      baseline+"&&"+xcuts[ic], procs, plt_types).Weight(wgt).Tag(cr_sample+"_shape_bcats")
      .RatioTitle("Bkg. nb","Bkg. 2b");
    pm.Push<Hist1D>(Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}),
      baseline+"&&"+xcuts[ic], procs_trub, plt_types).Weight(wgt).Tag(cr_sample+"_shape_trub");

    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
      baseline+"&&"+xcuts[ic], procs, plt_types).Weight(wgt).Tag(cr_sample+"_shape_bcats")
      .RatioTitle("Bkg. nb","Bkg. 2b");
    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
      baseline+"&&"+xcuts[ic], procs_trub, plt_types).Weight(wgt).Tag(cr_sample+"_shape_trub");

  }

  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    
      {"unblind", no_argument, 0, 'u'},   
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "js:tau", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'u':
      unblind = true;
      break;
    case 'a':
      do_allbkg = true;
      break;
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
