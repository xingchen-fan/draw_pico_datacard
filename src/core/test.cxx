#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 35.9;

  string base_path = "";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    base_path = "/net/cms2";
  }
  string bkg_dir = base_path+"/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/mc/skim_met150/";
  string sig_dir = base_path+"/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/unskimmed/";
  string data_dir = base_path+"/cms2r0/pico/NanoAODv5/higgsino_angeles/2016/data/skim_met150/";

  Palette colors("txt/colors.txt", "default");

  // Define 1D plot types of interest
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes).ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  // vector of the plot varieties to make for each Hist1D instance
  // vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
  //                                   log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};
  vector<PlotOpt> all_plot_types = {lin_shapes};

  // Define 2D plot types of interest, and define vectors of plot-type combinations to be passed to each Hist2D instance
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style2D().Stack(StackType::data_norm).Title(TitleType::preliminary)};
  vector<PlotOpt> bkg_pts = {style2D().Stack(StackType::lumi_shapes).Title(TitleType::info)};

  // define paths to various ntuples of interest
  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({bkg_dir+"*TTJets_*Lept*", bkg_dir+"*_TTZ*.root", bkg_dir+"*_TTW*.root",
                                   bkg_dir+"*_TTGJets*.root", bkg_dir+"*ttHTobb*.root",bkg_dir+"*_TTTT*.root"});
  mctags["vjets"]   = set<string>({bkg_dir+"*_ZJet*.root", bkg_dir+"*_WJetsToLNu*.root", bkg_dir+"*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({bkg_dir+"*_ST_*.root"});
  mctags["qcd"]     = set<string>({bkg_dir+"*QCD_HT100to200_Tune*", bkg_dir+"*QCD_HT200to300_Tune*",
                                   bkg_dir+"*QCD_HT300to500_Tune*", 
                                   bkg_dir+"*QCD_HT500to700_Tune*",
                                   bkg_dir+"*QCD_HT700to1000_Tune*", bkg_dir+"*QCD_HT1000to1500_Tune*", 
                                   bkg_dir+"*QCD_HT1500to2000_Tune*", bkg_dir+"*QCD_HT2000toInf_Tune*"});
  mctags["other"]   = set<string>({bkg_dir+"*_WH_HToBB*.root", bkg_dir+"*_ZH_HToBB*.root",
                                     bkg_dir+"*_WWTo*.root", bkg_dir+"*_WZ*.root", bkg_dir+"*_ZZ_*.root"});

  // Define processes of interest
  // Note: One can apply cuts at the time of loading a process
  NamedFunc base_func = "1";
  auto ttx = Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background, colors("tt_1l"),
      mctags["ttx"], base_func+"pass && stitch");
  auto vjets = Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1, 
    mctags["vjets"], base_func+"pass && stitch");
  auto singlet = Process::MakeShared<Baby_pico>("Single t", Process::Type::background, kViolet-7, 
    mctags["singlet"], base_func+"pass && stitch");
  auto qcd = Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),    
    mctags["qcd"], base_func+"pass && stitch");
  auto other = Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGreen-2, 
    mctags["other"], base_func+"pass && stitch");  

  auto tchi700 = Process::MakeShared<Baby_pico>("TChiHH(700,1)", Process::Type::signal, kRed+1,
    {sig_dir+"*_mChi-700_*.root"},"1");
  tchi700->SetMarkerStyle(21);
  tchi700->SetMarkerSize(0.9);

  auto tchi400 = Process::MakeShared<Baby_pico>("TChiHH(400,1)", Process::Type::signal, kRed+1,
    {sig_dir+"*_mChi-400_*.root"},"1");
  tchi400->SetLineStyle(2);
  tchi400->SetMarkerStyle(21);
  tchi400->SetMarkerSize(0.9);

  auto data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
    {data_dir+"*.root"},"pass && HLT_PFMET120_PFMHT120_IDTight");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.);

  // vector of processes to be passed to each Hist1D and/or Hist2D and/or table...
  vector<shared_ptr<Process> > procs = {tchi700, tchi400};
  // vector<shared_ptr<Process> > procs = {data, tchi700, tchi400, ttx, vjets, singlet, qcd , other};

  PlotMaker pm;

  // Standard 1D plots
  pm.Push<Hist1D>(Axis(20, 0, 1000., "met", "MET [GeV]", {200., 300., 450.}),
                  "1", procs, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0.5, 10.5, "njet", "N_{jets}", {3.5,5.5}),
                  "nvlep==0 && met>150", procs, all_plot_types)
    .Tag("changing_tags_and_weights").Weight("1.2345*weight").RatioTitle("Numerator","Denominator");

  // // Standard 2D plots
  // pm.Push<Hist2D>(Axis(48, 0., 1200., "hig_cand_am[0]", "#LTm#GT [GeV]", {250., 400.}),
  //                 Axis(25, 0., 700., "hig_cand_dm[0]", "#Deltam [GeV]", {140.}),
  //                 "nvlep==0 && met>150 && njet>=4 && nbt>=2",
  //                 procs, bkg_hist);
  // pm.Push<Hist2D>(Axis(48, 0., 1200., "hig_cand_am[0]", "#LTm#GT [GeV]", {250., 400.}),
  //                 Axis(25, 0., 700., "hig_cand_dm[0]", "#Deltam [GeV]", {140.}),
  //                 "nvlep==0 && met>150 && njet>=4 && nbt>=2",
  //                 vector<shared_ptr<Process> >{ttx, vjets, tchi400}, bkg_pts);

  // // Cutflow table
  // Table & cutflow = pm.Push<Table>("cutflow", vector<TableRow>{
  //     TableRow("Baseline"),
  //       TableRow("No Selection", "1"),
  //       TableRow("$0\\ell$, $E_{\\text{T}}^{\\text{miss}}>150$", "nvlep==0 && met>150"),
  //       TableRow("$N_{\\text{jets}}\\geq4$", "nvlep==0 && met>150&& njet>=4"),
  //       TableRow("$N_{b,T}\\geq2$", "nvlep==0 && met>150 && njet>=4 && nbt>=2")
  //       }, procs);

  // // Systematic variations
  // NamedFunc mm_wgt = Functions::MismeasurementWeight("txt/sys_weights.cfg", "off");
  // NamedFunc mm_crt = Functions::MismeasurementCorrection("txt/sys_weights.cfg", "off",
  //                                                        Functions::Variation::central);

  // // Event scan
  // pm.Push<EventScan>("scan", true, vector<NamedFunc>{"weight", "met", mm_wgt, mm_crt},
  //                    vector<shared_ptr<Process> >{ttx});

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);

  // vector<GammaParams> yields = cutflow.BackgroundYield(lumi);
  // for(const auto &yield: yields){
  //   cout << yield << endl;
  // }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
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
