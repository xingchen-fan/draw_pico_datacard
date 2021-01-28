///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 1;
  bool doSignal = false;
  bool resolved = true; 
  bool quick = false;
  string tag = "";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  string foldermc_base(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string foldermc_skim("mc/merged_higmc_preselect/");
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/");
  //set<int> years = {2016, 2017, 2018};
  set<int> years = {2016};

  map<string, set<string>> mctags; 
  mctags["!!DYJetsToLL"]   = set<string>({"*DYJetsToLL*.root"});
  mctags["!!QCD_HT"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*",
                                   "*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["!!ST_s-channel"]   = set<string>({"*ST_s-channel*.root"});
  mctags["!!ST_t-channel"]   = set<string>({"*ST_t-channel*.root"});
  mctags["!!ST_tW"]   = set<string>({"*ST_tW*.root"});
  mctags["!!TTGJets"]   = set<string>({"*TTGJets*.root"});
  mctags["!!TTJets_DiLept_Tune"]   = set<string>({"*TTJets_DiLept_Tune*.root"});
  mctags["!!TTJets_DiLept_genMET-150"]   = set<string>({"*TTJets_DiLept_genMET-150*.root"}); 
  mctags["!!TTJets_SingleLept_Tune"]   = set<string>({"*TTJets_SingleLeptFromTbar_Tune*.root","*TTJets_SingleLeptFromT_Tune*.root"});
  mctags["!!TTJets_SingleLept_genMET-150"]   = set<string>({"*TTJets_SingleLept*_genMET-150*.root"});
  mctags["!!TTTT"]   = set<string>({"*TTTT*.root"});
  mctags["!!TTWJetsToLNu"]   = set<string>({"*TTWJetsToLNu*.root"});
  mctags["!!TTWJetsToQQ"]   = set<string>({"*TTWJetsToQQ*.root"});
  mctags["!!TTZToLLNuNu_M-10"]   = set<string>({"*TTZToLLNuNu_M-10*.root"});
  mctags["!!TTZToQQ"]   = set<string>({"*TTZToQQ*.root"});
  mctags["!!WH_HToBB_WToLNu_M125"]   = set<string>({"*WH_HToBB_WToLNu_M125*.root"});
  mctags["!!WJetsToLNu_HT"]   = set<string>({"*WJetsToLNu_HT*.root"});
  mctags["!!WJetsToLNu_Tune"]   = set<string>({"*WJetsToLNu_Tune*.root"});
  mctags["!!WWTo2L2Nu"]   = set<string>({"*WWTo2L2Nu*.root"});
  mctags["!!WWToLNuQQ"]   = set<string>({"*WWToLNuQQ*.root"});
  mctags["!!WZTo1L1Nu2Q"]   = set<string>({"*WZTo1L1Nu2Q*.root"});
  mctags["!!WZTo1L3Nu"]   = set<string>({"*WZTo1L3Nu*.root"});
  mctags["!!WZTo2L2Q"]   = set<string>({"*WZTo2L2Q*.root"});
  mctags["!!WZTo3LNu"]   = set<string>({"*WZTo3LNu*.root"});
  mctags["!!ZH_HToBB_ZToNuNu"]   = set<string>({"*ZH_HToBB_ZToNuNu*.root"});
  mctags["!!ZJetsToNuNu_HT"]   = set<string>({"*ZJetsToNuNu_HT*.root"});
  mctags["!!ZZ_Tune"]   = set<string>({"*ZZ_Tune*.root"});
  mctags["!!ttHTobb"]   = set<string>({"*ttHTobb*.root"});

  //mctags["ttx"]     = set<string>({"*TTJets_*Lept*","*_TTZ*.root", "*_TTW*.root",
  //                                   "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  //mctags["vjets"]   = set<string>({"*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ZJet*.root"});
  //mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*",
  //                                 "*_QCD_HT500to700_*",
  //                                 "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
  //                                 "*_QCD_HT2000toInf_*"});
  //mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
  //                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*_ST_*.root"});

  //    Baseline definitions
  //-----------------------------------------
  NamedFunc wgt = "w_lumi*w_isr"*Higfuncs::eff_higtrig;
  if (years.size()==1 && *years.begin()==2016) wgt *= "137.";
  else wgt *= Higfuncs::w_years;
  string baseline = "stitch && nvlep==0 && ntk==0 && !low_dphi_met && met>150";
  string higtrim = "njet>=4&&njet<=5&&hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  if (resolved) {
    baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  }

  //    Define processes, including intersections
  //--------------------------------------------------
  Palette colors("txt/colors.txt", "default");

  NamedFunc base_filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  vector<shared_ptr<Process> > procs;
  //procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(foldermc_base,years,foldermc_skim, mctags["other"]), base_filters && baseline));
  //procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kGreen+1,
  //                attach_folder(foldermc_base,years,foldermc_skim,mctags["vjets"]), base_filters && baseline));
  //procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(foldermc_base,years,foldermc_skim, mctags["ttx"]), base_filters && baseline));
  //procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(foldermc_base,years,foldermc_skim, mctags["qcd"]), base_filters && baseline)); 
  // Loop over each sample
  for (auto it = mctags.begin(); it != mctags.end(); ++it) {
    procs.push_back(Process::MakeShared<Baby_pico>(it->first, Process::Type::background, kGray+2,
                    attach_folder(foldermc_base,years,foldermc_skim, it->second), base_filters && baseline));
  }

  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");
  vector<TString> vc_met;
  vc_met.push_back("met>150 && met<=200");
  vc_met.push_back("met>200 && met<=300");
  vc_met.push_back("met>300 && met<=400");
  vc_met.push_back("met>400");

  if (doSignal) {
    vector<string> sigm = {"175", "400", "650"};
    for (unsigned isig(0); isig<sigm.size(); isig++){
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        1, {foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}, base_filters && baseline));
      cout<<foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"<<endl;
    }
  }

  //        Define and fill tables
  //------------------------------------
  PlotMaker pm;
  vector<TableRow> tablerows;
  if (resolved) {
    for (auto &imet: vc_met) {
      for (auto &idrmax: vc_drmax) {
        NamedFunc res_reg = imet+"&&"+idrmax;
        tablerows.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+"$", res_reg, 0,0, wgt));
      }
    }
  }
  TString tabname = "table_"; tabname += (resolved ? "resolved" : "boosted");
  pm.Push<Table>(tabname.Data(), tablerows, procs, 0);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"no_signal", no_argument, 0, 'n'},    
      {"luminosity", required_argument, 0, 'l'},    
      {"tag", required_argument, 0, 't'},
      {"boo", no_argument, 0, 0},
      {"quick", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "nl:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      doSignal = false;
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "boo"){
        resolved = false;
      }else if(optname == "quick"){
        quick = true;
      } else {
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
