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
#include "higgsino/hig_utilities.hpp"

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 137;
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
  if(Contains(hostname, "cms") || Contains(hostname, "compute-") || Contains(hostname, "physics.ucsb.edu"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  string foldermc_base("/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string foldermc_skim("mc/merged_higmc_higloose/");
  string foldersig("/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_unskimmed/");
  set<int> years = {2016, 2017, 2018};

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*"});
  mctags["topx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"});

  //    Baseline definitions
  //-----------------------------------------
  NamedFunc wgt = "w_lumi*w_isr";//"weight"; 
  string baseline = "stitch && !lowDphiFix && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  if (resolved) {
    baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  } else {
    baseline += "&& nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600";
    baseline += " && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  }

  //    Define processes, including intersections
  //--------------------------------------------------
  Palette colors("txt/colors.txt", "default");

  NamedFunc base_filters = "1";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(foldermc_base,years,foldermc_skim, mctags["other"]), base_filters && baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(foldermc_base,years,foldermc_skim,mctags["zjets"]), base_filters && baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(foldermc_base,years,foldermc_skim,mctags["wjets"]), base_filters && baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc_base,years,foldermc_skim, mctags["tt"]), base_filters && baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("t/t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc_base,years,foldermc_skim, mctags["topx"]), base_filters && baseline));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(foldermc_base,years,foldermc_skim, mctags["qcd"]), base_filters && baseline)); 

  if (doSignal) {
    vector<int> sigm = {450, 700, 950};
    for (unsigned isig(0); isig<sigm.size(); isig++){
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+to_string(sigm[isig])+",1)", Process::Type::signal, 
        1, {foldersig+"*TChiHH_mChi-"+to_string(sigm[isig])+"_*.root"}, base_filters && baseline));
      cout<<"Adding: "<<foldersig+"*TChiHH_mChi-"+to_string(sigm[isig])+"_*.root"<<endl;
    }
  }


  //    Useful binning definitions
  //------------------------------------------
  string c_2b = "nbm==2";
  string c_3b = "nbm==3&&nbl==3";
  string c_4b = "nbm>=3&&nbl>=4";
  string hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
  string sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  vector<string> res_abcd = {c_2b+"&&"+sbd, c_3b+"&&"+sbd, c_4b+"&&"+sbd, 
                             c_2b+"&&"+hig, c_3b+"&&"+hig, c_4b+"&&"+hig};

  bool bin_drmax = false;
  bool bin_met = true;
  bool bin_dphi = true;

  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");
  if (!bin_drmax) vc_drmax = {"hig_cand_drmax[0]>0"};

  // assume common MET binning
  vector<TString> vc_met;
  if (resolved) {
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=400");
    vc_met.push_back("met>400");
  } else {
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=500");
    vc_met.push_back("met>500 && met<=700");
    vc_met.push_back("met>700");
    // vc_met.push_back("met>150");
  }   
  if (!bin_met) vc_met = {"met>150"};

  vector<pair<string, NamedFunc>> vc_dphi;
  vc_dphi.push_back({"0 \\leq \\Delta \\phi \\leq 0.5" ,min_jet_dphi>=0. && min_jet_dphi<=0.5});
  //vc_dphi.push_back({"0.5 < \\Delta \\phi \\leq 1", min_jet_dphi>0.5 && min_jet_dphi<=1.});
  //vc_dphi.push_back({"1 < \\Delta \\phi", min_jet_dphi>1.});
  vc_dphi.push_back({"0.5 < \\Delta \\phi", min_jet_dphi>0.5});
  if (!bin_dphi) vc_dphi = {{"0 \\leq \\Delta \\phi", min_jet_dphi>=0.}};

  //        Define and fill tables
  //------------------------------------
  PlotMaker pm;
  vector<TableRow> tablerows;
  unsigned nbins(0);
  if (resolved) {
    for (auto &imet: vc_met) {
      for (auto &idrmax: vc_drmax) {
        for (auto &idphi: vc_dphi) {
          // Label
          tablerows.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+","+idphi.first+"$"));
          for (auto &iabcd: res_abcd) {
            // Cut for row
            NamedFunc res_reg_ = imet+"&&"+idrmax+"&&"+iabcd&&idphi.second;
            // Add row to table
            tablerows.push_back(TableRow("$"+CodeToLatex(iabcd)+"$", res_reg_, 0,0, wgt));
            nbins++;
          }
        }
      }
    }
  } else { //boosted priority
    for (auto &imet: vc_met) {
      tablerows.push_back(TableRow("$"+CodeToLatex(imet.Data())+"$"));
      for (unsigned ireg(0); ireg<6; ireg++) {
        NamedFunc boo_reg_ = imet && "boostedRegionIdx=="+to_string(ireg);
        tablerows.push_back(TableRow(CodeToPlainText("boostedRegionIdx=="+to_string(ireg)), boo_reg_, 0,0, wgt));
        nbins++;
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
