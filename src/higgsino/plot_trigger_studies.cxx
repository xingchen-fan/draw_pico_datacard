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
#include "TFile.h"
#include "TH1.h"

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
#include "core/cross_sections.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018temp.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  //flags for which things to plot
  int year = 2018;

  //don't change this, automatically set below
  double lumi = 35.9;
  Palette colors("txt/colors.txt", "default");
  // Define 1D+2D plot types of interest
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow);
  vector<PlotOpt> all_plot_types = {lin_lumi};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).Overflow(OverflowType::overflow)};

  vector<shared_ptr<Process> > procs_data = {};
  vector<shared_ptr<Process> > procs_mc = {};
  shared_ptr<Process> pro_data;
  shared_ptr<Process> pro_mc;
  //_all processes are used only for plots that don't have a MET150 cut so we can speed up things with skim
  vector<shared_ptr<Process> > procs_data_all = {};
  shared_ptr<Process> pro_data_all;
  vector<shared_ptr<Process> > procs_mc_all = {};
  shared_ptr<Process> pro_mc_all;
  if (year == 2017) {
  	lumi = 41.5;
  	string data_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2017/data/";
  	string mc_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2017/mc/";
  	set<string> str_data({data_dir+"skim_met150/raw_pico_met150_SingleElectron*.root",data_dir+"skim_met150/raw_pico_met150_MET*.root",data_dir+"skim_met150/raw_pico_met150_SingleMuon*.root",data_dir+"skim_met150/raw_pico_met150_JetHT*.root"});
  	pro_data = Process::MakeShared<Baby_pico>("2017 Data", Process::Type::data, kBlack, str_data, "stitch");
  	set<string> str_mc({mc_dir+"skim_met150/pico_met150_TTJets_SingleLeptFromT*.root"});
  	pro_mc = Process::MakeShared<Baby_pico>("2017 MC", Process::Type::background, kBlack, str_mc, "stitch");
  	set<string> str_data_all({data_dir+"raw_pico/raw_pico_SingleElectron*.root",data_dir+"raw_pico/raw_pico_MET*.root",data_dir+"raw_pico/raw_pico_SingleMuon*.root",data_dir+"raw_pico/raw_pico_JetHT*.root"});
  	pro_data_all = Process::MakeShared<Baby_pico>("2017 Data", Process::Type::data, kBlack, str_data_all, "stitch");
  	set<string> str_mc_all({mc_dir+"unskimmed/pico_TTJets_SingleLeptFromT*.root"});
  	pro_mc_all = Process::MakeShared<Baby_pico>("2017 MC", Process::Type::background, kBlack, str_mc_all, "stitch");
  }
  else if (year == 2018) {
  	lumi = 60.0;
  	string data_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/data/";
  	string mc_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/mc/";
  	set<string> str_data({data_dir+"skim_met150/raw_pico_met150_EGamma*.root",data_dir+"skim_met150/raw_pico_met150_MET*.root",data_dir+"skim_met150/raw_pico_met150_SingleMuon*.root",data_dir+"skim_met150/raw_pico_met150_JetHT*.root"});
  	pro_data = Process::MakeShared<Baby_pico>("2018 Data", Process::Type::data, kBlack, str_data, "stitch");
  	set<string> str_mc({mc_dir+"skim_met150/pico_met150_TTJets_SingleLeptFromT*.root"});
  	pro_mc = Process::MakeShared<Baby_pico>("2018 MC", Process::Type::background, kBlack, str_mc, "stitch");
  	set<string> str_data_all({data_dir+"raw_pico/raw_pico_EGamma*.root",data_dir+"raw_pico/raw_pico_MET*.root",data_dir+"raw_pico/raw_pico_SingleMuon*.root",data_dir+"raw_pico/raw_pico_JetHT*.root"});
  	pro_data_all = Process::MakeShared<Baby_pico>("2018 Data", Process::Type::data, kBlack, str_data_all, "stitch");
  	set<string> str_mc_all({mc_dir+"unskimmed/pico_TTJets_SingleLeptFromT*.root"});
  	pro_mc_all = Process::MakeShared<Baby_pico>("2018 MC", Process::Type::background, kBlack, str_mc_all, "stitch");
  }
  else {
  	string data_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/data/";
  	string mc_dir = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_humboldt/2016/mc/";
	//default 2016
  	set<string> str_data({data_dir+"skim_met150/raw_pico_met150_SingleElectron*.root",data_dir+"skim_met150/raw_pico_met150_MET*.root",data_dir+"skim_met150/raw_pico_met150_SingleMuon*.root",data_dir+"skim_met150/raw_pico_met150_JetHT*.root"});
  	//set<string> str_data({data_dir+"raw_pico_MET__Run2016B_ver1__*.root",data_dir+"raw_pico_SingleElectron__Run2016C__*.root",data_dir+"raw_pico_SingleMuon__Run2016B_ver1__*.root",data_dir+"raw_pico_JetHT__Run2016B_ver1__*.root"});
  	pro_data = Process::MakeShared<Baby_pico>("2016 Data", Process::Type::data, kBlack, str_data, "stitch");
  	//set<string> str_mc({mc_dir+"/skim_met150/pico_met150*.root"});
  	set<string> str_mc({mc_dir+"skim_met150/pico_met150_TTJets_SingleLeptFromT*.root"});
  	pro_mc = Process::MakeShared<Baby_pico>("2016 MC", Process::Type::background, kBlack, str_mc, "stitch");
  	set<string> str_data_all({data_dir+"raw_pico/raw_pico_SingleElectron*.root",data_dir+"raw_pico/raw_pico_MET*.root",data_dir+"raw_pico/raw_pico_SingleMuon*.root",data_dir+"raw_pico/raw_pico_JetHT*.root"});
  	pro_data_all = Process::MakeShared<Baby_pico>("2016 Data", Process::Type::data, kBlack, str_data, "stitch");
  	set<string> str_mc_all({mc_dir+"unskimmed/pico_TTJets_SingleLeptFromT*.root"});
  	pro_mc_all = Process::MakeShared<Baby_pico>("2016 MC", Process::Type::background, kBlack, str_mc_all, "stitch");
  }
  procs_data.push_back(pro_data);
  procs_mc.push_back(pro_mc);
  procs_data_all.push_back(pro_data_all);
  procs_mc_all.push_back(pro_mc_all);

  // vector of processes to be passed to each Hist1D and/or Hist2D and/or table...
  //vector<shared_ptr<Process> > procs = {ttz, ttx, vjets, singlet, qcd , other};

  //named funcs
  //NamedFunc met_trigger = "HLT_PFMET110_PFMHT110_IDTight||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight||HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight||HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60";
  const NamedFunc met_trigger("met_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
		  bool r_met_trigger = b.HLT_PFMET110_PFMHT110_IDTight()||b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight()||b.HLT_PFMET120_PFMHT120_IDTight()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight();
		  return r_met_trigger;
  });

  const NamedFunc met_prescale_2017("met_prescale_2017", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_met_prescale_2017 = 1.0;
		  if ((b.HLT_PFMET110_PFMHT110_IDTight() || b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight()) && !(b.HLT_PFMET120_PFMHT120_IDTight() || b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight()))
		  	r_met_prescale_2017 = 1.0/230.0;
		  return r_met_prescale_2017;
  });

  const NamedFunc el_trigger("el_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
		  bool r_el_trigger = b.HLT_Ele27_WPTight_Gsf()||b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele115_CaloIdVT_GsfTrkIdT();
		  return r_el_trigger;
  });

  const NamedFunc mu_trigger("mu_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
		  bool r_mu_trigger = b.HLT_IsoMu24()||b.HLT_IsoMu27()||b.HLT_Mu50();
		  return r_mu_trigger;
  });

  const NamedFunc jet_trigger("jet_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
		  return b.HLT_PFJet500();
  });

  const NamedFunc high_pt_jet("high_pt_jet", [](const Baby &b) -> NamedFunc::ScalarType{
		  bool r_high_pt_jet = false;
		  if (b.njet()>0) {
		    if (b.jet_pt()->at(0) > 500) {
			    r_high_pt_jet = true;
		    }
		  }
		  return r_high_pt_jet;
  });

  const NamedFunc hig_cand_am_safe("hig_cand_am_safe", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_hig_cand_am_safe = -999;
		  if (b.hig_cand_am()->size() > 0) {
		    r_hig_cand_am_safe = b.hig_cand_am()->at(0);
		  }
		  return r_hig_cand_am_safe;
  });

  const NamedFunc hig_cand_drmax_safe("hig_cand_drmax_safe", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_hig_cand_drmax_safe = -999;
		  if (b.hig_cand_drmax()->size() > 0) {
		    r_hig_cand_drmax_safe = b.hig_cand_drmax()->at(0);
		  }
		  return r_hig_cand_drmax_safe;
  });

  const NamedFunc el_max_pt("el_max_pt", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_el_max_pt = -1;
		  for (unsigned el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
		    if (b.el_sig()->at(el_idx) && b.el_pt()->at(el_idx)>r_el_max_pt) {
		      r_el_max_pt = b.el_pt()->at(el_idx);
		    }
		  }
		  return r_el_max_pt;
  });

  const NamedFunc mu_max_pt("mu_max_pt", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_mu_max_pt = -1;
		  for (unsigned mu_idx = 0; mu_idx < b.mu_pt()->size(); mu_idx++) {
		    if (b.mu_sig()->at(mu_idx) && b.mu_pt()->at(mu_idx)>r_mu_max_pt) {
		      r_mu_max_pt = b.mu_pt()->at(mu_idx);
		    }
		  }
		  return r_mu_max_pt;
  });

  const NamedFunc nb_higgsino("nb_higgsino", [](const Baby &b) -> NamedFunc::ScalarType{
        	  int r_nb_higgsino;
        	  if (b.nbm() == 0) r_nb_higgsino = 0;
        	  else if (b.nbm() == 1) r_nb_higgsino = 1;
        	  else if (b.nbt() == 2 && b.nbm() == 2) r_nb_higgsino = 2;
        	  else if (b.nbt() >= 2 && b.nbm() == 3 && b.nbl() == 3) r_nb_higgsino = 3;
        	  else if (b.nbt() >= 2 && b.nbm() >= 3 && b.nbl() >= 4) r_nb_higgsino = 4;
		  else r_nb_higgsino = -1;
		  return r_nb_higgsino;
  });

  const NamedFunc nhig_cand("nhig_cand", [](const Baby &b) -> NamedFunc::ScalarType{
		  int r_nhig_cand = b.hig_cand_dm()->size();
		  return r_nhig_cand;
  });

  const NamedFunc trigger_efficiency("trigger_efficiency",[](const Baby &b) -> NamedFunc::ScalarType{
    float trig_eff = 1;
    if (b.SampleType()==2016) {
      trig_eff = Higfuncs::get_0l_trigeff2016.GetScalar(b);
    }
    else if (b.SampleType()==2017) {
      trig_eff = Higfuncs::get_0l_trigeff2017.GetScalar(b);
    }
    else if (b.SampleType()==2018) {
      trig_eff = Higfuncs::get_0l_trigeff2018temp.GetScalar(b);
    }
    return trig_eff;
  });

  PlotMaker pm;
  std::vector<double> true_met_bins{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350};
  std::vector<double> sys_met_bins{150,160,180,200,225,250,300,350};
  std::vector<double> fake_met_bins{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350,400,450,500,550};
  std::vector<double> ht_bins{0,200,600,800,1000,1200};
  std::vector<double> twodim_met_bins{0,110,120,130,140,150,160,170,180,190,200,210,250};
  std::vector<double> el_pt_bins{20,25,30,110,120,150};
  std::vector<double> mu_pt_bins{20,25,30,50,100};
  std::vector<double> twoel_pt_bins{40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,120};
  std::vector<double> twomu_pt_bins{40,45,50,60};

  NamedFunc baseline_nlep0 = "stitch&&pass&&nlep==0&&njet>=3";
  NamedFunc baseline_triggerregion = "stitch&&pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0";
  NamedFunc baseline_signalregion = "stitch&&pass&&nlep==0&&njet>=4&&njet<=5&&!low_dphi_met&&nel==1&&nmu==0";

  pm.Push<Hist2D>(Axis(50, 0., 300., "met", "MET [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               "stitch&&pass&&nlep==0&&njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_met_ambb").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 300., "met", "MET [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_met_ambb_baseline").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 800., "ht", "HT [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               "stitch&&pass&&nlep==0&&njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_ht_ambb").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 800., "ht", "HT [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_ht_ambb_baseline").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 300., "met", "MET [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               "stitch&&pass&&300<=ht&&ht<400&&nlep==0&&njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_met_ambb_ht300400").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 300., "met", "MET [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "njet>=4&&300<=ht&&ht<400", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_met_ambb_baseline_ht300400").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 800., "ht", "HT [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               "stitch&&pass&&300<=ht&&ht<400&&nlep==0&&njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_ht_ambb_ht300400").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 800., "ht", "HT [GeV]", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "njet>=4&&300<=ht&&ht<400", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_ht_ambb_baseline_ht300400").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 1., trigger_efficiency, "Trigger Efficiency", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "met>150&&njet>=4&&300<=ht&&ht<400", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_trigeff_ambb_baseline_ht300400").LuminosityTag("60 fb^{-1}");
  pm.Push<Hist2D>(Axis(50, 0., 1., trigger_efficiency, "Trigger Efficiency", {}),
  	               Axis(20, 0., 300., "hig_cand_am[0]", "<m> [GeV]", {}),
  	               baseline_triggerregion && "met>200&&njet>=4", procs_mc_all, twodim_plotopts).Weight("weight").Tag("FixName:higcandambbstudies_trigeff_ambb_baseline").LuminosityTag("60 fb^{-1}");

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  TFile* out_file = TFile::Open("ntuples/higcandamstudies.root","RECREATE");
  TH2D hist_met_ambb  = static_cast<Hist2D*>(pm.Figures()[0].get())->GetBkgHist(true);
  hist_met_ambb.SetName("hist_met_ambb");
  hist_met_ambb.Write();
  TH2D hist_met_ambb_baseline  = static_cast<Hist2D*>(pm.Figures()[1].get())->GetBkgHist(true);
  hist_met_ambb_baseline.SetName("hist_met_ambb_baseline");
  hist_met_ambb_baseline.Write();
  TH2D hist_ht_ambb  = static_cast<Hist2D*>(pm.Figures()[2].get())->GetBkgHist(true);
  hist_ht_ambb.SetName("hist_ht_ambb");
  hist_ht_ambb.Write();
  TH2D hist_ht_ambb_baseline  = static_cast<Hist2D*>(pm.Figures()[3].get())->GetBkgHist(true);
  hist_ht_ambb_baseline.SetName("hist_ht_ambb_baseline");
  hist_ht_ambb_baseline.Write();
  TH2D hist_met_ambb_ht300400  = static_cast<Hist2D*>(pm.Figures()[4].get())->GetBkgHist(true);
  hist_met_ambb_ht300400.SetName("hist_met_ambb_ht300400");
  hist_met_ambb_ht300400.Write();
  TH2D hist_met_ambb_baseline_ht300400  = static_cast<Hist2D*>(pm.Figures()[5].get())->GetBkgHist(true);
  hist_met_ambb_baseline_ht300400.SetName("hist_met_ambb_baseline_ht300400");
  hist_met_ambb_baseline_ht300400.Write();
  TH2D hist_ht_ambb_ht300400  = static_cast<Hist2D*>(pm.Figures()[6].get())->GetBkgHist(true);
  hist_ht_ambb_ht300400.SetName("hist_ht_ambb_ht300400");
  hist_ht_ambb_ht300400.Write();
  TH2D hist_ht_ambb_baseline_ht300400  = static_cast<Hist2D*>(pm.Figures()[7].get())->GetBkgHist(true);
  hist_ht_ambb_baseline_ht300400.SetName("hist_ht_ambb_baseline_ht300400");
  hist_ht_ambb_baseline_ht300400.Write();
  TH2D hist_trigeff_ambb_baseline_ht300400  = static_cast<Hist2D*>(pm.Figures()[8].get())->GetBkgHist(true);
  hist_trigeff_ambb_baseline_ht300400.SetName("hist_trigeff_ambb_baseline_ht300400");
  hist_trigeff_ambb_baseline_ht300400.Write();
  TH2D hist_trigeff_ambb_baseline  = static_cast<Hist2D*>(pm.Figures()[9].get())->GetBkgHist(true);
  hist_trigeff_ambb_baseline.SetName("hist_trigeff_ambb_baseline");
  hist_trigeff_ambb_baseline.Write();
  out_file->Close();

  time(&endtime);
  cout<<endl<<"Processing took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
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
