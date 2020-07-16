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

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

//helper function declaration
void generate_ratio_plots(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, const char* hist_title, string hist_name);
void generate_2d_efficiencies(PlotMaker* pm, int denominator_index, int numerator_index, const char* hist_title, string hist_name);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  //flags for which things to plot
  bool do_variables = true;
  bool do_systematics = false;
  bool do_efficiency = false;
  int year = 2016;
  bool do_controlregions = false; //for do_variables, false will omit plots of 1l and 2l CRs. This speeds up processing by about 10~20x

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
  }
  procs_data.push_back(pro_data);
  procs_mc.push_back(pro_mc);
  procs_data_all.push_back(pro_data_all);

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

  //plots of efficiency with respect to various analysis variables
  if (do_variables) {
	//eff vs MET (real MET)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && met_trigger, procs_data, all_plot_types);
	//eff vs MET (QCD MET)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&HLT_PFJet500&&low_dphi_met&&nvlep==0", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&HLT_PFJet500&&low_dphi_met&&nvlep==0" && met_trigger, procs_data, all_plot_types);
	//eff vs HT (low MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200" && met_trigger, procs_data, all_plot_types);
	//eff vs HT (high MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300" && met_trigger, procs_data, all_plot_types);
	//eff vs nb
  	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {-0.5,4.5}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {-0.5,4.5}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met" && met_trigger, procs_data, all_plot_types);
	//eff vs njet
  	pm.Push<Hist1D>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {2.5,7.5}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {2.5,7.5}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met" && met_trigger, procs_data, all_plot_types);
	//eff vs <m> higgs
	//temp: use only met(nomu)120 trigger
  	pm.Push<Hist1D>(Axis(20, 0., 250., "hig_cand_am[0]", "#LT m#RT [GeV]", {0.,250.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&200<met" && nhig_cand>0., procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(20, 0., 250., "hig_cand_am[0]", "#LT m#RT [GeV]", {0.,250.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&200<met&&(HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight)" && nhig_cand>0., procs_data, all_plot_types);
	//eff vs drmax
  	pm.Push<Hist1D>(Axis(20, 0., 4., "hig_cand_drmax[0]", "#Delta R_{max}", {0.,4.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met" && nhig_cand>0., procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(20, 0., 4., "hig_cand_drmax[0]", "#Delta R_{max}", {0.,4.}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&200<met" && nhig_cand>0. && met_trigger, procs_data, all_plot_types);
	if (do_controlregions) {
		//eff vs MET (1e region)
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nel==1", procs_data_all, all_plot_types);
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nel==1" && (met_trigger||el_trigger), procs_data_all, all_plot_types);
		//eff vs el_pt (1el region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., el_max_pt, "Offline Electron p_{T} [GeV]", {20.,170.}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nel==1&&150<met&&met<200", procs_data, all_plot_types);
  		pm.Push<Hist1D>(Axis(30, 20., 170., el_max_pt, "Offline Electron p_{T} [GeV]", {20.,170.}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nel==1&&150<met&&met<200" && (met_trigger||el_trigger), procs_data, all_plot_types);
		//eff vs MET (1mu region)
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nmu==1", procs_data_all, all_plot_types);
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nmu==1" && (met_trigger||mu_trigger), procs_data_all, all_plot_types);
		//eff vs mu_pt (1mu region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., mu_max_pt, "Offline Muon p_{T} [GeV]", {20.,170.}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nmu==1&&150<met&&met<200", procs_data, all_plot_types);
  		pm.Push<Hist1D>(Axis(30, 20., 170., mu_max_pt, "Offline Muon p_{T} [GeV]", {20.,170.}),
  		                "pass&&HLT_PFJet500&&njet>=2&&nmu==1&&150<met&&met<200" && (met_trigger||mu_trigger), procs_data, all_plot_types);
		//eff vs el_pt (2el region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., el_max_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                "pass&&njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger), procs_data_all, all_plot_types);
  		pm.Push<Hist1D>(Axis(30, 20., 170., el_max_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                "pass&&njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger) && el_trigger, procs_data_all, all_plot_types);
		//eff vs mu_pt (2mu region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., mu_max_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                "pass&&njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger), procs_data_all, all_plot_types);
  		pm.Push<Hist1D>(Axis(30, 20., 170., mu_max_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                "pass&&njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger) && mu_trigger, procs_data_all, all_plot_types);
	}
	//MET120 only (nominal)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0" && el_max_pt > 30, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0&&(HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight)" && el_max_pt > 30, procs_data, all_plot_types);
	//MET120 MC
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0)*weight" && el_max_pt > 30, procs_mc, all_plot_types);
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0&&(HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight))*weight" && el_max_pt>30, procs_mc, all_plot_types);
	//MET120 MC vs <m>
  	pm.Push<Hist1D>(Axis(20, 0., 250., "hig_cand_am[0]", "<m> [GeV]", {0,250.}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&met>=200&&!low_dphi_met&&nel==1&&nmu==0)*weight" && nhig_cand>0., procs_mc, all_plot_types);
  	pm.Push<Hist1D>(Axis(20, 0., 250., "hig_cand_am[0]", "<m> [GeV]", {0,250.}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&met>=200&&!low_dphi_met&&nel==1&&nmu==0&&(HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight))*weight" && nhig_cand>0., procs_mc, all_plot_types);
  }
  //0l systematics plots
  if (do_systematics) {
	//nj3 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && met_trigger, procs_data, all_plot_types);
	//nj4
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet==4&&!low_dphi_met&&nel==1", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet==4&&!low_dphi_met&&nel==1" && met_trigger, procs_data, all_plot_types);
	//nj5
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet==5&&!low_dphi_met&&nel==1", procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet==5&&!low_dphi_met&&nel==1" && met_trigger, procs_data, all_plot_types);
	//nb0
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino==0., procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino==0. && met_trigger, procs_data, all_plot_types);
	//nb1
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino==1., procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino==1. && met_trigger, procs_data, all_plot_types);
	//nb2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino>=2., procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nb_higgsino>=2. && met_trigger, procs_data, all_plot_types);
	//am<140
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_am_safe<=140, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_am_safe<=140 && met_trigger, procs_data, all_plot_types);
	//am>140
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_am_safe>140, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_am_safe>140 && met_trigger, procs_data, all_plot_types);
	//drmax < 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_drmax_safe<=2.2, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_drmax_safe<=2.2 && met_trigger, procs_data, all_plot_types);
	//drmax > 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_drmax_safe>2.2, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && nhig_cand>0. && hig_cand_drmax_safe>2.2 && met_trigger, procs_data, all_plot_types);
	//jet_pt>500 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet && met_trigger, procs_data, all_plot_types);
	//jet_pt>500, ref trigger jet500
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet && met_trigger, procs_data, all_plot_types);
	//MET120 only (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0" && el_max_pt > 30, procs_data, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "pass&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0&&(HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)" && el_max_pt > 30, procs_data, all_plot_types);
	//MET120 MC (AN cuts)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "(stitch&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&pass&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0)*weight" && el_max_pt > 30, procs_mc, all_plot_types);
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                "(stitch&&(HLT_Ele35_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&pass&&njet>=3&&!low_dphi_met&&nel==1&&(met/met_calo<5)&&nmu==0&&(HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60))*weight" && el_max_pt > 30, procs_mc, all_plot_types);
  }
  // trigger efficiency plots
  if (do_efficiency) {
  	pm.Push<Hist2D>(Axis(true_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1", procs_data, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(true_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                "pass&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1" && met_trigger, procs_data, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                "pass&&HLT_PFJet500&&low_dphi_met&&nvlep==0", procs_data, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                "pass&&HLT_PFJet500&&low_dphi_met&&nvlep==0" && met_trigger, procs_data, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  	                Axis(el_pt_bins, "el_pt[0]", "Electron pt", {}),
  	                "pass&&HLT_PFJet500&&njet>=2&&nel==1", procs_data_all, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  	                Axis(el_pt_bins, "el_pt[0]", "Electron pt", {}),
  	                "pass&&HLT_PFJet500&&njet>=2&&nel==1" && (met_trigger||el_trigger), procs_data_all, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  	                Axis(mu_pt_bins, "mu_pt[0]", "Muon pt", {}),
  	                "pass&&HLT_PFJet500&&njet>=2&&nmu==1", procs_data_all, twodim_plotopts);
  	pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  	                Axis(mu_pt_bins, "mu_pt[0]", "Muon pt", {}),
  	                "pass&&HLT_PFJet500&&njet>=2&&nmu==1" && (met_trigger||mu_trigger), procs_data_all, twodim_plotopts);
  	pm.Push<Hist1D>(Axis(twoel_pt_bins, el_max_pt, "Offline Max Electron p_{T} [GeV]", {}),
  	                "pass&&njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger), procs_data_all, all_plot_types);
  	pm.Push<Hist1D>(Axis(twoel_pt_bins, el_max_pt, "Offline Max Electron p_{T} [GeV]", {}),
  	                "pass&&njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger) && el_trigger, procs_data_all, all_plot_types);
  	pm.Push<Hist1D>(Axis(twomu_pt_bins, mu_max_pt, "Offline Max Muon p_{T} [GeV]", {}),
  	                "pass&&njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger), procs_data_all, all_plot_types);
  	pm.Push<Hist1D>(Axis(twomu_pt_bins, mu_max_pt, "Offline Max Muon p_{T} [GeV]", {}),
  	                "pass&&njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jet_trigger) && mu_trigger, procs_data_all, all_plot_types);
  }

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.print_2d_figures_ = false;
  pm.MakePlots(lumi);

  //save plots to root file
  TFile* out_file = TFile::Open("ntuples/triggereff.root","RECREATE");
  int pm_idx = 0;
  if (do_variables) {
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET#leq 200, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_htlowmet");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_hthighmet");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline N_{b}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline N_{j}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandam");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline #Delta R_{max}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcanddrmax");
	pm_idx += 2;
	if (do_controlregions) {
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elmet");
		pm_idx += 2;
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline Electron Pt [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elpt");
		pm_idx += 2;
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_mumet");
		pm_idx += 2;
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline Muon Pt [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_mupt");
		pm_idx += 2;
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{e}=2 80<m_{ee}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Electron Pt [GeV]; Efficiency [Ele(27_WPTight||35_WPTight||115)]","hist_elel_show");
		pm_idx += 2;
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{#mu}=2 80<m_{#mu#mu}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Muon Pt [GeV]; Efficiency [IsoMu(24||27)||Mu50]","hist_mumu_show");
		pm_idx += 2;
	}
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_datamet120");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_mc, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, Summer16 t#bar{t} 1l; Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_mcmet120");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_mc, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, Summer16 t#bar{t} 1l; Offline #LT m#RT [GeV]; Efficiency [MET120]","hist_mcmet120higcandam");
	pm_idx += 2;
  }
  if (do_systematics) {
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj3");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}= 4 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj4");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}= 5 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj5");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}=0 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb0");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}=1 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb1");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}#geq 1 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb2");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #LT m#RT #leq 140 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_amlow");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #LT m#RT > 140 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_amhigh");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #Delta R_{max}#leq 2.2 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_drmaxlow");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #Delta R_{max}>2.2 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_drmaxhigh");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 p_{Tjet}>500 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj3jet500");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_jetjet500");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_met120");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_mc, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_mc");
	pm_idx += 2;
  }
  if (do_efficiency) {
  	generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_realmetht");
	pm_idx += 2;
  	generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemetht");
	pm_idx += 2;
  	generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||Ele(27_WPTight||35_WPTight||115)] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], Electron pt [GeV]", "hist_metelpt");
	pm_idx += 2;
  	generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||IsoMu(24||27)||Mu50] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], Muon pt [GeV]", "hist_metmupt");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{e}=2 80<m_{ee}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Electron Pt [GeV]; Efficiency [Ele(27_WPTight||35_WPTight||115)]","hist_elel");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_data_all, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{#mu}=2 80<m_{#mu#mu}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Muon Pt [GeV]; Efficiency [IsoMu(24||27)||Mu50]","hist_mumu");
	pm_idx += 2;
  }
  out_file->Close();

  // vector<GammaParams> yields = cutflow.BackgroundYield(lumi);
  // for(const auto &yield: yields){
  //   cout << yield << endl;
  // }

  time(&endtime);
  cout<<endl<<"Processing took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void generate_ratio_plots(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, const char* hist_title, string hist_name) {
  Hist1D * den_h1d = static_cast<Hist1D*>(pm->Figures()[denominator_index].get());
  Hist1D::SingleHist1D* den_sh = static_cast<Hist1D::SingleHist1D*>(den_h1d->GetComponent(proc.get()));
  TH1D den_h = den_sh->scaled_hist_;
  Hist1D * num_h1d = static_cast<Hist1D*>(pm->Figures()[numerator_index].get());
  Hist1D::SingleHist1D* num_sh = static_cast<Hist1D::SingleHist1D*>(num_h1d->GetComponent(proc.get()));
  TH1D num_h = num_sh->scaled_hist_;
  TGraphAsymmErrors* ratio_h = new TGraphAsymmErrors(&num_h,&den_h,"cp");
  ratio_h->SetTitle(hist_title);
  ratio_h->GetYaxis()->SetRangeUser(0,1.4);
  ratio_h->SetMarkerStyle(kFullCircle);
  ratio_h->SetMarkerSize(1);
  ratio_h->GetXaxis()->SetTitleSize(0.04);
  ratio_h->GetYaxis()->SetTitleSize(0.04);
  ratio_h->GetXaxis()->SetTitleOffset(1.0);
  ratio_h->GetYaxis()->SetTitleOffset(1.0);
  den_h.Write((hist_name+"_denominator").c_str());
  num_h.Write((hist_name+"_numerator").c_str());
  ratio_h->Write((hist_name+"_ratio").c_str());
}

void generate_2d_efficiencies(PlotMaker* pm, int denominator_index, int numerator_index, const char* hist_title, string hist_name) {
  Hist2D * effhist_den = static_cast<Hist2D*>(pm->Figures()[denominator_index].get());
  TH2D eff_den_h = effhist_den->GetDataHist();
  Hist2D * effhist_num = static_cast<Hist2D*>(pm->Figures()[numerator_index].get());
  TH2D eff_num_h = effhist_num->GetDataHist();
  TH2D *eff_ratio_h = static_cast<TH2D*>(eff_num_h.Clone());
  eff_ratio_h->Divide(&eff_den_h);
  eff_ratio_h->SetTitle(hist_title);
  eff_den_h.Write((hist_name+"_denominator").c_str());
  eff_num_h.Write((hist_name+"_numerator").c_str());
  eff_ratio_h->Write((hist_name+"_ratio").c_str());
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
