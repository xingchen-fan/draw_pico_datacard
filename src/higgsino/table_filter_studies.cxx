#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TRandom.h"
#include "TError.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TVector2.h"
#include "TLorentzVector.h"

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
#include "higgsino/hig_functions.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"

using namespace std;
using namespace Higfuncs;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
  bool do_metsync = false;
  bool do_mc_efficiency = false;
  bool do_data_efficiency = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  float lumi = 137.;
  Palette colors("txt/colors.txt", "default");
  // Define 1D+2D plot types of interest
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow);
  vector<PlotOpt> all_plot_types = {lin_lumi};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).Overflow(OverflowType::overflow)};

  string production_dir = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  shared_ptr<Process> pro_met2016c = Process::MakeShared<Baby_pico>("MET\\_2016RunC", Process::Type::data, kBlack, {production_dir+"2016/data/raw_pico/raw_pico_MET__Run2016C*.root"}, "1");
  shared_ptr<Process> pro_met2018d = Process::MakeShared<Baby_pico>("MET\\_2018RunD", Process::Type::data, kBlack, {production_dir+"2018/data/raw_pico/raw_pico_MET__Run2018D*.root"}, "1");
  shared_ptr<Process> pro_egamma2018d = Process::MakeShared<Baby_pico>("EGamma\\_2018RunD", Process::Type::data, kBlack, {production_dir+"2018/data/skim_met150/raw_pico_met150_EGamma__Run2018D*.root"}, "1");
  vector<shared_ptr<Process> > procs_metsync = {pro_met2016c};

  string production2017_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/";
  shared_ptr<Process> pro_met2017c = Process::MakeShared<Baby_pico>("MET\\_2017RunC", Process::Type::data, kBlack, {production_dir+"2017/data/raw_pico/raw_pico_MET__Run2017C*.root"}, "1");
  vector<shared_ptr<Process> > procs_met2017c = {pro_met2017c};

  // Set MC 
  vector<shared_ptr<Process> > mc_procs;
  map<string, set<string>> mctags; 
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  string mc_skim_folder = "mc/skim_met150/";
  std::set<int> years = {2016};
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_SingleLept*","TTJets_DiLept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  // Set mc processes
  mc_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  mc_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  mc_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  mc_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  mc_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));
  //add signal pts
  vector<string> mchi_sigm = {"175","500","900"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  vector<string> mchi_sigm2d = {"250","350","450"}; 
  vector<string> mlsp_sigm2d = {"50","200","100"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue}; // need sigm.size() <= sig_colors.size()
  vector<int> sig_colors2d = {kOrange, kYellow, kCyan};
  string foldersig = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/skim_met150/";
  for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
    mc_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
    std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
  }
  for (unsigned isig(1); isig<mchi_sigm2d.size(); isig++) {
    mc_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
    std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
  }

  const NamedFunc mixed_model_weight("mixed_model_weight",[](const Baby &b) -> NamedFunc::ScalarType{
    set<int> mchi_sigm_int = {175, 500, 900}; 
    double trig_eff = 1.0;
    double lumi_scale = 1.0;
    if (b.SampleType()==2016) {
      trig_eff = Higfuncs::get_0l_trigeff2016.GetScalar(b);
      lumi_scale = 35.9;
    }
    else if (b.SampleType()==2017) {
      trig_eff = Higfuncs::get_0l_trigeff2017.GetScalar(b);
      lumi_scale = 41.5;
    }
    else if (b.SampleType()==2018) {
      trig_eff = Higfuncs::get_0l_trigeff2018.GetScalar(b);
      lumi_scale = 60.0;
    }
    if (b.mprod() == -999) {
      //not signal
      return b.weight()*trig_eff*lumi_scale;
    }
    if (mchi_sigm_int.count(b.mprod()) > 0) {
      //1d signal
      return b.weight()*trig_eff*lumi_scale;
    }
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return b.weight()/xsec1d*xsec2d*trig_eff*lumi_scale;
  });
  
  const NamedFunc boosted_mht_phi("boosted_mht_phi", [](const Baby &b) -> NamedFunc::ScalarType{
    //recalculate MHT using only jets with |eta|<5
    TLorentzVector mht_vec;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && fabs(b.jet_eta()->at(jet_idx)) < 5.) {
        TLorentzVector ijet_v4;
        ijet_v4.SetPtEtaPhiM(b.jet_pt()->at(jet_idx), b.jet_eta()->at(jet_idx), b.jet_phi()->at(jet_idx), b.jet_m()->at(jet_idx));
	mht_vec -= ijet_v4;
      }
    }
    return mht_vec.Phi();
  });

  const NamedFunc old_pass_low_neutral_jet("old_pass_low_neutral_jet", [](const Baby &b) -> NamedFunc::ScalarType{
    //emulate pass_low_neutral_jet used in humboldt production
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx)<=30 || fabs(b.jet_eta()->at(jet_idx))>2.4) continue;
      float dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
      if (b.jet_ne_emef()->at(jet_idx) <0.03 && dphi>(TMath::Pi()-0.4))
        return false;
      break; //only apply to leading jet
    }
    return true;
  });

  const NamedFunc pass_htratio_dphi_tight_fixed("pass_htratio_dphi_tight_fixed", [](const Baby &b) -> NamedFunc::ScalarType{
    bool r_pass_htratio_dphi_tight_fixed = true;
    float ht_ratio = b.ht5()/b.ht();
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) < 30 || fabs(b.jet_eta()->at(jet_idx)) > 2.4) continue;
      if (ht_ratio > 1.2 && fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi())) < (5.3*ht_ratio - 4.78)) 
        r_pass_htratio_dphi_tight_fixed = false;
      break;
    }
    return r_pass_htratio_dphi_tight_fixed;
  });

  const NamedFunc htratio("htratio", [](const Baby &b) -> NamedFunc::ScalarType{
    float r_ht_ratio = b.ht5()/b.ht();
    return r_ht_ratio;
  });

  const NamedFunc pass_boostedhemveto("pass_boostedhemveto", [](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned int el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
      if (b.el_miniso()->at(el_idx) < 0.1 && -3.0 < b.el_eta()->at(el_idx) && b.el_eta()->at(el_idx) < -1.4 && -1.57 < b.el_phi()->at(el_idx) && b.el_phi()->at(el_idx) < -0.87) {
        return false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && -3.2 < b.jet_eta()->at(jet_idx) && b.jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.jet_phi()->at(jet_idx) && b.jet_phi()->at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
	if (dphi < 0.5) {
	  return false;
	}
      }
    }
    return true;
  });

  NamedFunc baseline = "stitch&&150<=met&&nvlep==0&&ntk==0&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc oneel_baseline = "stitch&&150<=met&&nlep==1&&nel==1&&ntk==0&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";

  PlotMaker pm;

  pm.Push<Hist2D>(Axis(100, 0, 3.14, "jet_met_dphi", "#Delta#phi", {}), Axis(100, 0, 500, "jet_pt", "Jet p_{T} [GeV]", {}),
    "2.4<=jet_eta&&jet_eta<5.0",
    procs_met2017c, twodim_plotopts).Weight(mixed_model_weight).Tag("FixName:EcalNoiseJetFilter_plot").LuminosityTag("137 fb^{-1}");
  pm.Push<Hist2D>(Axis(100, 0, 3.14, "jet_met_dphi[0]", "Jet 1 p_{T} [GeV]", {}), Axis(100, 0, 3.0, htratio, "H_{T5}/H_{T}", {}), 
    "njet>=1",
    procs_met2017c, twodim_plotopts).Weight(mixed_model_weight).Tag("FixName:HTRatioDPhiTightFilter_plot").LuminosityTag("137 fb^{-1}");

  if (do_data_efficiency) {
  ////no cut table for 1l region
  //pm.Push<Table>("cutflow", vector<TableRow>{
  //TableRow("MET>150", 
  //  "1",0,0, mixed_model_weight),
  //TableRow("Trigger HLT_Ele35_WPTight_Gsf", 
  //  "HLT_Ele35_WPTight_Gsf",0,0, mixed_model_weight),
  //TableRow("N_{l}=1", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1",0,0, mixed_model_weight),
  //TableRow("High \\Delta \\phi", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met",0,0, mixed_model_weight),
  //TableRow("$N_{pv good}>0$", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0",0,0, mixed_model_weight),
  //TableRow("globalSuperTightHalo2016Filter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
  //TableRow("HBHENoiseFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
  //TableRow("HBHENoiseIsoFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
  //TableRow("EcalDeadCellTriggerPrimitiveFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
  //TableRow("BadPFMuonFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
  //TableRow("eeBadScFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
  //TableRow("EcalNoiseJetFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
  //TableRow("MuonJetFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
  //TableRow("$p_{T}^{miss}/p_{Tcalo}^{miss}<2.0$", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
  //TableRow("$p_{T}^{miss}/H_{T}^{miss}<2.0$", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
  //TableRow("LowNeutralJetFilter (Usynced)", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
  //TableRow("HTRatioDPhiTightFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed ,0,0, mixed_model_weight),
  //TableRow("HEMDPhiVetoFilter", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, mixed_model_weight),
  //TableRow("JetID", 
  //  "HLT_Ele35_WPTight_Gsf&&nlep==1&&!low_dphi_met&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, mixed_model_weight),
  }
  
  if (do_mc_efficiency) {
    //MET datasets for sync with boosted
    pm.Push<Table>("cutflow", vector<TableRow>{
    TableRow("Baseline", 
      baseline,0,0, mixed_model_weight),
    TableRow("$N_{pv good}>0$", 
      baseline && "npv_good>0",0,0, mixed_model_weight),
    TableRow("globalSuperTightHalo2016Filter", 
      baseline && "npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
    TableRow("HBHENoiseFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
    TableRow("HBHENoiseIsoFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
    TableRow("BadPFMuonFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
    TableRow("eeBadScFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
    TableRow("EcalNoiseJetFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("MuonJetFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_{T}^{miss}/p_{Tcalo}^{miss}<2.0$", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_{T}^{miss}/H_{T}^{miss}<2.0$", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("LowNeutralJetFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HTRatioDPhiTightFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed ,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, mixed_model_weight),
    TableRow("JetID", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_jets" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, mixed_model_weight),
    },mc_procs,false,true,false,false,true,false);
  }

  if (do_metsync) {
    //MET datasets for sync with boosted
    pm.Push<Table>("cutflow", vector<TableRow>{
    TableRow("Trigger HLT\\_PFMET120\\_PFMHT120\\_IDTight", 
      "HLT_PFMET120_PFMHT120_IDTight",0,0, "1"),
    TableRow("$N_{pv good}>0$", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0",0,0, "1"),
    TableRow("globalSuperTightHalo2016Filter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight",0,0, "1"),
    TableRow("HBHENoiseFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, "1"),
    TableRow("HBHENoiseIsoFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, "1"),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, "1"),
    TableRow("BadPFMuonFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, "1"),
    TableRow("eeBadScFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, "1"),
    TableRow("EcalNoiseJetFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, "1"),
    TableRow("MuonJetFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, "1"),
    TableRow("$p_{T}^{miss}/p_{Tcalo}^{miss}<2.0$", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, "1"),
    TableRow("$p_{T}^{miss}/H_{T}^{miss}<2.0$", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, "1"),
    TableRow("LowNeutralJetFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && old_pass_low_neutral_jet && pass_ecalnoisejet,0,0, "1"),
    TableRow("HTRatioDPhiTightFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && old_pass_low_neutral_jet && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed ,0,0, "1"),
    TableRow("HEMDPhiVetoFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && old_pass_low_neutral_jet && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, "1"),
    TableRow("JetID", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_jets" && old_pass_low_neutral_jet && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed && pass_boostedhemveto,0,0, "1"),
    },procs_metsync,false,true,false,false,true,false);
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Processing took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"do_metsync", no_argument, 0, 0},
      {"do_mceff", no_argument, 0, 0},
      {"do_dataeff", no_argument, 0, 0},
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
      if(optname == "do_metsync"){
        do_metsync = true;
      } else if(optname == "do_mceff"){
        do_mc_efficiency = true;
      } else if(optname == "do_dateff"){
        do_data_efficiency = true;
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
