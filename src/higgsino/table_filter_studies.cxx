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
  bool do_ele_efficiency = false;
  bool do_mu_efficiency = false;
  std::string year_string = "2016";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  float lumi = 137.;
  if (year_string == "2016") {
    lumi = 35.9;
  }
  else if (year_string == "2017") {
    lumi = 41.5;
  }
  else if (year_string == "2018") {
    lumi = 60.0;
  }
  Palette colors("txt/colors.txt", "default");
  // Define 1D+2D plot types of interest
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow);
  vector<PlotOpt> all_plot_types = {lin_lumi};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).YAxis(YAxisType::log).Overflow(OverflowType::overflow)};

  string production_dir = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string old_production_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  shared_ptr<Process> pro_met2016c = Process::MakeShared<Baby_pico>("MET\\_2016RunC", Process::Type::data, kBlack, {old_production_dir+"2016/data/raw_pico/raw_pico_MET__Run2016C*.root"}, "1");
  shared_ptr<Process> pro_met2018d = Process::MakeShared<Baby_pico>("MET\\_2018RunD", Process::Type::data, kBlack, {production_dir+"2018/data/raw_pico/raw_pico_MET__Run2018D*.root"}, "1");
  shared_ptr<Process> pro_egamma2018d = Process::MakeShared<Baby_pico>("EGamma\\_2018RunD", Process::Type::data, kBlack, {production_dir+"2018/data/skim_met150/raw_pico_met150_EGamma__Run2018D*.root"}, "1");
  vector<shared_ptr<Process> > procs_metsync = {pro_met2016c};

  string onelep_dir = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/"+year_string+"/data/skim_higlep1T/raw_pico_higlep1T_";
  vector<shared_ptr<Process> > procs_onelep = {};
  if (year_string == "2018") {
    shared_ptr<Process> pro_onelep = Process::MakeShared<Baby_pico>("1l data", Process::Type::data, kBlack, {onelep_dir+"MET*.root", onelep_dir+"EGamma*.root",onelep_dir+"SingleMuon*.root"}, "1");
    procs_onelep.push_back(pro_onelep);
  }
  else {
    shared_ptr<Process> pro_onelep = Process::MakeShared<Baby_pico>("1l data", Process::Type::data, kBlack, {onelep_dir+"MET*.root", onelep_dir+"SingleElectron*.root",onelep_dir+"SingleMuon*.root"}, "1");
    procs_onelep.push_back(pro_onelep);
  }

  //use Proccess::Type::background for 2d plots so draw_pico doesn't crash trying to plots millions of points on a pdf
  string production2017_dir = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/";
  shared_ptr<Process> pro_met2017c = Process::MakeShared<Baby_pico>("MET\\_2017RunC", Process::Type::background, kBlack, {production_dir+"2017/data/raw_pico/raw_pico_MET__Run2017C*.root"}, "1");
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
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));
  //add signal pts
  vector<string> mchi_sigm = {"175","500","900"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  vector<string> mchi_sigm2d = {"250","350","450"}; 
  vector<string> mlsp_sigm2d = {"50","200","100"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue}; // need sigm.size() <= sig_colors.size()
  vector<int> sig_colors2d = {kOrange, kYellow, kCyan};
  string foldersig = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/"+year_string+"/SMS-TChiHH_2D/skim_met150/";
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
      trig_eff = Higfuncs::get_0l_trigeff2016.GetVector(b)[0];
      lumi_scale = 35.9;
    }
    else if (b.SampleType()==2017) {
      trig_eff = Higfuncs::get_0l_trigeff2017.GetVector(b)[0];
      lumi_scale = 41.5;
    }
    else if (b.SampleType()==2018) {
      trig_eff = Higfuncs::get_0l_trigeff2018.GetVector(b)[0];
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
    //only apply for relevant runs
    bool pass_hem = true;
    if (year_string != "2018") return 1.0;
    if (b.type()/1000 == 0 && b.run() < 319077) return 1.0;
    for (unsigned int el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
      if (b.el_miniso()->at(el_idx) < 0.1 && -3.0 < b.el_eta()->at(el_idx) && b.el_eta()->at(el_idx) < -1.4 && -1.57 < b.el_phi()->at(el_idx) && b.el_phi()->at(el_idx) < -0.87) {
        pass_hem = false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && -3.2 < b.jet_eta()->at(jet_idx) && b.jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.jet_phi()->at(jet_idx) && b.jet_phi()->at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
	if (dphi < 0.5) {
	  pass_hem = false;
	}
      }
    }
    if (b.type()/1000 == 0) {
      //data
      if (pass_hem) return 1.0;
      return 0.0;
    }
    //mc
    if (pass_hem) return 1.0;
    //weight for amount of 2018 luminosity after HEM issue occured
    return 0.645;
  });

  const NamedFunc el_trigger("el_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_el_trigger = b.HLT_Ele25_WPTight_Gsf()||b.HLT_Ele27_WPTight_Gsf()||b.HLT_Ele28_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG()||b.HLT_Ele32_WPTight_Gsf()||b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele45_WPLoose_Gsf()||b.HLT_Ele105_CaloIdVT_GsfTrkIdT()||b.HLT_Ele115_CaloIdVT_GsfTrkIdT()||b.HLT_Ele135_CaloIdVT_GsfTrkIdT()||b.HLT_Ele145_CaloIdVT_GsfTrkIdT()||b.HLT_Ele25_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPLoose_Gsf()||b.HLT_Ele20_WPLoose_Gsf()||b.HLT_Ele20_eta2p1_WPLoose_Gsf()||b.HLT_Ele25_eta2p1_WPLoose_Gsf()||b.HLT_Ele15_IsoVVVL_PFHT350()||b.HLT_Ele15_IsoVVVL_PFHT400()||b.HLT_Ele15_IsoVVVL_PFHT450()||b.HLT_Ele15_IsoVVVL_PFHT600()||b.HLT_Ele50_IsoVVVL_PFHT450();
    return r_el_trigger;
  });

  const NamedFunc mu_trigger("mu_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_mu_trigger = b.HLT_IsoMu20()||b.HLT_IsoMu22()||b.HLT_IsoMu24()||b.HLT_IsoMu27()||b.HLT_IsoTkMu20()||b.HLT_IsoTkMu22()||b.HLT_IsoTkMu24()||b.HLT_Mu50()||b.HLT_Mu55()||b.HLT_TkMu50()||b.HLT_IsoMu22_eta2p1()||b.HLT_IsoMu24_eta2p1()||b.HLT_Mu45_eta2p1()||b.HLT_Mu15_IsoVVVL_PFHT350()||b.HLT_Mu15_IsoVVVL_PFHT400()||b.HLT_Mu15_IsoVVVL_PFHT450()||b.HLT_Mu15_IsoVVVL_PFHT600()||b.HLT_Mu50_IsoVVVL_PFHT400()||b.HLT_Mu50_IsoVVVL_PFHT450();
    return r_mu_trigger;
  });

  const NamedFunc jetid_modified_cuts("jetid_modified_cuts", [](const Baby &b) -> NamedFunc::ScalarType{
    //implements njet, nbt, and higcand cuts when jets failing Id are rejected
    unsigned int njet = 0, nbt = 0;
    //if no bad jets, just return normal cuts
    if (b.pass_jets()) {
      if (b.njet() < 4) return false;
      return (b.njet()>=4)&&(b.njet()<=5)&&(b.nbt()>=2)&&(b.hig_cand_dm()->at(0)<=40)&&(b.hig_cand_drmax()->at(0)<2.2)&&(b.hig_cand_am()->at(0)<200)&&(!b.low_dphi_met());
    }
    //otherwise, we have to recalculate higgs candidates
    std::vector<std::pair<unsigned int, float>> jet_idx_and_btagdisc;
    std::vector<std::pair<float, float>> jet_pt_and_dphi;
    bool recalculate_higgs = false;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) {
        if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
        continue;
      }
      njet += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.8953 && year_string == "2016") nbt += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.8001 && year_string == "2017") nbt += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.7527 && year_string == "2018") nbt += 1;
      jet_idx_and_btagdisc.push_back(std::make_pair(jet_idx, b.jet_deepcsv()->at(jet_idx)));
      jet_pt_and_dphi.push_back(std::make_pair(b.jet_pt()->at(jet_idx), b.jet_met_dphi()->at(jet_idx)));
    }
    //njet, nbt, dphi cuts, also if bad jet was the non-higgs cand jet, just return original cuts
    if (njet < 4 || njet > 5 || nbt < 2) return false;
    std::sort(jet_pt_and_dphi.begin(), jet_pt_and_dphi.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.first > bb.first;
      });
    if (jet_pt_and_dphi[0].second < 0.5 || jet_pt_and_dphi[1].second < 0.5 || jet_pt_and_dphi[2].second < 0.3 || jet_pt_and_dphi[3].second < 0.3) return false;
    if (!recalculate_higgs) return (b.hig_cand_dm()->at(0)<=40)&&(b.hig_cand_drmax()->at(0)<2.2)&&(b.hig_cand_am()->at(0)<200);
    //get the lowest Delta m Higgs candidate pairing with the first two indices corresponding to the higher pt higgs candidate
    std::sort(jet_idx_and_btagdisc.begin(), jet_idx_and_btagdisc.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.second > bb.second;
      });
    std::vector<TLorentzVector> b_jet_p;
    for (unsigned int bjet_idx = 0; bjet_idx < 4; bjet_idx++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_eta()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_phi()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_m()->at(jet_idx_and_btagdisc[bjet_idx].first));
      b_jet_p.push_back(lv);
    }
    std::vector<std::vector<unsigned int>> higgs_combinations = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    float min_delta_m = 999;
    float min_dm_drmax = 0;
    float min_dm_am = 0;
    bool first = true;
    for (unsigned int higgs_idx = 0; higgs_idx < 3; higgs_idx++) {
       TLorentzVector higgs1_p = b_jet_p[higgs_combinations[higgs_idx][0]]+b_jet_p[higgs_combinations[higgs_idx][1]];
       TLorentzVector higgs2_p = b_jet_p[higgs_combinations[higgs_idx][2]]+b_jet_p[higgs_combinations[higgs_idx][3]];
       float delta_m = TMath::Abs(higgs1_p.M()-higgs2_p.M());
       if ((delta_m < min_delta_m) || first) {
         first = false;
         min_delta_m = delta_m;
	 float dr_1 = deltaR(b_jet_p[higgs_combinations[higgs_idx][0]].Eta(), b_jet_p[higgs_combinations[higgs_idx][0]].Phi(), b_jet_p[higgs_combinations[higgs_idx][1]].Eta(),b_jet_p[higgs_combinations[higgs_idx][1]].Phi());
	 float dr_2 = deltaR(b_jet_p[higgs_combinations[higgs_idx][2]].Eta(), b_jet_p[higgs_combinations[higgs_idx][2]].Phi(), b_jet_p[higgs_combinations[higgs_idx][3]].Eta(),b_jet_p[higgs_combinations[higgs_idx][3]].Phi());
	 min_dm_drmax = dr_1 > dr_2 ? dr_1 : dr_2;
	 min_dm_am = (higgs1_p.M() + higgs2_p.M())/2.0;
       }
    }
    return (min_delta_m <= 40)&&(min_dm_drmax < 2.2)&&(min_dm_am<200);
  });

  NamedFunc baseline = "stitch&&150<=met&&nvlep==0&&ntk==0&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc baseline_jetid = "stitch&&150<=met&&nvlep==0&&ntk==0" && jetid_modified_cuts;
  NamedFunc singleele_baseline = el_trigger && "stitch&&150<=met&&nlep==1&&nel==1&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc singleele_baseline_jetid = el_trigger && "stitch&&150<=met&&nlep==1&&nel==1" && jetid_modified_cuts;
  NamedFunc singlemu_baseline = mu_trigger && "stitch&&150<=met&&nlep==1&&nmu==1&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc singlemu_baseline_jetid = mu_trigger && "stitch&&150<=met&&nlep==1&&nmu==1" && jetid_modified_cuts;

  PlotMaker pm;

  //pm.Push<Hist2D>(Axis(100, 0, 3.14, "jet_met_dphi", "#Delta#phi", {}), Axis(100, 30, 500, "jet_pt", "Jet p_{T} [GeV]", {}),
  //  "2.4<=jet_eta&&jet_eta<5.0&&met>=150",
  //  procs_met2017c, twodim_plotopts).Weight(mixed_model_weight).Tag("FixName:EcalNoiseJetFilter_plot").LuminosityTag("9.63");
  //pm.Push<Hist2D>(Axis(100, 0, 3.14, "jet_met_dphi[0]", "Jet 1 p_{T} [GeV]", {}), Axis(100, 1.0, 3.0, htratio, "H_{T5}/H_{T}", {}), 
  //  "njet>=1&&met>=150",
  //  procs_met2017c, twodim_plotopts).Weight(mixed_model_weight).Tag("FixName:HTRatioDPhiTightFilter_plot").LuminosityTag("9.63");

  if (do_ele_efficiency) {
    //1e region with baseline selection
    pm.Push<Table>("cutflow_oneele_baseline_"+year_string, vector<TableRow>{
    TableRow("Events passing single electron triggers", 
      el_trigger,0,0, mixed_model_weight),
    TableRow("Baseline, $p_\\text{T}^\\text{miss} \\geq 150$ GeV, $\\Delta \\phi$ cuts", 
      singleele_baseline,0,0, mixed_model_weight),
    TableRow("goodVerticesFilter", 
      singleele_baseline && "npv_good>0",0,0, mixed_model_weight),
    TableRow("globalSuperTightHalo2016Filter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
    TableRow("HBHENoiseFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
    TableRow("HBHENoiseIsoFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
    TableRow("BadPFMuonFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
    TableRow("eeBadScFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
    TableRow("EcalNoiseJetFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("MuonJetFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/p_\\text{Tcalo}^\\text{miss}<2.0$", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/H_\\text{T}^\\text{miss}<2.0$", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("LowNeutralJetFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HTRatioDPhiTightFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    TableRow("JetID", 
      singleele_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet ,0,0, mixed_model_weight*pass_boostedhemveto),
    },procs_onelep,false,true,false,false,false,false);

    pm.Push<Table>("cutflow_oneele_baseline_jetid_"+year_string, vector<TableRow>{
    TableRow("Events passing single electron triggers", 
      el_trigger,0,0, mixed_model_weight),
    TableRow("Baseline, $p_\\text{T}^\\text{miss} \\geq 150 GeV$, $\\Delta \\phi$ cuts, with Jet ID", 
      singleele_baseline_jetid,0,0, mixed_model_weight),
    TableRow("goodVerticesFilter", 
      singleele_baseline_jetid && "npv_good>0",0,0, mixed_model_weight),
    TableRow("globalSuperTightHalo2016Filter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
    TableRow("HBHENoiseFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
    TableRow("HBHENoiseIsoFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
    TableRow("BadPFMuonFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
    TableRow("eeBadScFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
    TableRow("EcalNoiseJetFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("MuonJetFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/p_\\text{Tcalo}^\\text{miss}<2.0$", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/H_\\text{T}^\\text{miss}<2.0$", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("LowNeutralJetFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HTRatioDPhiTightFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    TableRow("JetID", 
      singleele_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    },procs_onelep,false,true,false,false,false,false);
  }
  if (do_mu_efficiency) {
    //1mu region with baseline selection
    pm.Push<Table>("cutflow_onemu_baseline_"+year_string, vector<TableRow>{
    TableRow("Events passing single electron triggers", 
      mu_trigger,0,0, mixed_model_weight),
    TableRow("Baseline, $p_\\text{T}^\\text{miss} \\geq 150$ GeV, $\\Delta \\phi$ cuts", 
      singlemu_baseline,0,0, mixed_model_weight),
    TableRow("goodVerticesFilter", 
      singlemu_baseline && "npv_good>0",0,0, mixed_model_weight),
    TableRow("globalSuperTightHalo2016Filter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
    TableRow("HBHENoiseFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
    TableRow("HBHENoiseIsoFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
    TableRow("BadPFMuonFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
    TableRow("eeBadScFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
    TableRow("EcalNoiseJetFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("MuonJetFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/p_\\text{Tcalo}^\\text{miss}<2.0$", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/H_\\text{T}^\\text{miss}<2.0$", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("LowNeutralJetFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HTRatioDPhiTightFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    TableRow("JetID", 
      singlemu_baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    },procs_onelep,false,true,false,false,false,false);

    pm.Push<Table>("cutflow_onemu_baseline_jetid_"+year_string, vector<TableRow>{
    TableRow("Events passing single electron triggers", 
      mu_trigger,0,0, mixed_model_weight),
    TableRow("Baseline, $p_\\text{T}^\\text{miss} \\geq 150$ GeV, $\\Delta \\phi$ cuts, with Jet ID", 
      singlemu_baseline_jetid,0,0, mixed_model_weight),
    TableRow("goodVerticesFilter", 
      singlemu_baseline_jetid && "npv_good>0",0,0, mixed_model_weight),
    TableRow("globalSuperTightHalo2016Filter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight",0,0, mixed_model_weight),
    TableRow("HBHENoiseFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe",0,0, mixed_model_weight),
    TableRow("HBHENoiseIsoFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso",0,0, mixed_model_weight),
    TableRow("EcalDeadCellTriggerPrimitiveFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell",0,0, mixed_model_weight),
    TableRow("BadPFMuonFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu",0,0, mixed_model_weight),
    TableRow("eeBadScFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc",0,0, mixed_model_weight),
    TableRow("EcalNoiseJetFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("MuonJetFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/p_\\text{Tcalo}^\\text{miss}<2.0$", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}/H_\\text{T}^\\text{miss}<2.0$", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("LowNeutralJetFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HTRatioDPhiTightFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    TableRow("JetID", 
      singlemu_baseline_jetid && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    },procs_onelep,false,true,false,false,false,false);
  }
  
  if (do_mc_efficiency) {
    //signal/MC efficiency
    pm.Push<Table>("cutflow_filters_mc_"+year_string, vector<TableRow>{
    TableRow("Baseline", 
      baseline,0,0, mixed_model_weight),
    TableRow("goodVerticesFilter", 
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
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight),
    TableRow("HEMDPhiVetoFilter", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
    TableRow("JetID", 
      baseline && "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_jets" && pass_ecalnoisejet,0,0, mixed_model_weight*pass_boostedhemveto),
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
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet,0,0, "1"),
    TableRow("HTRatioDPhiTightFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed ,0,0, "1"),
    TableRow("HEMDPhiVetoFilter", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed,0,0, pass_boostedhemveto),
    TableRow("JetID", 
      "HLT_PFMET120_PFMHT120_IDTight&&npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_jets" && pass_ecalnoisejet && pass_htratio_dphi_tight_fixed,0,0, pass_boostedhemveto),
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
      {"year", required_argument, 0, 0},
      {"do_mceff", no_argument, 0, 0},
      {"do_ele_eff", no_argument, 0, 0},
      {"do_mu_eff", no_argument, 0, 0},
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
      } else if(optname == "do_ele_eff"){
        do_ele_efficiency = true;
      } else if(optname == "do_mu_eff"){
        do_mu_efficiency = true;
      } else if(optname == "year"){
        year_string = optarg;
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
