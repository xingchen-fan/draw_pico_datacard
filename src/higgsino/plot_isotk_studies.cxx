//This script generates the following supplementary plots and tables: pie charts, cutflows
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options: 
//ttbar_reco,signal_isotkeff,cr_isotkeff,pies

#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TPie.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/cross_sections.hpp"
#include "core/event_scan.hpp"
#include "core/functions.hpp"
#include "core/gamma_params.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(options.unblind);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(options.unblind);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();

  set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  int lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);

  std::vector<std::shared_ptr<Process>> ttbar_procs = 
      script_utilities::getall_processes(
      years, {}, "ttbar", options.unblind);

  std::vector<std::shared_ptr<Process>> zll_procs = 
      script_utilities::getall_processes(
      years, {}, "zll", options.unblind);

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  //mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mctags["other_and_single_t"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root","*_ST_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/";
  string mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string mc_zll_skim_folder = "mc/merged_higmc_higlep2T/";
  string sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";

  std::vector<std::shared_ptr<Process>> lep_procs;
  lep_procs.push_back(Process::MakeShared<Baby_pico>("0 leptons", 
      Process::Type::background,kBlue+3,
      attach_folder(mc_base_folder,years,mc_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==0)"));
  lep_procs.push_back(Process::MakeShared<Baby_pico>("1 lepton", 
      Process::Type::background,kGreen-3,
      attach_folder(mc_base_folder,years,mc_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==1)"));
  lep_procs.push_back(Process::MakeShared<Baby_pico>("2+ leptons (no hadronic taus)", 
      Process::Type::background,kRed+2,
      attach_folder(mc_base_folder,years,mc_skim_folder,mctags["all"]),"stitch&&(ntrulep>=2)&&(ntrutauh==0)"));
  lep_procs.push_back(Process::MakeShared<Baby_pico>("2+ leptons (hadronic taus)", 
      Process::Type::background,kOrange+1,
      attach_folder(mc_base_folder,years,mc_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh>=2)&&(ntrutauh>=1)"));

  std::vector<std::shared_ptr<Process>> zll_lep_procs;
  zll_lep_procs.push_back(Process::MakeShared<Baby_pico>("0 leptons", 
      Process::Type::background,kBlue+3,
      attach_folder(mc_base_folder,years,mc_zll_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==0)"));
  zll_lep_procs.push_back(Process::MakeShared<Baby_pico>("1 lepton", 
      Process::Type::background,kGreen-3,
      attach_folder(mc_base_folder,years,mc_zll_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==1)"));
  zll_lep_procs.push_back(Process::MakeShared<Baby_pico>("2 leptons", 
      Process::Type::background,kRed+2,
      attach_folder(mc_base_folder,years,mc_zll_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==2)"));
  zll_lep_procs.push_back(Process::MakeShared<Baby_pico>("3+ leptons", 
      Process::Type::background,kOrange+1,
      attach_folder(mc_base_folder,years,mc_zll_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh>=3)"));

  std::vector<std::shared_ptr<Process>> signal_procs;
  std::vector<float> signal_masses;
  int higgsino_mass = 125;
  while (higgsino_mass <= 1250) {
    int real_higgsino_mass = higgsino_mass;
    if (real_higgsino_mass==125) real_higgsino_mass+=2;
    signal_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G("+std::to_string(real_higgsino_mass)+",1)", 
        Process::Type::signal,kBlack,
        attach_folder(sig_base_folder,years,sig_skim_folder,{"*SMS-TChiHH_mChi-"+std::to_string(real_higgsino_mass)+"_mLSP-0_*.root"}),"stitch"));
    signal_masses.push_back(static_cast<float>(real_higgsino_mass));
    higgsino_mass += 25;
  }

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc filters = Higfuncs::final_pass_filters;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc sr_baseline = filters&&script_utilities::search_resolved;
  NamedFunc ttbar_resolved = filters&&script_utilities::ttbar_resolved;
  NamedFunc zll_resolved = filters&&script_utilities::zll_resolved;

  const NamedFunc ntruleptot("ntruleptot",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns number of truth leptons including hadronic taus
    return b.ntrulep()+b.ntrutauh();
  });

  const NamedFunc ttbar_reco("ttbar_reco",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        if (!(60 <= w.M() && w.M() <= 100)) continue;
        for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
          ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
          ROOT::Math::PtEtaPhiMVector top = w+bjet;
          if (150 <= top.M() && top.M() <= 200) return 1;
        }
      }
    }
    return 0;
  });

  const NamedFunc ttbar_reco_w("ttbar_reco_w",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    float min_x = 999.0;
    float min_x_wmass = -999.0;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
          ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
          ROOT::Math::PtEtaPhiMVector top = w+bjet;
          float x = TMath::Sqrt((w.M()-80.4)*(w.M()-80.4)/64.64+(top.M()-172.5)*(top.M()-172.5)/297.56);
          if (x < min_x) {
            min_x = x;
            min_x_wmass = w.M();
          }
        }
      }
    }
    return min_x_wmass;
  });

  const NamedFunc ttbar_reco_w_unbiased("ttbar_reco_w_unbiased",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    float min_topmass_diff = 999.0;
    float best_wmass = -999.0;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
          ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
          ROOT::Math::PtEtaPhiMVector top = w+bjet;
          float topmass_diff = TMath::Abs(top.M()-172.5);
          if (topmass_diff < min_topmass_diff) {
            min_topmass_diff = topmass_diff;
            best_wmass = w.M();
          }
        }
      }
    }
    return best_wmass;
  });

  const NamedFunc ttbar_reco_top_unbiased("ttbar_reco_top_unbiased",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    float min_wmass_diff = 999.0;
    float best_topmass = -999.0;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        float wmass_diff = TMath::Abs(w.M()-80.4);
        if (wmass_diff < min_wmass_diff) {
          min_wmass_diff = wmass_diff;
          float min_topmass_diff = 999.0;
          for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
            ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
            ROOT::Math::PtEtaPhiMVector top = w+bjet;
            float topmass_diff = TMath::Abs(top.M()-172.5);
            if (topmass_diff < min_topmass_diff) {
              min_topmass_diff = topmass_diff;
              best_topmass = top.M();
            }
          }
        } // /if wmass_diff < min_wmass_diff
      }
    }
    return best_topmass;
  });

  const NamedFunc ttbar_reco_top("ttbar_reco_top",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    float min_x = 999.0;
    float min_x_topmass = -999.0;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
          ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
          ROOT::Math::PtEtaPhiMVector top = w+bjet;
          float x = TMath::Sqrt((w.M()-80.4)*(w.M()-80.4)/64.64+(top.M()-172.5)*(top.M()-172.5)/297.56);
          if (x < min_x) {
            min_x = x;
            min_x_topmass = top.M();
          }
        }
      }
    }
    return min_x_topmass;
  });

  const NamedFunc ttbar_reco_x("ttbar_reco_x",[](const Baby &b) -> NamedFunc::ScalarType{
    //tries to reconstruct a top quark from 1 b-jet and two other jets in the event
    if (b.njet() < 4) return false;
    //first, determine 2 jets with highest b-tag score to be b candidates - should be 2 since nbt=2
    unsigned highest_b_idx[2] = {999,999};
    float highest_b_score[2] = {-1,-1};
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      for (unsigned highest_idx = 0; highest_idx < 2; highest_idx++) {
        if (b.jet_deepcsv()->at(ijet) > highest_b_score[highest_idx]) {
          highest_b_idx[highest_idx] = ijet;
          highest_b_score[highest_idx] = b.jet_deepcsv()->at(ijet);
          continue;
        }
      }
    }
    //then try some combinatorics
    float min_x = 999.0;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (ijet1 == highest_b_idx[0] || ijet1 == highest_b_idx[1]) continue;
      for (unsigned ijet2 = 0; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (ijet2 == highest_b_idx[0] || ijet2 == highest_b_idx[1]) continue;
        ROOT::Math::PtEtaPhiMVector wjet1(b.jet_pt()->at(ijet1),b.jet_eta()->at(ijet1),b.jet_phi()->at(ijet1),b.jet_m()->at(ijet1));
        ROOT::Math::PtEtaPhiMVector wjet2(b.jet_pt()->at(ijet2),b.jet_eta()->at(ijet2),b.jet_phi()->at(ijet2),b.jet_m()->at(ijet2));
        ROOT::Math::PtEtaPhiMVector w = wjet1+wjet2;
        for (unsigned ibjet = 0; ibjet < 2; ibjet++) {
          ROOT::Math::PtEtaPhiMVector bjet(b.jet_pt()->at(highest_b_idx[ibjet]),b.jet_eta()->at(highest_b_idx[ibjet]),b.jet_phi()->at(highest_b_idx[ibjet]),b.jet_m()->at(highest_b_idx[ibjet]));
          ROOT::Math::PtEtaPhiMVector top = w+bjet;
          float x = TMath::Sqrt((w.M()-80.4)*(w.M()-80.4)/64.64+(top.M()-172.5)*(top.M()-172.5)/297.56);
          if (x < min_x) {
            min_x = x;
          }
        }
      }
    }
    return min_x;
  });

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;
  if(HigUtilities::is_in_string_options(options.string_options, "pies")) {
    //nlep pie chart for basic selection
    pm.Push<Table>("FixName:isotk__pie__1lcrleps", vector<TableRow> ({TableRow("", ttbar_resolved, 0, 0, weight)}), lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__1lcr2b_leps", vector<TableRow> ({TableRow("", ttbar_resolved&&"nbm==2", 0, 0, weight)}), lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2lprocs", vector<TableRow> ({TableRow("", ttbar_resolved&&"(ntrulep+ntrutauh>=2)", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__1lprocs", vector<TableRow> ({TableRow("", ttbar_resolved&&"(ntrulep+ntrutauh==1)", 0, 0, weight)}), ttbar_procs, true, true, true);

    pm.Push<Table>("FixName:isotk__pie__2lcrleps", vector<TableRow> ({TableRow("", zll_resolved, 0, 0, weight)}), zll_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2l2bcrleps", vector<TableRow> ({TableRow("", zll_resolved&&"nbt==2", 0, 0, weight)}), zll_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2lcrprocs", vector<TableRow> ({TableRow("", zll_resolved, 0, 0, weight)}), zll_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2l2bcrprocs", vector<TableRow> ({TableRow("", zll_resolved&&"nbt==2", 0, 0, weight)}), zll_procs, true, true, true);
  }

  if(HigUtilities::is_in_string_options(options.string_options, "cr_isotkeff")) {
    pm.Push<Hist1D>(Axis(25, 0.0, 100.0, "tk_pt", "p_{T}^{tk} [GeV]", {}),
        zll_resolved,
        zll_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__2lcr")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0.0, 100.0, "tk_pt", "p_{T}^{tk} [GeV]", {}),
        zll_resolved&&"nbt==2",
        zll_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__2l2bcr")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Table>("zllcr_isotkeff_"+options.year_string, vector<TableRow>{
      TableRow("2l CR baseline", 
          zll_resolved,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          zll_resolved&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          zll_resolved&&"ntk==0",0,0,weight),
      TableRow("2l CR baseline + 2b", 
          zll_resolved&&"nbt==2",0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          zll_resolved&&"nbt==2&&ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          zll_resolved&&"nbt==2&&ntk==0",0,0,weight),
    },zll_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);

    pm.Push<Table>("ttbarcr_isotkeff_"+options.year_string, vector<TableRow>{
      TableRow("1l CR baseline", 
          ttbar_resolved,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          ttbar_resolved&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          ttbar_resolved&&"ntk==0",0,0,weight),
    },ttbar_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
  }
  
  if(HigUtilities::is_in_string_options(options.string_options, "signal_isotkeff")) {
    pm.Push<Table>("signal_isotkeff_"+options.year_string, vector<TableRow>{
      TableRow("baseline minus isotk", 
          filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40"
          ,0,0,weight),
      TableRow("baseline", 
          filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40"
          ,0,0,weight)
    },signal_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Tag("signal_isotkeff").Precision(1);
  }

  if(HigUtilities::is_in_string_options(options.string_options, "ttbar_reco")) {
    pm.Push<Hist1D>(Axis(50, 0.0, 10.0, ttbar_reco_x, "minimum X_{Wt}", {}),
        ttbar_resolved,
        ttbar_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__x__proctype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 0.0, 10.0, ttbar_reco_x, "minimum X_{Wt}", {}),
        ttbar_resolved,
        lep_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__x__leptype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 0.0, 150.0, ttbar_reco_w, "m_{W} (min x_{Wt})", {}),
        ttbar_resolved,
        ttbar_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mw__proctype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 0.0, 150.0, ttbar_reco_w, "m_{W} (min x_{Wt})", {}),
        ttbar_resolved,
        lep_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mw__leptype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 50.0, 300.0, ttbar_reco_top, "m_{top} (min x_{Wt})", {}),
        ttbar_resolved,
        ttbar_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mtop__proctype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 50.0, 300.0, ttbar_reco_top, "m_{top} (min x_{Wt})", {}),
        ttbar_resolved,
        lep_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mtop__leptype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 0.0, 150.0, ttbar_reco_w_unbiased, "m_{W} (best m_{top})", {}),
        ttbar_resolved,
        ttbar_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mw_lessbias__proctype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 0.0, 150.0, ttbar_reco_w_unbiased, "m_{W} (best m_{top})", {}),
        ttbar_resolved,
        lep_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mw_lessbias__leptype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 50.0, 300.0, ttbar_reco_top_unbiased, "m_{top} (best m_{W})", {}),
        ttbar_resolved,
        ttbar_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mtop_lessbias__proctype")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(50, 50.0, 300.0, ttbar_reco_top_unbiased, "m_{top} (best m_{W})", {}),
        ttbar_resolved,
        lep_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__mtop_lessbias__leptype")
        .LuminosityTag(total_luminosity_string);
  }
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  //------------------------------------------------------------------------------------
  //                                     post processing
  //------------------------------------------------------------------------------------

  if(HigUtilities::is_in_string_options(options.string_options, "signal_isotkeff")) {
    Table* isotkeff_table = static_cast<Table*>(pm.GetFigure("signal_isotkeff").get());
    TH1D signal_eff_num("signal_eff_num","",45,125,1250);
    TH1D signal_eff_den("signal_eff_num","",45,125,1250);
    //std::vector<float> signal_effs;
    for (unsigned signal_proc_idx = 0; signal_proc_idx < signal_procs.size(); signal_proc_idx++) {
      std::vector<GammaParams> proc_yield = isotkeff_table->Yield(signal_procs[signal_proc_idx].get(),1.0);
      signal_eff_num.SetBinContent(signal_proc_idx+1,proc_yield[1].Yield());
      signal_eff_den.SetBinContent(signal_proc_idx+1,proc_yield[0].Yield());
      signal_eff_num.SetBinError(signal_proc_idx+1,proc_yield[1].Uncertainty());
      signal_eff_den.SetBinError(signal_proc_idx+1,proc_yield[0].Uncertainty());
      //signal_effs.push_back(proc_yield[1].Yield()/proc_yield[0].Yield());
    }
    TCanvas can("can","",1000,1000);
    can.SetGrid();
    can.SetMargin(0.15,0.15,0.15,0.15);
    TGraphAsymmErrors signal_eff_graph(&signal_eff_num,&signal_eff_den,"cp");
    //TGraph signal_eff_graph(signal_procs.size(),&signal_masses[0],&signal_effs[0]);
    signal_eff_graph.SetMarkerStyle(kFullCircle);
    signal_eff_graph.SetMarkerSize(1);
    signal_eff_graph.SetTitle("Isolated Track Veto Efficiency, baseline selection excluding isolated track veto; m(#tilde{#chi}_{1}^{0}) [GeV]; Efficiency");
    signal_eff_graph.GetYaxis()->SetRangeUser(0,1.0);
    signal_eff_graph.Draw("AP");
    can.Print("plots/isotk__signaleff.pdf");
    std::cout << "Open plots/isotk__signaleff.pdf";
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

