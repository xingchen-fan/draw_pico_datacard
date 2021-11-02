//This script is used to perform various studies on isolated tracks targeted by a veto in the HH+MET analysis
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options: 
//ttbar_reco,signal_isotkeff,cr_isotkeff,ttbar2l,fastsimverify,fastsimisotk,tt2lan

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
#include "core/efficiency_plot.hpp"
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

  std::vector<std::shared_ptr<Process>> search_procs = 
      script_utilities::getall_processes(
      years, {}, "search", false);

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

  string new_base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath/";
  string mc_2l_skim_folder = "mc/skim_2l/";
  string data_2l_skim_folder = "data/skim_2l/";
  string mc_unskimmed_folder = "mc/unskimmed/";
  string fastsim_unskimmed_folder = "mc_FastSimJmeCorrection/unskimmed/";
  string mc_met150_folder = "mc/skim_met150/";
  string fastsim_met150_folder = "mc_FastSimJmeCorrection/skim_met150/";

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

  std::vector<std::shared_ptr<Process>> tt2l_procs;
  std::vector<std::shared_ptr<Process>> tt2l_lep_procs;

  if(HigUtilities::is_in_string_options(options.string_options, "ttbar2l")
    ||HigUtilities::is_in_string_options(options.string_options, "tt2lan")) {
    //globbing cms17 takes ages, so only do it if necessary
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                    attach_folder(new_base_folder, years, mc_2l_skim_folder, mctags["tt"]),"stitch&&(ntrutauh==0)"));
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                    attach_folder(new_base_folder, years, mc_2l_skim_folder, mctags["tt"]),"stitch&&(ntrutauh>0)"));
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                    attach_folder(new_base_folder, years, mc_2l_skim_folder,mctags["zjets"]),"stitch"));
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                    attach_folder(new_base_folder, years, mc_2l_skim_folder,mctags["wjets"]),"stitch"));
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                    attach_folder(new_base_folder, years, mc_2l_skim_folder, mctags["qcd"]),"stitch")); 
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background, colors("single_t"),
                    attach_folder(new_base_folder, years, mc_2l_skim_folder, mctags["single_t"]),"stitch")); 
    tt2l_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                    attach_folder(new_base_folder, years, mc_2l_skim_folder, mctags["other"]),"stitch"));
    if (options.unblind) {
      tt2l_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                      attach_folder(new_base_folder, years, data_2l_skim_folder, {"*.root"}), Higfuncs::el_trigger||Higfuncs::mu_trigger));
    }

    //tt2l_lep_procs.push_back(Process::MakeShared<Baby_pico>("0 leptons", 
    //    Process::Type::background,kBlue+3,
    //    attach_folder(new_base_folder,years,mc_2l_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==0)"));
    //tt2l_lep_procs.push_back(Process::MakeShared<Baby_pico>("1 lepton", 
    //    Process::Type::background,kGreen-3,
    //    attach_folder(new_base_folder,years,mc_2l_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==1)"));
    tt2l_lep_procs.push_back(Process::MakeShared<Baby_pico>("2 leptons", 
        Process::Type::background,kRed+2,
        attach_folder(new_base_folder,years,mc_2l_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh==2)"));
    tt2l_lep_procs.push_back(Process::MakeShared<Baby_pico>("3+ leptons", 
        Process::Type::background,kOrange+1,
        attach_folder(new_base_folder,years,mc_2l_skim_folder,mctags["all"]),"stitch&&(ntrulep+ntrutauh>=3)"));
  }

  //std::vector<std::shared_ptr<Process>> fastsimcheck_procs;
  //fastsimcheck_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FullSIM", Process::Type::background,colors("tt_1l"),
  //                attach_folder(new_base_folder, years, mc_unskimmed_folder, {"*TTJets_*Lept*"}),"1"));
  //fastsimcheck_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FastSIM", Process::Type::background,colors("tt_1l"),
  //                attach_folder(new_base_folder, years, fastsim_unskimmed_folder, {"*TTJets_*Lept*"}),"1"));

  std::vector<std::shared_ptr<Process>> fastsim_procs;
  fastsim_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FullSIM", Process::Type::signal,colors("tt_1l"),
                  attach_folder(new_base_folder, years, mc_met150_folder, {"*TTJets_SingleLeptFromT_Tune*","*TTJets_SingleLeptFromTbar_Tune*"}),"1"));
  fastsim_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FastSIM", Process::Type::signal, kBlack,
                  attach_folder(new_base_folder, years, fastsim_unskimmed_folder, {"*TTJets_*Lept*"}),"((met+met_tru)/2)>150"));
  //fastsim_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FullSIM", Process::Type::signal,colors("tt_1l"),
  //                attach_folder(new_base_folder, years, mc_unskimmed_folder, {"*TTJets_SingleLeptFromT_Tune*","*TTJets_SingleLeptFromTbar_Tune*"}),"1"));
  //fastsim_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} 1l FastSIM", Process::Type::signal, kBlack,
  //                attach_folder(new_base_folder, years, fastsim_unskimmed_folder, {"*TTJets_*Lept*"}),"1"));

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
  //currently no met cut, have to use SF leptons because no trigger efficiencies for opposite flavor
  NamedFunc ttbar2lcr = "met/met_calo<5&&met/mht<5&&weight<1.5&&"
                        "nlep==2&&(nel==2||nmu==2)&&nbt>=2&&nll>=1&&(ll_m[0]<80||ll_m[0]>100)"
                        &&filters&&Higfuncs::lead_signal_lepton_pt>30;
  NamedFunc ttbar2ljets = "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&"
                          "hig_cand_am[0]<200&&hig_cand_dm[0]<40";
  NamedFunc baseline_notkveto = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "!low_dphi_met&&nlep==1&&njet>=4&&njet<=5&&nbt>=2&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40"
                         &&filters; //also missing met cut since this differs between fastSIM and fullSIM


  const NamedFunc jet_has_isotk("jet_has_isotk",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> jet_has_isotk_;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      bool found_tk = false;
      for (unsigned itk = 0; itk < b.tk_pt()->size(); itk++) {
        float tk_jet_dr = deltaR(b.tk_eta()->at(itk),b.tk_phi()->at(itk),b.jet_eta()->at(ijet),b.jet_phi()->at(ijet));
        if (tk_jet_dr < 0.4) {
          found_tk = true;
        }
      }
      if (found_tk) {
        jet_has_isotk_.push_back(1);
      }
      else {
        jet_has_isotk_.push_back(0);
      }
    }
    return jet_has_isotk_;
  });

  const NamedFunc ntruleptot("ntruleptot",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns number of truth leptons including hadronic taus
    return b.ntrulep()+b.ntrutauh();
  });

  const NamedFunc isotkwgt_new("isotkwgt_new",[jet_has_isotk](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType() < 0) return 1.0; //data
    //new isotk sfs: 2016
    std::vector<double> jet_has_isotk_ = jet_has_isotk.GetVector(b);
    bool nonjet_isotk = false;
    if (b.ntk()>0) nonjet_isotk = true;
    float sf = 1.0;
    for (unsigned ijet = 0; ijet<b.jet_pt()->size(); ijet++) {
      if (!(b.jet_isgood()->at(ijet))) continue;
      if (jet_has_isotk_[ijet] > 0.5) nonjet_isotk = false;
      float dcsvmwp = 1.0;
      if (abs(b.SampleType())==2016) dcsvmwp = 0.6321;
      if (abs(b.SampleType())==2017) dcsvmwp = 0.4941;
      if (abs(b.SampleType())==2018) dcsvmwp = 0.4148;
      bool is_b = (b.jet_deepcsv()->at(ijet)>dcsvmwp);
      float pt = b.jet_pt()->at(ijet);
      if (abs(b.SampleType())==2016) {
        if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 0.987949;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 0.992924;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 0.991273;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 0.996162;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 0.99273;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.998874;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 0.993096;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 0.994477;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 0.991008;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 0.996098;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 0.991569;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.99908;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 1.20174;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 1.208;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 1.40904;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 1.26908;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 1.86818;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.27395;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 1.09692;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 1.15057;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 1.4038;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 1.26259;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 1.94433;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.19661;
      }
      else if (abs(b.SampleType())==2017) {
        if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 0.994873;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 1.0005;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 0.995164;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 0.997693;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 0.995211;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.998194;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 1.00284;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 1.00332;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 0.995773;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 0.997343;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 0.994794;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.9986;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 1.07672;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 0.987972;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 1.18702;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 1.13165;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 1.50718;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.31995;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 0.965735;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 0.928357;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 1.14894;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 1.1379;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 1.48881;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.16937;
      }
      else if (abs(b.SampleType())==2018) {
        if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 0.990637;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 0.996963;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 0.994828;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 0.992855;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 0.994328;
        else if ((jet_has_isotk_[ijet]<0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.998497;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 0.995294;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 0.999124;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 0.994703;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 0.992042;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 0.993308;
        else if ((jet_has_isotk_[ijet]<0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 0.997792;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=30&&pt<60)) sf *= 1.14531;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=60&&pt<90)) sf *= 1.07676;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=90&&pt<120)) sf *= 1.20435;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=120&&pt<150)) sf *= 1.41748;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=150&&pt<250)) sf *= 1.51097;
        else if ((jet_has_isotk_[ijet]>0.5)&&(!is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.29265;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=30&&pt<60)) sf *= 1.0601;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=60&&pt<90)) sf *= 1.0198;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=90&&pt<120)) sf *= 1.19043;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=120&&pt<150)) sf *= 1.42707;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=150&&pt<250)) sf *= 1.53061;
        else if ((jet_has_isotk_[ijet]>0.5)&&(is_b)&&(pt>=250&&pt<9999.0)) sf *= 1.33105;
      }
    }
    if (abs(b.SampleType())==2016) {
      if (b.njet() <= 2) {
        if (nonjet_isotk) {
          sf *= 1.38087;
        } else {
          sf *= 0.990393;
        }
      } else if (b.njet() == 3) {
        if (nonjet_isotk) {
          sf *= 1.39492;
        } else {
          sf *= 0.987977;
        }
      } else if (b.njet() == 4) {
        if (nonjet_isotk) {
          sf *= 1.33345;
        } else {
          sf *= 0.988068;
        }
      } else {
        if (nonjet_isotk) {
          sf *= 1.18028;
        } else {
          sf *= 0.992574;
        }
      }
    }
    else if (abs(b.SampleType())==2017) {
      if (b.njet() <= 2) {
        if (nonjet_isotk) {
          sf *= 1.3157;
        } else {
          sf *= 0.9922;
        }
      } else if (b.njet() == 3) {
        if (nonjet_isotk) {
          sf *= 1.32952;
        } else {
          sf *= 0.990666;
        }
      } else if (b.njet() == 4) {
        if (nonjet_isotk) {
          sf *= 1.162;
        } else {
          sf *= 0.994671;
        }
      } else {
        if (nonjet_isotk) {
          sf *= 1.16833;
        } else {
          sf *= 0.993896;
        }
      }
    }
    else if (abs(b.SampleType())==2018) {
      if (b.njet() <= 2) {
        if (nonjet_isotk) {
          sf *= 1.53728;
        } else {
          sf *= 0.987156;
        }
      } else if (b.njet() == 3) {
        if (nonjet_isotk) {
          sf *= 1.43281;
        } else {
          sf *= 0.987155;
        }
      } else if (b.njet() == 4) {
        if (nonjet_isotk) {
          sf *= 1.42807;
        } else {
          sf *= 0.985739;
        }
      } else {
        if (nonjet_isotk) {
          sf *= 1.23447;
        } else {
          sf *= 0.990622;
        }
      }
    }
    return sf;
  });

  const NamedFunc isotkwgt("isotkwgt",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns scale factors to make MC and data agree for isotk veto in tt2l region
    if (b.SampleType() < 0) return 1.0; //data
    if (b.SampleType() == 2017) {
      if (b.njet() <= 2) {
        if (b.ntk()>0)
          return 1.0321056;
        return 0.87757538/0.88138363; //~0.996
      }
      else if (b.njet() == 3) {
        if (b.ntk()>0)
          return 1.1245581;
        return 0.83468019/0.85299131; //~0.98
      }
      else if (b.njet() == 4) {
        if (b.ntk()>0)
          return 1.1161442;
        return 0.80132926/0.82200262; //~0.97
      }
      else {
        if (b.ntk()>0)
          return 1.1486036;
        return 0.76456767/0.79502734; //~0.96
      }
    }
    if (b.SampleType() == 2018) {
      if (b.njet() <= 2) {
        if (b.ntk()>0)
          return 1.1917092;
        return 0.86733216/0.88867432; //~0.98
      }
      else if (b.njet() == 3) {
        if (b.ntk()>0)
          return 1.2015485;
        return 0.82563291/0.85488135; //~0.97
      }
      else if (b.njet() == 4) {
        if (b.ntk()>0)
          return 1.2384423;
        return 0.78591616/0.82713459; //~0.95
      }
      else {
        if (b.ntk()>0)
          return 1.2741405;
        return 0.74345709/0.79865415; //~0.93
      }
    }
    else { //2016
      if (b.njet() <= 2) {
        if (b.ntk()>0)
          return 1.2493318;
        return 0.87450296/0.89954867; //~0.97
      }
      else if (b.njet() == 3) {
        if (b.ntk()>0)
          return 1.2393575;
        return 0.83343459/0.86560342; //~0.96
      }
      else if (b.njet() == 4) {
        if (b.ntk()>0)
          return 1.2534025;
        return 0.79139505/0.83356906; //~0.95
      }
      else {
        if (b.ntk()>0)
          return 1.2953800;
        return 0.75597554/0.81161940; //~0.93
      }
    }
    return 1.0;
  });

  NamedFunc weight_isotk = weight*isotkwgt;
  NamedFunc weight_isotk_new = weight*isotkwgt_new;

  const NamedFunc isotk_is_lepton("isotk_is_lepton",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns true if isolated track is close to a truth lepton
    if (b.ntk()==0) {
      return 0;
    }
    for (unsigned itk = 0; itk < b.tk_pt()->size(); itk++) {
      for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
        int pdgid = abs(b.mc_id()->at(imc));
        if(pdgid == 11 || pdgid == 13 || pdgid == 15) {
          float tk_lep_dr = deltaR(b.tk_eta()->at(itk),b.tk_phi()->at(itk),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
          if (tk_lep_dr < 0.4)
           return 1;
        }
      }
    }
    return 0;
  });

  const NamedFunc nonjet_isotk("nonjet_isotk",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns true if there is an isolated track not associated with a jet
    if (b.ntk()==0) {
      return false;
    }
    bool nonjet_isotk_ = true;
    for (unsigned itk = 0; itk < b.tk_pt()->size(); itk++) {
      for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
        if (b.jet_isgood()->at(ijet)) {
          float tk_jet_dr = deltaR(b.tk_eta()->at(itk),b.tk_phi()->at(itk),b.jet_eta()->at(ijet),b.jet_phi()->at(ijet));
          if (tk_jet_dr < 0.4)
            return nonjet_isotk_ = false;
        }
      }
    }
    return nonjet_isotk_;
  });

  const NamedFunc dr_tk_jet("dr_tk_jet",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns min delta-r between isolated track and jet in an event
    if (b.ntk()==0) {
      return 999.0;
    }
    float min_dr = 999.0;
    for (unsigned itk = 0; itk < b.tk_pt()->size(); itk++) {
      for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
        if (b.jet_isgood()->at(ijet)) {
          float tk_jet_dr = deltaR(b.tk_eta()->at(itk),b.tk_phi()->at(itk),b.jet_eta()->at(ijet),b.jet_phi()->at(ijet));
          if (tk_jet_dr < min_dr)
            min_dr = tk_jet_dr;
        }
      }
    }
    return min_dr;
  });

  const NamedFunc jet_is_bm("jet_is_bm",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> jet_is_bm_;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      float dcsvmwp = 1.0;
      if (abs(b.SampleType())==2016) dcsvmwp = 0.6321;
      if (abs(b.SampleType())==2017) dcsvmwp = 0.4941;
      if (abs(b.SampleType())==2018) dcsvmwp = 0.4148;
      if (b.jet_deepcsv()->at(ijet)>dcsvmwp) {
        jet_is_bm_.push_back(1);
      }
      else {
        jet_is_bm_.push_back(0);
      }
    }
    return jet_is_bm_;
  });

  const NamedFunc dr_tk_bjet("dr_tk_bjet",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns min delta-r between isolated track and b-jet in an event
    if (b.ntk()==0) {
      return 999.0;
    }
    float min_dr = 999.0;
    for (unsigned itk = 0; itk < b.tk_pt()->size(); itk++) {
      for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
        float dcsvmwp = 1.0;
        if (abs(b.SampleType())==2016) dcsvmwp = 0.6321;
        if (abs(b.SampleType())==2017) dcsvmwp = 0.4941;
        if (abs(b.SampleType())==2018) dcsvmwp = 0.4148;
        if (b.jet_isgood()->at(ijet) && b.jet_deepcsv()->at(ijet)>dcsvmwp) {
          float tk_jet_dr = deltaR(b.tk_eta()->at(itk),b.tk_phi()->at(itk),b.jet_eta()->at(ijet),b.jet_phi()->at(ijet));
          if (tk_jet_dr < min_dr)
            min_dr = tk_jet_dr;
        }
      }
    }
    return min_dr;
  });

  const NamedFunc min_jet_pt("avr_jet_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns lowest jet pt
    float min_jet_pt_ = 999.0;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (b.jet_isgood()->at(ijet)) {
        if (b.jet_pt()->at(ijet) < min_jet_pt_)
          min_jet_pt_ = b.jet_pt()->at(ijet);
      }
    }
    return min_jet_pt_;
  });

  const NamedFunc avr_jet_pt("avr_jet_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    //returns average jet pt
    float tot_jet_pt = 0;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (b.jet_isgood()->at(ijet)) {
        tot_jet_pt += b.jet_pt()->at(ijet);
      }
    }
    return tot_jet_pt/static_cast<float>(b.njet());
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

  const NamedFunc nh("nh",[](const Baby &b) -> NamedFunc::ScalarType{
    int nh_ = 0;
    for (unsigned ifjet = 0; ifjet < b.fjet_deep_md_hbb_btv()->size(); ifjet++) {
      if (ifjet > 1) break;
      if (b.fjet_deep_md_hbb_btv()->at(ifjet) > 0.7) nh_++;
    }
    return nh_;
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

  const std::vector<double> avjetpt_bins = {0,75,100,150,250,750};
  const std::vector<double> jetpt_bins = {30,60,90,120,150,250,500};

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;
  if(HigUtilities::is_in_string_options(options.string_options, "cr_isotkeff")) {
    //nlep pie chart for basic selection
    pm.Push<Table>("FixName:isotk__pie__1lcrleps_"+options.year_string, vector<TableRow> ({TableRow("", ttbar_resolved, 0, 0, weight)}), lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__1lcr2b_leps_"+options.year_string, vector<TableRow> ({TableRow("", ttbar_resolved&&"nbm==2", 0, 0, weight)}), lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2lprocs_"+options.year_string, vector<TableRow> ({TableRow("", ttbar_resolved&&"(ntrulep+ntrutauh>=2)", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__1lprocs_"+options.year_string, vector<TableRow> ({TableRow("", ttbar_resolved&&"(ntrulep+ntrutauh==1)", 0, 0, weight)}), ttbar_procs, true, true, true);

    pm.Push<Table>("FixName:isotk__pie__2lcrleps_"+options.year_string, vector<TableRow> ({TableRow("", zll_resolved, 0, 0, weight)}), zll_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2l2bcrleps_"+options.year_string, vector<TableRow> ({TableRow("", zll_resolved&&"nbt==2", 0, 0, weight)}), zll_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2lcrprocs_"+options.year_string, vector<TableRow> ({TableRow("", zll_resolved, 0, 0, weight)}), zll_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__2l2bcrprocs_"+options.year_string, vector<TableRow> ({TableRow("", zll_resolved&&"nbt==2", 0, 0, weight)}), zll_procs, true, true, true);
    //vars
    pm.Push<EfficiencyPlot>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {}),
        zll_resolved && "njet>=4",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_nb_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 40.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        zll_resolved && "njet>=4",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 0.0, 40.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        zll_resolved && "njet>=4",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 0.0, 2.2, "hig_cand_drmax[0]", "#DeltaR_{max}", {}),
        zll_resolved && "njet>=4",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 400, "met", "p_{T}^{miss} [GeV]", {}),
        zll_resolved && "njet>=4",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_met_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(3, -0.5, 2.5, nh, "N_{H}", {}),
        zll_resolved,
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_nh_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 200, "fjet_msoftdrop[0]", "m_{J1} [GeV]", {}),
        zll_resolved&&"nfjet>0",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_mj1_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 200, "fjet_msoftdrop[1]", "m_{J2} [GeV]", {}),
        zll_resolved&&"nfjet>1",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__2lcr_mj2_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    //1l & 2l CR(DY) plots
    pm.Push<EfficiencyPlot>(Axis({0,20,40,60,80,100,150,200,400}, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        zll_resolved,
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__2lcr__tkeff").YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis({0,20,40,60,80,100,150,200,400}, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        zll_resolved&&"nbt==2",
        "ntk==0",
        zll_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__2l2bcr__tkeff").YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis({0,20,40,60,80,100,150,200,400}, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar_resolved,
        "ntk==0",
        ttbar_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__1lcr__tkeff").YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
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
    if (options.unblind == false) {
      pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_jet, "min #Delta R(tk,jet)", {0.4}),
          ttbar_resolved&&"ntrutauh>0&&ntk>0",
          ttbar_procs, plt_lin).Weight(weight)
          .Tag("FixName:isotk__1lcr_mindr_"+options.year_string)
          .LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_bjet, "min #Delta R(tk,b-jet)", {0.4}),
          ttbar_resolved&&"ntrutauh>0&&ntk>0",
          ttbar_procs, plt_lin).Weight(weight)
          .Tag("FixName:isotk__1lcr_mindrb_"+options.year_string)
          .LuminosityTag(total_luminosity_string);
    }
  }

  if(HigUtilities::is_in_string_options(options.string_options, "sr_isotkeff")) {
    //some MC plots from the search region
    pm.Push<Hist1D>(Axis(10, 5.0, 30.0, "tk_pt[0]", "Track p_{T} [GeV]", {}),
        baseline_notkveto&&isotk_is_lepton&&"met>150&&ntrutah==0&&ntk==1",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_leptkpt_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0.0, 0.3, "tk_reliso_chg[0]", "Track I_{rel}", {}),
        baseline_notkveto&&isotk_is_lepton&&"met>150&&ntrutah==0&&ntk==1",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_leptkreliso_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0.0, 0.3, "tk_miniso_chg[0]", "Track I_{mini}", {}),
        baseline_notkveto&&isotk_is_lepton&&"met>150&&ntrutah==0&&ntk==1",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_leptkminiso_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_jet, "min #Delta R(tk,jet)", {0.4}),
        baseline_notkveto&&"met>150&&ntrutauh>0&&ntk>0",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_mindr_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_bjet, "min #Delta R(tk,b-jet)", {0.4}),
        baseline_notkveto&&"met>150&&ntrutauh>0&&ntk>0",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_mindrb_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_jet, "min #Delta R(tk,jet)", {0.4}),
        baseline_notkveto&&isotk_is_lepton&&"met>150&&ntk>0",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_mindr_leptks_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_jet, "min #Delta R(tk,jet)", {0.4}),
        baseline_notkveto&&!isotk_is_lepton&&"met>150&&ntk>0",
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__sr_mindr_nonleptks_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Table>("isotk__sr_isotktable_"+options.year_string, vector<TableRow>{
      TableRow("SR baseline except track veto - events with had. $\\tau$s", 
          baseline_notkveto&&"met>150&&ntrutauh>0",0,0,weight),
      TableRow("SR baseline except track veto - events with had. $\\tau$s and track", 
          baseline_notkveto&&"met>150&&ntrutauh>0&&ntk>0",0,0,weight),
      TableRow("SR baseline except track veto - events without had. $\\tau$s", 
          baseline_notkveto&&"met>150&&ntrutauh==0",0,0,weight),
      TableRow("SR baseline except track veto - events without had. $\\tau$s and with track", 
          baseline_notkveto&&"met>150&&ntrutauh==0&&ntk>0",0,0,weight),
      TableRow("SR baseline except track veto - events with lep. $\\tau$s", 
          baseline_notkveto&&"met>150&&ntrutauh==0&&ntrutaul>0",0,0,weight),
      TableRow("SR baseline except track veto - events with lep. $\\tau$s and with track", 
          baseline_notkveto&&"met>150&&ntrutauh==0&&ntrutaul>0&&ntk>0",0,0,weight),
      TableRow("SR baseline except track veto - events with track", 
          baseline_notkveto&&"met>150&&ntk>0",0,0,weight),
      TableRow("SR baseline except track veto - events with track that is lepton", 
          baseline_notkveto&&isotk_is_lepton&&"met>150&&ntk>0",0,0,weight),
      TableRow("SR baseline except track veto - events with track that is lepton (had tau)", 
          baseline_notkveto&&isotk_is_lepton&&"met>150&&ntk>0&&ntrutauh>0",0,0,weight),
      TableRow("SR baseline except track veto - events with track that is lepton (lep tau)", 
          baseline_notkveto&&isotk_is_lepton&&"met>150&&ntk>0&&ntrutaul>0",0,0,weight),
      TableRow("SR baseline except track veto - events with track that is lepton (not tau)", 
          baseline_notkveto&&isotk_is_lepton&&"met>150&&ntk>0&&ntrutauh==0&&ntrutaul==0",0,0,weight),
    },search_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
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

  if(HigUtilities::is_in_string_options(options.string_options, "fastsimisotk")) {
    pm.Push<Table>("isotk__fastfullsim_table_"+options.year_string, vector<TableRow>{
      TableRow("tt 2l CR baseline", 
          baseline_notkveto,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          baseline_notkveto&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          baseline_notkveto&&"ntk==0",0,0,weight),
    },fastsim_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, min_jet_pt, "Minimum jet p_{T} [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_minpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, min_jet_pt, "Minimum jet p_{T} [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_minpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_avrpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_avrpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(6, 150.0, 300.0, "met", "p_{T}^{miss} [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_met_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(10, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&!low_dphi_met&&nvlep==0&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(10, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&!low_dphi_met&&nvlep==0&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(10, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&!low_dphi_met&&nvlep==0&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_dm[0]<40&&hig_cand_am[0]<200",
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        filters&&"met/mht<2 && met/met_calo<2&&weight<1.5&&!low_dphi_met&&nvlep==0&&njet>=4&&njet<=5&&nbt>=2&&hig_cand_dm[0]<40&&hig_cand_am[0]<200",
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(2, 3.5, 5.5, "njet", "N_{j}", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__fastfullsim__tkeff_njet_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(2, 3.5, 5.5, "njet", "N_{j}", {}),
        baseline_notkveto,
        "ntk==0",
        fastsim_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__fastfullsim__tkeff_corr_njet_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
  }

  if(HigUtilities::is_in_string_options(options.string_options, "fastsimverify")) {
    pm.Push<Table>("isotk__fastsimverify_"+options.year_string, vector<TableRow>{
      TableRow("filters and $p_\\mathrm{T}^\\mathrm{miss}>150$ GeV", 
          filters&&"met/met_calo<5&&met/mht<5&&weight<1.5",0,0,weight),
      TableRow("$N_\\mathrm{l}=1$", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1",0,0,weight),
      TableRow("$4\\leq N_\\mathrm{j}\\leq 5$", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5",0,0,weight),
      TableRow("$N_\\mathrm{b}\\geq 2$", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5&&nbt>=2",0,0,weight),
      TableRow("$M_\\mathrm{T}<100$ GeV", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5&&nbt>=2&&mt<=100",0,0,weight),
      TableRow("$\\Delta R_\\mathrm{max}<2.2$", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5&&nbt>=2&&mt<=100&&hig_cand_drmax[0]<2.2",0,0,weight),
      TableRow("$\\langle m_\\mathrm{bb}\\rangle<200$ GeV", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5&&nbt>=2&&mt<=100&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",0,0,weight),
      TableRow("$\\Delta m_\\mathrm{bb}<40$ GeV", 
          filters&&(Higfuncs::lead_signal_lepton_pt>30)&&"met/met_calo<5&&met/mht<5&&weight<1.5&&nlep==1&&njet>=4&&njet<=5&&nbt>=2&&mt<=100&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40",0,0,weight),
    },fastsim_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
  }

  if(HigUtilities::is_in_string_options(options.string_options, "tt2lan")) {
    //numbers - tables
    pm.Push<Table>("tt2lcr_isotkeff_"+options.year_string, vector<TableRow>{
      TableRow("tt 2l CR baseline", 
          ttbar2lcr,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          ttbar2lcr&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          ttbar2lcr&&"ntk==0",0,0,weight),
      TableRow("tt 2l CR baseline + jets", 
          ttbar2lcr&&ttbar2ljets,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          ttbar2lcr&&ttbar2ljets&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          ttbar2lcr&&ttbar2ljets&&"ntk==0",0,0,weight),
    },tt2l_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
    //composition pie charts
    pm.Push<Table>("FixName:isotk__pie__tt2lcrleps_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr, 0, 0, weight)}), tt2l_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2ljetscrleps_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr&&ttbar2ljets, 0, 0, weight)}), tt2l_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2lcrprocs_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr, 0, 0, weight)}), tt2l_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2ljetscrprocs_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr&&ttbar2ljets, 0, 0, weight)}), tt2l_procs, true, true, true);
    //variables
    pm.Push<EfficiencyPlot>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {}),
        ttbar2lcr && "njet>=4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_nb_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 40.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr && "njet>=4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        ttbar2lcr && "njet>=4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0.0, 4.4, "hig_cand_drmax[0]", "#DeltaR_{max}", {}),
        ttbar2lcr && "njet>=4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 400, "met", "p_{T}^{miss} [GeV]", {}),
        ttbar2lcr && "njet>=4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_met_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(3, -0.5, 2.5, nh, "N_{H}", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_nh_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 200, "fjet_msoftdrop[0]", "m_{J1} [GeV]", {}),
        ttbar2lcr&&"nfjet>0",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_mj1_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0, 200, "fjet_msoftdrop[1]", "m_{J2} [GeV]", {}),
        ttbar2lcr&&"nfjet>1",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__vars__tt2lcr_mj2_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
  }

  if(HigUtilities::is_in_string_options(options.string_options, "ttbar2l")) {
    //Keith scheme - efficiency per jet and for non-jet tracks
    pm.Push<EfficiencyPlot>(Axis(jetpt_bins, "jet_pt", "Jet p_{T} [GeV]", {}),
        ttbar2lcr&&"jet_isgood",
        !jet_has_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__tkeff_jetpt_"+options.year_string).YTitle("No Iso. Track").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(jetpt_bins, "jet_pt", "b-Jet p_{T} [GeV]", {}),
        ttbar2lcr&&jet_is_bm&&"jet_isgood",
        !jet_has_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__tkeff_bjetpt_"+options.year_string).YTitle("No Iso. Track").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 1.5, 5.5, "njet", "N_{j}", {}),
        ttbar2lcr,
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__nonjettkeff_njet_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0, 400.0, avr_jet_pt, "Average jet p_{T}", {}),
        ttbar2lcr,
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__nonjettkeff_avrpt_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 50.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__nonjettkeff_am_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__nonjettkeff_dm_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__nonjettkeff_drmax_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    //verification plots for Keith's scheme
    pm.Push<EfficiencyPlot>(Axis(jetpt_bins, "jet_pt", "Jet p_{T} [GeV]", {}),
        ttbar2lcr&&"jet_isgood",
        !jet_has_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__tkeff_corrnew_jetpt_"+options.year_string).YTitle("No Iso. Track").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(jetpt_bins, "jet_pt", "b-Jet p_{T} [GeV]", {}),
        ttbar2lcr&&jet_is_bm&&"jet_isgood",
        !jet_has_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__tkeff_corrnew_bjetpt_"+options.year_string).YTitle("No Iso. Track").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 1.5, 5.5, "njet", "N_{j}", {}),
        ttbar2lcr,
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__nonjettkeff_corrnew_njet_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0, 400.0, avr_jet_pt, "Average jet p_{T}", {}),
        ttbar2lcr,
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__nonjettkeff_corrnew_avrpt_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 50.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__nonjettkeff_corrnew_am_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__nonjettkeff_corrnew_dm_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        ttbar2lcr&&"njet>=4",
        !nonjet_isotk,
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__nonjettkeff_corrnew_drmax_"+options.year_string).YTitle("N_{njtk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2lcr__tkeff_corrnew_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 50.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2ljetscr__tkeff_corrnew_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {40.0}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2ljetscr__tkeff_corrnew_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2ljetscr__tkeff_corrnew_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 1.5, 5.5, "njet", "N_{j}", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk_new).Tag("FixName:isotk__tt2ljetscr__tkeff_corrnew_njet_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    //previous nj scheme
    //pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, min_jet_pt, "Minimum jet p_{T} [GeV]", {}),
    //    ttbar2lcr,
    //    "ntk==0",
    //    tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__tkeff_minpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, min_jet_pt, "Minimum jet p_{T} [GeV]", {}),
    //    ttbar2lcr,
    //    "ntk==0",
    //    tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2lcr__tkeff_corr_minpt_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2lcr__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(16, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2lcr__tkeff_corr_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(8, 0.0, 400.0, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr&&ttbar2ljets,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(avjetpt_bins, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr&&"njet==2",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr_nj2__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(avjetpt_bins, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr&&"njet==3",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr_nj3__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(avjetpt_bins, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr&&"njet==4",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr_nj4__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(avjetpt_bins, avr_jet_pt, "Average jet p_{T} [GeV]", {}),
        ttbar2lcr&&"njet==5",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr_nj5__tkeff_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 50.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr__tkeff_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 50.0, 200.0, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2ljetscr__tkeff_corr_am_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {40.0}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr__tkeff_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 100.0, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {40.0}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2ljetscr__tkeff_corr_dm_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr__tkeff_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 0.0, 4.4, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        ttbar2lcr && "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40",
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2ljetscr__tkeff_corr_drmax_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 1.5, 5.5, "njet", "N_{j}", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight).Tag("FixName:isotk__tt2ljetscr__tkeff_njet_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(4, 1.5, 5.5, "njet", "N_{j}", {}),
        ttbar2lcr,
        "ntk==0",
        tt2l_procs,true,plt_lin).Weight(weight_isotk).Tag("FixName:isotk__tt2ljetscr__tkeff_corr_njet_"+options.year_string).YTitle("N_{tk}=0").LuminosityTag(total_luminosity_string);
    pm.Push<Table>("FixName:isotk__pie__tt2lcrleps_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr, 0, 0, weight)}), tt2l_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2ljetscrleps_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr&&ttbar2ljets, 0, 0, weight)}), tt2l_lep_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2lcrprocs_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr, 0, 0, weight)}), tt2l_procs, true, true, true);
    pm.Push<Table>("FixName:isotk__pie__tt2ljetscrprocs_"+options.year_string, vector<TableRow> (
        {TableRow("", ttbar2lcr&&ttbar2ljets, 0, 0, weight)}), tt2l_procs, true, true, true);
    pm.Push<Table>("tt2lcr_isotkeff_"+options.year_string, vector<TableRow>{
      TableRow("tt 2l CR baseline", 
          ttbar2lcr,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          ttbar2lcr&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          ttbar2lcr&&"ntk==0",0,0,weight),
      TableRow("tt 2l CR baseline + jets", 
          ttbar2lcr&&ttbar2ljets,0,0,weight),
      TableRow("$\\geq 1$ Iso.Tk.", 
          ttbar2lcr&&ttbar2ljets&&"ntk>=1",0,0,weight),
      TableRow("0 Iso.Tk.", 
          ttbar2lcr&&ttbar2ljets&&"ntk==0",0,0,weight),
    },tt2l_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Precision(1);
    pm.Push<Hist1D>(Axis(25, 0.0, 100.0, "tk_pt", "p_{T}^{tk} [GeV]", {}),
        ttbar2lcr,
        tt2l_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__tt2lcr_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0.0, 100.0, "tk_pt", "p_{T}^{tk} [GeV]", {}),
        ttbar2lcr&&ttbar2ljets,
        tt2l_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__tt2ljetscr_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_jet, "min #Delta R(tk,jet)", {0.4}),
        ttbar2lcr&&"ntk>0",
        tt2l_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__tt2lcr_mindr_"+options.year_string)
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 0.0, 3.0, dr_tk_bjet, "min #Delta R(tk,b-jet)", {0.4}),
        ttbar2lcr&&"ntk>0",
        tt2l_procs, plt_lin).Weight(weight)
        .Tag("FixName:isotk__tkpt__tt2lcr_mindrb_"+options.year_string)
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

  if(HigUtilities::is_in_string_options(options.string_options, "ttbar2l")) {
    std::string eff_plot_name = "isotk__tt2ljetscr__tkeff_njet_"+options.year_string;
    //EfficiencyPlot * eff_isotk_tt2l = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //TGraphAsymmErrors* eff_isotk_tt2l_data = static_cast<TGraphAsymmErrors*>((eff_isotk_tt2l->data_ratio_plots_[0]).get()->Clone());
    //TGraphAsymmErrors* eff_isotk_tt2l_mc = static_cast<TGraphAsymmErrors*>((eff_isotk_tt2l->background_ratio_plot_).get()->Clone());

    eff_plot_name = "isotk__tt2lcr__tkeff_jetpt_"+options.year_string;
    EfficiencyPlot * eff_isotk_jetpt = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    TGraphAsymmErrors* eff_isotk_jetpt_data = static_cast<TGraphAsymmErrors*>((eff_isotk_jetpt->data_ratio_plots_[0]).get()->Clone());
    TGraphAsymmErrors* eff_isotk_jetpt_mc = static_cast<TGraphAsymmErrors*>((eff_isotk_jetpt->background_ratio_plot_).get()->Clone());

    eff_plot_name = "isotk__tt2lcr__tkeff_bjetpt_"+options.year_string;
    EfficiencyPlot * eff_isotk_bjetpt = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    TGraphAsymmErrors* eff_isotk_bjetpt_data = static_cast<TGraphAsymmErrors*>((eff_isotk_bjetpt->data_ratio_plots_[0]).get()->Clone());
    TGraphAsymmErrors* eff_isotk_bjetpt_mc = static_cast<TGraphAsymmErrors*>((eff_isotk_bjetpt->background_ratio_plot_).get()->Clone());

    eff_plot_name = "isotk__tt2lcr__nonjettkeff_njet_"+options.year_string;
    EfficiencyPlot * eff_isotk_njet = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    TGraphAsymmErrors* eff_isotk_njet_data = static_cast<TGraphAsymmErrors*>((eff_isotk_njet->data_ratio_plots_[0]).get()->Clone());
    TGraphAsymmErrors* eff_isotk_njet_mc = static_cast<TGraphAsymmErrors*>((eff_isotk_njet->background_ratio_plot_).get()->Clone());

    //std::vector<TGraphAsymmErrors*> data_effplots;
    //std::vector<TGraphAsymmErrors*> mc_effplots;
    //eff_plot_name = "isotk__tt2ljetscr_nj2__tkeff_"+options.year_string;
    //EfficiencyPlot * eff_isotk_tt2lnj2 = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //data_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj2->data_ratio_plots_[0]).get()->Clone()));
    //mc_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj2->background_ratio_plot_).get()->Clone()));

    //eff_plot_name = "isotk__tt2ljetscr_nj3__tkeff_"+options.year_string;
    //EfficiencyPlot * eff_isotk_tt2lnj3 = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //data_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj3->data_ratio_plots_[0]).get()->Clone()));
    //mc_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj3->background_ratio_plot_).get()->Clone()));

    //eff_plot_name = "isotk__tt2ljetscr_nj4__tkeff_"+options.year_string;
    //EfficiencyPlot * eff_isotk_tt2lnj4 = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //data_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj4->data_ratio_plots_[0]).get()->Clone()));
    //mc_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj4->background_ratio_plot_).get()->Clone()));

    //eff_plot_name = "isotk__tt2ljetscr_nj5__tkeff_"+options.year_string;
    //EfficiencyPlot * eff_isotk_tt2lnj5 = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //data_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj5->data_ratio_plots_[0]).get()->Clone()));
    //mc_effplots.push_back(static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lnj5->background_ratio_plot_).get()->Clone()));

    TFile * isotkeff_file = TFile::Open(("tables/isotk__effnew_"+options.year_string+".root").c_str(),"RECREATE");
    eff_isotk_jetpt_data->Write("eff_isotk_jetpt_data");
    eff_isotk_jetpt_mc->Write("eff_isotk_jetpt_mc");
    eff_isotk_bjetpt_data->Write("eff_isotk_bjetpt_data");
    eff_isotk_bjetpt_mc->Write("eff_isotk_bjetpt_mc");
    eff_isotk_njet_data->Write("eff_isotk_njet_data");
    eff_isotk_njet_mc->Write("eff_isotk_njet_mc");

    std::cout << "Printing code for scale factors" << std::endl;
    std::cout << "  std::vector<double> jet_has_isotk_ = jet_has_isotk.GetVector(b);" << std::endl;
    std::cout << "  bool nonjet_isotk = false;" << std::endl;
    std::cout << "  if (b.ntk()>0) nonjet_isotk = true;" << std::endl;
    std::cout << "  float sf = 1.0;" << std::endl;
    std::cout << "  for (unsigned ijet = 0; ijet<b.jet_pt()->size(); ijet++) {" << std::endl;
    std::cout << "    if (!(b.jet_isgood()->at(ijet))) continue;" << std::endl;
    std::cout << "    if (jet_has_isotk_[ijet] > 0.5) nonjet_isotk = false;" << std::endl;
    std::cout << "    float dcsvmwp = 1.0;" << std::endl;
    std::cout << "    if (abs(b.SampleType())==2016) dcsvmwp = 0.6321;" << std::endl;
    std::cout << "    if (abs(b.SampleType())==2017) dcsvmwp = 0.4941;" << std::endl;
    std::cout << "    if (abs(b.SampleType())==2018) dcsvmwp = 0.4148;" << std::endl;
    std::cout << "    bool is_b = (b.jet_deepcsv()->at(ijet)>dcsvmwp);" << std::endl;
    std::cout << "    float pt = b.jet_pt()->at(ijet);" << std::endl;
    bool first = true;
    double * nonb_mc_frac = eff_isotk_jetpt_mc->GetY();
    double * nonb_data_frac = eff_isotk_jetpt_data->GetY();
    double * b_mc_frac = eff_isotk_bjetpt_mc->GetY();
    double * b_data_frac = eff_isotk_bjetpt_data->GetY();
    double * njet_mc_frac = eff_isotk_njet_mc->GetY();
    double * njet_data_frac = eff_isotk_njet_data->GetY();
    for (unsigned ijet_tk = 0; ijet_tk < 2; ijet_tk++) { //0 means no track, 1 means track
      for (unsigned ijet_b = 0; ijet_b < 2; ijet_b++) { //0 means non-b, 1 means b-jet
        for (unsigned ijet_pt = 0; ijet_pt < (jetpt_bins.size()-1); ijet_pt++) {
          float sf = 1.0;
          if (ijet_b==0) {
            if (ijet_tk==0) {
              //scale down
              sf *= nonb_data_frac[ijet_pt]/nonb_mc_frac[ijet_pt];
            }
            else {
              //scale up
              sf *= nonb_mc_frac[ijet_pt]/(1.0-nonb_mc_frac[ijet_pt])*(1.0-nonb_data_frac[ijet_pt]/nonb_mc_frac[ijet_pt])+1.0;
            }
          }
          else {
            if (ijet_tk==0) {
              //scale down
              sf *= b_data_frac[ijet_pt]/b_mc_frac[ijet_pt];
            }
            else {
              //scale up
              sf *= b_mc_frac[ijet_pt]/(1.0-b_mc_frac[ijet_pt])*(1.0-b_data_frac[ijet_pt]/b_mc_frac[ijet_pt])+1.0;
            }
          }
          if (first) {
            std::cout << "    if (";
            first = false;
          }
          else {
            std::cout << "    else if (";
          }
          if (ijet_tk==0) std::cout << "(jet_has_isotk_[ijet]<0.5)";
          else std::cout << "(jet_has_isotk_[ijet]>0.5)";
          if (ijet_b==0) std::cout << "&&(!is_b)";
          else std::cout << "&&(is_b)";
          std::cout << "&&(pt>=" << jetpt_bins[ijet_pt];
          if (ijet_pt != (jetpt_bins.size()-2)) std::cout << "&&pt<" << jetpt_bins[ijet_pt+1] << ")) ";
          else std::cout << "&&pt<9999.0)) ";
          std::cout << "sf *= " << sf << ";" << std::endl;
        }
      }
    }
    for (unsigned inj = 2; inj <= 5; inj++) {
      if (inj==2) std::cout << "  if (b.njet() <= 2) {" << std::endl;
      else if (inj==3) std::cout << "  } else if (b.njet() == 3) {" << std::endl;
      else if (inj==4) std::cout << "  } else if (b.njet() == 4) {" << std::endl;
      else std::cout << "  } else {" << std::endl;
      std::cout << "    if (nonjet_isotk) {" << std::endl;
      std::cout << "      sf *= " << (njet_mc_frac[inj-2]/(1.0-njet_mc_frac[inj-2])*(1.0-njet_data_frac[inj-2]/njet_mc_frac[inj-2])+1.0) << ";" << std::endl;
      std::cout << "    } else {" << std::endl;
      std::cout << "      sf *= " << njet_data_frac[inj-2]/njet_mc_frac[inj-2] << ";" << std::endl;
      std::cout << "    }" << std::endl;;
      
    }
    std::cout << "  }" << std::endl;
    std::cout << "  return sf;" << std::endl;

    //eff_isotk_tt2l_data->Write("isotkeff_tt2lcr_data");
    //eff_isotk_tt2l_mc->Write("isotkeff_tt2lcr_mc");

    //std::cout << "Printing code for scale factors" << std::endl;
    //for (unsigned inj = 0; inj < 4; inj++) {
    //  if (inj==0)
    //    std::cout << "      if (b.njet() <= 2) {" << std::endl;
    //  else if (inj==1)
    //    std::cout << "      else if (b.njet() == 3) {" << std::endl;
    //  else if (inj==2)
    //    std::cout << "      else if (b.njet() == 4) {" << std::endl;
    //  else
    //    std::cout << "      else {" << std::endl;
    //  double * mc_frac = mc_effplots[inj]->GetY();
    //  double * data_frac = data_effplots[inj]->GetY();
    //  mc_effplots[inj]->Write(("isotkeff_tt2lcr_mc_nj"+std::to_string(inj+2)).c_str());
    //  data_effplots[inj]->Write(("isotkeff_tt2lcr_data_nj"+std::to_string(inj+2)).c_str());
    //  for (unsigned iavpt = 0; iavpt < (avjetpt_bins.size()-1); iavpt++) {
    //    float pass_sf = data_frac[iavpt]/mc_frac[iavpt];
    //    float fail_sf = mc_frac[iavpt]/(1.0-mc_frac[iavpt])*(1.0-pass_sf)+1.0;
    //    if (iavpt == 0)
    //      std::cout << "        if (avr_jet_pt >= " << avjetpt_bins[iavpt] << " && avr_jet_pt < " << avjetpt_bins[iavpt+1] << ") {" << std::endl;
    //    else if (iavpt != (avjetpt_bins.size()-1))
    //      std::cout << "        else if (avr_jet_pt >= " << avjetpt_bins[iavpt] << " && avr_jet_pt < " << avjetpt_bins[iavpt+1] << ") {" << std::endl;
    //    else 
    //      std::cout << "        else {" << std::endl;
    //    std::cout << "          if (b.ntk()>0)" << std::endl;
    //    std::cout << "            return " << fail_sf << ";" << std::endl;
    //    std::cout << "          return " << pass_sf << ";" << std::endl;
    //    std::cout << "        }" << std::endl;
    //  }
    //  std::cout << "      }" << std::endl;
    //}
    //EfficiencyPlot * eff_isotk_tt2lj = static_cast<EfficiencyPlot*>(pm.GetFigure(eff_plot_name).get());
    //TGraphAsymmErrors* eff_isotk_tt2lj_data = static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lj->data_ratio_plots_[0]).get()->Clone());
    //TGraphAsymmErrors* eff_isotk_tt2lj_mc = static_cast<TGraphAsymmErrors*>((eff_isotk_tt2lj->background_ratio_plot_).get()->Clone());
    //eff_isotk_tt2lj_data->Write("isotkeff_tt2ljcr_data");
    //eff_isotk_tt2lj_mc->Write("isotkeff_tt2ljcr_mc");
    isotkeff_file->Write();

  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

