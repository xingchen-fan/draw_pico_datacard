//This script generates di-Higgs invariant mass plots
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options (default indicates it is an AN plot): 

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
#include "TLorentzVector.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

using namespace std;

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "plot_variables,plot_nminus1,make_cutflow,make_pies");
  std::vector<PlotOpt> plt_lin = script_utilities::plot_lin(options.unblind);
  std::vector<PlotOpt> plt_log = script_utilities::plot_log(options.unblind);
  std::vector<PlotOpt> plt_shapes = script_utilities::plot_shapes();
  std::vector<PlotOpt> plt_log_shapes = script_utilities::plot_log_shapes();
  set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string);

  std::vector<std::shared_ptr<Process>> search_procs = 
      script_utilities::getall_processes(
      years, {"TChiHH(175,0)","TChiHH(500,0)","TChiHH(900,0)",
      "TChiHH(350,200)","TChiHH(450,100)"}, "search", options.unblind);

  //// Set MC
  //map<string, set<string>> mctags;
  //// Set base tags
  //mctags["tt"]     = set<string>({"*TTJets_*Lept*",
  //                                "*_TTZ*.root", "*_TTW*.root",
  //                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  //mctags["single_t"] = set<string>({"*_ST_*.root"});
  //mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  //mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  //mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  //mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
  //                                 "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
  //                                 "*_QCD_HT2000toInf_*"});
  //mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
  //                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  //// Combine all tags
  //mctags["all"] = set<string>({"*TTJets_SingleLept*",
  //                             "*TTJets_DiLept*",
  //                             "*_TTZ*.root", "*_TTW*.root",
  //                             "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
  //                             "*_WJetsToLNu*.root", "*_ZJet*.root",
  //                             "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
  //                             "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
  //                             "*_QCD_HT2000toInf_*",
  //                             "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
  //                             "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  //});

  //vector<shared_ptr<Process> > search_procs;
  //// Set mc processes
  ////search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  ////                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  ////search_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  ////                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, data_skim_folder, {"*.root"}),"stitch"));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  NamedFunc dihiggs_mass("dihiggs_mass",[](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<TLorentzVector> jet_momentum; 
      for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
        if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) {
          TLorentzVector jet_p;
          jet_p.SetPtEtaPhiM(b.jet_pt()->at(jet_idx),b.jet_eta()->at(jet_idx),b.jet_phi()->at(jet_idx),
              b.jet_m()->at(jet_idx));
          jet_momentum.push_back(jet_p);
        }
      }
      if (jet_momentum.size() != 4) return 0;
      TLorentzVector hh_p = jet_momentum[0] + jet_momentum[1] + jet_momentum[2] + jet_momentum[3];
      return hh_p.M();
  });
  
  NamedFunc filters = Higfuncs::final_pass_filters;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc sr_baseline = filters&&script_utilities::search_resolved;
  NamedFunc sr_baseline_sr = filters&&script_utilities::search_resolved &&
    "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3";
  NamedFunc mixed_model_weight = script_utilities::mixed_model_weight;
  NamedFunc hig_nb = script_utilities::hig_nb;

  //------------------------------------------------------------------------------------
  //                                     make plots
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  // 2b3b4b 
  //pm.Push<Table>("FixName:selection__search_pies__2b3b4b_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved, 0, 0, weight)}), search_procs, true, true, true);
  
  //if (HigUtilities::is_in_string_options(options.string_options,"plot_variables")) {
  pm.Push<Hist1D>(Axis(70, 100, 1500, dihiggs_mass, "m_{HH} [GeV]", {}),
    sr_baseline,
    search_procs, plt_lin).Weight(mixed_model_weight)
    .Tag("FixName:hhmass__basline_"+options.year_string)
    .LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(70, 100, 1500, dihiggs_mass, "m_{HH} [GeV]", {}),
    sr_baseline_sr,
    search_procs, plt_lin).Weight(mixed_model_weight)
    .Tag("FixName:hhmass__sr_"+options.year_string)
    .LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(70, 100, 1500, dihiggs_mass, "m_{HH} [GeV]", {}),
    sr_baseline_sr && "hig_cand_drmax[0]<1.1",
    search_procs, plt_lin).Weight(mixed_model_weight)
    .Tag("FixName:hhmass__sr_lowdrmax_"+options.year_string)
    .LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(7, 100, 1500, dihiggs_mass, "m_{HH} [GeV]", {}),
    sr_baseline_sr && "nbm==3&&nbl==3&&met>=300&&met<400&&hig_cand_drmax[0]<1.1",
    search_procs, plt_lin).Weight(mixed_model_weight)
    .Tag("FixName:hhmass__sr_bin11_"+options.year_string)
    .LuminosityTag(total_luminosity_string);

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

