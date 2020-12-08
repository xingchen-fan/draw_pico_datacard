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
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

const NamedFunc old_lepton_trig("old_lepton_trig",  [](const Baby &b) -> NamedFunc::ScalarType{
  return b.HLT_Ele27_WPTight_Gsf() || b.HLT_Ele35_WPTight_Gsf() || b.HLT_Ele115_CaloIdVT_GsfTrkIdT() || b.HLT_IsoMu24() || b.HLT_IsoMu27() || b.HLT_Mu50();
});

namespace{
  bool single_thread = false;
  bool do_prefire = false;
  string year_string = "2016";
  bool run_f = false;
  bool unblind = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  //------------------------------------------------------------------------------------
  //                                 plot opts
  //------------------------------------------------------------------------------------

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_norm_nooverflow = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Overflow(OverflowType::none);
  PlotOpt lin_norm_nooverflow_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Overflow(OverflowType::none).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_nooverflow = {lin_norm_nooverflow};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_lin_nooverflow = {lin_norm_nooverflow_data};
  if (unblind) plt_log = {log_norm_data};
  vector<PlotOpt> plt_lin_mc = {lin_norm};
  vector<PlotOpt> plt_log_mc = {log_norm};

  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).YAxis(YAxisType::linear).Overflow(OverflowType::overflow)};
  vector<PlotOpt> twodim_plotopts_log = {style2D().Title(TitleType::info).YAxis(YAxisType::log).Overflow(OverflowType::overflow)};

  //------------------------------------------------------------------------------------
  //                                 samples
  //------------------------------------------------------------------------------------

  // Set options
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string mc_unskimmed_folder = "mc/unskimmed/";
  //string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string ttbar_mc_loose_skim_folder = "mc/skim_1l2j/";
  //string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  //string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";
  string ttbar_mc_skim_folder = "mc/skim_higlep1T/";
  string zll_mc_skim_folder = "mc/skim_higlep2T/";
  string qcd_mc_skim_folder = "mc/skim_met150/";
  string met150_mc_skim_folder = "mc/skim_met150/";

  string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  //string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  string ttbar_data_loose_skim_folder = "data/skim_1l2j/";
  //string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  //string qcd_data_skim_folder = "data/merged_higdata_higqcd/";
  string ttbar_data_skim_folder = "data/skim_higlep1T/";
  string zll_data_skim_folder = "data/skim_higlep2T/";
  string qcd_data_skim_folder = "data/skim_met150/";
  string met150_data_skim_folder = "data/skim_met150/";

  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string search_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
  //string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";
  string ttbar_sig_loose_skim_folder = "SMS-TChiHH_2D/skim_1l2j/";
  //string zll_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep2T/";
  //string qcd_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higqcd/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D/skim_higlep1T/";
  string zll_sig_skim_folder = "SMS-TChiHH_2D/skim_higlep2T/";
  string qcd_sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";
  string met150_sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  //years = {2016, 2017, 2018};
  //years = {2016};
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
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

  // set qcd procs
  vector<shared_ptr<Process> > qcd_procs;
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["zjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["wjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["single_t"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["qcd"]),"stitch")); 
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    qcd_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
		    (met_trigger)));
  }

  vector<shared_ptr<Process> > qcd_procs_data_as_mc;
  if (!run_f) {
    qcd_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                      attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}), (met_trigger)));
  }
  else {
    qcd_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                      attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}), (met_trigger && "run>=305044")));
  }

  // set ttbar procs
  vector<shared_ptr<Process> > ttbar_procs;
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
		    (el_trigger || mu_trigger || met_trigger)));
  }

  vector<shared_ptr<Process> > ttbar_procs_loose;
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder, mctags["tt"]),"stitch"));
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder,mctags["zjets"]),"stitch"));
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder,mctags["wjets"]),"stitch"));
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder, mctags["single_t"]),"stitch"));
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder, mctags["qcd"]),"stitch")); 
  ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_unskimmed_folder, mctags["other"]),"stitch"));
  if (unblind) {
    ttbar_procs_loose.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_loose_skim_folder, {"*.root"}),
		    (el_trigger || mu_trigger || met_trigger)));
  }

  vector<shared_ptr<Process> > ttbar_procs_data_as_mc;
  if (!run_f) {
    ttbar_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
		    (el_trigger || mu_trigger || met_trigger)));
  }
  else {
    ttbar_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
		    (el_trigger || mu_trigger || met_trigger) && "run>=305044"));
  }
  
  // set zll procs
  vector<shared_ptr<Process> > zll_procs;
  zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["zjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["wjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["single_t"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["qcd"]),"stitch")); 
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    zll_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    (el_trigger || mu_trigger)));
  }

  // set sr procs
  vector<shared_ptr<Process> > sr_procs;
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder,mctags["zjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder,mctags["wjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["single_t"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["qcd"]),"stitch")); 
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["other"]),"stitch"));

  vector<shared_ptr<Process> > sr_procs_data_as_mc;
  if (!run_f) {
    sr_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, met150_data_skim_folder, {"*.root"}),
		    (met_trigger)));
  }
  else {
    sr_procs_data_as_mc.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, met150_data_skim_folder, {"*.root"}),
		    (met_trigger) && "run>=305044"));
  }

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  const NamedFunc weight_nopref("weight_nopref",[](const Baby &b) -> NamedFunc::ScalarType{
    //data
    if (b.SampleType()<0) return Higfuncs::final_weight.GetScalar(b);
    //mc
    if (b.w_prefire()==0.) {
      return b.w_lumi()*Higfuncs::eff_higtrig_run2.GetScalar(b)*Higfuncs::w_years.GetScalar(b);
    }
    else {
      return (Higfuncs::final_weight.GetScalar(b))/(b.w_prefire());
    }
  });

  const NamedFunc weight_eraef("weight_eraef",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return Higfuncs::final_weight.GetScalar(b);
    return Higfuncs::final_weight.GetScalar(b)*(22.8/41.5);
  });

  const NamedFunc weight_nopref_eraef("weight_nopref_eraef",[](const Baby &b) -> NamedFunc::ScalarType{
    //data
    if (b.SampleType()<0) return Higfuncs::final_weight.GetScalar(b);
    //mc
    if (b.w_prefire()==0.) {
      return b.w_lumi()*Higfuncs::eff_higtrig_run2.GetScalar(b)*Higfuncs::w_years.GetScalar(b)*(22.8/41.5);
    }
    else {
      return (Higfuncs::final_weight.GetScalar(b))/(b.w_prefire())*(22.8/41.5);
    }
  });

  const NamedFunc pass_nohem("pass_nohem",[](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //was modified
    if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    return true;
  });

  const NamedFunc pass_noecalnoisejet("pass_noecalnoisejet",[](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //was modified
    //if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if (!Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  const NamedFunc pass_nohtratio("pass_nohtratio",[](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    //if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //was modified
    if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if (!Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  const NamedFunc pass_nohtratio_noecalnoisejet("pass_nohtratio_noecalnoisejet",[](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    //if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //was modified
    //if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if (!Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  const NamedFunc ht_ratio("ht_ratio",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ht5()/b.ht();
  });
 

  NamedFunc sr_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && !jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc badmet_baseline = pass_filters && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb<=2) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc badmet_baseline_noecalnoisejet = pass_noecalnoisejet && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb<=2) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc badmet_baseline_nohtratio_noecalnoisejet = pass_nohtratio_noecalnoisejet && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb<=2) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc badmet_baseline_nohtratio = pass_nohtratio && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb<=2) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline_nohig = pass_filters && "met>150&&njet>=1&&nlep==1&&nbm>=2" &&
    (jetid_njet>=3) && (lead_signal_lepton_pt>30);

  NamedFunc ttbar_baseline_nohem = pass_nohem && "nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline_nohig_noecalnoisejet = pass_noecalnoisejet && "met>150&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30);

  NamedFunc ttbar_baseline_nohig_nohtratio = pass_nohtratio && "met>150&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30);

  NamedFunc ttbar_baseline_nohig_nohtratio_noecalnoisejet = pass_nohtratio_noecalnoisejet && "met>150&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30);

  NamedFunc zll_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2&&ll_m[0]>80&&ll_m[0]<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);
  
  NamedFunc qcd_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);
  
  NamedFunc qcd_baseline_nohig = pass_filters && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met;

  NamedFunc qcd_baseline_nohem = pass_nohem && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc qcd_baseline_nohig_noecalnoisejet = pass_noecalnoisejet && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met;

  NamedFunc qcd_baseline_nohig_nohtratio = pass_nohtratio && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met;

  NamedFunc qcd_baseline_nohig_nohtratio_noecalnoisejet = pass_nohtratio_noecalnoisejet && "met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met;


  //------------------------------------------------------------------------------------
  //                                     plots
  //------------------------------------------------------------------------------------

  PlotMaker pm;

  //prefire weights plots
  if (year_string != "2018" && do_prefire) {
    pm.Push<Hist1D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}),
      ttbar_baseline_nohig && "jet_pt[0]>100&&jet_isgood[0]",
      ttbar_procs_loose, plt_lin).Weight(weight_nopref).Tag("FixName:dataquality_jeteta_noprefireweights_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}),
      ttbar_baseline_nohig && "jet_pt[0]>100&&jet_isgood[0]",
      ttbar_procs_loose, plt_lin).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_jeteta_"+year_string).LuminosityTag(total_luminosity_string);
    if (year_string == "2017") {
      pm.Push<Hist1D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}),
        ttbar_baseline_nohig && "jet_pt[0]>100&&jet_isgood[0]&&(run>303825||type/1000>=1)",
        ttbar_procs_loose, plt_lin).Weight(weight_nopref_eraef).Tag("FixName:dataquality_jeteta_noprefireweights_2017ef_"+year_string).LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}),
        ttbar_baseline_nohig && "jet_pt[0]>100&&jet_isgood[0]&&(run>303825||type/1000>=1)",
        ttbar_procs_loose, plt_lin).Weight(weight_eraef).Tag("FixName:dataquality_jeteta_2017ef_"+year_string).LuminosityTag(total_luminosity_string);
    }
  }
  //HEM veto plots
  if (year_string == "2018") {
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}), Axis(40,-3.14,3.14,"jet_phi[0]","#phi_{j1}",{}),
      ttbar_baseline && "jet_isgood[0]",
      ttbar_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_1lcr_jet_phieta_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}), Axis(40,-3.14,3.14,"jet_phi[0]","#phi_{j1}",{}),
      ttbar_baseline_nohem && "jet_isgood[0]",
      ttbar_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_1lcr_jet_phieta_nohem_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "el_eta[0]", "#eta_{e}", {}), Axis(40,-3.14,3.14,"el_phi[0]","#phi_{e}",{}),
      ttbar_baseline && "nel==1&&el_sig[0]",
      ttbar_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_1lcr_el_phieta_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "el_eta[0]", "#eta_{e}", {}), Axis(40,-3.14,3.14,"el_phi[0]","#phi_{e}",{}),
      ttbar_baseline_nohem && "nel==1&&el_sig[0]",
      ttbar_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_1lcr_el_phieta_nohem_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}), Axis(40,-3.14,3.14,"jet_phi[0]","#phi_{j1}",{}),
      qcd_baseline && "jet_isgood[0]",
      qcd_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_0lcr_jet_phieta_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(Axis(40, -2.5, 2.5, "jet_eta[0]", "#eta_{j1}", {}), Axis(40,-3.14,3.14,"jet_phi[0]","#phi_{j1}",{}),
      qcd_baseline_nohem && "jet_isgood[0]",
      qcd_procs_data_as_mc, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_0lcr_jet_phieta_nohem_"+year_string).LuminosityTag(total_luminosity_string);
  }
  //ecal noise jet filter plots
  std::vector<NamedFunc> ttbar_filter_combinations = {ttbar_baseline_nohig, ttbar_baseline_nohig_noecalnoisejet, ttbar_baseline_nohig_nohtratio, ttbar_baseline_nohig_nohtratio_noecalnoisejet};
  std::vector<NamedFunc> qcd_filter_combinations = {qcd_baseline_nohig, qcd_baseline_nohig_noecalnoisejet, qcd_baseline_nohig_nohtratio, qcd_baseline_nohig_nohtratio_noecalnoisejet};
  std::vector<NamedFunc> lowb_filter_combinations = {badmet_baseline, badmet_baseline_noecalnoisejet, badmet_baseline_nohtratio, badmet_baseline_nohtratio_noecalnoisejet};
  std::vector<std::string> filter_combinations_names = {"nohig","nohig_noecalnoise","nohig_nohtratio","nohig_nohtratio_noecalnoise"};
  std::string run_string = "";
  if (run_f) run_string = "RunF";
  for (unsigned int region_idx = 0; region_idx < 3; region_idx++) {
    vector<shared_ptr<Process>> *procs = &ttbar_procs_data_as_mc;
    std::string procs_name = "1lcr";
    if (region_idx==1) {
      procs = &qcd_procs_data_as_mc;
      procs_name = "0lcr";
    }
    if (region_idx==2) {
      procs = &sr_procs_data_as_mc;
      procs_name = "0bcr";
    }
    if (region_idx!=2) continue; //temp only fake met region

    for (unsigned int endcap_idx = 0; endcap_idx < 2; endcap_idx++) {
      NamedFunc endcap_condition = "1";
      std::string endcap_name = "fulldetector";
      if (endcap_idx==1) {
        endcap_condition = "((jet_eta[0]>2.4&&jet_eta[0]<5.0)||(jet_eta[0]<-2.4&&jet_eta[0]>-5.0))";
        endcap_name = "endcap";
      }

      for (unsigned int filter_idx = 0; filter_idx < ttbar_filter_combinations.size(); filter_idx++) {
	NamedFunc filter_combinations = ttbar_filter_combinations[filter_idx];
        if (region_idx==1) filter_combinations = qcd_filter_combinations[filter_idx];
        if (region_idx==2) filter_combinations = lowb_filter_combinations[filter_idx];
	std::string filter_name = filter_combinations_names[filter_idx];

        pm.Push<Hist2D>(Axis(50, 0, 3.14, "jet_met_dphi[0]", "#Delta #phi", {}), Axis(50,30.,900.,"jet_pt[0]","p_{T} [GeV]",{}),
          filter_combinations && "jet_pt[0]>30" && endcap_condition,
          *procs, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_"+procs_name+"_"+endcap_name+"_dphipt_"+filter_name+"_"+year_string+run_string).LuminosityTag(total_luminosity_string);

        pm.Push<Hist2D>(Axis(50,0.5,4.,ht_ratio,"H_{T5}/H_{T}",{}), Axis(50, 0, 3.14, "jet_met_dphi[0]", "#Delta #phi_{j1}", {}), 
          filter_combinations && "jet_pt[0]>30" && endcap_condition,
          *procs, twodim_plotopts_log).Weight(Higfuncs::final_weight).Tag("FixName:dataquality_"+procs_name+"_"+endcap_name+"_htratiodphi_"+filter_name+"_"+year_string+run_string).LuminosityTag(total_luminosity_string);

      }
    }
  }


  
  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"year", required_argument, 0, 0},
      {"unblind", no_argument, 0, 0},
      {"doprefire", no_argument, 0, 0},
      {"runf", no_argument, 0, 0},
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
      if(optname == "year"){
        year_string = optarg;
      } else if (optname == "unblind") {
        unblind = true;
      } else if (optname == "runf") {
        run_f = true;
      } else if (optname == "doprefire") {
        do_prefire = true;
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
