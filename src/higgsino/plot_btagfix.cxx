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
//

#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TFile.h"
#include "TLatex.h"
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
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt lin_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  //PlotOpt lin_norm_paper_ratio = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt log_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt log_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::info);
  PlotOpt log_shapes_paper = log_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::info);
  PlotOpt lin_shapes_info_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_paper_data = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::supplementary);

  vector<PlotOpt> plt_lin = {lin_norm_paper};
  vector<PlotOpt> plt_shapes = {lin_shapes_paper};
  vector<PlotOpt> plt_log_shapes = {log_shapes_paper};
  vector<PlotOpt> plt_lin_ratio = {lin_norm_paper_data};
  vector<PlotOpt> plt_log = {log_norm_paper};

  set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  int lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);

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
  string klamath_v2_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v2/";
  string klamath_v3_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string sig_unskimmed_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";
  string mc_klamath_v2_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath/";
  string mc_klamath_v3_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath_v3/";

  std::vector<std::shared_ptr<Process>> search_procs;
  std::shared_ptr<Process> proc_ttbar = Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch");
  std::shared_ptr<Process> proc_zjets = Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch");
  std::shared_ptr<Process> proc_wjets = Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch");
  std::shared_ptr<Process> proc_qcd = Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch"); 
  std::shared_ptr<Process> proc_other = Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other_and_single_t"]),"stitch");

  std::vector<std::shared_ptr<Process>> bkg_procs;
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (2b) (pre-fix)", Process::Type::signal, kGreen,
              attach_folder(mc_klamath_v2_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==2)&&"stitch"));
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (3b) (pre-fix)", Process::Type::signal, kRed,
              attach_folder(mc_klamath_v2_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (4b) (pre-fix)", Process::Type::signal, kBlue,
              attach_folder(mc_klamath_v2_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==4)&&"stitch"));
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (2b) (post-fix)", Process::Type::signal, kYellow+1,
              attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==2)&&"stitch"));
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (3b) (post-fix)", Process::Type::signal, kViolet,
              attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  bkg_procs.push_back(Process::MakeShared<Baby_pico>("Background MC (4b) (post-fix)", Process::Type::signal, kCyan,
              attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["all"]),(Higfuncs::hig_bcat==4)&&"stitch"));
  for (auto &proc : bkg_procs) {
    proc->SetLineWidth(1);
  }

  std::vector<std::shared_ptr<Process>> threeb_procs;
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (3b) (post-fix)", Process::Type::signal, colors("tt_htau"),
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["tt"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("Single t (3b) (post-fix)", Process::Type::signal, colors("single_t"),
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["single_t"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets (3b) (post-fix)", Process::Type::signal, colors("zjets"),
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["zjets"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("W+jets (3b) (post-fix)", Process::Type::signal, colors("wjets"),
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["wjets"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("QCD (3b) (post-fix)", Process::Type::signal, colors("other"),
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["qcd"]),(Higfuncs::hig_bcat==3)&&"stitch"));
  threeb_procs.push_back(Process::MakeShared<Baby_pico>("Other (3b) (post-fix)", Process::Type::signal, kGray+2,
           attach_folder(mc_klamath_v3_folder, years, mc_skim_folder, mctags["other"]),(Higfuncs::hig_bcat==3)&&"stitch"));


  search_procs.push_back(proc_ttbar);
  search_procs.push_back(proc_zjets);
  search_procs.push_back(proc_wjets);
  search_procs.push_back(proc_qcd);
  search_procs.push_back(proc_other);

  std::vector<std::string> proc_names = {"t#bar{t}+X","Z+jets","W+jets","QCD","Other"};
  std::vector<Int_t> proc_colors = {TColor::GetColor(102,102,255),kOrange+1,kGreen+1,TColor::GetColor(255,200,0),kGray+2};

  //std::vector<std::shared_ptr<Process>> signal450_procs;
  //std::shared_ptr<Process> proc_other = Process::MakeShared<Baby_pico>("TChiHH-G(450,1)", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-450_mLSP-0*.root"}),"stitch");
  
  std::vector<std::shared_ptr<Process>> cutflow_procs;
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(175,1)", Process::Type::signal, kGray+2,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-175_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)", Process::Type::signal, kGray+2,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(900,1)", Process::Type::signal, kGray+2,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-900_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(350,200)", Process::Type::signal, kGray+2,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-350_mLSP-200*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(450,100)", Process::Type::signal, kGray+2,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-450_mLSP-100*.root"}),"stitch"));

  std::vector<std::string> cutflow_proc_names = {"TChiHH-G(175,1)","TChiHH-G(500,1)","TChiHH-G(900,1)","TChiHH(350,200)","TChiHH(450,100)"};
  std::vector<std::string> cutflow_proc_outnames = {"175_1","500_1","900_1","350_200","450_100"};
  std::vector<std::string> cutflow_proc_filenames = {"*TChiHH_mChi-175_mLSP-0*.root","*TChiHH_mChi-500_mLSP-0*.root","*TChiHH_mChi-900_mLSP-0*.root","*TChiHH_mChi-350_mLSP-200*.root","*TChiHH_mChi-450_mLSP-100*.root"};

  std::vector<std::vector<std::shared_ptr<Process>>> pair_compare_procs;
  for (unsigned signal_idx = 0; signal_idx < cutflow_proc_names.size(); signal_idx++) {
    pair_compare_procs.push_back(std::vector<std::shared_ptr<Process>>());
    pair_compare_procs[signal_idx].push_back(Process::MakeShared<Baby_pico>(cutflow_proc_names[signal_idx]+"(pre-fix)", Process::Type::signal, kRed,
        attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {cutflow_proc_filenames[signal_idx]}),"stitch"));
    pair_compare_procs[signal_idx].push_back(Process::MakeShared<Baby_pico>(cutflow_proc_names[signal_idx]+"(post-fix)", Process::Type::signal, kBlue,
        attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {cutflow_proc_filenames[signal_idx]}),"stitch"));
  }

  std::vector<std::shared_ptr<Process>> bcat_procs;
  bcat_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)(2b)", Process::Type::signal, kRed,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),(Higfuncs::hig_bcat==2)));
  bcat_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)(3b)", Process::Type::signal, kGreen,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),(Higfuncs::hig_bcat==3)));
  bcat_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)(4b)", Process::Type::signal, kBlue,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),(Higfuncs::hig_bcat==4)));

  std::vector<std::shared_ptr<Process>> bdist_procs;
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(175,1)(pre-fix)", Process::Type::signal, kGreen,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-175_mLSP-0*.root"}),"stitch"));
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)(pre-fix)", Process::Type::signal, kRed,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),"stitch"));
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(900,1)(pre-fix)", Process::Type::signal, kBlue,
                  attach_folder(klamath_v2_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-900_mLSP-0*.root"}),"stitch"));
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(175,1)(post-fix)", Process::Type::signal, kYellow+1,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-175_mLSP-0*.root"}),"stitch"));
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)(post-fix)", Process::Type::signal, kViolet,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),"stitch"));
  bdist_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(900,1)(post-fix)", Process::Type::signal, kCyan,
                  attach_folder(klamath_v3_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-900_mLSP-0*.root"}),"stitch"));


  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc filters = Higfuncs::final_pass_filters;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc sr_baseline = filters&&script_utilities::search_resolved;

  NamedFunc metnbbin("metnbbin",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbm()<3) {
      if (b.met() < 200) return 1;
      else return 2;
    }
    else if (b.nbm()==3 && b.nbl()<4) {
      if (b.met() < 200) return 3;
      else return 4;
    }
    if (b.met() < 200) return 5;
    return 6;
  });

  NamedFunc drmaxnbbin("drmaxnbbin",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbm()<3) {
      if (b.hig_cand_drmax()->at(0) > 1.1) return 1;
      else return 2;
    }
    else if (b.nbm()==3 && b.nbl()<4) {
      if (b.hig_cand_drmax()->at(0) > 1.1) return 3;
      else return 4;
    }
    if (b.hig_cand_drmax()->at(0) > 1.1) return 5;
    return 6;
  });

  NamedFunc effplot_bin_x("effplot_bin_x",[](const Baby &b) -> NamedFunc::ScalarType{
      double hig_cand_am = b.hig_cand_am()->at(0);
      bool hig_window = (hig_cand_am > 100 && hig_cand_am < 140);
      if (b.nbm()<3) {
        if (!hig_window) return 2;
        return 3;
      }
      else if (b.nbm()==3 && b.nbl()<4) {
        if (!hig_window) return 4;
        return 5;
      }
      if (!hig_window) return 6;
      return 7;
  });

  NamedFunc effplot_bin_y("effplot_bin_y",[](const Baby &b) -> NamedFunc::ScalarType{
      double bin = 4;
      if (b.met() < 200) {
        bin = 1;
      }
      else if (b.met() < 300) {
        bin = 2;
      }
      else if (b.met() < 400) {
        bin = 3;
      }
      if (b.hig_cand_drmax()->at(0) < 1.1) {
        bin += 4;
      }
      return bin;
  });

  //NamedFunc one("one",[](const Baby &b) -> NamedFunc::ScalarType{
  //    return 1;
  //});

  NamedFunc mixed_model_weight = script_utilities::mixed_model_weight;

  NamedFunc nbm_exclusive = "nbm-nbt";
  NamedFunc nbl_exclusive = "nbl-nbm";


  NamedFunc mixed_model_weight_lumionly("mixed_model_weight_lumionly",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type() / 1000 != 106) return b.w_lumi()*Higfuncs::w_years_search.GetScalar(b); //not TChiHH
    if (b.mlsp() <= 1) return b.w_lumi()*Higfuncs::w_years_search.GetScalar(b); //N1N2
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return b.w_lumi()*Higfuncs::w_years_search.GetScalar(b)/xsec1d*xsec2d;
  });

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;
  //if(HigUtilities::is_in_string_options(options.string_options, "pie")) {

  //for (unsigned signal_idx = 0; signal_idx < cutflow_proc_names.size(); signal_idx++) {
  //  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {}),
  //      "1",
  //      pair_compare_procs[signal_idx], plt_lin_ratio).Weight(weight)
  //      .Tag(("FixName:bfix__nb_comp_"+cutflow_proc_outnames[signal_idx]).c_str())
  //      .LuminosityTag(total_luminosity_string).RatioTitle("post-fix","pre-fix");
  //  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {}),
  //      sr_baseline,
  //      pair_compare_procs[signal_idx], plt_lin_ratio).Weight(weight)
  //      .Tag(("FixName:bfix__nb_comp_baseline_"+cutflow_proc_outnames[signal_idx]).c_str())
  //      .LuminosityTag(total_luminosity_string).RatioTitle("post-fix","pre-fix");
  //}
  //pm.Push<Hist1D>(Axis(4, -0.5, 3.5, nbm_exclusive, "N_{bm}", {}),
  //    sr_baseline,
  //    bcat_procs, plt_lin).Weight("1")
  //    .Tag("FixName:bfix__nbm_bcat")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(4, -0.5, 3.5, nbl_exclusive, "N_{bl}", {}),
  //    sr_baseline,
  //    bcat_procs, plt_lin).Weight("1")
  //    .Tag("FixName:bfix__nbl_bcat")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(14, 150, 850, "met", "p_{T}^{miss} [GeV]", {}),
  //    sr_baseline,
  //    bdist_procs, plt_log_shapes).Weight(weight)
  //    .Tag("FixName:bfix__dist_met")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 100, "hig_cand_dm[0]", "#Delta m_{bb} [GeV]", {}),
  //    filters && "met/mht<2 && met/met_calo<2&&weight<1.5&&"
  //    "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&"
  //    "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200",
  //    bdist_procs, plt_shapes).Weight(weight)
  //    .Tag("FixName:bfix__dist_dm")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {}),
  //    filters && "met/mht<2 && met/met_calo<2&&weight<1.5&&"
  //    "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&nbt>=2&&"
  //    "hig_cand_am[0]<200&&hig_cand_dm[0]<40",
  //    bdist_procs, plt_log_shapes).Weight(weight)
  //    .Tag("FixName:bfix__dist_drmax")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {}),
  //    sr_baseline,
  //    bdist_procs, plt_log_shapes).Weight(weight)
  //    .Tag("FixName:bfix__dist_nb")
  //    .LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {}),
  //    sr_baseline,
  //    bdist_procs, plt_shapes).Weight(weight)
  //    .Tag("FixName:bfix__dist_am")
  //    .LuminosityTag(total_luminosity_string);

  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100,140}),
      sr_baseline,
      bkg_procs, plt_shapes).Weight(weight)
      .Tag("FixName:bfix__dist_am_bkg")
      .LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100,140}),
      sr_baseline,
      threeb_procs, plt_shapes).Weight(weight)
      .Tag("FixName:bfix__dist_am_3b")
      .LuminosityTag(total_luminosity_string);

  //  pm.Push<Table>("sr_cutflow_"+options.year_string, vector<TableRow>{
  //    TableRow("No cuts(lumi and cross section $\\PH\\to\\PQb\\PAQb$ only)", 
  //        "1",0,0,mixed_model_weight_lumionly),
  //    TableRow("Filters, $p_\\text{T}^\\text{miss}>150 \\text{ GeV}$", 
  //        Higfuncs::pass_filters && "met>150",0,0,mixed_model_weight),
  //    TableRow("$N_\\text{vl}=N_\\text{tk}=0$", 
  //        Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0",0,0,mixed_model_weight),
  //    TableRow("$4 \\leq N_\\text{jet} \\leq 5$", 
  //        Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5",
  //        0,0,mixed_model_weight),
  //    TableRow("$N_\\text{b}\\geq 2$", 
  //        Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2",
  //        0,0,mixed_model_weight),
  //    TableRow("$p_\\text{T}^\\text{miss}$ quality and $\\Delta\\phi$ cuts",        
  //        Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met",
  //        0,0,mixed_model_weight),
  //    TableRow("$\\Delta m_\\text{bb}<40 \\text{ GeV}$, $\\langle m_\\text{bb} \\rangle<200 \\text{ GeV}$",        
  //        Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40&&hig_cand_am[0]<200",
  //        0,0,mixed_model_weight),
  //    TableRow("$\\Delta R_\\text{max}<2.2$",        
  //        sr_baseline,0,0,mixed_model_weight),
  //    TableRow("\\hline\n $100 < \\langle m\\rangle < 140 \\text{ GeV}$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140",0,0,mixed_model_weight),
  //    TableRow("$N_\\text{b}\\geq 3$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3",0,0,mixed_model_weight),
  //    TableRow("$N_\\text{b}=4$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4",
  //        0,0,mixed_model_weight),
  //    TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>200",
  //        0,0,mixed_model_weight),
  //    TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>300",
  //        0,0,mixed_model_weight),
  //    TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
  //        sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>400",
  //        0,0,mixed_model_weight)
  //  },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Tag("cutflow").Precision(1);

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

