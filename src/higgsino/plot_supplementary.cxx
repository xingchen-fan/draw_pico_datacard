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
//pie,cutflow

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
  PlotOpt lin_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt lin_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt log_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt log_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_info_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_paper_data = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::supplementary);

  vector<PlotOpt> plt_lin = {lin_norm_paper};
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

  string mc_base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath_v3/";
  //string data_base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string sig_unskimmed_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";

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
                  attach_folder(sig_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-175_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(500,1)", Process::Type::signal, kGray+2,
                  attach_folder(sig_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-500_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH-G(900,1)", Process::Type::signal, kGray+2,
                  attach_folder(sig_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-900_mLSP-0*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(350,200)", Process::Type::signal, kGray+2,
                  attach_folder(sig_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-350_mLSP-200*.root"}),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(450,100)", Process::Type::signal, kGray+2,
                  attach_folder(sig_base_folder, years, sig_unskimmed_folder, {"*TChiHH_mChi-450_mLSP-100*.root"}),"stitch"));

  std::vector<std::string> cutflow_proc_names = {"TChiHH-G(175,1)","TChiHH-G(500,1)","TChiHH-G(900,1)","TChiHH(350,200)","TChiHH(450,100)"};

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
  if(HigUtilities::is_in_string_options(options.string_options, "pie")) {
    pm.Push<Hist1D>(Axis(6, 0.5, 6.5, metnbbin, "MET-Nb category", {}),
        sr_baseline,
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:supppies__metnbcat")
        .LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, 0.5, 6.5, drmaxnbbin, "Drmax-Nb category", {}),
        sr_baseline,
        search_procs, plt_lin).Weight(weight)
        .Tag("FixName:supppies__drmaxnbcat")
        .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist2D>(Axis(7, 0.5, 7.5, effplot_bin_x, "ABCD region", {}),
    //    Axis(9, 0.5, 9.5, effplot_bin_y, "ABCD plane", {}),
    //    sr_baseline,
    //    signal450_procs, plt_lin).Weight(weight)
    //    .Tag("FixName:suppeff_part1")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist2D>(Axis(7, 0.5, 7.5, effplot_bin_x, "ABCD region", {}),
    //    Axis(9, 0.5, 9.5, one, "ABCD plane", {}),
    //    sr_baseline,
    //    signal450_procs, plt_lin).Weight(weight)
    //    .Tag("FixName:suppeff_part2")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist2D>(Axis(7, 0.5, 7.5, one, "ABCD region", {}),
    //    Axis(9, 0.5, 9.5, effplot_bin_y, "ABCD plane", {}),
    //    sr_baseline,
    //    signal450_procs, plt_lin).Weight(weight)
    //    .Tag("FixName:suppeff_part3")
    //    .LuminosityTag(total_luminosity_string);
    //pm.Push<Hist2D>(Axis(7, 0.5, 7.5, one, "ABCD region", {}),
    //    Axis(9, 0.5, 9.5, effplot_bin_y, "ABCD plane", {}),
    //    sr_baseline,
    //    signal450_procs, plt_lin).Weight(weight)
    //    .Tag("FixName:suppeff_part4")
    //    .LuminosityTag(total_luminosity_string);
  }
  if(HigUtilities::is_in_string_options(options.string_options, "cutflow")) {
    pm.Push<Table>("sr_cutflow_"+options.year_string, vector<TableRow>{
      TableRow("No cuts(lumi and cross section $\\PH\\to\\PQb\\PAQb$ only)", 
          "1",0,0,mixed_model_weight_lumionly),
      TableRow("Filters, $p_\\text{T}^\\text{miss}>150 \\text{ GeV}$", 
          Higfuncs::pass_filters && "met>150",0,0,mixed_model_weight),
      TableRow("$N_\\text{vl}=N_\\text{tk}=0$", 
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0",0,0,mixed_model_weight),
      TableRow("$4 \\leq N_\\text{jet} \\leq 5$", 
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5",
          0,0,mixed_model_weight),
      TableRow("$N_\\text{b}\\geq 2$", 
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}$ quality and $\\Delta\\phi$ cuts",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met",
          0,0,mixed_model_weight),
      TableRow("$\\Delta m_\\text{bb}<40 \\text{ GeV}$, $\\langle m_\\text{bb} \\rangle<200 \\text{ GeV}$",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40&&hig_cand_am[0]<200",
          0,0,mixed_model_weight),
      TableRow("$\\Delta R_\\text{max}<2.2$",        
          sr_baseline,0,0,mixed_model_weight),
      TableRow("\\hline\n $100 < \\langle m\\rangle < 140 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140",0,0,mixed_model_weight),
      TableRow("$N_\\text{b}\\geq 3$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3",0,0,mixed_model_weight),
      TableRow("$N_\\text{b}=4$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>200",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>300",
          0,0,mixed_model_weight),
      TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
          sr_baseline && "hig_cand_am[0]>100&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>400",
          0,0,mixed_model_weight)
    },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Tag("cutflow").Precision(1);
    pm.Push<Table>("lep_veto_cutflow_"+options.year_string, vector<TableRow>{
      TableRow("No cuts (weighted)", 
          "1",0,0,mixed_model_weight),
      TableRow("Lepton veto", 
          "nvlep==0",0,0,mixed_model_weight),
      TableRow("Track veto", 
          "ntk==0",0,0,mixed_model_weight),
      TableRow("<2b", 
          "nbt<2",0,0,mixed_model_weight),
      TableRow("2b", 
          "nbt>=2&&nbm==2",0,0,mixed_model_weight),
      TableRow("3b", 
          "nbt>=2&&nbm>=3&&nbl==3",0,0,mixed_model_weight),
      TableRow("4b", 
          "nbt>=2&&nbm>=3&&nbl>=4",0,0,mixed_model_weight),
      TableRow("$N_\\mathrm{j}\\geq 4$", 
          "njet>=4",0,0,mixed_model_weight),
      TableRow("$N_\\mathrm{j}\\leq 5$", 
          "njet<=5",0,0,mixed_model_weight),
      TableRow("$p_\\mathrm{T}^\\mathrm{miss}>150$ GeV", 
          "met>=150",0,0,mixed_model_weight),
      TableRow("$N_\\mathrm{j}\\geq 4$ and $\\Delta m_\\mathrm{bb}<40$ GeV", 
          "njet>=4&&hig_cand_dm[0]<40",0,0,mixed_model_weight),
      TableRow("$N_\\mathrm{j}\\geq 4$ and $\\langle m_\\mathrm{bb}\\rangle<200$ GeV", 
          "njet>=4&&hig_cand_am[0]<200",0,0,mixed_model_weight),
      TableRow("$N_\\mathrm{j}\\geq 4$ and $\\Delta R_\\mathrm{max}<2.2$", 
          "njet>=4&&hig_cand_drmax[0]<2.2",0,0,mixed_model_weight),


      TableRow("Baseline except with 4b",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && nbt>=2 && nbm>=3 && nbl>=4 && njet >= 4 && njet <= 5 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except with 3b",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && nbt>=2 && nbm>=3 && nbl==3 && njet >= 4 && njet <= 5 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except with 2b",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && nbt>=2 && nbm==2 && njet >= 4 && njet <= 5 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except with <2b",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && nbt<2 && njet >= 4 && njet <= 5 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except 2b cut",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except $p_\\mathrm{T}^\\mathrm{miss}>150$ GeV cut",        
          Higfuncs::pass_filters && "nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except $N_\\mathrm{j}\\leq 5$ cut",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except $\\Delta R_\\mathrm{max}<2.2$ cut",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200",
          0,0,mixed_model_weight),
      TableRow("Baseline except $\\langle m_\\mathrm{bb}\\rangle<200$ GeV cut",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except $\\Delta m_\\mathrm{bb}<40$ GeV cut",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except track veto",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline except lepton veto",        
          Higfuncs::pass_filters && "met>150 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
      TableRow("Baseline",        
          Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0 && njet >= 4 && njet <= 5 && nbt >= 2 && met/mht<2 && met/met_calo<2 && !low_dphi_met && hig_cand_dm[0]<40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<2.2",
          0,0,mixed_model_weight),
    },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string).Tag("lep_veto_cutflow").Precision(1);
  }
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  std::vector<std::string> table_cut_strings = {
      "No cuts(lumi and cross section $\\PH\\to\\PQb\\PAQb$ only)", 
      "Filters, $p_\\text{T}^\\text{miss}>150 \\text{ GeV}$", 
      "$N_\\text{vl}=N_\\text{tk}=0$", 
      "$4 \\leq N_\\text{jet} \\leq 5$", 
      "$N_\\text{b}\\geq 2$", 
      "$p_\\text{T}^\\text{miss}$ quality and $\\Delta\\phi$ cuts",        
      "$\\Delta m_\\text{bb}<40 \\text{ GeV}$, $\\langle m_\\text{bb} \\rangle<200 \\text{ GeV}$",        
      "$\\Delta R_\\text{max}<2.2$",        
      "\\hline\n $100 < \\langle m\\rangle < 140 \\text{ GeV}$",        
      "$N_\\text{b}\\geq 3$",        
      "$N_\\text{b}=4$",        
      "$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
      "$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
      "$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$"};

  if(HigUtilities::is_in_string_options(options.string_options, "pie")) {
    TCanvas pie_canvas("pie_canvas","pie_canvas",700,1000);
    Hist1D* metnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("supppies__metnbcat").get());
    std::vector<TH1D*> metnb_proc_hist;
    Hist1D* drmaxnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("supppies__drmaxnbcat").get());
    std::vector<TH1D*> drmaxnb_proc_hist;
    for (auto proc : search_procs) {
      metnb_proc_hist.push_back(new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(proc.get())))->scaled_hist_));
      drmaxnb_proc_hist.push_back(new TH1D((static_cast<Hist1D::SingleHist1D*>(drmaxnb_hist1d->GetComponent(proc.get())))->scaled_hist_));
    }

    pie_canvas.cd();
    TPad cms_title_pad("cms_title_pad","cms_title_pad",0,0.9,0.571,1.0);
    cms_title_pad.Draw();
    cms_title_pad.cd();
    TLatex cms_label;
    cms_label.SetTextSize(0.3);
    cms_label.SetNDC(kTRUE);
    cms_label.SetTextAlign(22);
    cms_label.DrawLatex(0.5, 0.5,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}");

    pie_canvas.cd();
    TPad lowmet_2b_pad("lowmet_2b_pad","lowmet_2b_pad",0,0.7,0.286,0.9);
    lowmet_2b_pad.Draw();
    lowmet_2b_pad.cd();
    std::vector<double> lowmet_2b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowmet_2b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(1));
    }
    TPie lowmet_2b_pie("lowmet_2b_pie","2b, 150 < p_{T}^{miss} < 200 GeV",5,&lowmet_2b_proc_comp[0],&proc_colors[0]);
    lowmet_2b_pie.SetLabelFormat("%perc");
    lowmet_2b_pie.SetCircle(0.5, 0.48, 0.35);
    lowmet_2b_pie.Draw();

    pie_canvas.cd();
    TPad lowmet_3b_pad("lowmet_3b_pad","lowmet_3b_pad",0.286,0.7,0.571,0.9);
    lowmet_3b_pad.Draw();
    lowmet_3b_pad.cd();
    std::vector<double> lowmet_3b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowmet_3b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(3));
    }
    TPie lowmet_3b_pie("lowmet_3b_pie","3b, 150 < p_{T}^{miss} < 200 GeV",5,&lowmet_3b_proc_comp[0],&proc_colors[0]);
    lowmet_3b_pie.SetLabelFormat("%perc");
    lowmet_3b_pie.SetCircle(0.5, 0.48, 0.35);
    lowmet_3b_pie.Draw();

    pie_canvas.cd();
    TPad lowmet_4b_pad("lowmet_4b_pad","lowmet_4b_pad",0.571,0.7,0.857,0.9);
    lowmet_4b_pad.Draw();
    lowmet_4b_pad.cd();
    std::vector<double> lowmet_4b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowmet_4b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(5));
    }
    TPie lowmet_4b_pie("lowmet_4b_pie","4b, 150 < p_{T}^{miss} < 200 GeV",5,&lowmet_4b_proc_comp[0],&proc_colors[0]);
    lowmet_4b_pie.SetLabelFormat("%perc");
    lowmet_4b_pie.SetCircle(0.5, 0.48, 0.35);
    lowmet_4b_pie.Draw();

    pie_canvas.cd();
    TPad legend_pad("legend_pad","legend_pad",0.857,0.7,1.0,0.9);
    legend_pad.Draw();
    legend_pad.cd();
    TLegend leg(0.0, 0.0, 1.0, 1.0);
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      leg.AddEntry(metnb_proc_hist[proc_idx], proc_names[proc_idx].c_str(), "f");
    }
    leg.Draw();

    pie_canvas.cd();
    TPad highmet_2b_pad("highmet_2b_pad","highmet_2b_pad",0,0.5,0.286,0.7);
    highmet_2b_pad.Draw();
    highmet_2b_pad.cd();
    std::vector<double> highmet_2b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highmet_2b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(2));
    }
    TPie highmet_2b_pie("highmet_2b_pie","2b, p_{T}^{miss} > 200 GeV",5,&highmet_2b_proc_comp[0],&proc_colors[0]);
    highmet_2b_pie.SetLabelFormat("%perc");
    highmet_2b_pie.SetCircle(0.5, 0.48, 0.35);
    highmet_2b_pie.Draw();

    pie_canvas.cd();
    TPad highmet_3b_pad("highmet_3b_pad","highmet_3b_pad",0.286,0.5,0.571,0.7);
    highmet_3b_pad.Draw();
    highmet_3b_pad.cd();
    std::vector<double> highmet_3b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highmet_3b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(4));
    }
    TPie highmet_3b_pie("highmet_3b_pie","3b, p_{T}^{miss} > 200 GeV",5,&highmet_3b_proc_comp[0],&proc_colors[0]);
    highmet_3b_pie.SetLabelFormat("%perc");
    highmet_3b_pie.SetCircle(0.5, 0.48, 0.35);
    highmet_3b_pie.Draw();

    pie_canvas.cd();
    TPad highmet_4b_pad("highmet_4b_pad","highmet_4b_pad",0.571,0.5,0.857,0.7);
    highmet_4b_pad.Draw();
    highmet_4b_pad.cd();
    std::vector<double> highmet_4b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highmet_4b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(6));
    }
    TPie highmet_4b_pie("highmet_4b_pie","4b, p_{T}^{miss} > 200 GeV",5,&highmet_4b_proc_comp[0],&proc_colors[0]);
    highmet_4b_pie.SetLabelFormat("%perc");
    highmet_4b_pie.SetCircle(0.5, 0.48, 0.35);
    highmet_4b_pie.Draw();

    pie_canvas.cd();
    TPad highdrmax_2b_pad("highdrmax_2b_pad","highdrmax_2b_pad",0,0.2,0.286,0.4);
    highdrmax_2b_pad.Draw();
    highdrmax_2b_pad.cd();
    std::vector<double> highdrmax_2b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highdrmax_2b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(1));
    }
    TPie highdrmax_2b_pie("highdrmax_2b_pie","2b, 1.1 < #DeltaR_{max} < 2.2",5,&highdrmax_2b_proc_comp[0],&proc_colors[0]);
    highdrmax_2b_pie.SetLabelFormat("%perc");
    highdrmax_2b_pie.SetCircle(0.5, 0.48, 0.35);
    highdrmax_2b_pie.Draw();

    pie_canvas.cd();
    TPad highdrmax_3b_pad("highdrmax_3b_pad","highdrmax_3b_pad",0.286,0.2,0.571,0.4);
    highdrmax_3b_pad.Draw();
    highdrmax_3b_pad.cd();
    std::vector<double> highdrmax_3b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highdrmax_3b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(3));
    }
    TPie highdrmax_3b_pie("highdrmax_3b_pie","3b, 1.1 < #DeltaR_{max} < 2.2",5,&highdrmax_3b_proc_comp[0],&proc_colors[0]);
    highdrmax_3b_pie.SetLabelFormat("%perc");
    highdrmax_3b_pie.SetCircle(0.5, 0.48, 0.35);
    highdrmax_3b_pie.Draw();

    pie_canvas.cd();
    TPad highdrmax_4b_pad("highdrmax_4b_pad","highdrmax_4b_pad",0.571,0.2,0.857,0.4);
    highdrmax_4b_pad.Draw();
    highdrmax_4b_pad.cd();
    std::vector<double> highdrmax_4b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      highdrmax_4b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(5));
    }
    TPie highdrmax_4b_pie("highdrmax_4b_pie","4b, 1.1 < #DeltaR_{max} < 2.2",5,&highdrmax_4b_proc_comp[0],&proc_colors[0]);
    highdrmax_4b_pie.SetLabelFormat("%perc");
    highdrmax_4b_pie.SetCircle(0.5, 0.48, 0.35);
    highdrmax_4b_pie.Draw();

    pie_canvas.cd();
    TPad lowdrmax_2b_pad("lowdrmax_2b_pad","lowdrmax_2b_pad",0,0.0,0.286,0.2);
    lowdrmax_2b_pad.Draw();
    lowdrmax_2b_pad.cd();
    std::vector<double> lowdrmax_2b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowdrmax_2b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(2));
    }
    TPie lowdrmax_2b_pie("lowdrmax_2b_pie","2b, #DeltaR_{max} < 1.1",5,&lowdrmax_2b_proc_comp[0],&proc_colors[0]);
    lowdrmax_2b_pie.SetLabelFormat("%perc");
    lowdrmax_2b_pie.SetCircle(0.5, 0.48, 0.35);
    lowdrmax_2b_pie.Draw();

    pie_canvas.cd();
    TPad lowdrmax_3b_pad("lowdrmax_3b_pad","lowdrmax_3b_pad",0.286,0.0,0.571,0.2);
    lowdrmax_3b_pad.Draw();
    lowdrmax_3b_pad.cd();
    std::vector<double> lowdrmax_3b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowdrmax_3b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(4));
    }
    TPie lowdrmax_3b_pie("lowdrmax_3b_pie","3b, #DeltaR_{max} < 1.1",5,&lowdrmax_3b_proc_comp[0],&proc_colors[0]);
    lowdrmax_3b_pie.SetLabelFormat("%perc");
    lowdrmax_3b_pie.SetCircle(0.5, 0.48, 0.35);
    lowdrmax_3b_pie.Draw();

    pie_canvas.cd();
    TPad lowdrmax_4b_pad("lowdrmax_4b_pad","lowdrmax_4b_pad",0.571,0.0,0.857,0.2);
    lowdrmax_4b_pad.Draw();
    lowdrmax_4b_pad.cd();
    std::vector<double> lowdrmax_4b_proc_comp;
    for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
      lowdrmax_4b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(6));
    }
    TPie lowdrmax_4b_pie("lowdrmax_4b_pie","4b, #DeltaR_{max} < 1.1",5,&lowdrmax_4b_proc_comp[0],&proc_colors[0]);
    lowdrmax_4b_pie.SetLabelFormat("%perc");
    lowdrmax_4b_pie.SetCircle(0.5, 0.48, 0.35);
    lowdrmax_4b_pie.Draw();

    pie_canvas.Print("plots/supp_res_pies.pdf");
    std::cout << "Opening plots/supp_res_pies.pdf\n";
  }

  if(HigUtilities::is_in_string_options(options.string_options, "cutflow")) {
    Table* cutflow_table = static_cast<Table*>(pm.GetFigure("cutflow").get());
    TH2D cutflow_hist("signal_cutflow","",5,0.5,5.5,14,0.5,14.5);
    //float lumi_value = HigUtilities::getLuminosity(options.year_string);
    int proc_idx = 0;
    for (std::shared_ptr<Process> proc : cutflow_procs) {
      cutflow_hist.GetXaxis()->SetBinLabel(proc_idx+1, cutflow_proc_names[proc_idx].c_str());
      std::vector<GammaParams> proc_yield = cutflow_table->Yield(proc.get(),1.0);
      for (unsigned table_row = 0; table_row < proc_yield.size(); table_row++) {
        cutflow_hist.SetBinContent(proc_idx+1,table_row+1,proc_yield[table_row].Yield());
        cutflow_hist.SetBinError(proc_idx+1,table_row+1,proc_yield[table_row].Uncertainty());
        //std::cout << cutflow_proc_names[proc_idx] << ", " << table_cut_strings[table_row] << ": " << proc_yield[table_row].Yield() << std::endl;
      }
      proc_idx++;
    }
    for (unsigned table_row = 0; table_row < table_cut_strings.size(); table_row++) {
      cutflow_hist.GetYaxis()->SetBinLabel(table_row+1, table_cut_strings[table_row].c_str());
    }
    std::cout << "open tables/CMS-SUS-20-004_aux_Table_005.root" << std::endl;
    TFile * cutflow_file = TFile::Open("tables/CMS-SUS-20-004_aux_Table_005.root","RECREATE");
    cutflow_hist.Write();
    cutflow_file->Close();
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

