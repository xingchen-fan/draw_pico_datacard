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
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "2016";
  bool unblind = false;
  //possible string_options: paper_style
  string string_options = "";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

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

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_data = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};
  if (HigUtilities::is_in_string_options(string_options, "paper_style")) {
     plt_shapes = {lin_shapes_paper};
     plt_shapes_data = {lin_shapes_paper_data};
     if (!unblind) {
       plt_lin = {lin_norm_paper};
       plt_log = {log_norm_paper};
     }
     else {
       plt_lin = {lin_norm_paper_data};
       plt_log = {log_norm_paper_data};
     }
  }

  // Set options
  //string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  //string mc_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado";
  string mc_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";

  //string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  string data_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  string zll_data_skim_folder = "data/merged_higdata_higlep2T/";

  //string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  //string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string search_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
  string zll_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
  string qcd_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";

  //years = {2016, 2017, 2018};
  //years = {2016};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  int lumi_precision = 1;
  if (HigUtilities::is_in_string_options(string_options, "paper_style"))
    lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string, lumi_precision);

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years*Functions::w_pileup;
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc triggers_data = lepton_triggers;
  //NamedFunc triggers_data = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50"
  //                    "||HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";

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

  vector<shared_ptr<Process> > zll_procs;
  if (HigUtilities::is_in_string_options(string_options, "paper_style")) {
    // Set mc processes
    zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch"));
    //zll_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
    //                attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["vjets"]),"stitch"));
    zll_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["zjets"]),"stitch"));
    zll_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["wjets"]),"stitch"));
    zll_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["qcd"]),"stitch")); 
    zll_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["other_and_single_t"]),"stitch"));
  }
  else {
    // Set mc processes
    //zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
    //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch"));
    zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
    zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
    //zll_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
    //                attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["vjets"]),"stitch"));
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
  }

  //vector<string> sigm = {"175", "500", "950"};
  ////vector<string> sigm = {"400"};
  ////vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  //vector<int> sig_colors = {kSpring, kPink, kAzure, kOrange}; // need sigm.size() >= sig_colors.size()
  //for (unsigned isig(0); isig<sigm.size(); isig++){
  //  zll_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
  //    sig_colors[isig], attach_folder(sig_base_folder, years, zll_sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), 
  //    "stitch"));
  //}

  if (unblind) {
    zll_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    triggers_data
                    ));
  }


  // Set processes according to btag
  // Getting colors
  // TColor * color; float red, green, blue;
  // color = gROOT->GetColor(kAzure+1); color->GetRGB(red,green,blue); cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl;
  vector<shared_ptr<Process> > zll_procs_btag;
  zll_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  zll_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  //zll_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==2)"));

  // Set processes according to true number of b
  vector<shared_ptr<Process> > zll_procs_trueB;
  zll_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  zll_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  zll_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));
  zll_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, colors("true_3b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),  "stitch"&& Functions::ntrub==3));
  zll_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("4 B-hadrons", Process::Type::background, colors("true_4b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==4));

  // Set processes according to true number of b simple version
  vector<shared_ptr<Process> > zll_procs_trueB012;
  zll_procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  zll_procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  zll_procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));

  vector<shared_ptr<Process> > zll_data_procs_btag;
  if (unblind) {
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("0b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    "(nbm==0)" && triggers_data
                    ));
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("1b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    "(nbm==1)" && triggers_data
                    ));
    //zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("2b Data", Process::Type::data, kBlack,
    //                attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
    //                "stitch&&(nbm==2)" && triggers_data
    //                ));
  } else {
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
    //zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
    //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==2)"));
  }


  vector<shared_ptr<Process> > zll_procs_nisr;
  zll_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 0", Process::Type::background,colors("nisr_0"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&nisr==0"));
  zll_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 1", Process::Type::background,colors("nisr_1"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&nisr==1"));
  zll_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 2", Process::Type::background,colors("nisr_2"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&nisr==2"));

  NamedFunc base_filters = Higfuncs::final_zll_pass_filters;
  //NamedFunc base_filters = Functions::hem_veto && "pass && met/met_calo<5";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  // Could add additional generation cuts

  // resolved cuts
  //NamedFunc base_resolved = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         &&Higfuncs::lead_signal_lepton_pt>30;
  NamedFunc zll_resolved =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  PlotMaker pm;

  // According to true b
  // (0b, 1b, 2b tag)  low drmax 
  pm.Push<Table>("FixName:syst__zll_pies_trueb__0b_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies_trueb__1b_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies_trueb__2b_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);
  // (2b, 3b, 4b tag)  high drmax 
  pm.Push<Table>("FixName:syst__zll_pies_trueb__0b_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies_trueb__1b_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies_trueb__2b_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs_trueB, true, true, true);

  // 0b (pll[0]: 75, 150, 200, 300) low drmax
  pm.Push<Table>("FixName:syst__zll_pies__0b_met0_lowdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0               &&ll_pt[0]<=75  &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met75_lowdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met150_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met200_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met300_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>300                 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  // 1b (pll[0]: 75, 150, 200, 300) low drmax
  pm.Push<Table>("FixName:syst__zll_pies__1b_met0_lowdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1               &&ll_pt[0]<=75  &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met75_lowdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met150_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met200_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met300_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>300                 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  // 2b (pll[0]: 75, 150, 200, 300) low drmax
  pm.Push<Table>("FixName:syst__zll_pies__2b_met0_lowdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2               &&ll_pt[0]<=75  &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met75_lowdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met150_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met200_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met300_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>300                 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), zll_procs, true, true, true);
  // 0b (pll[0]: 75, 150, 200, 300) high drmax
  pm.Push<Table>("FixName:syst__zll_pies__0b_met0_highdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0               &&ll_pt[0]<=75  &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met75_highdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met150_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met200_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__0b_met300_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==0&&ll_pt[0]>300                 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  // 1b (pll[0]: 75, 150, 200, 300) high drmax
  pm.Push<Table>("FixName:syst__zll_pies__1b_met0_highdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1               &&ll_pt[0]<=75  &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met75_highdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met150_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met200_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__1b_met300_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==1&&ll_pt[0]>300                 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  // 2b (pll[0]: 75, 150, 200, 300) high drmax
  pm.Push<Table>("FixName:syst__zll_pies__2b_met0_highdrmax"    , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2               &&ll_pt[0]<=75  &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met75_highdrmax"   , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>75  &&ll_pt[0]<=150 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met150_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>150 &&ll_pt[0]<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met200_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>200 &&ll_pt[0]<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);
  pm.Push<Table>("FixName:syst__zll_pies__2b_met300_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm==2&&ll_pt[0]>300                 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), zll_procs, true, true, true);

  // 0b1b
  pm.Push<Table>("FixName:syst__zll_pies__0b1b"  , vector<TableRow> ({TableRow("", base_filters&&zll_resolved&&"nbm<=1", 0, 0, weight)}), zll_procs, true, true, true);

  // comparison of data with mc
  // [drmax;no drmax] (0b, 1b, 0b+1b)
  pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
    "(nbm==0)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_drmax_0b").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
    "(nbm==1)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_drmax_1b").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
    "(nbm<=1)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_drmax_0b1b").LuminosityTag(total_luminosity_string);

  // [<m>;no <m>, low_drmax] (0b, 1b)
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "hig_cand_drmax[0]<1.1&&"
    "(nbm==0)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_0b_lowdrmax").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "hig_cand_drmax[0]<1.1&&"
    "(nbm==1)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_1b_lowdrmax").LuminosityTag(total_luminosity_string);
  // [<m>;no <m>, high_drmax] (0b, 1b)
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "hig_cand_drmax[0]>=1.1&&"
    "(nbm==0)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_0b_highdrmax").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "hig_cand_drmax[0]>=1.1&&"
    "(nbm==1)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_1b_highdrmax").LuminosityTag(total_luminosity_string);


  // [<m>;no <m>] (0b): Extra plot to compare with AN-2016
  pm.Push<Hist1D>(Axis(25, 0, 250, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "(nbm==0)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_0b").LuminosityTag(total_luminosity_string);
  // [<m>;no <m>] (1b): Extra plot to compare with AN-2016
  pm.Push<Hist1D>(Axis(25, 0, 250, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==2&&njet>=4&&njet<=5&&met<50&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
    "(nbm==1)", 
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_amjj_1b").LuminosityTag(total_luminosity_string);

  // [p_ll] (0b+1b)
  pm.Push<Hist1D>(Axis(16, 0, 400, "ll_pt[0]", "p_{T}^{ll} [GeV]", {}),
    base_filters&&zll_resolved&&"nbm<=1",
    zll_procs, plt_lin).Weight(weight).Tag("FixName:syst__zll_pt_ll").LuminosityTag(total_luminosity_string);

  // [<m> shapes (true 0b, 1b, 2b);(low, high) drmax]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]<1.1",
    zll_procs_trueB012, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_trueb__lowdrmax").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]>=1.1",
    zll_procs_trueB012, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_trueb__highdrmax").LuminosityTag(total_luminosity_string);
  // [<m> shapes (0b, 1b, 2b);(low, high) drmax]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]<1.1",
    zll_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_btag__lowdrmax").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]>=1.1",
    zll_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_btag__highdrmax").LuminosityTag(total_luminosity_string);

  // [<m> shapes (0b data, 1b data); (low, high) drmax]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]<1.1",
    zll_data_procs_btag, plt_shapes_data).Weight(weight).Tag("FixName:syst__zll_amjj_btag__lowdrmax__data").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]>=1.1",
    zll_data_procs_btag, plt_shapes_data).Weight(weight).Tag("FixName:syst__zll_amjj_btag__highdrmax__data").LuminosityTag(total_luminosity_string);

  // [<m> shapes (0b mc, 1b mc); (low, high) drmax]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]<1.1",
    zll_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_btag__lowdrmax__mc").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&zll_resolved&&"hig_cand_drmax[0]>=1.1",
    zll_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__zll_amjj_btag__highdrmax__mc").LuminosityTag(total_luminosity_string);

  // kappa plot
  if (unblind) {
    (void)!system(("./run/higgsino/plot_kappas.exe --sample zll --scen data --unblind -o paper_style --year "+year_string).c_str());
  } else {
    (void)!system(("./run/higgsino/plot_kappas.exe --sample zll --scen mc_as_data --year "+year_string).c_str());
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
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "so:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "year"){
        year_string = optarg;
      } else if(optname == "unblind"){
        unblind = true;
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
