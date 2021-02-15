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

namespace{
  bool single_thread = false;
  string year_string = "2018";
  bool unblind = true;
  const std::vector<double> PT_BINS = {25,50,75,100,125,150,175,200,250,300,400,500,600};
}

std::vector<double> get_pt_sfs(PlotMaker &pm, unsigned int pm_index);

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
  PlotOpt lin_norm_nooverflow = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Overflow(OverflowType::none);
  PlotOpt lin_norm_nonorm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Stack(StackType::lumi_shapes);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_norm_data_nooverflow = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).Overflow(OverflowType::none);
  PlotOpt lin_norm_data_nonorm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).Stack(StackType::lumi_shapes);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_nooverflow = {lin_norm_nooverflow};
  vector<PlotOpt> plt_lin_nonorm = {lin_norm_nonorm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) {
    plt_lin = {lin_norm_data};
    plt_log = {log_norm_data};
    plt_lin_nooverflow = {lin_norm_data_nooverflow};
    plt_lin_nonorm = {lin_norm_data_nonorm};
  }

  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).YAxis(YAxisType::linear).Overflow(OverflowType::overflow)};
  vector<PlotOpt> twodim_plotopts_log = {style2D().Title(TitleType::info).YAxis(YAxisType::log).Overflow(OverflowType::overflow)};

  string higgsino_version = "";

  // Set options
  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt"+higgsino_version+"/";
  //string mc_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";

  string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt"+higgsino_version+"/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  //string ttbar_data_skim_folder = "data/skim_higlep1T/";

  string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt"+higgsino_version+"/";
  //string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";

  //years = {2016, 2017, 2018};
  //years = {2016};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  NamedFunc weight = Higfuncs::final_weight;
  //if (years.size()==1 && *years.begin()==2016) weight *= "137.";
  //else weight *= w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;

  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  //NamedFunc lepton_triggers = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50 || HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";
  //NamedFunc met_triggers = "(HLT_PFMET110_PFMHT110_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)";
  NamedFunc triggers_data = lepton_triggers || met_triggers;

  NamedFunc base_filters = Higfuncs::final_ttbar_pass_filters;
  //NamedFunc base_filters = Functions::hem_veto && "pass && met/met_calo<5 && weight<1.5";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/met_calo<5 && weight<10"; //since pass_fsjets is not quite usable...

  // resolved cuts
  //NamedFunc search_resolved = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         &&Higfuncs::lead_signal_lepton_pt>30;


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

  vector<shared_ptr<Process> > ttbar_procs;
  // Set mc processes
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["vjets"]),"stitch"));
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
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data));
  }

  vector<shared_ptr<Process> > nonorm_procs;
  nonorm_procs.push_back(Process::MakeShared<Baby_pico>("MC", Process::Type::background, colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch"));
  nonorm_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                  attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data));

  vector<shared_ptr<Process> > ttbar_procs_nodata;
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
  ttbar_procs_nodata.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));

  vector<shared_ptr<Process> > ttbar_procs_data;
  ttbar_procs_data.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                  attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}), triggers_data));

  //// Set signal mc
  //vector<string> signal_mass = {"200", "450", "700", "950"}; 
  //vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need signal_mass.size() >= sig_colors.size()
  //for (unsigned isig(0); isig<signal_mass.size(); isig++){
  //  procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+signal_mass[isig]+",1)", Process::Type::signal, 
  //    sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+signal_mass[isig]+"_mLSP-0*.root"}), "1"));
  //}

  // Set processes according to btag
  // Getting colors
  // TColor * color;; float red, green, blue;
  // color = gROOT->GetColor(kAzure+1); color->GetRGB(red,green,blue); cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl;
  vector<shared_ptr<Process> > ttbar_procs_btag;
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));

  // Set processes according to true number of b
  vector<shared_ptr<Process> > ttbar_procs_trueB;
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, colors("true_3b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),  "stitch"&& Functions::ntrub==3));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("4 B-hadrons", Process::Type::background, colors("true_4b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==4));

  vector<shared_ptr<Process> > ttbar_procs_nisr;
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 0", Process::Type::background,colors("nisr_0"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==0"));
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 1", Process::Type::background,colors("nisr_1"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==1"));
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 2", Process::Type::background,colors("nisr_2"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==2"));

  vector<shared_ptr<Process> > ttbar_data_procs_btag;
  if (unblind) {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("2b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt==2&&nbm==2)" && triggers_data
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("3b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm==3&&nbl==3)" && triggers_data
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("4b Data", Process::Type::data, kGreen,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm>=3&&nbl>=4)" && triggers_data
                    ));
  } else {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  }

  vector<shared_ptr<Process> > ttbar_mc_procs_btag;
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  
  NamedFunc onelep_baseline = base_filters && (Higfuncs::lead_signal_lepton_pt>30) && "nlep==1&&mt<=100&&njet>=4&&njet<=5&&hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&nbt>=2";

  const NamedFunc dphi_min("dphi_min",[](const Baby &b) -> NamedFunc::ScalarType {
    double min_dphi = 3.1415;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_isgood()->at(jet_idx)) {
        if (fabs(b.jet_met_dphi()->at(jet_idx)) < min_dphi) 
          min_dphi = fabs(b.jet_met_dphi()->at(jet_idx));
      }
    }
    return min_dphi;
  });

  PlotMaker pm;
  unsigned int pm_idx = 0;
  unsigned int pt_scale_hist_idx = 0;

  //checking defecit/excess in drmax-met space
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_lowdrmax_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1&&met<125",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_lowdrmaxlowmet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1&&met>125&&met<175",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_lowdrmaxmidmet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1&&met>175",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_lowdrmaxhimet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]>1.1&&met<125",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_hidrmaxlowmet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]>1.1&&met>125&&met<175",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_hidrmaxmidmet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]>1.1&&met>175",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__am_hidrmaxhimet_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  //conclusion: no met dependence

  //omit low amm region
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1",
    nonorm_procs, plt_lin_nonorm).Weight(weight).Tag("FixName:defecit__am_nonorm_lowdrmax_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline,
    nonorm_procs, plt_lin_nonorm).Weight(weight).Tag("FixName:defecit__am_nonorm_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  //conclusion: excess in low <mbb> not defecit in high
  
  //low level plots
  pm.Push<Hist1D>(Axis(25, 0, 600, "jet_pt", "Jet p_{T} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__jetpt_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 600, "jet_pt", "Jet p_{T} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__jetpt_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline,
    ttbar_procs_nodata, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_mc_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline,
    ttbar_procs_data, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_data_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs_nodata, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_mc_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs_data, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_data_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "jet_pt[0]>200",
    ttbar_procs_nodata, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_mc_highjetpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "jet_pt[0]>200",
    ttbar_procs_data, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_data_highjetpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "hig_cand_am[0]>110&&jet_pt[0]>200",
    ttbar_procs_nodata, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_mc_higamhighjetpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist2D>(Axis(25, -2.5, 2.5, "jet_eta", "Jet #eta", {}),
    Axis(25, -3.14, 3.14, "jet_phi", "Jet #phi", {}),
    onelep_baseline && "hig_cand_am[0]>110&&jet_pt[0]>200",
    ttbar_procs_data, twodim_plotopts).Weight(weight).Tag("FixName:defecit__jetetaphi_data_highamhighjetpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  //no obvious detector failures
  
  //plots split by mbb
  pm.Push<Hist1D>(Axis(25, 0, 3.0, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1,2.2}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__drmax_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 3.0, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1,2.2}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__drmax_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 40.0, "hig_cand_dm[0]", "#Delta m", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__dm_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 40.0, "hig_cand_dm[0]", "#Delta m", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__dm_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 600.0, "ht", "H_{T} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__ht_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 600.0, "ht", "H_{T} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__ht_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 400, "met", "p_{T}^{miss} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__met_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 400, "met", "p_{T}^{miss} [GeV]", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__met_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::jetid_nb, "N_{b}", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nb_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::jetid_nb, "N_{b}", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nb_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "nmu", "N_{#mu}", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nmu_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(2, -0.5, 1.5, "nmu", "N_{#mu}", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nmu_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(2, 3.5, 5.5, "njet", "N_{j}", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nj_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(2, 3.5, 5.5, "njet", "N_{j}", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__nj_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 3.14, dphi_min, "#Delta #phi_{min}", {}),
    onelep_baseline && "hig_cand_am[0]>110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__dphimet_higham_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(25, 0, 3.14, dphi_min, "#Delta #phi_{min}", {}),
    onelep_baseline && "hig_cand_am[0]<110",
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__dphimet_lowam_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm_idx++;
  
  pm.Push<Hist1D>(Axis(PT_BINS, "jet_pt", "Jet p_{T} [GeV]", {}),
    onelep_baseline,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:defecit__scalejetpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pt_scale_hist_idx = pm_idx;
  pm_idx++;

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  std::vector<double> pt_sfs = get_pt_sfs(pm, pt_scale_hist_idx);

  const NamedFunc w_pt("w_pt",[pt_sfs](const Baby &b) -> NamedFunc::ScalarType {
    double w_pt_ = 1.0;
    if (b.SampleType()<0) return 1.0; //data
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_isgood()->at(jet_idx)) {
        for (unsigned int pt_bin = 0; pt_bin < PT_BINS.size()-1; pt_bin++) {
          double max_value = PT_BINS[pt_bin+1];
          if (pt_bin == PT_BINS.size()-2) max_value = 9999.0;
          if (b.jet_pt()->at(jet_idx)>=PT_BINS[pt_bin] && b.jet_pt()->at(jet_idx) < max_value) {
            w_pt_ = w_pt_*pt_sfs[pt_bin];
            break;
          }
        }
      }
    }
    return w_pt_;
  });

  PlotMaker pm2;
  pm2.Push<Hist1D>(Axis(PT_BINS, "jet_pt", "Jet p_{T} [GeV]", {}),
    onelep_baseline,
    ttbar_procs, plt_lin).Weight(weight*w_pt).Tag("FixName:defecit__jetptscale_wpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm2.Push<Hist1D>(Axis(25, 0, 600, "jet_pt", "Jet p_{T} [GeV]", {}),
    onelep_baseline,
    ttbar_procs, plt_lin).Weight(weight*w_pt).Tag("FixName:defecit__jetpt_wpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm2.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline,
    ttbar_procs, plt_lin).Weight(weight*w_pt).Tag("FixName:defecit__am_wpt_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm2.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    onelep_baseline && "hig_cand_drmax[0]<1.1",
    ttbar_procs, plt_lin).Weight(weight*w_pt).Tag("FixName:defecit__am_wpt_lowdrmax_1lcr_"+year_string).LuminosityTag(total_luminosity_string);
  pm2.multithreaded_ = !single_thread;
  pm2.min_print_ = true;
  pm2.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

std::vector<double> get_pt_sfs(PlotMaker &pm, unsigned int pm_index) {
  double bot_min = 0, bot_max = 0;
  Hist1D *variations_hist = static_cast<Hist1D*>(pm.Figures()[pm_index].get());
  TH1D variations_bottom_hist = variations_hist->GetBottomPlots(bot_min, bot_max).at(0);
  std::vector<double> pt_sfs;
  for (unsigned int pt_bin = 0; pt_bin < PT_BINS.size()-1; pt_bin++) {
    pt_sfs.push_back(variations_bottom_hist.GetBinContent(static_cast<int>(pt_bin)+1));
    std::cout << "This bin content: " << variations_bottom_hist.GetBinContent(static_cast<int>(pt_bin)) << std::endl;
  }
  return pt_sfs;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"unblind", no_argument, 0, 0},
      {"year", required_argument, 0, 0},
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
