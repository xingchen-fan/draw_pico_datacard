#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "higgsino/ordered_dict.hpp"

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

namespace{
  bool single_thread = false;
  bool unblind = false;
  bool unblind_signalregion = false;
  // sample can be search, ttbar, zll, qcd
  string sample_name = "search";
  // year can be 2016, 2017, 2018 or 2016,2017,2018 or run2
  string year_string = "run2";
  // lepton type can be electron or muon. Used for ttbar and zll sample.
  string lepton_type = "";
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  // Options: plot_additional_variables,plot_ht_correlation,plot_eta_vs_phi,plot_in_bins, cut_3b4b
  string string_options = "";
}

const NamedFunc goodJet_phi("goodJet_phi", [](const Baby &b) -> NamedFunc::VectorType{
  vector<double> jet_phi;
  jet_phi.reserve((*b.jet_isgood()).size());
  for (unsigned iJet = 0; iJet < (*b.jet_isgood()).size(); ++iJet) {
    bool jet_isgood = (*b.jet_isgood())[iJet];
    if (jet_isgood) {
      jet_phi.push_back((*b.jet_phi())[iJet]);
    }
  }
  return jet_phi;
});
const NamedFunc goodJet_eta("goodJet_eta", [](const Baby &b) -> NamedFunc::VectorType{
  vector<double> jet_eta;
  jet_eta.reserve((*b.jet_isgood()).size());
  for (unsigned iJet = 0; iJet < (*b.jet_isgood()).size(); ++iJet) {
    bool jet_isgood = (*b.jet_isgood())[iJet];
    if (jet_isgood) {
      jet_eta.push_back((*b.jet_eta())[iJet]);
    }
  }
  return jet_eta;
});
const NamedFunc goodJet0_eta("goodJet0_eta", [](const Baby &b) -> NamedFunc::VectorType{
  vector<double> jet_eta;
  jet_eta.reserve((*b.jet_isgood()).size());
  for (unsigned iJet = 0; iJet < (*b.jet_isgood()).size(); ++iJet) {
    bool jet_isgood = (*b.jet_isgood())[iJet];
    if (jet_isgood) {
      jet_eta.push_back((*b.jet_eta())[iJet]);
      break;
    }
  }
  return jet_eta;
});

const NamedFunc weight_ht("weight_ht", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>100 && b.ht()<=150) weight = 1.424;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.115;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.030;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.016;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.994;
  else if (b.ht()>350 && b.ht()<=400) weight = 1.037;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.901;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.919;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.914;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.874;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.812;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.954;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.825;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.911;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.700;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.860;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.705;
  else if (b.ht()>950) weight = 0.874;
  return weight;
});

void combine_bins(vector<pair<string, NamedFunc> > & combined_bins, vector<pair<string, NamedFunc> > const & bins_a,  vector<pair<string, NamedFunc> > const & bins_b) {
  for (auto const & bin_a : bins_a) {
    for (auto const & bin_b : bins_b) {
      combined_bins.push_back({bin_a.first+"_"+bin_b.first,bin_a.second && bin_b.second});
      //cout<<bin_a.first+"_"+bin_b.first<<" "<<(bin_a.second && bin_b.second).Name()<<endl;
    }
  }
}
void combine_bins(vector<pair<string, NamedFunc> > & combined_bins, vector<pair<string, NamedFunc> > const & bins_a,  vector<pair<string, NamedFunc> > const & bins_b, vector<pair<string, NamedFunc> > const & bins_c) {
  for (auto const & bin_a : bins_a) {
    for (auto const & bin_b : bins_b) {
      for (auto const & bin_c : bins_c) {
        combined_bins.push_back({bin_a.first+"_"+bin_b.first+"_"+bin_c.first, bin_a.second && bin_b.second && bin_c.second});
        //cout<<bin_a.first+"_"+bin_b.first+"_"+bin_c.first<<" "<<(bin_a.second && bin_b.second && bin_c.second).Name()<<endl;
      }
    }
  }
}

int main(int argc, char *argv[]){
  // ./run/higgsino/plot_n_minus_1.exe --sample search --year 2016
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    //.LegendColumns(3);
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_print = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).PrintVals(true);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_norm_data_print = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_print = {lin_norm_print};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_lin_print = {lin_norm_data_print};
  if (unblind) plt_log = {log_norm_data};

  string higgsino_version = "";
  bool plot_additional_variables = HigUtilities::is_in_string_options(string_options, "plot_additional_variables");

  // Set options
  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath";
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  //string search_mc_skim_folder = "mc/merged_higmc_higloose/";
  string search_mc_skim_folder = "mc/merged_higmc_preselect/";
  // higlep1T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
  //   nlep==1 && 
  //   (Max$(el_pt*el_sig)>40 || Max$(mu_pt*mu_sig)>40) // pass_1l_trig40 (sig is signal lepton)
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  // higlep2T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nlep==2 && 
  //   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  // higqcd with met150:
  //   met>150  // Since applied to met150 skim
  string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string search_data_skim_folder = "data/merged_higdata_higloose/";
  string search_data_skim_folder = "data/merged_higdata_preselect/";
  string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  //string ttbar_data_skim_folder = "data/skim_higlep1T/";
  string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  //string zll_data_skim_folder = "data/skim_higlep2T/";
  string qcd_data_skim_folder = "data/merged_higdata_higqcd/";

  string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string search_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
  string search_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
  string zll_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
  string qcd_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";

  //years = {2016, 2017, 2018};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  NamedFunc weight = Higfuncs::final_weight;

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

  string mc_skim_folder;
  if (sample_name == "ttbar") mc_skim_folder = ttbar_mc_skim_folder;
  else if (sample_name == "zll") mc_skim_folder = zll_mc_skim_folder;
  else if (sample_name == "qcd") mc_skim_folder = qcd_mc_skim_folder;
  else mc_skim_folder = search_mc_skim_folder;
  string data_skim_folder;
  if (sample_name == "ttbar") data_skim_folder = ttbar_data_skim_folder;
  else if (sample_name == "zll") data_skim_folder = zll_data_skim_folder;
  else if (sample_name == "qcd") data_skim_folder = qcd_data_skim_folder;
  else data_skim_folder = search_data_skim_folder;
  string sig_skim_folder;
  if (sample_name == "ttbar") sig_skim_folder = ttbar_sig_skim_folder;
  else if (sample_name == "zll") sig_skim_folder = zll_sig_skim_folder;
  else if (sample_name == "qcd") sig_skim_folder = qcd_sig_skim_folder;
  else sig_skim_folder = search_sig_skim_folder;

  NamedFunc triggers_data = "1";
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;
  if (sample_name == "zll") triggers_data = lepton_triggers;
  else if (sample_name == "ttbar") triggers_data = lepton_triggers || met_triggers;
  else if (sample_name == "qcd") triggers_data = met_triggers;
  else  triggers_data = met_triggers;

  NamedFunc base_resolved = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  //NamedFunc ttbar_resolved_cuts = 
  //                       "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  //NamedFunc zll_resolved_cuts =
  //                       "nlep==2&&njet>=4&&njet<=5&&met<50&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  //NamedFunc qcd_resolved_cuts =
  //                       "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  torch::OrderedDict<string, NamedFunc> search_resolved_cuts;
  search_resolved_cuts.insert("ntk", "ntk==0");
  search_resolved_cuts.insert("low_dphi_met", "!low_dphi_met");
  search_resolved_cuts.insert("nvlep", "nvlep==0");
  search_resolved_cuts.insert("met", "met>150");
  search_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  search_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  search_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  search_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  search_resolved_cuts.insert("btags", "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))");
  search_resolved_cuts.insert("mht_filter", "met/mht<2");
  search_resolved_cuts.insert("met_calo_filter", "met/met_calo<2");
  search_resolved_cuts.insert("weight", "weight<1.5");
  torch::OrderedDict<string, NamedFunc> ttbar_resolved_cuts;
  if (lepton_type == "electron") ttbar_resolved_cuts.insert("nel", "nel==1");
  else if (lepton_type == "muon") ttbar_resolved_cuts.insert("nmu", "nmu==1");
  else ttbar_resolved_cuts.insert("nlep", "nlep==1");
  //ttbar_resolved_cuts.insert("low_dphi_met", "!low_dphi_met");
  ttbar_resolved_cuts.insert("lep_pt", Higfuncs::lead_signal_lepton_pt>30);
  ttbar_resolved_cuts.insert("mt", "mt<=100");
  ttbar_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  ttbar_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  ttbar_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  ttbar_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  ttbar_resolved_cuts.insert("btags", "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))");
  ttbar_resolved_cuts.insert("met_calo_filter", "met/met_calo<5");
  ttbar_resolved_cuts.insert("weight", "weight<1.5");
  torch::OrderedDict<string, NamedFunc> zll_resolved_cuts;
  if (lepton_type == "electron") zll_resolved_cuts.insert("nel", "nel==2");
  else if (lepton_type == "muon") zll_resolved_cuts.insert("nmu", "nmu==2");
  else zll_resolved_cuts.insert("nlep", "nlep==2");
  zll_resolved_cuts.insert("lep_pt", Higfuncs::lead_signal_lepton_pt>30);
  zll_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  zll_resolved_cuts.insert("met", "met<50");
  zll_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  zll_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  zll_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  zll_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");
  zll_resolved_cuts.insert("met_calo_filter", "met/met_calo<5");
  zll_resolved_cuts.insert("weight", "weight<1.5");
  torch::OrderedDict<string, NamedFunc> qcd_resolved_cuts;
  qcd_resolved_cuts.insert("low_dphi_met", "low_dphi_met");
  qcd_resolved_cuts.insert("nvlep", "nvlep==0");
  qcd_resolved_cuts.insert("met", "met>150");
  qcd_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  qcd_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  qcd_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  qcd_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  qcd_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");
  qcd_resolved_cuts.insert("mht_filter", "met/mht<2");
  qcd_resolved_cuts.insert("met_calo_filter", "met/met_calo<2");
  //qcd_resolved_cuts.insert("weight", "weight<1.5");

  vector<shared_ptr<Process> > procs;
  // Set mc processes
  //procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  vector<string> sigm = {"175", "500", "950"};
  //vector<string> sigm = {"400"};
  //vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  vector<int> sig_colors = {kSpring, kPink, kAzure, kOrange}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sigm.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
      sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), 
      "stitch"));
  }

  vector<shared_ptr<Process> > procs_data;
  if (unblind) {
    string dataFiles = "*.root";
    procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, data_skim_folder, {dataFiles}),triggers_data));

    procs_data.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, data_skim_folder, {dataFiles}),triggers_data));
  }

  NamedFunc base_filters = "1";
  if (sample_name=="search") base_filters = Higfuncs::final_pass_filters;
  else if (sample_name=="ttbar") base_filters = Higfuncs::final_ttbar_pass_filters;
  else if (sample_name=="zll") base_filters = Higfuncs::final_zll_pass_filters;
  else if (sample_name=="qcd") base_filters = Higfuncs::final_qcd_pass_filters;

  PlotMaker pm;

  torch::OrderedDict<string, Axis> axis_dict;
  axis_dict.insert("ntk",Axis(10, -0.5, 9.5, "ntk", "Number of iso tk", {0.5}));
  axis_dict.insert("npv",Axis(81, -0.5, 80.5, "npv", "Number of reconstructed primary vertices", {}));
  axis_dict.insert("npu_tru_mean",Axis(81, -0.5, 80.5, "npu_tru_mean", "Mean true interaction vertices", {}));
  axis_dict.insert("low_dphi_met",Axis(2, -0.5, 1.5, "low_dphi_met", "Low #Delta #phi", {0.5}));
  axis_dict.insert("nvlep",Axis(6, -0.5, 5.5, "nvlep", "N_{veto leps}", {0.5}));
  axis_dict.insert("nlep",Axis(5, 0.5, 5.5, "nlep", "N_{leps}", {}));
  axis_dict.insert("nmu",Axis(5, 0.5, 5.5, "nmu", "N_{muons}", {}));
  axis_dict.insert("nel",Axis(5, 0.5, 5.5, "nel", "N_{electrons}", {}));
  axis_dict.insert("met",Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}));
  axis_dict.insert("met_detail",Axis(28, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}));
  axis_dict.insert("mht_detail",Axis(28, 150, 850., "mht", "MHT [GeV]", {}));
  axis_dict.insert("ht",Axis(20, 0, 1000., "ht", "HT [GeV]", {}));
  axis_dict.insert("ht_wide",Axis(30, 0, 1500., "ht", "HT [GeV]", {}));
  axis_dict.insert("mht",Axis(20, 0, 500., "mht", "MHT [GeV]", {}));
  axis_dict.insert("mht_phi",Axis(31, -3.14, 3.14, "mht_phi", "MHT_Phi [rad]", {}));
  axis_dict.insert("met_phi",Axis(31, -3.14, 3.14, "met_phi", "MET_Phi [rad]", {}));
  axis_dict.insert("jet_met_dphi",Axis(63, 0, 3.14, "jet_met_dphi", "DeltaR(jet,MET)", {}));
  axis_dict.insert("jet_met_dphi0",Axis(63, 0, 3.14, "jet_met_dphi[0]", "DeltaR(jet,MET)", {}));
  axis_dict.insert("jet_met_dphi1",Axis(63, 0, 3.14, "jet_met_dphi[1]", "DeltaR(jet,MET)", {}));
  axis_dict.insert("jet_met_dphi2",Axis(63, 0, 3.14, "jet_met_dphi[2]", "DeltaR(jet,MET)", {}));
  axis_dict.insert("jet_met_dphi3",Axis(63, 0, 3.14, "jet_met_dphi[3]", "DeltaR(jet,MET)", {}));
  axis_dict.insert("h1b1_jetid",Axis(20, -0.5, 19.5, h1b1_jetid, "Jet id lead higgs high-tag b jet [GeV]", {0.5}));
  axis_dict.insert("h1b2_jetid",Axis(20, -0.5, 19.5, h1b2_jetid, "Jet id lead higgs low-tag b jet [GeV]", {0.5}));
  axis_dict.insert("h2b1_jetid",Axis(20, -0.5, 19.5, h2b1_jetid, "Jet id sublead higgs high-tag b jet [GeV]", {0.5}));
  axis_dict.insert("h2b2_jetid",Axis(20, -0.5, 19.5, h2b2_jetid, "Jet id sublead higgs low-tag b jet [GeV]", {0.5}));
  axis_dict.insert("h1b1_pt",Axis(16, 0, 480., h1b1_pt, "p_{T} lead higgs high-tag b jet [GeV]", {}));
  axis_dict.insert("h1b2_pt",Axis(16, 0, 480., h1b2_pt, "p_{T} lead higgs low-tag b jet [GeV]", {}));
  axis_dict.insert("h2b1_pt",Axis(16, 0, 480., h2b1_pt, "p_{T} sublead higgs high-tag b jet [GeV]", {}));
  axis_dict.insert("h2b2_pt",Axis(16, 0, 480., h2b2_pt, "p_{T} sublead higgs low-tag b jet [GeV]", {}));
  axis_dict.insert("h1_dr",Axis(20,0,4,h1_dr, "Lead higgs #DeltaR", {1.1, 2.2}));
  axis_dict.insert("h2_dr",Axis(20,0,4,h2_dr, "Sublead higgs #DeltaR", {1.1, 2.2}));
  axis_dict.insert("h1_mass",Axis(10, 0, 200, h1_mass, "Lead m_{bb} [GeV]", {100, 140}));
  axis_dict.insert("h2_mass",Axis(10, 0, 200, h2_mass, "Sublead m_{bb} [GeV]", {100, 140}));
  axis_dict.insert("jet_pt[0]",Axis(20, 0, 600., "jet_pt[0]", "p_{T} jet [0] [GeV]", {}));
  axis_dict.insert("jet_pt[1]",Axis(16, 0, 480., "jet_pt[1]", "p_{T} jet [1] [GeV]", {}));
  axis_dict.insert("jet_pt[2]",Axis(16, 0, 480., "jet_pt[2]", "p_{T} jet [2] [GeV]", {}));
  axis_dict.insert("jet_pt[3]",Axis(16, 0, 480., "jet_pt[3]", "p_{T} jet [3] [GeV]", {}));
  axis_dict.insert("jet_pt[4]",Axis(16, 0, 480., "jet_pt[4]", "p_{T} jet [4] [GeV]", {}));
  axis_dict.insert("jet_phi",Axis(31, -3.14, 3.14, "jet_phi", "jet phi [rad]", {}));
  axis_dict.insert("jet_phi[0]",Axis(31, -3.14, 3.14, "jet_phi[0]", "jet [0] phi [rad]", {}));
  axis_dict.insert("jet_phi[1]",Axis(31, -3.14, 3.14, "jet_phi[1]", "jet [1] phi [rad]", {}));
  axis_dict.insert("jet_phi[2]",Axis(31, -3.14, 3.14, "jet_phi[2]", "jet [2] phi [rad]", {}));
  axis_dict.insert("jet_phi[3]",Axis(31, -3.14, 3.14, "jet_phi[3]", "jet [3] phi [rad]", {}));
  axis_dict.insert("jet_phi[4]",Axis(31, -3.14, 3.14, "jet_phi[4]", "jet [4] phi [rad]", {}));
  axis_dict.insert("jet_eta",Axis(50, -5, 5, "jet_eta", "jet eta", {}));
  axis_dict.insert("jet_eta[0]",Axis(50, -5, 5, "jet_eta[0]", "jet [0] eta", {}));
  axis_dict.insert("jet_eta[1]",Axis(50, -5, 5, "jet_eta[1]", "jet [1] eta", {}));
  axis_dict.insert("jet_eta[2]",Axis(50, -5, 5, "jet_eta[2]", "jet [2] eta", {}));
  axis_dict.insert("jet_eta[3]",Axis(50, -5, 5, "jet_eta[3]", "jet [3] eta", {}));
  axis_dict.insert("jet_eta[4]",Axis(50, -5, 5, "jet_eta[4]", "jet [4] eta", {}));
  axis_dict.insert("goodJet_phi",Axis(31, -3.14, 3.14, goodJet_phi, "jet phi [rad]", {}));
  axis_dict.insert("goodJet_eta",Axis(50, -5, 5, goodJet_eta, "jet eta", {}));
  axis_dict.insert("goodJet0_eta",Axis(50, -5, 5, goodJet0_eta, "jet_{0} eta", {}));
  axis_dict.insert("met_zll",Axis(20, 0, 150., "met", "p_{T}^{miss} [GeV]", {30}));
  axis_dict.insert("met_ttbar",Axis(16, 0, 800., "met", "p_{T}^{miss} [GeV]", {}));
  axis_dict.insert("njet",Axis(12, -0.5, 11.5, "njet", "N_{jets}", {3.5, 5.5}));
  axis_dict.insert("hig_cand_drmax",Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}));
  axis_dict.insert("hig_cand_am",Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}));
  axis_dict.insert("hig_cand_dm",Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}));
  axis_dict.insert("btags",Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}));
  //axis_dict.insert("lep_pt",Axis(10, 0, 300., "lep_pt[0]", "p_{l} [GeV]", {30}));
  axis_dict.insert("lep_pt",Axis(10, 0, 300., Higfuncs::lead_signal_lepton_pt, "p_{t}^{l} [GeV]", {30}));
  axis_dict.insert("mt",Axis(10, 0, 200., "mt", "m_{T} [GeV]", {100}));
  axis_dict.insert("dbtags",Axis(5, -0.5, 4.5, "nbm", "N_{b medium}", {}));
  axis_dict.insert("ll_pt",Axis(16, 0, 400., "ll_pt[0]", "p_{T}^{ll} [GeV]", {75., 150., 200., 300}));
  axis_dict.insert("ll_m",Axis(16, 0, 400., "ll_m", "m^{ll} [GeV]", {80., 100.}));
  axis_dict.insert("mht_filter",Axis(10, 0, 10, "met/mht", "MET/MHT", {2}));
  axis_dict.insert("met_calo_filter",Axis(10, 0, 10, "met/met_calo", "MET/MET_calo", {2,5}));
  axis_dict.insert("weight",Axis(80, 0, 20, "weight", "weight", {1.5}));
  vector<string> log_plots = {"met", "btags", "weight", "h1b1_jetid", "h1b2_jetid", "h2b1_jetid", "h1b2_jetid"};


  // Draw n-1
  // Selection according to sample
  torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  if (sample_name == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  else if (sample_name == "zll") map_resolved_cuts = zll_resolved_cuts;
  else if (sample_name == "qcd") map_resolved_cuts = qcd_resolved_cuts;
  else map_resolved_cuts = search_resolved_cuts;

  // Variables to draw
  set<string> target_variables;
  // Add all cut variables
  for (auto const & target_var : map_resolved_cuts) {
    target_variables.insert(target_var.key());
  }
  // Add additional variables
  if (sample_name == "qcd") target_variables.insert("ht_wide");
  else target_variables.insert("ht");
  target_variables.insert("npv");
  if (sample_name=="zll") {
    target_variables.insert("ll_pt");
    //target_variables.insert("ll_m");
  }
  if (sample_name=="ttbar") {
    target_variables.insert("met");
    if (plot_additional_variables) {
      target_variables.insert("jet_pt[0]");
      target_variables.insert("jet_pt[1]");
      target_variables.insert("jet_pt[2]");
      target_variables.insert("jet_pt[3]");
      target_variables.insert("jet_pt[4]");
      target_variables.insert("jet_eta[0]");
      target_variables.insert("jet_eta[1]");
      target_variables.insert("jet_eta[2]");
      target_variables.insert("jet_eta[3]");
      target_variables.insert("jet_eta[4]");
      target_variables.insert("jet_phi[0]");
      target_variables.insert("jet_phi[1]");
      target_variables.insert("jet_phi[2]");
      target_variables.insert("jet_phi[3]");
      target_variables.insert("jet_phi[4]");
    }
  }
  if (plot_additional_variables) {
    target_variables.insert("npu_tru_mean");
    target_variables.insert("h1b1_pt");
    target_variables.insert("h1b2_pt");
    target_variables.insert("h2b1_pt");
    target_variables.insert("h2b2_pt");
    target_variables.insert("h1b1_jetid");
    target_variables.insert("h1b2_jetid");
    target_variables.insert("h2b1_jetid");
    target_variables.insert("h2b2_jetid");
    target_variables.insert("h1_dr");
    target_variables.insert("h2_dr");
    target_variables.insert("h1_mass");
    target_variables.insert("h2_mass");
    target_variables.insert("mht_filter");
    target_variables.insert("jet_pt[0]");
  }

  NamedFunc basic_cut = "1";
  if (sample_name == "search") basic_cut = "(nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)";
  if (HigUtilities::is_in_string_options(string_options, "cut_3b4b")) basic_cut = "(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)";
  //basic_cut = "njet>=4&&hig_cand_am[0]>100&&hig_cand_am[0]<=140&&hig_cand_drmax[0]<=1.1";

  // Loop over variables
  for (auto const & target_var : target_variables) {

    // Make n-1 cut for target_var
    NamedFunc n_minus_1_cut = basic_cut;
    if (sample_name=="search" && unblind && !unblind_signalregion) n_minus_1_cut = "njet>=4&&njet<=5&&!(((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&(hig_cand_am[0]>100&&hig_cand_am[0]<=140))";
    //n_minus_1_cut = n_minus_1_cut && base_resolved;
    for (auto const & cut_item : map_resolved_cuts) {
      // Remove plotting variable
      if (cut_item.key() == target_var) continue;
      // Below are special conditions. Ignore certain cuts.
      if (target_var == "njet") {
        if (cut_item.key() == "hig_cand_drmax") continue;
        if (cut_item.key() == "hig_cand_am") continue;
        if (cut_item.key() == "hig_cand_dm") continue;
      }
      if (target_var == "nlep") {
        if (cut_item.key() == "mt") continue;
        if (cut_item.key() == "lep_pt") continue;
      }
      if (target_var == "sig_el_pt") {
        n_minus_1_cut = n_minus_1_cut && "nel>=1";
      }
      if (target_var == "sig_mu_pt") {
        n_minus_1_cut = n_minus_1_cut && "nmu>=1";
      }
      n_minus_1_cut = n_minus_1_cut && cut_item.value();
    }
    
    // Draw target_var
    // For special case
    if (target_var == "met" && sample_name == "zll") {
      pm.Push<Hist1D>(axis_dict[target_var+"_zll"],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    } else if (target_var == "met" && sample_name == "ttbar") {
      pm.Push<Hist1D>(axis_dict[target_var+"_ttbar"],
        base_filters&&n_minus_1_cut,
        procs, plt_log).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    } else if (target_var == "jet_pt[4]"||target_var == "jet_phi[4]"||target_var == "jet_eta[4]") {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut&&"njet>=5",
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    } else if (target_var == "ht") {
      // print weights
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut,
        procs, plt_lin_print).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // For log plots
    } else if (std::find(log_plots.begin(), log_plots.end(), target_var) != log_plots.end()) {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut,
        procs, plt_log).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")+"_log").LuminosityTag(total_luminosity_string);
    // Normal case
    } else {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }
  }

  NamedFunc resolved_cuts = basic_cut;
  for (auto & item : map_resolved_cuts) {
    resolved_cuts = resolved_cuts && item.value();
  }

  bool plot_ht_correlation = HigUtilities::is_in_string_options(string_options, "plot_ht_correlation");
  if (plot_ht_correlation) {
    // Draw 2D ht vs met (MC)
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["met_detail"],
      base_filters&&resolved_cuts, procs, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_met__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

    if (unblind) {
      // Draw 2D ht vs met (Data)
      pm.Push<Hist2D>(axis_dict["ht"], axis_dict["met_detail"],
        base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_met__data__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

      // Draw 2D ht vs mht (Data)
      pm.Push<Hist2D>(axis_dict["ht"], axis_dict["mht_detail"],
        base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_mht__data__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }

    // Draw 2D ht vs met for signal
    for (unsigned isig(0); isig<sigm.size(); isig++){
      vector<shared_ptr<Process> > procs_signal;
      procs_signal.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::background, 
        sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), 
        "stitch"));
      pm.Push<Hist2D>(axis_dict["ht"], axis_dict["met_detail"],
        base_filters&&resolved_cuts, procs_signal, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_sig"+sigm[isig]+"_ht_vs_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }

    // met vs eta
    pm.Push<Hist2D>(axis_dict["met"], axis_dict["goodJet_eta"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_met_vs_goodJetEta__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

    // ht vs mht
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["mht"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_mht__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs mht_phi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["mht_phi"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_mht_phi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs met_phi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["met_phi"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_met_phi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs hig_cand_drmax
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["hig_cand_drmax"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_hig_cand_drmax__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs jet_met_dphi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_met_dphi"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jet_met_dphi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs jet_met_dphi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_met_dphi0"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jet_met_dphi0__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_met_dphi1"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jet_met_dphi1__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_met_dphi2"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jet_met_dphi2__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_met_dphi3"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jet_met_dphi3__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  }

  bool plot_eta_vs_phi = HigUtilities::is_in_string_options(string_options, "plot_eta_vs_phi");
  if (plot_eta_vs_phi) {
     //Draw 2D jet eta vs phi 
    pm.Push<Hist2D>(axis_dict["jet_phi"], axis_dict["jet_eta"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jetEta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs eta
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_eta"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jetEta__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs phi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["jet_phi"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // Draw 2D jet eta vs phi 
    pm.Push<Hist2D>(axis_dict["goodJet_phi"], axis_dict["goodJet_eta"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_goodJetEta_vs_goodJetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs eta
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["goodJet0_eta"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_goodJetEta__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs phi
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["goodJet_phi"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_goodJetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // ht vs njet
    pm.Push<Hist2D>(axis_dict["ht"], axis_dict["njet"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_njet__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["jet_phi[0]"], axis_dict["jet_eta[0]"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jet0Eta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["jet_phi[1]"], axis_dict["jet_eta[1]"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jet1Eta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["jet_phi[2]"], axis_dict["jet_eta[2]"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jet2Eta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["jet_phi[3]"], axis_dict["jet_eta[3]"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jet3Eta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(axis_dict["jet_phi[4]"], axis_dict["jet_eta[4]"],
      base_filters&&resolved_cuts, procs_data, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_jet4Eta_vs_jetPhi__mc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  }

  bool plot_in_bins = HigUtilities::is_in_string_options(string_options, "plot_in_bins");
  if (plot_in_bins) {
    // Plot in bins
    vector<pair<string, NamedFunc> > nb_bins = {{"2b","nbt==2&&nbm==2"},{"3b","nbt>=2&&nbm==3&&nbl==3"},{"4b","nbt>=2&&nbm>=3&&nbl>=4"}};
    vector<pair<string, NamedFunc> > mass_bins = {{"SDB","hig_cand_am[0]<=100||hig_cand_am[0]>140"},{"HIG","hig_cand_am[0]>100&&hig_cand_am[0]<=140"}};
    vector<pair<string, NamedFunc> > met_bins = {{"MET75", "met<=75"},{"75MET150","met>75&&met<=150"},{"150MET200","met>150&&met<=200"},{"200MET", "met>200"}};
    vector<pair<string, NamedFunc> > drmax_bins = {{"low-drmax","hig_cand_drmax[0]<=1.1"}, {"high-drmax","hig_cand_drmax[0]>1.1"}};
    // combine bins
    vector<pair<string, NamedFunc> > nb_mass_bins;
    combine_bins(nb_mass_bins, nb_bins, mass_bins);
    vector<pair<string, NamedFunc> > met_nb_mass_bins;
    combine_bins(met_nb_mass_bins, met_bins, nb_bins, mass_bins);
    vector<pair<string, NamedFunc> > met_drmax_bins;
    combine_bins(met_drmax_bins, met_bins, drmax_bins);
    vector<pair<string, NamedFunc> > met_drmax_nb_mass_bins;
    combine_bins(met_drmax_nb_mass_bins, met_drmax_bins, nb_bins, mass_bins);
    NamedFunc resolved_nodrmax_cuts = basic_cut;
    for (auto & cut_item : map_resolved_cuts) {
      if (cut_item.key() == "hig_cand_drmax") continue;
      resolved_nodrmax_cuts = resolved_nodrmax_cuts && cut_item.value();
    }
    // Plot drmax in each bin
    for (auto const & bin : met_nb_mass_bins) {
      pm.Push<Hist1D>(axis_dict["hig_cand_drmax"],base_filters&&resolved_nodrmax_cuts&&bin.second, procs, plt_lin).Weight(weight).Tag("FixName:n-1_drmax_"+bin.first+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }
  }

  //// Example
  //NamedFunc resolved_cuts = "1";
  //torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  //if (sample_name == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  //else if (sample_name == "zll") map_resolved_cuts = zll_resolved_cuts;
  //else if (sample_name == "qcd") map_resolved_cuts = qcd_resolved_cuts;
  //else map_resolved_cuts = search_resolved_cuts;
  //for (auto & item : map_resolved_cuts) {
  //  resolved_cuts = resolved_cuts && item.value();
  //}
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&resolved_cuts&&"hig_cand_drmax[0]<1.1",
  //  procs, plt_lin).Weight(weight).Tag("FixName:n-1_hig_cand_am_2016");

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
      {"unblind_signalregion", no_argument, 0, 'a'},
      {"unblind", no_argument, 0, 'u'},
      {"sample", required_argument, 0, 0},
      {"year", required_argument, 0, 'y'},
      {"string_options", required_argument, 0, 'o'},
      {"lepton_type", required_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "so:uy:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 'u':
      unblind = true;
      break;
    case 'a':
      unblind_signalregion = true;
      break;
    case 'y':
      year_string = optarg;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "sample"){
        sample_name = optarg;
      } else if (optname == "lepton_type") {
        lepton_type = optarg;
      } else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

