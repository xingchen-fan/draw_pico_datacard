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
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "run2";
}

vector<pair<int, int> > parseMassPoints(string const & mass_points_string) {
  vector<string> v_mass_point;
  HigUtilities::stringToVectorString(mass_points_string, v_mass_point, ",");
  vector<pair<int, int> > mass_points;
  for (auto const & mass_point : v_mass_point) {
    vector<string> v_mass;
    HigUtilities::stringToVectorString(mass_point, v_mass, "_");
    int nlsp_mass = stoi(v_mass[0]); int lsp_mass = stoi(v_mass[1]);
    mass_points.push_back({nlsp_mass, lsp_mass});
  }
  return mass_points;
}

// preselect:
//   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
//   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
// higloose: 
//   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
//   met>150 && nvlep==0
// higlep1T:
//   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
//   nlep==1 && 
//   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) 
// higlep2T:
//   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
//   nlep==2 && 
//   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
//   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
// higqcd with met150:
//   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
//   nvlep==0 && ntk==0 && low_dphi_met &&
//   met>150  // Since applied to met150 skim
void setProcsDict(string const & production, string const & nanoAodFolder, string const & year_string, string const & sample_name, map<string, vector<shared_ptr<Process> > > & procsDict) {
  set<int> years;
  //years_a = {2016};
  HigUtilities::parseYears(year_string, years);

  // Data cuts
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  NamedFunc triggers_data = "1";
  if (sample_name == "zll") triggers_data = lepton_triggers;
  else if (sample_name == "ttbar") triggers_data = lepton_triggers || met_triggers;
  else if (sample_name == "qcd") triggers_data = met_triggers;
  else  triggers_data = met_triggers;

  // Set folders
  string mc_production_folder = nanoAodFolder+"/"+production;
  string data_production_folder = nanoAodFolder+"/"+production;
  string signal_production_folder = nanoAodFolder+"/"+production;
  string mc_skim_folder, data_skim_folder, signal_skim_folder;
  if (sample_name == "search") {
    mc_skim_folder = "mc/merged_higmc_higloose/";
    data_skim_folder = "data/merged_higdata_higloose/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
  } else if (sample_name == "ttbar") {
    mc_skim_folder = "mc/merged_higmc_higlep1T/";
    data_skim_folder = "data/merged_higdata_higlep1T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
  } else if (sample_name == "zll") {
    mc_skim_folder = "mc/merged_higmc_higlep2T/";
    data_skim_folder = "data/merged_higdata_higlep2T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
  } else if (sample_name == "qcd") {
    mc_skim_folder = "mc/merged_higmc_higqcd/";
    data_skim_folder = "data/merged_higdata_higqcd/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";
  }

  // Set mc file names
  map<string, set<string>> mc_filenames; 
  mc_filenames["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mc_filenames["single_t"] = set<string>({"*_ST_*.root"});
  mc_filenames["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mc_filenames["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mc_filenames["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mc_filenames["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mc_filenames["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  // Set signal file names
  string mass_points_string = "175_0, 500_0, 950_0";
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // Requires mass_points.size() >= sig_colors.size()
  // signal_filenames["TChiHH(nlsp,lsp)"] = {filename}
  torch::OrderedDict<string, set<string> > signal_filenames;
  // mass_points = [[nlsp, lsp]]
  vector<pair<int, int> > mass_points = parseMassPoints(mass_points_string);
  for (auto const & mass_point : mass_points) {
    signal_filenames.insert("TChiHH("+to_string(mass_point.first)+","+to_string(mass_point.second)+")", {"*TChiHH_mChi-"+to_string(mass_point.first)+"_mLSP-"+to_string(mass_point.second)+"_*.root"});
  }

  // Set mc background+signal procs
  Palette colors("txt/colors.txt", "default");
  procsDict["mc"];
  //procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["single_t"]),"stitch"));
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["other"]),"stitch"));

  procsDict["mcTtbar"];
  //procsDict["mcTtbar"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  procsDict["mcTtbar"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  procsDict["mcTtbar"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));

  procsDict["mc_and_sig"];
  //procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["single_t"]),"stitch"));
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["other"]),"stitch"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key(), Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "1"));
  }

  // Set mc btag procs
  procsDict["mc_btag"];
  procsDict["mc_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt==2&&nbm==2)"));
  procsDict["mc_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  procsDict["mc_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));

  // Set mc isr procs
  procsDict["mc_nisr"];
  procsDict["mc_nisr"].push_back(Process::MakeShared<Baby_pico>("N_{ISR jets}=0", Process::Type::background,kBlack,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&nisr==0"));
  procsDict["mc_nisr"].push_back(Process::MakeShared<Baby_pico>("N_{ISR jets}=1", Process::Type::background,kPink+5,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&nisr==1"));
  procsDict["mc_nisr"].push_back(Process::MakeShared<Baby_pico>("N_{ISR jets}=2", Process::Type::background,kPink+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&nisr==2"));

  // Set mc met procs
  procsDict["mcTtbar_met"];
  procsDict["mcTtbar_met"].push_back(Process::MakeShared<Baby_pico>("150<p_{T}^{miss} #leq 200", Process::Type::background,kGreen,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>150&&met<=200"));
  procsDict["mcTtbar_met"].push_back(Process::MakeShared<Baby_pico>("200<p_{T}^{miss} #leq 300", Process::Type::background,kBlue,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>200&&met<=300"));
  procsDict["mcTtbar_met"].push_back(Process::MakeShared<Baby_pico>("300<p_{T}^{miss} #leq 400", Process::Type::background,kOrange,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>300&&met<=400"));
  procsDict["mcTtbar_met"].push_back(Process::MakeShared<Baby_pico>("400<p_{T}^{miss}", Process::Type::background,kRed,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>400"));

  procsDict["mcTtbar_lowmet"];
  procsDict["mcTtbar_lowmet"].push_back(Process::MakeShared<Baby_pico>("0<p_{T}^{miss} #leq 75", Process::Type::background,kGreen,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>0&&met<=75"));
  procsDict["mcTtbar_lowmet"].push_back(Process::MakeShared<Baby_pico>("75<p_{T}^{miss} #leq 150", Process::Type::background,kBlue,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>75&&met<=150"));
  procsDict["mcTtbar_lowmet"].push_back(Process::MakeShared<Baby_pico>("150<p_{T}^{miss} #leq 200", Process::Type::background,kOrange,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>150&&met<=200"));
  procsDict["mcTtbar_lowmet"].push_back(Process::MakeShared<Baby_pico>("200<p_{T}^{miss}", Process::Type::background,kRed,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&met>200"));
  
  // Set mc and data
  procsDict["mc_and_data"];
  //procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["single_t"]),"stitch"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["other"]),"stitch"));
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, 
                   attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}), triggers_data));

  // Set data btag procs
  procsDict["data_btag"];
  procsDict["data_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt==2&&nbm==2)"));
  procsDict["data_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt>=2&&nbm==3&&nbl==3)"));
  procsDict["data_btag"].push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt>=2&&nbm>=3&&nbl>=4)"));

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
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  // Set folders according
  // production, nanoAODFolder, sample_name, year_string, 
  // Set procs
  // Set baseline, filter according: sample_name
  //string production = "higgsino_klamath"; 
  //string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7";
  string production = "higgsino_klamath"; 
  string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7";

  // [mc, mcTtbar, mc_and_sig, mc_btag, mc_nisr, mcTtbar_met, mcTtbar_lowmet, mc_and_data, data_btag]
  map<string, vector<shared_ptr<Process> > > procs_search;
  setProcsDict(production, nanoAodFolder, year_string, "search", procs_search); 
  map<string, vector<shared_ptr<Process> > > procs_ttbar;
  setProcsDict(production, nanoAodFolder, year_string, "ttbar", procs_ttbar); 

  NamedFunc weight = Higfuncs::final_weight; //"weight"*eff_higtrig_run2*w_years*Functions::w_pileup;
  NamedFunc weight_notrgeff = Higfuncs::final_weight_notrgeff; //"weight"*w_years*Functions::w_pileup

  // Filters for each sample
  NamedFunc search_filters = Higfuncs::final_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2&&weight<1.5"
  NamedFunc ttbar_filters = Higfuncs::final_ttbar_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc zll_filters = Higfuncs::final_zll_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc qcd_filters = Higfuncs::final_qcd_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2"

  // resolved cuts
  NamedFunc search_resolved_cuts = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc trk_resolved_cuts = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk>=1&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved_cuts = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         && Higfuncs::lead_signal_lepton_pt>30;
  NamedFunc zll_resolved_cuts =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<=50&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved_cuts =
                         "met/mht<2 && met/met_calo<2&&"
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<=2.2&&hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  vector<NamedFunc> bin_cuts_search;
  bin_cuts_search.push_back("met>150&&met<=200 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>400           && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]>1.1");
  bin_cuts_search.push_back("met>400           && hig_cand_drmax[0]>1.1");

  vector<NamedFunc> bin_cuts_ttbar;
  bin_cuts_ttbar.push_back("met>0&&met<=75 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>75&&met<=150 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>150&&met<=200 && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>200           && hig_cand_drmax[0]<=1.1");
  bin_cuts_ttbar.push_back("met>0&&met<=75 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>75&&met<=150 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");
  bin_cuts_ttbar.push_back("met>200           && hig_cand_drmax[0]>1.1");


  PlotMaker pm;

  // Show explicitly why kappa is not 1
  bool plot_explicit = true;
  if (plot_explicit) {
    bool explicit_detail = true;
    if (explicit_detail) {
      // average mass; 2b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_lowdrmax__search").LuminosityTag(total_luminosity_string);
      // average mass; 3b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_lowdrmax__search").LuminosityTag(total_luminosity_string);
      // average mass; 4b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_lowdrmax__search").LuminosityTag(total_luminosity_string);

      // average mass; 2b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_highdrmax__search").LuminosityTag(total_luminosity_string);
      // average mass; 3b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_highdrmax__search").LuminosityTag(total_luminosity_string);
      // average mass; 4b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_highdrmax__search").LuminosityTag(total_luminosity_string);
    }

    // average mass (mc 2bvs3bvs4b); low-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax__search").LuminosityTag(total_luminosity_string);
    // average mass (mc 2bvs3bvs4b); high-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax__search").LuminosityTag(total_luminosity_string);

    if (explicit_detail) {
      // average mass; 2b,low-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]<=1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_lowdrmax__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 3b,low-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]<=1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_lowdrmax__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 4b,low-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]<=1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_lowdrmax__ttbar").LuminosityTag(total_luminosity_string);

      // average mass; 2b,high-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]>1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_highdrmax__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 3b,high-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]>1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_highdrmax__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 4b,high-drmax [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]>1.1", procs_ttbar["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_highdrmax__ttbar").LuminosityTag(total_luminosity_string);
    }

    // average mass (mc 2bvs3bvs4b); low-drmax [ttbar] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax__ttbar").LuminosityTag(total_luminosity_string);
    // average mass (mc 2bvs3bvs4b); high-drmax [ttbar] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax__ttbar").LuminosityTag(total_luminosity_string);

    // average mass (data 2bvs3bvs4b); low-drmax [ttbar] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__lowdrmax__ttbar").LuminosityTag(total_luminosity_string);
    // average mass (data 2bvs3bvs4b); high-drmax [ttbar] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__highdrmax__ttbar").LuminosityTag(total_luminosity_string);

    if (explicit_detail) {
      // average mass; 2b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_lowdrmax__trk").LuminosityTag(total_luminosity_string);
      // average mass; 3b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_lowdrmax__trk").LuminosityTag(total_luminosity_string);
      // average mass; 4b,low-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]<=1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_lowdrmax__trk").LuminosityTag(total_luminosity_string);

      // average mass; 2b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt==2&&nbm==2)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_highdrmax__trk").LuminosityTag(total_luminosity_string);
      // average mass; 3b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_highdrmax__trk").LuminosityTag(total_luminosity_string);
      // average mass; 4b,high-drmax [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&trk_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_drmax[0]>1.1", procs_search["mc_and_sig"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_highdrmax__trk").LuminosityTag(total_luminosity_string);
    }

    // average mass (mc 2bvs3bvs4b); low-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&trk_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax__trk").LuminosityTag(total_luminosity_string);
    // average mass (mc 2bvs3bvs4b); high-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&trk_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax__trk").LuminosityTag(total_luminosity_string);

    // average mass (data 2bvs3bvs4b); low-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&trk_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_search["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__lowdrmax__trk").LuminosityTag(total_luminosity_string);
    // average mass (data 2bvs3bvs4b); high-drmax [search] 
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&trk_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_search["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__highdrmax__trk").LuminosityTag(total_luminosity_string);

  }

  // Show kappa relation with ISR
  bool plotKappaISRRelation = false;
  if (plotKappaISRRelation) {
    // # isr jets (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
      search_filters&&search_resolved_cuts, procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__search").LuminosityTag(total_luminosity_string);
    // average mass (mc 0isrvs1isrvs2isr) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts, procs_search["mc_nisr"], plt_shapes).Weight(weight).Tag("FixName:amjj_0isr1isr2isr__search").LuminosityTag(total_luminosity_string);
    // average mass (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts, procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__search").LuminosityTag(total_luminosity_string);

    //// Isn't the reason low-drmax and high-drmax mbb is different
    //// # isr jets (mc 2bvs3bvs4b);low-drmax [search]
    //pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
    //  search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__lowdrmax__search").LuminosityTag(total_luminosity_string);
    //// # isr jets (mc 2bvs3bvs4b);high-drmax [search]
    //pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
    //  search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__highdrmax__search").LuminosityTag(total_luminosity_string);

    // # isr jets (mc 2bvs3bvs4b) [ttbar]
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
      ttbar_filters&&ttbar_resolved_cuts, procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__ttbar").LuminosityTag(total_luminosity_string);
    // average mass (mc 0isrvs1isrvs2isr) [ttbar]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts, procs_ttbar["mc_nisr"], plt_shapes).Weight(weight).Tag("FixName:amjj_0isr1isr2isr__ttbar").LuminosityTag(total_luminosity_string);
    // average mass (mc 2bvs3bvs4b) [ttbar]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_filters&&ttbar_resolved_cuts, procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__ttbar").LuminosityTag(total_luminosity_string);

    //// Isn't the reason low-drmax and high-drmax mbb is different
    //// # isr jets (mc 2bvs3bvs4b);low-drmax [ttbar]
    //pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
    //  ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__lowdrmax__ttbar").LuminosityTag(total_luminosity_string);
    //// # isr jets (mc 2bvs3bvs4b);high-drmax [ttbar]
    //pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
    //  ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:nisr_2b3b4b__highdrmax__ttbar").LuminosityTag(total_luminosity_string);
  }

  bool plotMetRelation = false;
  if (plotMetRelation) {
    bool metRelationIsrDetail = false;
    if (metRelationIsrDetail) {
      // Show kappa relation with MET
      // # isr jets; 150<met<=200 [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"met>150&&met<=200", procs_search["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_150met200__search").LuminosityTag(total_luminosity_string);
      // # isr jets; 200<met<=300 [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"met>200&&met<=300", procs_search["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_200met300__search").LuminosityTag(total_luminosity_string);
      // # isr jets; 300<met<=400 [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"met>300&&met<=400", procs_search["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_300met400__search").LuminosityTag(total_luminosity_string);
      // # isr jets; 400<met      [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"met>400", procs_search["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_400met__search").LuminosityTag(total_luminosity_string);

      // # isr jets; 0<met<=75 [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>0&&met<=75", procs_ttbar["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_0met75__ttbar").LuminosityTag(total_luminosity_string);
      // # isr jets; 75<met<=150 [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>75&&met<=150", procs_ttbar["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_75met150__ttbar").LuminosityTag(total_luminosity_string);
      // # isr jets; 150<met<=200 [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>150&&met<=200", procs_ttbar["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_150met200__ttbar").LuminosityTag(total_luminosity_string);
      // # isr jets; 200<met [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>200", procs_ttbar["mcTtbar"], plt_lin).Weight(weight).Tag("FixName:nisr_200met__ttbar").LuminosityTag(total_luminosity_string);

      // # isr jets (mc met) [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts, procs_search["mcTtbar_met"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_met__search").LuminosityTag(total_luminosity_string);

      // # isr jets (mc met); low-drmax [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_search["mcTtbar_met"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_met__lowdrmax__search").LuminosityTag(total_luminosity_string);

      // # isr jets (mc met); high-drmax [search]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_search["mcTtbar_met"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_met__highdrmax__search").LuminosityTag(total_luminosity_string);

      // # isr jets (mc lowmet) [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts, procs_ttbar["mcTtbar_lowmet"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_lowmet__ttbar").LuminosityTag(total_luminosity_string);

      // # isr jets (mc lowmet); low-drmax [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1", procs_ttbar["mcTtbar_lowmet"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_lowmet__lowdrmax__ttbar").LuminosityTag(total_luminosity_string);

      // # isr jets (mc lowmet); high-drmax [ttbar]
      pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nisr", "N_{ISR jets}", {}),
        ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1", procs_ttbar["mcTtbar_lowmet"], plt_shapes_info).Weight(weight).Tag("FixName:nisr_lowmet__highdrmax__ttbar").LuminosityTag(total_luminosity_string);
    }

    // average mass; 150<met<=200 (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"met>150&&met<=200", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__150met200__search").LuminosityTag(total_luminosity_string);
    // average mass; 200<met<=300 (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"met>200&&met<=300", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__200met300__search").LuminosityTag(total_luminosity_string);
    // average mass; 300<met<=400 (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"met>300&&met<=400", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__300met400__search").LuminosityTag(total_luminosity_string);
    // average mass; 400<met      (mc 2bvs3bvs4b) [search]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&"met>400"          , procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__400met__search").LuminosityTag(total_luminosity_string);

    bool metRelationDetail = false;
    if (metRelationDetail) {
      // average mass; low-drmax,150<met<=200 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1&&met>150&&met<=200", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax_150met200__search").LuminosityTag(total_luminosity_string);
      // average mass; low-drmax,200<met<=300 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1&&met>200&&met<=300", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax_200met300__search").LuminosityTag(total_luminosity_string);
      // average mass; low-drmax,300<met<=400 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1&&met>300&&met<=400", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax_300met400__search").LuminosityTag(total_luminosity_string);
      // average mass; low-drmax,400<met      (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]<=1.1&&met>400"          , procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__lowdrmax_400met__search").LuminosityTag(total_luminosity_string);

      // average mass; high-drmax,150<met<=200 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1&&met>150&&met<=200", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax_150met200__search").LuminosityTag(total_luminosity_string);
      // average mass; high-drmax,200<met<=300 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1&&met>200&&met<=300", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax_200met300__search").LuminosityTag(total_luminosity_string);
      // average mass; high-drmax,300<met<=400 (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1&&met>300&&met<=400", procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax_300met400__search").LuminosityTag(total_luminosity_string);
      // average mass; high-drmax,400<met      (mc 2bvs3bvs4b) [search]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&"hig_cand_drmax[0]>1.1&&met>400"          , procs_search["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__highdrmax_400met__search").LuminosityTag(total_luminosity_string);

      for (unsigned iBin = 0; iBin < bin_cuts_search.size(); ++iBin) {
        // average mass; 2b,bin_cut [search]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          search_filters&&search_resolved_cuts&&"(nbt==2&&nbm==2)"&&bin_cuts_search[iBin], procs_search["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_bin"+to_string(iBin)+"__search").LuminosityTag(total_luminosity_string);
        // average mass; 3b,bin_cut [search]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)"&&bin_cuts_search[iBin], procs_search["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_bin"+to_string(iBin)+"__search").LuminosityTag(total_luminosity_string);
        // average mass; 4b,bin_cut [search]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          search_filters&&search_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)"&&bin_cuts_search[iBin], procs_search["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_bin"+to_string(iBin)+"__search").LuminosityTag(total_luminosity_string);
      }

      // average mass; 0<met<=75   (mc 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>0&&met<=75", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__0met75__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 75<met<=150 (mc 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>75&&met<=150", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__75met150__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 150<met<=200 (mc 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>150&&met<=200", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__150met200__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 200<met     (mc 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>200", procs_ttbar["mc_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_2b3b4b__200met__ttbar").LuminosityTag(total_luminosity_string);

      // average mass; 0<met<=75   (data 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>0&&met<=75", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__0met75__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 75<met<=150 (data 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>75&&met<=150", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__75met150__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 150<met<=200 (data 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>150&&met<=200", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__150met200__ttbar").LuminosityTag(total_luminosity_string);
      // average mass; 200<met     (data 2bvs3bvs4b) [ttbar]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        ttbar_filters&&ttbar_resolved_cuts&&"met>200", procs_ttbar["data_btag"], plt_shapes).Weight(weight).Tag("FixName:amjj_data2b3b4b__200met__ttbar").LuminosityTag(total_luminosity_string);

      for (unsigned iBin = 0; iBin < bin_cuts_ttbar.size(); ++iBin) {
        // average mass; 2b,bin_cut [ttbar]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          ttbar_filters&&ttbar_resolved_cuts&&"(nbt==2&&nbm==2)"&&bin_cuts_ttbar[iBin], procs_ttbar["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__2b_bin"+to_string(iBin)+"__ttbar").LuminosityTag(total_luminosity_string);
        // average mass; 3b,bin_cut [ttbar]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm==3&&nbl==3)"&&bin_cuts_ttbar[iBin], procs_ttbar["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__3b_bin"+to_string(iBin)+"__ttbar").LuminosityTag(total_luminosity_string);
        // average mass; 4b,bin_cut [ttbar]
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          ttbar_filters&&ttbar_resolved_cuts&&"(nbt>=2&&nbm>=3&&nbl>=4)"&&bin_cuts_ttbar[iBin], procs_ttbar["mc"], plt_lin).Weight(weight).Tag("FixName:amjj__4b_bin"+to_string(iBin)+"__ttbar").LuminosityTag(total_luminosity_string);
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
