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
  string year_string = "run2";
  string string_options = "plot_kappa";
  // possible string_options
  // plot_notrgeff - plot studies of am with no trigger efficiencies applied
  // plot_kappa - plot MC kappas for the search region
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

  // Set options
  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";

  string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string data_skim_folder = "data/merged_higdata_higloose/";

  string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrections/merged_higmc_higloose/";

  //years = {2016, 2017, 2018};
  //years = {2016};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc weight_notrgeff = Higfuncs::final_weight_notrgeff;

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

  vector<shared_ptr<Process> > procs;
  // Set mc processes
  procs.push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
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

  //// Set signal mc
  vector<string> signal_mass = {"200", "450", "700", "950"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need signal_mass.size() >= sig_colors.size()
  for (unsigned isig(0); isig<signal_mass.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+signal_mass[isig]+",1)", Process::Type::signal, 
      sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+signal_mass[isig]+"_mLSP-0*.root"}), "1"));
  }

  // Set processes according to btag
  // Getting colors
  // TColor * color;; float red, green, blue;
  // color = gROOT->GetColor(kAzure+1); color->GetRGB(red,green,blue); cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl;
  vector<shared_ptr<Process> > procs_btag;
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));

  // Set processes according to true number of b
  vector<shared_ptr<Process> > procs_trueB;
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, colors("true_3b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),  "stitch"&& Functions::ntrub==3));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("4 B-hadrons", Process::Type::background, colors("true_4b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==4));


  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...
  NamedFunc base_filters = Higfuncs::final_pass_filters;

  // resolved cuts
  NamedFunc base_resolved = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";

  PlotMaker pm;

  
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    base_filters&&base_resolved, procs_btag, plt_shapes).Weight(weight).Tag("FixName:bkgest__amjj_shapes").LuminosityTag(total_luminosity_string);

  pm.Push<Table>("FixName:pie_bkgest__trueb_pies_onemet__2b", vector<TableRow> ({TableRow("", base_filters&&base_resolved&&"(nbt==2&&nbm==2)", 0, 0, weight)}), procs_trueB, true, true, true);
  pm.Push<Table>("FixName:pie_bkgest__trueb_pies_onemet__3b", vector<TableRow> ({TableRow("", base_filters&&base_resolved&&"(nbt>=2&&nbm==3&&nbl==3)", 0, 0, weight)}), procs_trueB, true, true, true);
  pm.Push<Table>("FixName:pie_bkgest__trueb_pies_onemet__4b", vector<TableRow> ({TableRow("", base_filters&&base_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)", 0, 0, weight)}), procs_trueB, true, true, true);

  if (HigUtilities::is_in_string_options(string_options,"plot_notrgeff")) {
    // No trigger study
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&base_resolved, procs_btag, plt_shapes).Weight(weight_notrgeff).Tag("FixName:bkgest__amjj_shapes_notrgeff").LuminosityTag(total_luminosity_string);
    // Plane MET vs drmax
    vector<NamedFunc> plane_cuts;
    plane_cuts.push_back("met>150&&met<=200 && hig_cand_drmax[0]<=1.1");
    plane_cuts.push_back("met>200&&met<=300 && hig_cand_drmax[0]<=1.1");
    plane_cuts.push_back("met>300&&met<=400 && hig_cand_drmax[0]<=1.1");
    plane_cuts.push_back("met>400           && hig_cand_drmax[0]<=1.1");
    plane_cuts.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");
    plane_cuts.push_back("met>200&&met<=300 && hig_cand_drmax[0]>1.1");
    plane_cuts.push_back("met>300&&met<=400 && hig_cand_drmax[0]>1.1");
    plane_cuts.push_back("met>400           && hig_cand_drmax[0]>1.1");
    for (unsigned iPlane = 0; iPlane < plane_cuts.size(); ++iPlane) {
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&base_resolved&&plane_cuts[iPlane], procs_btag, plt_shapes).Weight(weight).Tag("FixName:bkgest__amjj_shapes_plane"+to_string(iPlane)).LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&base_resolved&&plane_cuts[iPlane], procs_btag, plt_shapes).Weight(weight_notrgeff).Tag("FixName:bkgest__amjj_shapes_plane"+to_string(iPlane)+"_notrgeff").LuminosityTag(total_luminosity_string);
    }
  }

  //pm.Push<Table>("syst__ttbar_pies__2b_met300", vector<TableRow> ({TableRow("", base_filters&&base_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), procs, true, true, true);

  // kappa plot
  if (HigUtilities::is_in_string_options(string_options,"plot_kappa")) {
    system(("./run/higgsino/plot_kappas.exe --sample search --scen mc -o paper_style --year "+year_string).c_str());
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
