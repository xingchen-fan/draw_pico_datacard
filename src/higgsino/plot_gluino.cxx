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
  bool unblind = false;
  bool unblind_signalregion = false;
}

const NamedFunc deepcsv_2016("deepcsv_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  //for (auto const & process : processes_) {
  //  if (Contains(process.name_, "bbbb")) {
  //  }
  //}
  if (Higfuncs::hig_bcat.GetScalar(b)==2) {
    weight = 2 - (2-1.34)/2000 * Higfuncs::hig_pt1.GetScalar(b);
    if (weight<1) weight = 1;
  }
  //for (unsigned iJet = 0; iJet < b.jet_deepcsv()->size(); ++iJet) {
  //  if      ((*b.jet_deepcsv())[iJet]>0   && (*b.jet_deepcsv())[iJet]<=0.1) weight *= 170.204/131.958;
  //  else if ((*b.jet_deepcsv())[iJet]>0.1 && (*b.jet_deepcsv())[iJet]<=0.2) weight *= 48.3374/37.9629;
  //  else if ((*b.jet_deepcsv())[iJet]>0.2 && (*b.jet_deepcsv())[iJet]<=0.3) weight *= 25.2594/20.08;
  //  else if ((*b.jet_deepcsv())[iJet]>0.3 && (*b.jet_deepcsv())[iJet]<=0.4) weight *= 20.7668/17.3701;
  //  else if ((*b.jet_deepcsv())[iJet]>0.4 && (*b.jet_deepcsv())[iJet]<=0.5) weight *= 15.1035/12.3414;
  //  else if ((*b.jet_deepcsv())[iJet]>0.5 && (*b.jet_deepcsv())[iJet]<=0.6) weight *= 17.8889/15.447;
  //  else if ((*b.jet_deepcsv())[iJet]>0.6 && (*b.jet_deepcsv())[iJet]<=0.7) weight *= 16.513/14.9131;
  //  else if ((*b.jet_deepcsv())[iJet]>0.7 && (*b.jet_deepcsv())[iJet]<=0.8) weight *= 19.9071/18.6731;
  //  else if ((*b.jet_deepcsv())[iJet]>0.8 && (*b.jet_deepcsv())[iJet]<=0.9) weight *= 33.2177/32.1704;
  //  else if ((*b.jet_deepcsv())[iJet]>0.9 && (*b.jet_deepcsv())[iJet]<=1.0) weight *= 328.241/281.812;
  //}
  //weight = pow(weight, 1./b.jet_deepcsv()->size());
  return weight;
});

const NamedFunc only_hbb("only_hbb",[](const Baby &b) -> NamedFunc::ScalarType{
  bool hbb = true;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_mom()->at(i)!=25) continue;
    if (fabs(b.mc_id()->at(i))!=5) hbb=false;
  }
  return hbb;
});

// 0: bbbb, 1: bbww, 2: bbgluongluon, 3: bbtautau, 4: bbcc, 5: bbzz, 6: other
const NamedFunc hh_decay("hh_decay",[](const Baby &b) -> NamedFunc::ScalarType{
  // Collect hh daughters
  vector<int> hh_daughters;
  hh_daughters.reserve(4);
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_mom()->at(i)!=25) continue;
    hh_daughters.push_back(b.mc_id()->at(i));
  }
  // Count
  int nb = 0, nw = 0, ngluon = 0, ntau = 0, nc = 0, nz = 0;
  for (unsigned iDaughter = 0; iDaughter < hh_daughters.size(); ++iDaughter) {
    if (fabs(hh_daughters[iDaughter]) == 5) nb++;
    else if (fabs(hh_daughters[iDaughter]) == 24) nw++;
    else if (fabs(hh_daughters[iDaughter]) == 21) ngluon++;
    else if (fabs(hh_daughters[iDaughter]) == 15) ntau++;
    else if (fabs(hh_daughters[iDaughter]) == 4) nc++;
    else if (fabs(hh_daughters[iDaughter]) == 23) nz++;
  }
  // Tag
  int decay = 6;
  if (nb==4) decay = 0;
  else if (nb==2 && nw==2) decay = 1;
  else if (nb==2 && ngluon ==2) decay = 2;
  else if (nb==2 && ntau == 2) decay = 3;
  else if (nb==2 && nc == 2) decay = 4;
  else if (nb==2 && nz == 2) decay = 5;
  return decay;
});

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
  string gluino_production_folder = nanoAodFolder+"/"+production;
  string mc_skim_folder, data_skim_folder, signal_skim_folder, gluino_fast_skim_folder, gluino_full_skim_folder;
  if (sample_name == "search") {
    //// search
    //mc_skim_folder = "mc/merged_higmc_higloose/";
    //data_skim_folder = "data/merged_higdata_higloose/";
    //signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
    //gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higloose/";
    //gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_higloose/";
    // preselect
    mc_skim_folder = "mc/merged_higmc_preselect/";
    data_skim_folder = "data/merged_higdata_preselect/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
    gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_preselect/";
    gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_preselect/";
    //// met150
    //mc_skim_folder = "mc/merged_higmc_met150/";
    //data_skim_folder = "data/merged_higdata_met150/";
    //signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_met150/";
    //gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_met150/";
    //gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_met150/";
  } else if (sample_name == "ttbar") {
    mc_skim_folder = "mc/merged_higmc_higlep1T/";
    data_skim_folder = "data/merged_higdata_higlep1T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/";
    gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep1T/";
    gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_higlep1T/";
  } else if (sample_name == "zll") {
    mc_skim_folder = "mc/merged_higmc_higlep2T/";
    data_skim_folder = "data/merged_higdata_higlep2T/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/";
    gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep2T/";
    gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_higlep2T/";
  } else if (sample_name == "qcd") {
    mc_skim_folder = "mc/merged_higmc_higqcd/";
    data_skim_folder = "data/merged_higdata_higqcd/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/";
    gluino_fast_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higqcd/";
    gluino_full_skim_folder = "SMS-T5qqqqZH_FullSim/merged_higmc_higqcd/";
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
  string mass_points_string = "175_0, 500_0, 950_0"; // Original
  //string mass_points_string = "225_0, 700_0, 950_0"; // Below looks good 
  //string mass_points_string = "225_0, 400_0"; // 2016 AN
  //string mass_points_string = "600_0, 800_0, 1000_0"; // Boosted
  //if (unblind && !unblind_signalregion) mass_points_string = "";
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange, kPink, kCyan}; // Requires mass_points.size() >= sig_colors.size()
  // signal_filenames["TChiHH(nlsp,lsp)"] = {filename}
  torch::OrderedDict<string, set<string> > signal_filenames;
  // mass_points = [[nlsp, lsp]]
  vector<pair<int, int> > mass_points = parseMassPoints(mass_points_string);
  for (auto const & mass_point : mass_points) {
    signal_filenames.insert("TChiHH("+to_string(mass_point.first)+","+to_string(mass_point.second)+")", {"*TChiHH_mChi-"+to_string(mass_point.first)+"_mLSP-"+to_string(mass_point.second)+"_*.root"});
  }

  //string gluino_points_string = "1600_1, 2000_1, 2200_1"; // Boosted
  //string gluino_points_string = "1000_1, 1800_1, 2500_1"; // Boosted
  string gluino_points_string = "1000_1, 1100_1, 1200_1"; // Boosted
  // gluino_filenames["Gluino(nlsp,lsp)"] = {filename}
  torch::OrderedDict<string, set<string> > gluino_fast_filenames;
  torch::OrderedDict<string, set<string> > gluino_full_filenames;
  vector<pair<int, int> > gluino_points = parseMassPoints(gluino_points_string);
  for (auto const & gluino_point : gluino_points) {
    gluino_fast_filenames.insert("Gluino fast("+to_string(gluino_point.first)+","+to_string(gluino_point.second)+")", {"*T5qqqqZH*mGluino-"+to_string(gluino_point.first)+"*mChi-"+to_string(gluino_point.first-50)+"*mLSP-"+to_string(gluino_point.second)+"_*.root"});
    gluino_full_filenames.insert("Gluino full("+to_string(gluino_point.first)+","+to_string(gluino_point.second)+")", {"*T5qqqqZH*mGluino-"+to_string(gluino_point.first)+"*mChi-"+to_string(gluino_point.first-50)+"*mLSP-"+to_string(gluino_point.second)+"_*.root"});
  }

  // Set mc background+signal procs
  Palette colors("txt/colors.txt", "default");
  procsDict["mc"];
  procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  //procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  //procsDict["mc"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
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
  procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  //procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  //procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
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
  vector<int> gluino_colors = {kOrange, kOrange+5, kViolet, kViolet+5, kAzure, kAzure+5, kGreen}; // Requires mass_points.size() >= sig_colors.size()
  for (unsigned iMassPoint = 0; iMassPoint < gluino_fast_filenames.size(); ++iMassPoint) {
    //procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>(gluino_fast_filenames[iMassPoint].key(), Process::Type::signal, 
    //  sig_colors[iMassPoint], attach_folder(gluino_production_folder, years, gluino_fast_skim_folder, gluino_fast_filenames[iMassPoint].value() ), "1"));
    procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key(), Process::Type::signal, 
      gluino_colors[2*iMassPoint+1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), "1"));
    procsDict["mc_and_sig"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key()+" bbbb", Process::Type::signal, 
      gluino_colors[2*iMassPoint], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), only_hbb));
  }


  procsDict["signal_shapes"];
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["signal_shapes"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key(), Process::Type::background, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "1"));
  }
  for (unsigned iMassPoint = 0; iMassPoint < 1; ++iMassPoint) {
    procsDict["signal_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key(), Process::Type::background, 
      gluino_colors[2*iMassPoint+1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), "1"));
  }

  procsDict["gluino"];
  for (unsigned iMassPoint = 0; iMassPoint < gluino_fast_filenames.size(); ++iMassPoint) {
    //procsDict["sig"].push_back(Process::MakeShared<Baby_pico>(gluino_fast_filenames[iMassPoint].key(), Process::Type::signal, 
    //  gluino_colors[2*iMassPoint], attach_folder(gluino_production_folder, years, gluino_fast_skim_folder, gluino_fast_filenames[iMassPoint].value() ), "1"));
    procsDict["gluino"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key()+" bbbb", Process::Type::signal, 
      gluino_colors[2*iMassPoint], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), only_hbb));
    procsDict["gluino"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key(), Process::Type::signal, 
      gluino_colors[2*iMassPoint+1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), "1"));
  }

  procsDict["gluino_shapes"];
  for (unsigned iMassPoint = 0; iMassPoint < gluino_fast_filenames.size(); ++iMassPoint) {
    //procsDict["sig_as_mc"].push_back(Process::MakeShared<Baby_pico>(gluino_fast_filenames[iMassPoint].key(), Process::Type::background, 
    //  gluino_colors[2*iMassPoint], attach_folder(gluino_production_folder, years, gluino_fast_skim_folder, gluino_fast_filenames[iMassPoint].value() ), "1"));
    procsDict["gluino_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key()+" bbbb", Process::Type::background, 
      gluino_colors[2*iMassPoint], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), only_hbb));
    procsDict["gluino_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[iMassPoint].key(), Process::Type::background, 
      gluino_colors[2*iMassPoint+1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[iMassPoint].value() ), "1"));
  }

  int gluino_index = 2;

  procsDict["gluino_HvsHbb_shapes"];
  procsDict["gluino_HvsHbb_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbbb", Process::Type::background, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), only_hbb));
  procsDict["gluino_HvsHbb_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key(), Process::Type::background, 
      gluino_colors[1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), "1"));

  procsDict["gluino_HH_all"];
  procsDict["gluino_HH_all"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key(), Process::Type::signal, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), "1"));

  procsDict["gluino_HH_bbbb"];
  procsDict["gluino_HH_bbbb"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbbb", Process::Type::signal, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==0.));

  procsDict["gluino_HH"];
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbbb", Process::Type::signal, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==0.));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbww", Process::Type::signal, 
      gluino_colors[1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==1));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbgluongluon", Process::Type::signal, 
      gluino_colors[2], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==2));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbtautatu", Process::Type::signal, 
      gluino_colors[3], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==3));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbcc", Process::Type::signal, 
      gluino_colors[4], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==4));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbzz", Process::Type::signal, 
      gluino_colors[5], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==5));
  procsDict["gluino_HH"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" other", Process::Type::signal, 
      gluino_colors[6], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==6));

  procsDict["gluino_HH_shapes"];
  procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbbb", Process::Type::background, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==0.));
  procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbww", Process::Type::background, 
      gluino_colors[1], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==1));
  procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbgluongluon", Process::Type::background, 
      gluino_colors[2], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==2));
  //procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbtautatu", Process::Type::background, 
  //    gluino_colors[3], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==3));
  procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbcc", Process::Type::background, 
      gluino_colors[4], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==4));
  //procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" bbzz", Process::Type::background, 
  //    gluino_colors[5], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==5));
  //procsDict["gluino_HH_shapes"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[0].key()+" other", Process::Type::background, 
  //    gluino_colors[6], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[0].value() ), hh_decay==6));
  
  procsDict["gluino_HH_bbvsall"];
  procsDict["gluino_HH_bbvsall"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" bbbb", Process::Type::background, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), hh_decay==0.));
  procsDict["gluino_HH_bbvsall"].push_back(Process::MakeShared<Baby_pico>(gluino_full_filenames[gluino_index].key()+" all", Process::Type::data, 
      gluino_colors[0], attach_folder(gluino_production_folder, years, gluino_full_skim_folder, gluino_full_filenames[gluino_index].value() ), "1"));

  

  procsDict["mc_and_sig_and_data"];
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["zjets"]),"stitch"));
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_production_folder, years, mc_skim_folder,mc_filenames["wjets"]),"stitch"));
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["single_t"]),"stitch"));
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["qcd"]),"stitch")); 
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["other"]),"stitch"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key(), Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "1"));
  }
  procsDict["mc_and_sig_and_data"].push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, 
                   attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}), triggers_data));

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
  procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch"));
  //procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh>0"));
  //procsDict["mc_and_data"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["tt"]),"stitch&&ntrutauh==0"));
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

  procsDict["mc_3btag"];
  procsDict["mc_3btag"].push_back(Process::MakeShared<Baby_pico>("Norm. MC 2b", Process::Type::background,kAzure+10,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt==2&&nbm==2)"));
  procsDict["mc_3btag"].back()->SetFillColor(kWhite);
  procsDict["mc_3btag"].back()->SetLineColor(kAzure+10);
  procsDict["mc_3btag"].back()->SetLineWidth(2);
  procsDict["mc_3btag"].push_back(Process::MakeShared<Baby_pico>("MC 3b", Process::Type::data,kBlack,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["mc_3btag"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key()+" 3b", Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "(nbt>=2&&nbm==3&&nbl==3)"));
    procsDict["mc_3btag"].back()->SetLineWidth(2);
  }

  procsDict["mc_4btag"];
  procsDict["mc_4btag"].push_back(Process::MakeShared<Baby_pico>("Norm. MC 2b", Process::Type::background,kAzure+10,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt==2&&nbm==2)"));
  procsDict["mc_4btag"].back()->SetFillColor(kWhite);
  procsDict["mc_4btag"].back()->SetLineColor(kAzure+10);
  procsDict["mc_4btag"].back()->SetLineWidth(2);
  procsDict["mc_4btag"].push_back(Process::MakeShared<Baby_pico>("MC 4b", Process::Type::data,kBlack,
                  attach_folder(mc_production_folder, years, mc_skim_folder, mc_filenames["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["mc_4btag"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key()+" 4b", Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "(nbt>=2&&nbm>=3&&nbl>=4)"));
    procsDict["mc_4btag"].back()->SetLineWidth(2);
  }

  // Set data btag procs
  procsDict["data_3btag"];
  procsDict["data_3btag"].push_back(Process::MakeShared<Baby_pico>("Norm. Data 2b", Process::Type::background,kAzure+10,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt==2&&nbm==2)"));
  procsDict["data_3btag"].back()->SetFillColor(kWhite);
  procsDict["data_3btag"].back()->SetLineColor(kAzure+10);
  procsDict["data_3btag"].back()->SetLineWidth(2);
  procsDict["data_3btag"].push_back(Process::MakeShared<Baby_pico>("Data 3b", Process::Type::data,kBlack,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt>=2&&nbm==3&&nbl==3)"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["data_3btag"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key()+" 3b", Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "(nbt>=2&&nbm==3&&nbl==3)"));
    procsDict["data_3btag"].back()->SetLineWidth(2);
  }

  procsDict["data_4btag"];
  procsDict["data_4btag"].push_back(Process::MakeShared<Baby_pico>("Norm. Data 2b", Process::Type::background,kAzure+10,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt==2&&nbm==2)"));
  procsDict["data_4btag"].back()->SetFillColor(kWhite);
  procsDict["data_4btag"].back()->SetLineColor(kAzure+10);
  procsDict["data_4btag"].back()->SetLineWidth(2);
  procsDict["data_4btag"].push_back(Process::MakeShared<Baby_pico>("Data 4b", Process::Type::data,kBlack,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt>=2&&nbm>=3&&nbl>=4)"));
  for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
    procsDict["data_4btag"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key()+" 4b", Process::Type::signal, 
      sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "(nbt>=2&&nbm>=3&&nbl>=4)"));
    procsDict["data_4btag"].back()->SetLineWidth(2);
  }

  procsDict["data_3b4btag"];
  procsDict["data_3b4btag"].push_back(Process::MakeShared<Baby_pico>("Norm. Data 2b", Process::Type::background,kAzure+10,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt==2&&nbm==2)"));
  procsDict["data_3b4btag"].back()->SetFillColor(kWhite);
  procsDict["data_3b4btag"].back()->SetLineColor(kAzure+10);
  procsDict["data_3b4btag"].back()->SetLineWidth(2);
  procsDict["data_3b4btag"].push_back(Process::MakeShared<Baby_pico>("Data 3b|4b", Process::Type::data,kBlack,
                  attach_folder(data_production_folder, years, data_skim_folder, {"*.root"}),triggers_data&&"(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)"));
  //for (unsigned iMassPoint = 0; iMassPoint < signal_filenames.size(); ++iMassPoint) {
  //  procsDict["data_3b4btag"].push_back(Process::MakeShared<Baby_pico>(signal_filenames[iMassPoint].key()+" 3b|4b", Process::Type::signal, 
  //    sig_colors[iMassPoint], attach_folder(signal_production_folder, years, signal_skim_folder, signal_filenames[iMassPoint].value() ), "(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)"));
  //  procsDict["data_3b4btag"].back()->SetLineWidth(2);
  //}

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
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(false);
  PlotOpt lin_lumi = lin_norm_info.Title(TitleType::data).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::data_norm).RatioMaximum(1.86).PrintVals(false);

  PlotOpt lin_test = lin_norm_info().Title(TitleType::info).YAxis(YAxisType::linear).Stack(StackType::lumi_shapes).Bottom(BottomType::ratio).RatioMaximum(3).RatioMinimum(0.9);
  vector<PlotOpt> plt_test = {lin_test};

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_lin_lumi = {lin_lumi};
  vector<PlotOpt> plt_lin_data = {lin_norm_data};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  // Set folders according
  // production, nanoAODFolder, sample_name, year_string, 
  // Set procs
  // Set baseline, filter according: sample_name
  //string production = "higgsino_inyo"; 
  //string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7";
  string production = "higgsino_klamath"; 
  string nanoAodFolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7";

  // [mc, mcTtbar, mc_and_sig, mc_and_sig_and_data, mc_btag, mc_nisr, mcTtbar_met, mcTtbar_lowmet, mc_and_data, data_3btag]
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
  bin_cuts_search.push_back("met>150&&met<=200 && hig_cand_drmax[0]>1.1");

  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>200&&met<=300 && hig_cand_drmax[0]>1.1");

  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]<=1.1");
  bin_cuts_search.push_back("met>300&&met<=400 && hig_cand_drmax[0]>1.1");

  bin_cuts_search.push_back("met>400           && hig_cand_drmax[0]<=1.1");
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

  // boosted table
  bool table_boosted = false;
  if (table_boosted) {
    vector<pair<string, NamedFunc> > cutflow_cuts;
    cutflow_cuts.push_back({"$\\rm{met>150,ntk=0,nlep=0,! low dphi}$", "met>150&&ntk==0&&nlep==0&&!low_dphi_met"});
    cutflow_cuts.push_back({"filters",search_filters});
    cutflow_cuts.push_back({"$\\rm{met>300,ht>600}$","met>300&&ht>600"});
    cutflow_cuts.push_back({"$\\rm{nfjet>=2}$","nfjet>=2"});
    cutflow_cuts.push_back({"$\\rm{fjet01 pt>300}$","fjet_pt[0]>300&&fjet_pt[1]>300"});
    cutflow_cuts.push_back({"$\\rm{fjet01 m = [60,260]}$","fjet_msoftdrop[0]>=60&&fjet_msoftdrop[0]<=260&&fjet_msoftdrop[1]>=60&&fjet_msoftdrop[1]<=260"});

    vector<TableRow> tablerows;
    NamedFunc cutflow_cut = "1";
    for (unsigned iCut=0; iCut<cutflow_cuts.size(); ++iCut) {
      cutflow_cut = cutflow_cut && cutflow_cuts[iCut].second;
      tablerows.push_back(TableRow(cutflow_cuts[iCut].first, cutflow_cut, 0,0, weight));
    }
    pm.Push<Table>("FixName:table_boosted_gluino", tablerows, procs_search["mc_and_sig"], 0, 1, 0, 1, 0, /*print_uncertainty*/ 1).LuminosityTag(total_luminosity_string);
  }

  // resolved table
  bool table_resolved = false;
  if (table_resolved) {
    vector<TableRow> tablerows;
    NamedFunc baseline = search_resolved_cuts&&search_filters;
    NamedFunc hig = "hig_cand_am[0]>100&&hig_cand_am[0]<=140";
    NamedFunc sdb = "hig_cand_am[0]<=100||hig_cand_am[0]>140";
    NamedFunc b2 = "nbt==2&&nbm==2";
    NamedFunc b3 = "nbt>=2&&nbm==3&&nbl==3";
    NamedFunc b4 = "nbt>=2&&nbm>=3&&nbl>=4";
    for (unsigned iPlane = 0; iPlane < bin_cuts_search.size(); ++iPlane) {
      tablerows.push_back(TableRow("$"+CodeToLatex(bin_cuts_search[iPlane].Name())+"$"));
      tablerows.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&bin_cuts_search[iPlane], 1,0, weight));
      tablerows.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&bin_cuts_search[iPlane], 0,0, weight));
      tablerows.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&bin_cuts_search[iPlane], 0,0, weight));
      tablerows.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&bin_cuts_search[iPlane], 0,0, weight));
      tablerows.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&bin_cuts_search[iPlane], 0,0, weight));
      tablerows.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&bin_cuts_search[iPlane], 0,1, weight));
    }
    pm.Push<Table>("FixName:table_resolved_gluino", tablerows, procs_search["mc_and_sig"], 0, 1, 0, 1, 0, /*print_uncertainty*/ 1).LuminosityTag(total_luminosity_string);

    vector<TableRow> tablerows_combined;
    tablerows_combined.push_back(TableRow("All"));
    tablerows_combined.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2, 1,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 2b", baseline&&hig&&b2, 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3, 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 3b", baseline&&hig&&b3, 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4, 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 4b", baseline&&hig&&b4, 0,1, weight));
    tablerows_combined.push_back(TableRow("low drmax"));
    tablerows_combined.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&"hig_cand_drmax[0]<=1.1", 1,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&"hig_cand_drmax[0]<=1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&"hig_cand_drmax[0]<=1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&"hig_cand_drmax[0]<=1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&"hig_cand_drmax[0]<=1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&"hig_cand_drmax[0]<=1.1", 0,1, weight));
    tablerows_combined.push_back(TableRow("high drmax"));
    tablerows_combined.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&"hig_cand_drmax[0]>1.1", 1,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&"hig_cand_drmax[0]>1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&"hig_cand_drmax[0]>1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&"hig_cand_drmax[0]>1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&"hig_cand_drmax[0]>1.1", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&"hig_cand_drmax[0]>1.1", 0,1, weight));
    tablerows_combined.push_back(TableRow("met$<=$300"));
    tablerows_combined.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&"met<=300", 1,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&"met<=300", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&"met<=300", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&"met<=300", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&"met<=300", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&"met<=300", 0,1, weight));
    tablerows_combined.push_back(TableRow("met$>$300"));
    tablerows_combined.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&"met>300", 1,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&"met>300", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&"met>300", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&"met>300", 0,0, weight));
    tablerows_combined.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&"met>300", 0,0, weight));
    tablerows_combined.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&"met>300", 0,1, weight));
    pm.Push<Table>("FixName:table_resolved_gluino_combined", tablerows_combined, procs_search["mc_and_sig"], 0, 1, 0, 1, 0, /*print_uncertainty*/ 1).LuminosityTag(total_luminosity_string);
  }

  // plot resolved
  bool plot_resolved = false;
  if (plot_resolved) {
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:amjj__fastfull").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
      search_filters&&search_resolved_cuts, procs_search["gluino"], plt_log).Weight(weight)
      .Tag("FixName:nb__fastfull").LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:amjj_2b__fastfull").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:amjj_3b__fastfull").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:amjj_4b__fastfull").LuminosityTag(total_luminosity_string);

    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:h1_pt__fastfull").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino"], plt_lin).Weight(weight)
      .Tag("FixName:h2_pt__fastfull").LuminosityTag(total_luminosity_string);

    // plot according to decay channel
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:amjj__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:amjj__2b__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:amjj__3b__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:amjj__4b__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:nb__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbt", "N_{b}^{t}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_log).Weight(weight)
      .Tag("FixName:nbt__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbm", "N_{b}^{m}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_log).Weight(weight)
      .Tag("FixName:nbm__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbl", "N_{b}^{l}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_log).Weight(weight)
      .Tag("FixName:nbl__gluino_HH").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 1, "jet_deepcsv", "deepcsv", {}),
        search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_lin).Weight(weight)
        .Tag("FixName:deepcsv__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_log).Weight(weight)
      .Tag("FixName:met__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_log).Weight(weight)
      .Tag("FixName:drmax__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:h1_pt__gluino_HH").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH"], plt_lin).Weight(weight)
      .Tag("FixName:h2_pt__gluino_HH").LuminosityTag(total_luminosity_string);

    // plot according to decay channel with normalization
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:amjj__gluino_HH_shapes").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:nb__gluino_HH_shapes").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(15, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:met__gluino_HH_shapes").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:drmax__gluino_HH_shapes").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:h1_pt__gluino_HH_shapes").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts, procs_search["gluino_HH_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:h2_pt__gluino_HH_shapes").LuminosityTag(total_luminosity_string);

  }

  bool plot_by_hig_pt = false;
  if (plot_by_hig_pt) {
    // Divide by lead H pt
    vector<NamedFunc> higPts= {Higfuncs::hig_pt1>0. &&Higfuncs::hig_pt1<=450,
                         Higfuncs::hig_pt1>450 &&Higfuncs::hig_pt1<=550,
                         Higfuncs::hig_pt1>550};


    // Compare between gluino mass points
    for (unsigned iHigPt=0; iHigPt < higPts.size(); ++iHigPt) {
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&higPts[iHigPt], procs_search["gluino_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_2b_hig_pt"+to_string(iHigPt)+"__gluino_mass_points").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3&&higPts[iHigPt], procs_search["gluino_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_3b_hig_pt"+to_string(iHigPt)+"__gluino_mass_points").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4&&higPts[iHigPt], procs_search["gluino_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_4b_hig_pt"+to_string(iHigPt)+"__gluino_mass_points").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
        search_filters&&search_resolved_cuts&&higPts[iHigPt], procs_search["gluino_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:nb_hig_pt"+to_string(iHigPt)+"__gluino_mass_points").LuminosityTag(total_luminosity_string);
    }

    // Compare between H vs bbbb
    for (unsigned iHigPt=0; iHigPt < higPts.size(); ++iHigPt) {
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&higPts[iHigPt], procs_search["gluino_HvsHbb_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_2b_hig_pt"+to_string(iHigPt)+"__HvsHbb").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3&&higPts[iHigPt], procs_search["gluino_HvsHbb_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_3b_hig_pt"+to_string(iHigPt)+"__HvsHbb").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4&&higPts[iHigPt], procs_search["gluino_HvsHbb_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:amjj_4b_hig_pt"+to_string(iHigPt)+"__HvsHbb").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 1, "jet_deepcsv", "deepcsv", {}),
        search_filters&&search_resolved_cuts&&higPts[iHigPt], procs_search["gluino_HvsHbb_shapes"], plt_shapes).Weight(weight)
        .Tag("FixName:deepcsv_hig_pt"+to_string(iHigPt)+"__HvsHbb_shapes").LuminosityTag(total_luminosity_string);
    }

    //higPts= {Higfuncs::hig_pt1>0. &&Higfuncs::hig_pt1<=200,
    //                     Higfuncs::hig_pt1>200 &&Higfuncs::hig_pt1<=300,
    //                     Higfuncs::hig_pt1>300 &&Higfuncs::hig_pt1<=400,
    //                     Higfuncs::hig_pt1>400 &&Higfuncs::hig_pt1<=500,
    //                     Higfuncs::hig_pt1>500 &&Higfuncs::hig_pt1<=600,
    //                     Higfuncs::hig_pt1>600 &&Higfuncs::hig_pt1<=700,
    //                     Higfuncs::hig_pt1>700};

    for (unsigned iHigPt=0; iHigPt < higPts.size(); ++iHigPt) {
      pm.Push<Hist1D>(Axis(10, 0, 1, "jet_deepcsv", "deepcsv", {}),
        search_filters&&search_resolved_cuts&&higPts[iHigPt], procs_search["gluino_HH"], plt_lin).Weight(weight)
        .Tag("FixName:deepcsv_hig_pt"+to_string(iHigPt)+"__HvsHbb").LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
        search_filters&&search_resolved_cuts&&higPts[iHigPt], procs_search["gluino_HH"], plt_lin).Weight(weight)
        .Tag("FixName:nb_hig_pt"+to_string(iHigPt)+"__HvsHbb").LuminosityTag(total_luminosity_string);
    }
  }

  bool plot_corrections = false;
  if (plot_corrections) {
    // Correction historgrams
    // 2b, 3b, 4b for lead H pT
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h1_pt_2b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h1_pt_3b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h1_pt_4b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");

    // 2b, 3b, 4b for sublead H pT
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h2_pt_2b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==3, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h2_pt_3b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt2, "sublead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==4, procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h2_pt_4b__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");

   // Large drmax
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&"hig_cand_drmax[0]>1.1", procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h1_pt_2b_highdrmax__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");
   // Small met
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&"met<200", procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
      .Tag("FixName:h1_pt_2b_lowmet__gluino_HH").LuminosityTag(total_luminosity_string).RatioTitle("All","bbbb");

    // Per signals
    pm.Push<Hist1D>(Axis(10, 0, 2000, Higfuncs::hig_pt1, "lead H p_{T}", {}),
      search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2, procs_search["signal_shapes"], plt_shapes).Weight(weight)
      .Tag("FixName:h1_pt_2b__signals").LuminosityTag(total_luminosity_string);
  }
  
  bool apply_correction = true;
  if (apply_correction) {
    // Check <m_bb> in planes for 2b
    // Before
    for (unsigned iPlane = 0; iPlane<bin_cuts_search.size(); ++iPlane) {
      pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&bin_cuts_search[iPlane], procs_search["gluino_HH_bbvsall"], plt_test).Weight(weight)
        .Tag("FixName:amjj_plane_"+to_string(iPlane)+"___before_correction__HH_bbbb").LuminosityTag(total_luminosity_string).RatioTitle("All", "bbbb");
    }
    // After
    for (unsigned iPlane = 0; iPlane<bin_cuts_search.size(); ++iPlane) {
      pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&bin_cuts_search[iPlane], procs_search["gluino_HH_all"], plt_lin).Weight(weight)
        .Tag("FixName:amjj_plane_"+to_string(iPlane)+"___after_correction__HH_all").LuminosityTag(total_luminosity_string);
    }
    for (unsigned iPlane = 0; iPlane<bin_cuts_search.size(); ++iPlane) {
      pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        search_filters&&search_resolved_cuts&&Higfuncs::hig_bcat==2&&bin_cuts_search[iPlane], procs_search["gluino_HH_bbbb"], plt_lin).Weight(weight*deepcsv_2016)
        .Tag("FixName:amjj_plane_"+to_string(iPlane)+"___after_correction__HH_bbbb").LuminosityTag(total_luminosity_string);
    }

    bool table_corrected = true;
    if (table_corrected) {
      vector<TableRow> tablerows;
      NamedFunc baseline = search_resolved_cuts&&search_filters;
      NamedFunc hig = "hig_cand_am[0]>100&&hig_cand_am[0]<=140";
      NamedFunc sdb = "hig_cand_am[0]<=100||hig_cand_am[0]>140";
      NamedFunc b2 = "nbt==2&&nbm==2";
      NamedFunc b3 = "nbt>=2&&nbm==3&&nbl==3";
      NamedFunc b4 = "nbt>=2&&nbm>=3&&nbl>=4";
      for (unsigned iPlane = 0; iPlane < bin_cuts_search.size(); ++iPlane) {
        tablerows.push_back(TableRow("$"+CodeToLatex(bin_cuts_search[iPlane].Name())+"$"));
        tablerows.push_back(TableRow("SBD, 2b", baseline&&sdb&&b2&&bin_cuts_search[iPlane], 1,0, weight*deepcsv_2016));
        tablerows.push_back(TableRow("HIG, 2b", baseline&&hig&&b2&&bin_cuts_search[iPlane], 0,0, weight*deepcsv_2016));
        tablerows.push_back(TableRow("SBD, 3b", baseline&&sdb&&b3&&bin_cuts_search[iPlane], 0,0, weight));
        tablerows.push_back(TableRow("HIG, 3b", baseline&&hig&&b3&&bin_cuts_search[iPlane], 0,0, weight));
        tablerows.push_back(TableRow("SBD, 4b", baseline&&sdb&&b4&&bin_cuts_search[iPlane], 0,0, weight));
        tablerows.push_back(TableRow("HIG, 4b", baseline&&hig&&b4&&bin_cuts_search[iPlane], 0,1, weight));
      }
      pm.Push<Table>("FixName:table_resolved_gluino_corrected", tablerows, procs_search["mc_and_sig"], 0, 1, 0, 1, 0, /*print_uncertainty*/ 1).LuminosityTag(total_luminosity_string);
    }
  }

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  // Get figure
  //for(auto & figure: pm.Figures()) {
  //  Hist1D * figure_hists = static_cast<Hist1D*>(figure.get());
  //  if (!(Contains(figure_hists->Name(), "nb_hig_pt") && Contains(figure_hists->Name(), "__HvsHbb"))) continue;
  //  cout<<"Figure: "<<figure_hists->tag_<<endl;
  //  for (auto & process : figure_hists->GetProcesses()) {
  //    cout<<"  Process: "<<process->name_<<endl;
  //    TH1D histogram = static_cast<Hist1D::SingleHist1D*>(figure_hists->GetComponent(process))->scaled_hist_;
  //    histogram.Print("all");
  //  }
  //}
  //for(auto & figure: pm.Figures()) {
  //  Hist1D * figure_hists = static_cast<Hist1D*>(figure.get());
  //  if (!(Contains(figure_hists->Name(), "deepcsv_hig_pt") && Contains(figure_hists->Name(), "__HvsHbb"))) continue;
  //  cout<<"Figure: "<<figure_hists->tag_<<endl;
  //  for (auto & process : figure_hists->GetProcesses()) {
  //    cout<<"  Process: "<<process->name_<<endl;
  //    TH1D histogram = static_cast<Hist1D::SingleHist1D*>(figure_hists->GetComponent(process))->scaled_hist_;
  //    histogram.Print("all");
  //  }
  //}
  // Collect histograms. Assume one process per histogram
  map<string, TH1D *> histograms;
  for(auto & figure: pm.Figures()) {
    Hist1D * figure_hists = static_cast<Hist1D*>(figure.get());
    if (!(Contains(figure_hists->Name(), "correction__HH"))) continue;
    cout<<"Figure: "<<figure_hists->tag_<<endl;
    for (auto & process : figure_hists->GetProcesses()) {
      cout<<"  Process: "<<process->name_<<endl;
      histograms[figure_hists->tag_] = &static_cast<Hist1D::SingleHist1D*>(figure_hists->GetComponent(process))->scaled_hist_;
    }
  }
  // Plot histograms together
  if (histograms.size() !=0) {
    for (unsigned iPlane = 0; iPlane < bin_cuts_search.size(); ++iPlane) {
      TCanvas c1("c1","c1", 500, 500);
      histograms["FixName:amjj_plane_"+to_string(iPlane)+"___after_correction__HH_bbbb"]->SetLineColor(kRed);
      histograms["FixName:amjj_plane_"+to_string(iPlane)+"___after_correction__HH_bbbb"]->Draw("hist");
      histograms["FixName:amjj_plane_"+to_string(iPlane)+"___after_correction__HH_all"]->Draw("same");
      c1.SaveAs(("plots/correction_plane_"+to_string(iPlane)+".pdf").c_str());
      cout<<"open plots/correction_plane_"<<iPlane<<".pdf"<<endl;
    }
  }


  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"unblind_signalregion", no_argument, 0, 0},
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
      } else if (optname == "unblind_signalregion") {
        unblind_signalregion = true;
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
