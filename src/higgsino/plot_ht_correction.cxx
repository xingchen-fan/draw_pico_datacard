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
#include "core/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "run2";
  bool unblind = false;
  bool unblind_signalregion = false;
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
void setProcsDict(string const & production, string const & nanoAodFolder, string const & local_year_string, string const & sample_name, map<string, vector<shared_ptr<Process> > > & procsDict) {
  set<int> years;
  //years_a = {2016};
  HigUtilities::parseYears(local_year_string, years);

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
    //mc_skim_folder = "mc/merged_higmc_higloose/";
    //data_skim_folder = "data/merged_higdata_higloose/";
    //signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/";
    mc_skim_folder = "mc/merged_higmc_preselect/";
    data_skim_folder = "data/merged_higdata_preselect/";
    signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";
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
  string mass_points_string = "175_0, 500_0, 950_0"; // Original
  //string mass_points_string = "225_0, 700_0, 950_0"; // Below looks good 
  //string mass_points_string = "225_0, 400_0"; // 2016 AN
  //if (unblind && !unblind_signalregion) mass_points_string = "";
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

  // Let's make a kappa corrected version

}

const NamedFunc weight_ht_2016("weight_ht_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.424;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.115;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.03;
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
  else if (b.ht()>800 && b.ht()<=850) weight = 0.7;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.86;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.705;
  else if (b.ht()>950) weight = 0.874;
  return weight;
});

const NamedFunc weight_ht_lowdrmax_2016("weight_ht_lowdrmax_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 2.159;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.216;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.07;
  else if (b.ht()>250 && b.ht()<=300) weight = 0.993;
  else if (b.ht()>300 && b.ht()<=350) weight = 1.054;
  else if (b.ht()>350 && b.ht()<=400) weight = 1.13;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.991;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.9;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.897;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.863;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.77;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.802;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.677;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.865;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.419;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.933;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.608;
  else if (b.ht()>950) weight = 0.811;
  return weight;
});

const NamedFunc weight_ht_highdrmax_2016("weight_ht_highdrmax_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.37;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.111;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.029;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.019;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.988;
  else if (b.ht()>350 && b.ht()<=400) weight = 1.026;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.886;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.919;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.913;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.869;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.815;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.999;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.882;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.913;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.9;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.803;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.76;
  else if (b.ht()>950) weight = 0.889;
  return weight;
});

const NamedFunc weight_ht_2017("weight_ht_2017", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.53;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.138;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.112;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.055;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.994;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.982;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.798;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.837;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.803;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.911;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.769;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.766;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.697;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.767;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.6;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.674;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.557;
  else if (b.ht()>950) weight = 0.782;
  return weight;
});

const NamedFunc weight_ht_lowdrmax_2017("weight_ht_lowdrmax_2017", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.484;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.254;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.181;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.064;
  else if (b.ht()>300 && b.ht()<=350) weight = 1.067;
  else if (b.ht()>350 && b.ht()<=400) weight = 1.093;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.937;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.786;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.805;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.97;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.703;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.644;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.803;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.597;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.498;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.608;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.434;
  else if (b.ht()>950) weight = 0.751;
  return weight;
});

const NamedFunc weight_ht_highdrmax_2017("weight_ht_highdrmax_2017", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.538;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.133;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.109;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.055;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.986;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.969;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.78;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.843;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.796;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.888;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.78;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.8;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.636;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.844;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.655;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.692;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.633;
  else if (b.ht()>950) weight = 0.766;
  return weight;
});


const NamedFunc weight_ht_2018("weight_ht_2018", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.658;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.256;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.133;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.063;
  else if (b.ht()>300 && b.ht()<=350) weight = 1.005;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.951;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.789;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.746;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.773;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.783;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.767;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.695;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.679;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.721;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.69;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.649;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.699;
  else if (b.ht()>950) weight = 0.612;
  return weight;
});
const NamedFunc weight_ht_lowdrmax_2018("weight_ht_lowdrmax_2018", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 2.681;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.386;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.18;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.11;
  else if (b.ht()>300 && b.ht()<=350) weight = 1.095;
  else if (b.ht()>350 && b.ht()<=400) weight = 1.081;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.918;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.769;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.748;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.684;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.833;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.811;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.65;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.701;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.776;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.611;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.867;
  else if (b.ht()>950) weight = 0.539;
  return weight;
});
const NamedFunc weight_ht_highdrmax_2018("weight_ht_highdrmax_2018", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.595;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.25;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.131;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.06;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.996;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.937;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.771;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.739;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.774;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.809;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.74;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.646;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.681;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.718;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.637;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.658;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.576;
  else if (b.ht()>950) weight = 0.663;
  return weight;
});

const NamedFunc weight_ht_sideband_2016("weight_ht_sideband_2016", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.492;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.111;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.044;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.019;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.976;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.982;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.889;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.836;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.859;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.723;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.871;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.802;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.747;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.791;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.852;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.802;
  else if (b.ht()>900 && b.ht()<=950) weight = 1.012;
  else if (b.ht()>950) weight = 0.971;
  return weight;
});

const NamedFunc weight_ht_sideband_2017("weight_ht_sideband_2017", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.5;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.164;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.098;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.051;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.977;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.897;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.775;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.779;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.741;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.697;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.717;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.657;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.661;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.718;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.603;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.8;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.426;
  else if (b.ht()>950) weight = 0.764;
  return weight;
});

const NamedFunc weight_ht_sideband_2018("weight_ht_sideband_2018", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;
  float weight = 1;
  if (b.ht()>0 && b.ht()<=50) weight = 0.0;
  else if (b.ht()>50 && b.ht()<=100) weight = 0.0;
  else if (b.ht()>100 && b.ht()<=150) weight = 1.629;
  else if (b.ht()>150 && b.ht()<=200) weight = 1.251;
  else if (b.ht()>200 && b.ht()<=250) weight = 1.136;
  else if (b.ht()>250 && b.ht()<=300) weight = 1.039;
  else if (b.ht()>300 && b.ht()<=350) weight = 0.956;
  else if (b.ht()>350 && b.ht()<=400) weight = 0.887;
  else if (b.ht()>400 && b.ht()<=450) weight = 0.738;
  else if (b.ht()>450 && b.ht()<=500) weight = 0.704;
  else if (b.ht()>500 && b.ht()<=550) weight = 0.687;
  else if (b.ht()>550 && b.ht()<=600) weight = 0.716;
  else if (b.ht()>600 && b.ht()<=650) weight = 0.686;
  else if (b.ht()>650 && b.ht()<=700) weight = 0.681;
  else if (b.ht()>700 && b.ht()<=750) weight = 0.667;
  else if (b.ht()>750 && b.ht()<=800) weight = 0.61;
  else if (b.ht()>800 && b.ht()<=850) weight = 0.82;
  else if (b.ht()>850 && b.ht()<=900) weight = 0.723;
  else if (b.ht()>900 && b.ht()<=950) weight = 0.47;
  else if (b.ht()>950) weight = 0.624;
  return weight;
});

const NamedFunc weight_ht("weight_ht", [](const Baby &b) -> NamedFunc::ScalarType {
  float weight = 1.;
  if (b.SampleType()==2016) weight = weight_ht_2016.GetScalar(b);
  else if (b.SampleType()==2017) weight = weight_ht_2017.GetScalar(b);
  else if (b.SampleType()==2018) weight = weight_ht_2018.GetScalar(b);
  return weight;
});

const NamedFunc weight_ht_lowdrmax("weight_ht_lowdrmax", [](const Baby &b) -> NamedFunc::ScalarType {
  float weight = 1.;
  if (b.SampleType()==2016) weight = weight_ht_lowdrmax_2016.GetScalar(b);
  else if (b.SampleType()==2017) weight = weight_ht_lowdrmax_2017.GetScalar(b);
  else if (b.SampleType()==2018) weight = weight_ht_lowdrmax_2018.GetScalar(b);
  return weight;
});

const NamedFunc weight_ht_highdrmax("weight_ht_highdrmax", [](const Baby &b) -> NamedFunc::ScalarType {
  float weight = 1.;
  if (b.SampleType()==2016) weight = weight_ht_highdrmax_2016.GetScalar(b);
  else if (b.SampleType()==2017) weight = weight_ht_highdrmax_2017.GetScalar(b);
  else if (b.SampleType()==2018) weight = weight_ht_highdrmax_2018.GetScalar(b);
  return weight;
});

const NamedFunc weight_ht_sideband("weight_ht_sideband", [](const Baby &b) -> NamedFunc::ScalarType {
  float weight = 1.;
  if (b.SampleType()==2016) weight = weight_ht_sideband_2016.GetScalar(b);
  else if (b.SampleType()==2017) weight = weight_ht_sideband_2017.GetScalar(b);
  else if (b.SampleType()==2018) weight = weight_ht_sideband_2018.GetScalar(b);
  return weight;
});

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
    //.LegendColumns(3);
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_print = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).PrintVals(true);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_norm_data_print = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_lumi = lin_norm_info.Title(TitleType::data).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::data_norm).RatioMaximum(1.86);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_print = {lin_norm_print};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_lin_lumi = {lin_lumi};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_lin_print = {lin_norm_data_print};
  if (unblind) plt_log = {log_norm_data};

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

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

  string proc_name = "mc_and_sig";
  if (unblind) proc_name += "_and_data";

  NamedFunc extra_cut = "1";
  NamedFunc signalReject = "njet>=4&&njet<=5&&!(((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&(hig_cand_am[0]>100&&hig_cand_am[0]<=140))";
  if (unblind && !unblind_signalregion) extra_cut = signalReject;

  // ht [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&ttbar_resolved_cuts,
    procs_ttbar[proc_name], plt_lin_print).Weight(weight).Tag("FixName:ht__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // ht; low-drmax [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1",
    procs_ttbar[proc_name], plt_lin_print).Weight(weight).Tag("FixName:ht__lowdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // ht; high-drmax [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1",
    procs_ttbar[proc_name], plt_lin_print).Weight(weight).Tag("FixName:ht__highdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // ht; sideband [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
    "(hig_cand_drmax[0]>2.2&&hig_cand_drmax[0]<=3)&&"
    "hig_cand_dm[0]<=40&&"
    "hig_cand_am[0]<=200&&"
    "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
    && Higfuncs::lead_signal_lepton_pt>30
    ,procs_ttbar[proc_name], plt_lin_print).Weight(weight).Tag("FixName:ht__sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am-ht; [ttbar]
  pm.Push<Hist2D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}), Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts, procs_ttbar["mc"], plt_2D).Weight(weight).Tag("FixName:ht_vs_am__ttbarmc__"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // drmax; no-drmax cut [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
    ttbar_filters&&
    "met/met_calo<5&&weight<1.5&&"
    "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
    "hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
    "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
    && Higfuncs::lead_signal_lepton_pt>30
    , procs_ttbar[proc_name], plt_lin).Weight(weight).Tag("FixName:drmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am; low-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__lowdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am; high-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__highdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

  // ht;ht_weight [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&ttbar_resolved_cuts,
    procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht).Tag("FixName:ht__htwgt__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight, low-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht).Tag("FixName:amjj__lowdrmax_htwgt__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight high-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht).Tag("FixName:amjj__highdrmax_htwgt__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(low-drmax), low-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_lowdrmax).Tag("FixName:amjj__lowdrmax_htwgt_lowdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(high-drmax), high-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_highdrmax).Tag("FixName:amjj__highdrmax_htwgt_lowdrmax__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

  // ht;ht_weight(sideband) [ttbar]
  pm.Push<Hist1D>(Axis(20, 0, 1000., "ht", "HT [GeV]", {}),
    ttbar_filters&&ttbar_resolved_cuts,
    procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:ht__htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(sideband), low-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__lowdrmax_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(sideband), high-drmax [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__highdrmax_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);


  // am;ht_weight(sideband), low-drmax,2b [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    "&&(nbt==2&&nbm==2)"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__lowdrmax_2b_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(sideband), low-drmax,3b|4b [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]<=1.1"
    "&&((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__lowdrmax_3b4b_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(sideband), high-drmax,2b [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    "&&(nbt==2&&nbm==2)"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__highdrmax_2b_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  // am;ht_weight(sideband), high-drmax,3b|4b [ttbar]
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
    ttbar_filters&&ttbar_resolved_cuts&&"hig_cand_drmax[0]>1.1"
    "&&((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
    , procs_ttbar[proc_name], plt_lin).Weight(weight*weight_ht_sideband).Tag("FixName:amjj__highdrmax_3b4b_htwgt_sideband__ttbar_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);



  //// ht [search]
  //// drmax; 2b, no-drmax cut [search]
  //pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
  //  search_filters&&
  //  "met/mht<2 && met/met_calo<2&&weight<1.5&&"
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
  //  "(nbt==2&&nbm==2)"&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:drmax__2b__search").LuminosityTag(total_luminosity_string);
  //// drmax; 3b, no-drmax cut [search]
  //pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
  //  search_filters&&
  //  "met/mht<2 && met/met_calo<2&&weight<1.5&&"
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<=200&&hig_cand_dm[0]<=40&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:drmax__3b4b__search").LuminosityTag(total_luminosity_string);

  //// average mass; 2b, low-drmax[search]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  search_filters&&search_resolved_cuts&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__2b_lowdrmax__search").LuminosityTag(total_luminosity_string);
  //// average mass; 3b|4b, low-drmax, no-average mass cut [search]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  search_filters&&search_resolved_cuts&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__3b4b_lowdrmax__search").LuminosityTag(total_luminosity_string);
  //// average mass; 2b, high-drmax, no-average mass cut [search]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  search_filters&&search_resolved_cuts&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__2b_highdrmax__search").LuminosityTag(total_luminosity_string);
  //// average mass; 3b|4b, high-drmax, no-average mass cut [search]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  search_filters&&search_resolved_cuts&&extra_cut
  //  , procs_search[proc_name], plt_lin).Weight(weight).Tag("FixName:amjj__3b4b_highdrmax__search").LuminosityTag(total_luminosity_string);


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
