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
#include "TLorentzVector.h"

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

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  //if (b.type()==106000) {
  //  return 35.9;
  //}
  if (b.SampleType()==2016){
    return weight*35.9;
  } else if (b.SampleType()==2017){
    return weight*41.5;
  } else {
    return weight*59.6;
  }
});

namespace{
  bool single_thread = false;
  string year_string = "2016";
  bool unblind = false;
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

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  //if (unblind) plt_lin = {lin_norm_data};
  //if (unblind) plt_log = {log_norm_data};

  string higgsino_version = "";

  // Set options
  //string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt"+higgsino_version+"/";
  //string mc_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado";
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";

  string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo"+higgsino_version+"/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  //string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  string ttbar_data_skim_folder = "data/skim_higlep1T/";

  //string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo"+higgsino_version+"/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";

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

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years*Functions::w_pileup;
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2_old*w_years;
  //if (years.size()==1 && *years.begin()==2016) weight *= "137.";
  //else weight *= w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;
  NamedFunc triggers_data = "1";
  NamedFunc lepton_triggers = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50 || HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";
  NamedFunc met_triggers = "(HLT_PFMET110_PFMHT110_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)";
  triggers_data = lepton_triggers || met_triggers;

  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/met_calo<5 && weight<10"; //since pass_fsjets is not quite usable...
  NamedFunc base_filters = Functions::hem_veto && "pass && met/met_calo<5 && weight<1.5";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  // resolved cuts
  //NamedFunc search_resolved = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved = 
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

  vector<shared_ptr<Process> > ttbar_procs;
  vector<shared_ptr<Process> > ttbar_procs_ratio;
  // Set mc processes
  ////ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  ////                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  ////ttbar_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  ////                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["vjets"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));

  if (unblind) {
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data-JetID event veto", Process::Type::data, kRed,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data&&"pass_jets"));
    ttbar_procs_ratio.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::background, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data));
    ttbar_procs_ratio.push_back(Process::MakeShared<Baby_pico>("Data-JetID event veto", Process::Type::data, kRed,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data&&"pass_jets"));
  }

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

  //----------------------------------------------------------------------------------
  //                                  namedfuncs
  //----------------------------------------------------------------------------------

  const NamedFunc pass_boostedhemveto("pass_boostedhemveto", [](const Baby &b) -> NamedFunc::ScalarType{
    //only apply for relevant runs
    bool pass_hem = true;
    if (year_string != "2018") return 1.0;
    if (b.type()/1000 == 0 && b.run() < 319077) return 1.0;
    for (unsigned int el_idx = 0; el_idx < b.el_pt()->size(); el_idx++) {
      if (b.el_miniso()->at(el_idx) < 0.1 && -3.0 < b.el_eta()->at(el_idx) && b.el_eta()->at(el_idx) < -1.4 && -1.57 < b.el_phi()->at(el_idx) && b.el_phi()->at(el_idx) < -0.87) {
        pass_hem = false;
      }
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_pt()->at(jet_idx) > 30. && -3.2 < b.jet_eta()->at(jet_idx) && b.jet_eta()->at(jet_idx) < -1.2 && -1.77 < b.jet_phi()->at(jet_idx) && b.jet_phi()->at(jet_idx) < -0.67) {
        double dphi = fabs(TVector2::Phi_mpi_pi(b.jet_phi()->at(jet_idx)-b.met_phi()));
	if (dphi < 0.5) {
	  pass_hem = false;
	}
      }
    }
    if (b.type()/1000 == 0) {
      //data
      if (pass_hem) return 1.0;
      return 0.0;
    }
    //mc
    if (pass_hem) return 1.0;
    //weight for amount of 2018 luminosity after HEM issue occured
    return 0.645;
  });

  const NamedFunc jetidmod_njet("jetidmod_njet", [](const Baby &b) -> NamedFunc::ScalarType{
    unsigned int njet = 0;
    if (b.pass_jets()) return b.njet();
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) continue;
      njet += 1;
    }
    return njet;
  });

  const NamedFunc nb("nb", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbt() < 2) return b.nbt();
    if (b.nbm() < 3) return 2;
    if (b.nbl()==3) return 3;
    return 4;
  });

  const NamedFunc jetidmod_nb("jetidmod_nb", [](const Baby &b) -> NamedFunc::ScalarType{
    unsigned int nbt = 0, nbm = 0, nbl = 0;
    if (b.pass_jets()) {
      if (b.nbt() < 2) return b.nbt();
      if (b.nbm() < 3) return 2;
      if (b.nbl()==3) return 3;
      return 4;
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) continue;
      if (b.jet_deepcsv()->at(jet_idx) > 0.8953 && year_string == "2016") nbt += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.8001 && year_string == "2017") nbt += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.7527 && year_string == "2018") nbt += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.6321 && year_string == "2016") nbm += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.4941 && year_string == "2017") nbm += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.4184 && year_string == "2018") nbm += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.2217 && year_string == "2016") nbl += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.1522 && year_string == "2017") nbl += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.0494 && year_string == "2018") nbl += 1;
    }
    if (nbt < 2) return nbt;
    if (nbm < 3) return 2;
    if (nbl==3) return 3;
    return 4;
  });

  const NamedFunc jetidmod_low_dphi_met("jetidmod_low_dphi_met", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.pass_jets()) return b.low_dphi_met();
    std::vector<std::pair<float, float>> jet_pt_and_dphi;
    unsigned int njet = 0;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) continue;
      njet += 1;
      jet_pt_and_dphi.push_back(std::make_pair(b.jet_pt()->at(jet_idx), b.jet_met_dphi()->at(jet_idx)));
    }
    std::sort(jet_pt_and_dphi.begin(), jet_pt_and_dphi.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.first > bb.first;
      });
    if (njet > 0) {
      if (jet_pt_and_dphi[0].second < 0.5) return true;
    }
    if (njet > 1) {
      if (jet_pt_and_dphi[1].second < 0.5) return true;
    }
    if (njet > 2) {
      if (jet_pt_and_dphi[2].second < 0.3) return true;
    }
    if (njet > 3) {
      if (jet_pt_and_dphi[3].second < 0.3) return true;
    }
    return false;
  });

  const NamedFunc jetidmod_hig_cand_am("jetidmod_hig_cand_am", [](const Baby &b) -> NamedFunc::ScalarType{
    //if no bad jets, just return normal cuts
    if (b.pass_jets()) {
      if (b.njet() < 4) return -999.;
      return b.hig_cand_am()->at(0);
    }
    //otherwise, we have to recalculate higgs candidates
    unsigned int njet = 0;
    std::vector<std::pair<unsigned int, float>> jet_idx_and_btagdisc;
    bool recalculate_higgs = false;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) {
        if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
        continue;
      }
      njet += 1;
      jet_idx_and_btagdisc.push_back(std::make_pair(jet_idx, b.jet_deepcsv()->at(jet_idx)));
    }
    //njet, nbt, dphi cuts, also if bad jet was the non-higgs cand jet, just return original cuts
    if (njet < 4) return -999.;
    if (b.njet()<4) std::cout << "ERROR: njet<4\n";
    if (b.hig_cand_am()->size()==0) std::cout << "ERROR: njet>=4 but no hig cand\n";
    if (!recalculate_higgs) return b.hig_cand_am()->at(0);
    //get the lowest Delta m Higgs candidate pairing with the first two indices corresponding to the higher pt higgs candidate
    std::sort(jet_idx_and_btagdisc.begin(), jet_idx_and_btagdisc.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.second > bb.second;
      });
    std::vector<TLorentzVector> b_jet_p;
    for (unsigned int bjet_idx = 0; bjet_idx < 4; bjet_idx++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_eta()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_phi()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_m()->at(jet_idx_and_btagdisc[bjet_idx].first));
      b_jet_p.push_back(lv);
    }
    std::vector<std::vector<unsigned int>> higgs_combinations = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    float min_delta_m = 999;
    float min_dm_am = 0;
    bool first = true;
    for (unsigned int higgs_idx = 0; higgs_idx < 3; higgs_idx++) {
       TLorentzVector higgs1_p = b_jet_p[higgs_combinations[higgs_idx][0]]+b_jet_p[higgs_combinations[higgs_idx][1]];
       TLorentzVector higgs2_p = b_jet_p[higgs_combinations[higgs_idx][2]]+b_jet_p[higgs_combinations[higgs_idx][3]];
       float delta_m = TMath::Abs(higgs1_p.M()-higgs2_p.M());
       if ((delta_m < min_delta_m) || first) {
         first = false;
         min_delta_m = delta_m;
	 min_dm_am = (higgs1_p.M() + higgs2_p.M())/2.0;
       }
    }
    return min_dm_am;
  });

  const NamedFunc jetidmod_hig_cand_dm("jetidmod_hig_cand_dm", [](const Baby &b) -> NamedFunc::ScalarType{
    //if no bad jets, just return normal cuts
    if (b.pass_jets()) {
      if (b.njet() < 4) return 999.;
      return b.hig_cand_dm()->at(0);
    }
    //otherwise, we have to recalculate higgs candidates
    unsigned int njet = 0;
    std::vector<std::pair<unsigned int, float>> jet_idx_and_btagdisc;
    bool recalculate_higgs = false;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) {
        if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
        continue;
      }
      njet += 1;
      jet_idx_and_btagdisc.push_back(std::make_pair(jet_idx, b.jet_deepcsv()->at(jet_idx)));
    }
    //njet, nbt, dphi cuts, also if bad jet was the non-higgs cand jet, just return original cuts
    if (njet < 4) return 999.;
    if (b.njet()<4) std::cout << "ERROR: njet<4\n";
    if (b.hig_cand_am()->size()==0) std::cout << "ERROR: njet>=4 but no hig cand\n";
    if (!recalculate_higgs) return b.hig_cand_dm()->at(0);
    //get the lowest Delta m Higgs candidate pairing with the first two indices corresponding to the higher pt higgs candidate
    std::sort(jet_idx_and_btagdisc.begin(), jet_idx_and_btagdisc.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.second > bb.second;
      });
    std::vector<TLorentzVector> b_jet_p;
    for (unsigned int bjet_idx = 0; bjet_idx < 4; bjet_idx++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_eta()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_phi()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_m()->at(jet_idx_and_btagdisc[bjet_idx].first));
      b_jet_p.push_back(lv);
    }
    std::vector<std::vector<unsigned int>> higgs_combinations = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    float min_delta_m = 999;
    bool first = true;
    for (unsigned int higgs_idx = 0; higgs_idx < 3; higgs_idx++) {
       TLorentzVector higgs1_p = b_jet_p[higgs_combinations[higgs_idx][0]]+b_jet_p[higgs_combinations[higgs_idx][1]];
       TLorentzVector higgs2_p = b_jet_p[higgs_combinations[higgs_idx][2]]+b_jet_p[higgs_combinations[higgs_idx][3]];
       float delta_m = TMath::Abs(higgs1_p.M()-higgs2_p.M());
       if ((delta_m < min_delta_m) || first) {
         first = false;
         min_delta_m = delta_m;
       }
    }
    return min_delta_m;
  });

  const NamedFunc jetidmod_hig_cand_drmax("jetidmod_hig_cand_drmax", [](const Baby &b) -> NamedFunc::ScalarType{
    //if no bad jets, just return normal cuts
    if (b.pass_jets()) {
      if (b.njet() < 4) return 999.;
      return b.hig_cand_drmax()->at(0);
    }
    //otherwise, we have to recalculate higgs candidates
    unsigned int njet = 0;
    std::vector<std::pair<unsigned int, float>> jet_idx_and_btagdisc;
    bool recalculate_higgs = false;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) {
        if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) recalculate_higgs = true;
        continue;
      }
      njet += 1;
      jet_idx_and_btagdisc.push_back(std::make_pair(jet_idx, b.jet_deepcsv()->at(jet_idx)));
    }
    //njet, nbt, dphi cuts, also if bad jet was the non-higgs cand jet, just return original cuts
    if (njet < 4) return 999.;
    if (b.njet()<4) std::cout << "ERROR: njet<4\n";
    if (b.hig_cand_am()->size()==0) std::cout << "ERROR: njet>=4 but no hig cand\n";
    if (!recalculate_higgs) return b.hig_cand_drmax()->at(0);
    //get the lowest Delta m Higgs candidate pairing with the first two indices corresponding to the higher pt higgs candidate
    std::sort(jet_idx_and_btagdisc.begin(), jet_idx_and_btagdisc.end(), 
      [](const std::pair<int, float> &a, const std::pair<int, float> &bb) -> bool {
        return a.second > bb.second;
      });
    std::vector<TLorentzVector> b_jet_p;
    for (unsigned int bjet_idx = 0; bjet_idx < 4; bjet_idx++) {
      TLorentzVector lv;
      lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_eta()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_phi()->at(jet_idx_and_btagdisc[bjet_idx].first), 
        b.jet_m()->at(jet_idx_and_btagdisc[bjet_idx].first));
      b_jet_p.push_back(lv);
    }
    std::vector<std::vector<unsigned int>> higgs_combinations = {{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    float min_delta_m = 999;
    float min_dm_drmax = 0;
    bool first = true;
    for (unsigned int higgs_idx = 0; higgs_idx < 3; higgs_idx++) {
       TLorentzVector higgs1_p = b_jet_p[higgs_combinations[higgs_idx][0]]+b_jet_p[higgs_combinations[higgs_idx][1]];
       TLorentzVector higgs2_p = b_jet_p[higgs_combinations[higgs_idx][2]]+b_jet_p[higgs_combinations[higgs_idx][3]];
       float delta_m = TMath::Abs(higgs1_p.M()-higgs2_p.M());
       if ((delta_m < min_delta_m) || first) {
         first = false;
         min_delta_m = delta_m;
	 float dr_1 = deltaR(b_jet_p[higgs_combinations[higgs_idx][0]].Eta(), b_jet_p[higgs_combinations[higgs_idx][0]].Phi(), b_jet_p[higgs_combinations[higgs_idx][1]].Eta(),b_jet_p[higgs_combinations[higgs_idx][1]].Phi());
	 float dr_2 = deltaR(b_jet_p[higgs_combinations[higgs_idx][2]].Eta(), b_jet_p[higgs_combinations[higgs_idx][2]].Phi(), b_jet_p[higgs_combinations[higgs_idx][3]].Eta(),b_jet_p[higgs_combinations[higgs_idx][3]].Phi());
	 min_dm_drmax = dr_1 > dr_2 ? dr_1 : dr_2;
       }
    }
    return min_dm_drmax;
  });

  const NamedFunc pass_filters = "npv_good>0&&pass_cschalo_tight&&pass_hbhe&&pass_hbheiso&&pass_ecaldeadcell&&pass_badpfmu&&pass_eebadsc&&pass_muon_jet&&(met/met_calo<2.0)&&(met/mht<2.0)&&pass_low_neutral_jet&&pass_htratio_dphi_tight" && pass_ecalnoisejet;

  //----------------------------------------------------------------------------------
  //                                plots n stuff
  //----------------------------------------------------------------------------------

  PlotMaker pm;

  //// 2b, (met 0, 75, 150, 200, 300), low drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

  //// 3b, (met 0, 75, 150, 200, 300), low drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

  //// 4b, (met 0, 75, 150, 200, 300), low drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

  //// 2b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__2b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

  //// 3b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__3b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

  //// 4b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies__4b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
  
  //various plots with all jets
  //met
  pm.Push<Hist1D>(Axis(10, 0, 600, "met", "p_{T}^{miss} [GeV]", {150,200,300,400}),
    pass_filters &&
    "!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_met_alljets").LuminosityTag(total_luminosity_string);
  //ht
  pm.Push<Hist1D>(Axis(10, 0, 800, "ht", "H_{T} [GeV]", {}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_ht_alljets").LuminosityTag(total_luminosity_string);
  //nb
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, nb, "N_{b}", {2.5}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_nb_alljets").LuminosityTag(total_luminosity_string);
  //ambb
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_alljets").LuminosityTag(total_luminosity_string);
  //dm
  pm.Push<Hist1D>(Axis(10, 0, 100, "hig_cand_dm[0]", "#Delta m [GeV]", {40}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_dm_alljets").LuminosityTag(total_luminosity_string);
  //drmax
  pm.Push<Hist1D>(Axis(16, 0, 4, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1,2.2}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_drmax_alljets").LuminosityTag(total_luminosity_string);
  //ambb - nb2
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&nbt==2&&nbm==2&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_nb2_alljets").LuminosityTag(total_luminosity_string);
  //ambb - nb>=3
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&nbm>=3&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_nb34_alljets").LuminosityTag(total_luminosity_string);

  //various plots with only ID'd jets
  //met
  pm.Push<Hist1D>(Axis(10, 0, 600, "met", "p_{T}^{miss} [GeV]", {150,200,300,400}),
    pass_filters &&
    !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_met_idjets").LuminosityTag(total_luminosity_string);
  //ht
  pm.Push<Hist1D>(Axis(10, 0, 800, "ht", "H_{T} [GeV]", {}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_ht_idjets").LuminosityTag(total_luminosity_string);
  //nb
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, jetidmod_nb, "N_{b}", {2.5}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_log).Weight(weight).Tag("FixName:jetid_studies_nb_idjets").LuminosityTag(total_luminosity_string);
  //ambb
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_idjets").LuminosityTag(total_luminosity_string);
  //dm
  pm.Push<Hist1D>(Axis(10, 0, 100, jetidmod_hig_cand_dm, "#Delta m [GeV]", {40}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_dm_idjets").LuminosityTag(total_luminosity_string);
  //drmax
  pm.Push<Hist1D>(Axis(16, 0, 4, jetidmod_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_drmax_idjets").LuminosityTag(total_luminosity_string);
  //ambb - nb2
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2==jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_nb2_idjets").LuminosityTag(total_luminosity_string);
  //ambb - nb>=3
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs, plt_lin).Weight(weight).Tag("FixName:jetid_studies_ambb_nb34_idjets").LuminosityTag(total_luminosity_string);

  //various ratio plots with all jets
  //met
  pm.Push<Hist1D>(Axis(10, 0, 600, "met", "p_{T}^{miss} [GeV]", {150,200,300,400}),
    pass_filters &&
    "!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_met_alljets").LuminosityTag(total_luminosity_string);
  //ht
  pm.Push<Hist1D>(Axis(10, 0, 800, "ht", "H_{T} [GeV]", {}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ht_alljets").LuminosityTag(total_luminosity_string);
  //nb
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, nb, "N_{b}", {2.5}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_nb_alljets").LuminosityTag(total_luminosity_string);
  //ambb
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_alljets").LuminosityTag(total_luminosity_string);
  //dm
  pm.Push<Hist1D>(Axis(10, 0, 100, "hig_cand_dm[0]", "#Delta m [GeV]", {40}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_dm_alljets").LuminosityTag(total_luminosity_string);
  //drmax
  pm.Push<Hist1D>(Axis(16, 0, 4, "hig_cand_drmax[0]", "#Delta R_{max}", {1.1,2.2}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&"
    "hig_cand_dm[0]<=40&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_drmax_alljets").LuminosityTag(total_luminosity_string);
  //ambb - nb2
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&nbt==2&&nbm==2&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_nb2_alljets").LuminosityTag(total_luminosity_string);
  //ambb - nb>=3
  pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met&&!low_dphi_met&&"
    "nlep==1&&mt<=100&&4<=njet&&njet<=5&&2<=nbt&&nbm>=3&&"
    "hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200"
    && Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_nb34_alljets").LuminosityTag(total_luminosity_string);

  //various ratio plots with only ID'd jets
  //met
  pm.Push<Hist1D>(Axis(10, 0, 600, "met", "p_{T}^{miss} [GeV]", {150,200,300,400}),
    pass_filters &&
    !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_met_idjets").LuminosityTag(total_luminosity_string);
  //ht
  pm.Push<Hist1D>(Axis(10, 0, 800, "ht", "H_{T} [GeV]", {}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ht_idjets").LuminosityTag(total_luminosity_string);
  //nb
  pm.Push<Hist1D>(Axis(3, 1.5, 4.5, jetidmod_nb, "N_{b}", {2.5}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_nb_idjets").LuminosityTag(total_luminosity_string);
  //ambb
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_idjets").LuminosityTag(total_luminosity_string);
  //dm
  pm.Push<Hist1D>(Axis(10, 0, 100, jetidmod_hig_cand_dm, "#Delta m [GeV]", {40}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_dm_idjets").LuminosityTag(total_luminosity_string);
  //drmax
  pm.Push<Hist1D>(Axis(16, 0, 4, jetidmod_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<=jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_drmax_idjets").LuminosityTag(total_luminosity_string);
  //ambb - nb2
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2==jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_nb2_idjets").LuminosityTag(total_luminosity_string);
  //ambb - nb>=3
  pm.Push<Hist1D>(Axis(10, 0, 200, jetidmod_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
    pass_filters &&
    "150<=met" && !jetidmod_low_dphi_met &&
    "nlep==1&&mt<=100" && (4<=jetidmod_njet) && (jetidmod_njet<=5) && (2<jetidmod_nb) &&
    (jetidmod_hig_cand_dm<=40) && (jetidmod_hig_cand_drmax<2.2) && (jetidmod_hig_cand_am <= 200) &&
    Higfuncs::lead_signal_lepton_pt>30,
    ttbar_procs_ratio, plt_shapes).Weight(weight).Tag("FixName:jetid_studies_shapes_ambb_nb34_idjets").LuminosityTag(total_luminosity_string);

  //// comparison of data with mc
  //// [drmax;no drmax] (2b, 3b+4b)
  //pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "(nbt==2&&nbm==2)"
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_drmax_2b").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_drmax_3b4b").LuminosityTag(total_luminosity_string);

  //// [<m>;no <m>, low_drmax] (2b, 3b+4b)
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "hig_cand_drmax[0]<1.1&&"
  //  "(nbt==2&&nbm==2)" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "hig_cand_drmax[0]<1.1&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax").LuminosityTag(total_luminosity_string);
  //// [<m>;no <m>, high_drmax] (2b, 3b+4b)
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "hig_cand_drmax[0]>=1.1&&"
  //  "(nbt==2&&nbm==2)" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_highdrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "hig_cand_drmax[0]>=1.1&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_highdrmax").LuminosityTag(total_luminosity_string);

  //// Extra plots [<m>;no <m>] (2b, 3b+4b)
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "(nbt==2&&nbm==2)" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  //"nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
  //  &&Higfuncs::lead_signal_lepton_pt>30,
  //  ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b").LuminosityTag(total_luminosity_string);


  //// [<m> (0b data, 1b data); (low, high) drmax]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
  //  ttbar_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax__data").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
  //  ttbar_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax__data").LuminosityTag(total_luminosity_string);

  //// [nisr shapes (2b,3b,4b);(low, high) drmax]
  //pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_nisr__lowdrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_nisr__highdrmax").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (nisr=0,1,2);(low, high) drmax]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
  //  ttbar_procs_nisr, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_nisr__lowdrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
  //  ttbar_procs_nisr, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_nisr__highdrmax").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);(low, high) drmax]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);very high drmax]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&
  //  "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
  //  "hig_cand_drmax[0]>=2.2",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__veryhighdrmax").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);low drmax, met (0,75,150,200,300)]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>0&&met<=75",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_0").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>75&&met<=150",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_75").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>150&&met<=200",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_150").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>200&&met<=300",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_200").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>300",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_300").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);high drmax, met (0,75,150,200,300)]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>0&&met<=75",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_0").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>75&&met<=150",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_75").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>150&&met<=200",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_150").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>200&&met<=300",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_200").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>300",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_300").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);low drmax, nisr (0,1,2)]
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==0",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_0").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_1").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==2",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_2").LuminosityTag(total_luminosity_string);
  //// [<m> shapes (2b, 3b, 4b);high drmax, nisr (0,1,2)]
  //// [TODO] Strange that drmax[1.1,2.2] shows a correlation. Doesn't agree with nisr explanation.
  //// Could study with xy, xy. Is b[pt high=x, high low=y] close to t1b, t1w, t2b, t2b, none==isr? => bb = 16 combinations, bbbb = 256 combinations
  //// Could study with xy, xy. Is b[pt high=x, high low=y] close to t1(b,w), t2(b,w), none==isr? => bb = 9 combinations, bbbb = 81 combinations
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==0",
  //  base_filters&&
  //  "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
  //  "hig_cand_drmax[0]>=1.1&&nisr==0",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_0").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==1",
  //  base_filters&&
  //  "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
  //  "hig_cand_drmax[0]>=1.1&&nisr==1",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_1").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==2",
  //  base_filters&&
  //  "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
  //  "hig_cand_drmax[0]>=1.1&&nisr==2",
  //  ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_2").LuminosityTag(total_luminosity_string);

  //// kappa plot
  //system(("./run/higgsino/plot_kappas.exe --sample ttbar --scen mc_as_data --year "+year_string).c_str());
  // Data example: ./run/higgsino/plot_kappas.exe --year 2016 --unblind --sample ttbar --scen data

  //// According to true b
  //// 2b, (met 0, 75, 150, 200, 300), low drmax 
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

  //// 3b, (met 0, 75, 150, 200, 300), low drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

  //// 4b, (met 0, 75, 150, 200, 300), low drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

  //// 2b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

  //// 3b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

  //// 4b, (met 0, 75, 150, 200, 300), high drmax
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  //pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

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
