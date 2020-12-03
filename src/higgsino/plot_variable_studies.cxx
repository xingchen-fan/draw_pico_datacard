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
#include "TMath.h"
#include "TRandomGen.h"
#include "TLorentzVector.h"

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
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "2016_all";
  bool unblind = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  PlotOpt plot_type("txt/plot_styles.txt", "CMSPaper");
  vector<PlotOpt> plt_lin = {plot_type.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::linear).Stack(StackType::data_norm).LegendColumns(3)};
  vector<PlotOpt> plt_lin_nooverflow = {plot_type.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::linear).Overflow(OverflowType::none).Stack(StackType::data_norm).LegendColumns(3)};
  vector<PlotOpt> plt_log = {plot_type.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::log).LogMinimum(.001).Stack(StackType::data_norm).LegendColumns(3)};
  vector<PlotOpt> plt_shapes = {plot_type.Title(TitleType::info).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::shapes).LegendColumns(3)};
  vector<PlotOpt> plt_lin_shapes = {plot_type.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::linear).Stack(StackType::data_norm).LegendColumns(3), plot_type.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::linear).Stack(StackType::shapes).LegendColumns(3)};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).YAxis(YAxisType::log).Overflow(OverflowType::overflow)};

  // Set options
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string mc_skim_folder = "mc/skim_met150/";
  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";
  string foldersig = mc_base_folder+year_string+"/SMS-TChiHH_2D/unskimmed/";
  if (year_string == "2016_all") foldersig = mc_base_folder+"/2016/SMS-TChiHH_2D/unskimmed/";

  set<int> years;
  //HigUtilities::parseYears(year_string, years);
  //years = {2016, 2017, 2018};
  years = {2016};
  if (year_string=="2017") years = {2017};
  if (year_string=="2018") years = {2018};
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  if (year_string == "2016_all") total_luminosity = 137.4;
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  NamedFunc weight_notrgeff = "w_lumi*w_isr"*w_years;

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_SingleLept*","TTJets_DiLept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
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

  vector<shared_ptr<Process> > search_signal_procs;
  vector<shared_ptr<Process> > search_procs;
  vector<shared_ptr<Process> > abcd_ttbar_procs;
  vector<shared_ptr<Process> > ttbar_procs;
  //signal settings
  //If you modify mchi_sigm, modify mchi_sigm_int in the NamedFunc below
  vector<string> mchi_sigm = {"175","500","900"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  vector<string> mchi_sigm2d = {"250","350","450"}; 
  vector<string> mlsp_sigm2d = {"50","200","100"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue}; // need sigm.size() <= sig_colors.size()
  vector<int> sig_colors2d = {kOrange, kYellow, kCyan};
  // Set mc processes
  search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  search_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));
  // Set ttbar
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  // Set ttbar for kappa plots
  abcd_ttbar_procs.push_back(Process::MakeShared<Baby_pico>("All bkgs 2b", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nbt==2&&nbm==2"));
  abcd_ttbar_procs.push_back(Process::MakeShared<Baby_pico>("All bkgds 3b", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nbt>=2&&nbm==3&&nbl==3"));
  abcd_ttbar_procs.push_back(Process::MakeShared<Baby_pico>("All bkgds 4b", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nbt>=2&&nbm>=3&&nbl>=4"));
  for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
    search_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
    search_signal_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
      std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
  }
  for (unsigned isig(1); isig<mchi_sigm2d.size(); isig++) {
    search_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
    search_signal_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
      std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
  }

  NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  //give the right weight for 2D Higgsino scans, 1D scans, and other stuff
  //also applies trigger efficiencies
  const NamedFunc mixed_model_weight("mixed_model_weight",[](const Baby &b) -> NamedFunc::ScalarType{
    set<int> mchi_sigm_int = {175, 500, 900}; 
    if (b.mprod() == -999 || mchi_sigm_int.count(b.mprod()) > 0) {
      //not signal
      return Higfuncs::final_weight.GetScalar(b);
    }
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return Higfuncs::final_weight.GetScalar(b)/xsec1d*xsec2d;
  });

  const NamedFunc nb77("nb77",[](const Baby &b) -> NamedFunc::ScalarType{
    //namedfunc simulating ATLAS WP77 b-tagging
    unsigned int r_nb77 = 0;
    TRandomMixMax rndgen;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_isgood()->size(); jet_idx++) {
      if (b.jet_isgood()->at(jet_idx)) {
        switch (b.jet_pflavor()->at(jet_idx)) {
          case -5:
          case 5: //b-quark
            if (rndgen.Uniform() < 0.77) r_nb77++;
            break;
          case -4:
          case 4: //c-quark
            if (rndgen.Uniform() < 0.17) r_nb77++;
            break;
          case -15:
          case 15: //hadronic tau
            if (rndgen.Uniform() < 0.045) r_nb77++;
            break;
          default: //usdg
            if (rndgen.Uniform() < 0.0075) r_nb77++;
        }
      }
    }
    //if (b.mprod() == -999) {
    //  if (r_nb77 >= 4) cout << "Found " << r_nb77 << "b signal event!" << endl;
    //}
    return r_nb77;
  });

  const NamedFunc nb70("nb70",[](const Baby &b) -> NamedFunc::ScalarType{
    //namedfunc simulating ATLAS WP70 b-tagging
    unsigned int r_nb70 = 0;
    TRandom3 rndgen;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_isgood()->size(); jet_idx++) {
      if (b.jet_isgood()->at(jet_idx)) {
        switch (b.jet_pflavor()->at(jet_idx)) {
          case -5:
          case 5: //b-quark
            if (rndgen.Uniform() < 0.70) r_nb70++;
            break;
          case -4:
          case 4: //c-quark
            if (rndgen.Uniform() < 0.083) r_nb70++;
            break;
          case -15:
          case 15: //hadronic tau
            if (rndgen.Uniform() < 0.018) r_nb70++;
            break;
          default: //usdg
            if (rndgen.Uniform() < 0.0026) r_nb70++;
        }
      }
    }
    return r_nb70;
  });

  const NamedFunc analysis_nb("analysis_nb",[](const Baby &b) -> NamedFunc::ScalarType{
    unsigned int r_analysis_nb =  0;
    if (b.nbt() < 2) {
      r_analysis_nb =  b.nbt();
    }
    else if (b.nbm() < 3) {
      r_analysis_nb = 2;
    }
    else if (b.nbl() < 4) {
      r_analysis_nb = 3;
    }
    else {
      r_analysis_nb = 4;
    }
    return r_analysis_nb;
  });

  const NamedFunc meff("meff",[](const Baby &b) -> NamedFunc::ScalarType{
    float ht_hig = 0;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_isgood()->size(); jet_idx++) {
      if (b.jet_h1d()->at(jet_idx) || b.jet_h2d()->at(jet_idx)) {
        ht_hig += b.jet_pt()->at(jet_idx);
      }
    }
    return b.met() + ht_hig;
  });

  const NamedFunc mtmin("mtmin",[](const Baby &b) -> NamedFunc::ScalarType{
    float r_mtmin = 9999;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_isgood()->size(); jet_idx++) {
      //if medium DeepCSV, calculate mT
      if (b.jet_isgood() && (b.jet_deepcsv()->at(jet_idx) > 0.6321)) {
        float jet_px = TMath::Cos(b.jet_phi()->at(jet_idx))*(b.jet_pt()->at(jet_idx)), jet_py = TMath::Sin(b.jet_phi()->at(jet_idx))*(b.jet_pt()->at(jet_idx)), jet_pt = (b.jet_pt()->at(jet_idx));
        float met_px = TMath::Cos(b.met_phi())*(b.met()), met_py = TMath::Sin(b.met_phi())*(b.met()), met_pt = (b.met());
        float mt = TMath::Sqrt((met_pt+jet_pt)*(met_pt+jet_pt)-(jet_px+met_px)*(jet_px+met_px)-(jet_py+met_py)*(jet_py+met_py));
	if (mt < r_mtmin) r_mtmin = mt;
      }
    }
    return r_mtmin;
  });

  const NamedFunc mt2("mt2",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector higgs_candidate1, higgs_candidate2;
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (b.jet_h1d()->at(jet_idx)) {
        TLorentzVector lv;
	lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx),b.jet_eta()->at(jet_idx),b.jet_phi()->at(jet_idx),b.jet_m()->at(jet_idx));
	higgs_candidate1 += lv;
      }
      else if (b.jet_h2d()->at(jet_idx)) {
        TLorentzVector lv;
	lv.SetPtEtaPhiM(b.jet_pt()->at(jet_idx),b.jet_eta()->at(jet_idx),b.jet_phi()->at(jet_idx),b.jet_m()->at(jet_idx));
	higgs_candidate2 += lv;
      }
    }
    float pt1 = higgs_candidate1.Pt();
    float pt2 = higgs_candidate2.Pt();
    float phi1 = higgs_candidate1.Phi();
    float phi2 = higgs_candidate2.Phi();
    float m1 = higgs_candidate1.M();
    float m2 = higgs_candidate2.M();
    float mlsp = 150;
    float et1 = TMath::Sqrt(pt1*pt1+m1*m1);
    float et2 = TMath::Sqrt(pt2*pt2+m2*m2);
    float etlsp = TMath::Sqrt(b.met()*b.met()+mlsp*mlsp);
    float higgs1_mt = TMath::Sqrt(m1*m1+mlsp*mlsp+2.0*(et1*etlsp-pt1*b.met()*TMath::Cos(phi1-b.met_phi())));
    float higgs2_mt = TMath::Sqrt(m2*m2+mlsp*mlsp+2.0*(et2*etlsp-pt2*b.met()*TMath::Cos(phi2-b.met_phi())));
    //logic below assumes mlsp = 0, otherwise, mt would need to be calculated for ptchi1 = -pt1-pt2-ptchi2
    if (higgs1_mt>=(m2+mlsp)) {
      if (higgs2_mt>=(m1+mlsp)) {
        //balanced solution
        float at = (et1*et2+pt1*pt2*TMath::Cos(phi1-phi2));
        return TMath::Sqrt(mlsp*mlsp+at+TMath::Sqrt((1.0+4.0*mlsp*mlsp/(2.0*at-m1*m1-m2*m2))*(at*at-m1*m1*m2*m2)));
      }
      else {
        //unbalanced solution
        return m1+mlsp;
      }
    }
    else {
      //unbalanced solution
      return m2+mlsp;
    }
    return 0.;
  });

  const NamedFunc atlas_hm_bin("atlas_hm_bin",[nb77, mtmin, meff](const Baby &b) -> NamedFunc::ScalarType{
    float event_meff = meff.GetScalar(b);
    unsigned int event_nb77 = nb77.GetScalar(b);
    unsigned int event_njet = b.njet();
    float event_drmax = b.hig_cand_drmax()->at(0);
    float event_mtmin = mtmin.GetScalar(b);
    if (event_meff >= 1100) {
      if (event_mtmin>130 && event_njet<=5 && event_nb77>=3) {
        return 3; //high meff
      }
    }
    else if (event_meff >= 850) {
      if (event_nb77 >= 4) {
        if (event_drmax<1.4) {
	  if (event_njet <= 6) {
	    return 6; //mid meff, 4b, low drmax
	  }
	}
	else if (event_drmax < 2.4) {
	  if (event_njet <= 6) {
	    return 7; //mid meff, 4b, high drmax
	  }
	}
      }
      else if (event_nb77 == 3) {
	if (event_mtmin>150 && event_njet<=5) {
          return 2; //mid meff, 3b
	}
      }
    }
    else if (event_meff >= 600) {
      if (event_nb77 >= 4) {
        if (event_drmax<1.4) {
	  if (event_njet <= 5) {
	    return 4; //low meff, 4b, low drmax
	  }
	}
	else if (event_drmax < 2.4) {
	  if (event_njet <= 5) {
	    return 5; //low meff, 4b, high drmax
	  }
	}
      }
      else if (event_nb77 == 3) {
	if (event_mtmin>150 && event_njet<=5) {
          return 1; //low meff, 3b
	}
      }
    }
    return 0;
  });

  const NamedFunc mh1("mh1",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.hig_cand_am()->at(0)+(b.hig_cand_dm()->at(0))/2.0;
  });

  const NamedFunc mh2("mh2",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.hig_cand_am()->at(0)-(b.hig_cand_dm()->at(0))/2.0;
  });

  const NamedFunc max_top_tag("max_top_tag",[](const Baby &b) -> NamedFunc::ScalarType{
    float r_max_top_tag = 0;
    for (unsigned int fat_jet_idx = 0; fat_jet_idx < b.fjet_deep_tvsqcd()->size(); fat_jet_idx++) {
      r_max_top_tag = r_max_top_tag > b.fjet_deep_tvsqcd()->at(fat_jet_idx) ? r_max_top_tag : b.fjet_deep_tvsqcd()->at(fat_jet_idx);
    }
    return r_max_top_tag;
  });

  std::vector<double> cms_met_bins{150,200,300,400,500};

  NamedFunc cms_baseline = Higfuncs::pass_filters&&"stitch&&150<=met&&nvlep==0&&ntk==0&&4<=njet&&njet<=5&&2<=nbt&&!low_dphi_met&&(met/met_calo)<2&&(met/mht)<2&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc cms_baseline_nob = "stitch&&150<=met&&nvlep==0&&ntk==0&&4<=njet&&njet<=5&&!low_dphi_met&&(met/met_calo)<2&&(met/mht)<2&&hig_cand_dm[0]<=40&&hig_cand_drmax[0]<2.2&&hig_cand_am[0]<=200";
  NamedFunc atlas_hm_baseline_nob = "stitch&&200<=met&&nvlep==0&&ntk==0&&4<=njet&&jet_met_dphi[0]>=0.4&&jet_met_dphi[1]>=0.4&&jet_met_dphi[2]>=0.4&&jet_met_dphi[3]>=0.4&&(met/met_calo)<2&&(met/mht)<2&&hig_cand_drmax[0]<=2.4" && (110<=mh1) && (mh1<150) && (90<=mh2) && (mh2<140);

  PlotMaker pm;

  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_nb3dr11").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&hig_cand_drmax[0]>=1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_nb3dr22").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_nb4dr11").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]>=1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_nb4dr22").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline_nob && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && nb77 == 3 && "hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_atlasbtag_nb3dr11").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline_nob && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && nb77 == 3 && "hig_cand_drmax[0]>=1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_atlasbtag_nb3dr22").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline_nob && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && nb77 >= 4 && "hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_atlasbtag_nb4dr11").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline_nob && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && nb77 >= 4 && "hig_cand_drmax[0]>=1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_atlasbtag_nb4dr22").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(cms_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
  //    cms_baseline_nob && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && nb77 >= 4 && "hig_cand_drmax[0]>=1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_cms_atlasbtag_nb4dr22").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(5, -0.5, 4.5, analysis_nb, "N_{bCMS}", {}),
  //    "stitch",
  //    search_signal_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_nbcms").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(5, -0.5, 4.5, nb77, "N_{b77}", {}),
  //    "stitch",
  //    search_signal_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:atlascomparison_nb77").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(7, 0.5, 7.5, atlas_hm_bin, "ATLAS bin no.", {}),
  //    atlas_hm_baseline_nob,
  //    search_procs, plt_lin_nooverflow).Weight(mixed_model_weight).Tag("FixName:atlascomparison_atlas_bins").LuminosityTag(total_luminosity_string);
  //mt min distributions
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&met<200&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_log).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met150nb3").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met<200&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met150nb4").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&met>=200&&met<300&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_log).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met200nb3").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>=200&&met<300&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met200nb4").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&met>=300&&met<400&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met300nb3").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>=300&&met<400&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met300nb4").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm==3&&nbl==3&&met>=400&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met400nb3").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 400, mtmin, "m_{Tmin} [GeV]", {}),
  //    cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&nbm>=3&&nbl>=4&&met>=400&&hig_cand_drmax[0]<1.1",
  //    search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:mtmin_lowdrmax_met400nb4").LuminosityTag(total_luminosity_string);
  ////mt min affect on ABCD
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met<200&&hig_cand_drmax[0]<1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met150lodrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met<200&&hig_cand_drmax[0]>=1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met150hidrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=200&&met<300&&hig_cand_drmax[0]<1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met200lodrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=200&&met<300&&hig_cand_drmax[0]>=1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met200hidrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=300&&met<400&&hig_cand_drmax[0]<1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met300lodrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=300&&met<400&&hig_cand_drmax[0]>=1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met300hidrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=400&&hig_cand_drmax[0]<1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met400lodrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=400&&hig_cand_drmax[0]>=1.1",
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met400hidrmax").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met<200&&hig_cand_drmax[0]<1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met150lodrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met<200&&hig_cand_drmax[0]>=1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met150hidrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=200&&met<300&&hig_cand_drmax[0]<1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met200lodrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=200&&met<300&&hig_cand_drmax[0]>=1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met200hidrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=300&&met<400&&hig_cand_drmax[0]<1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met300lodrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=300&&met<400&&hig_cand_drmax[0]>=1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met300hidrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=400&&hig_cand_drmax[0]<1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met400lodrmax_mtmin").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "#LT m_{bb} #GT", {}),
  //    cms_baseline && "met>=400&&hig_cand_drmax[0]>=1.1" && mtmin>150,
  //    abcd_ttbar_procs, plt_shapes).Weight(mixed_model_weight).Tag("FixName:abcd_met400hidrmax_mtmin").LuminosityTag(total_luminosity_string);
  ////2D plots
  //pm.Push<Hist2D>(Axis(30, 0., 300., mtmin, "m_{Tmin} [GeV]", {}),
  //    Axis(40, 0., 200., "hig_cand_am[0]", "<m> [GeV]", {}),
  //    cms_baseline, ttbar_procs, twodim_plotopts).Weight("weight").Tag("FixName:mtmin_higcandam").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist2D>(Axis(30, 0., 300., mtmin, "m_{Tmin} [GeV]", {}),
  //    Axis(3, 1.5, 4.5, analysis_nb, "N_{b}", {}),
  //    cms_baseline, ttbar_procs, twodim_plotopts).Weight("weight").Tag("FixName:mtmin_nb").LuminosityTag(total_luminosity_string);
  
  for (int met_bin_idx = 0; met_bin_idx < 4; met_bin_idx++) {
    NamedFunc met_bin_cuts = "";
    std::string met_filename = "";
    switch (met_bin_idx) {
      case 0:
        met_bin_cuts = "150<=met&&met<200";
	met_filename = "met150to200_";
	break;
      case 1:
        met_bin_cuts = "200<=met&&met<300";
	met_filename = "FixName:met200to300_";
	break;
      case 2:
        met_bin_cuts = "300<=met&&met<400";
	met_filename = "FixName:met300to400_";
	break;
      default:
        met_bin_cuts = "400<=met";
	met_filename = "FixName:met400toInf_";
    }
    for (int drmax_bin_idx = 0; drmax_bin_idx < 2; drmax_bin_idx++) {
      NamedFunc drmax_bin_cuts = "";
      std::string drmax_filename = "";
      switch (drmax_bin_idx) {
        case 0:
          drmax_bin_cuts = "hig_cand_drmax[0]<1.1";
	  drmax_filename = "lowdrmax_";
	  break;
        default:
          drmax_bin_cuts = "hig_cand_drmax[0]>=1.1";
	  drmax_filename = "highdrmax_";
      }
      for (int nb_bin_idx = 3; nb_bin_idx <= 4; nb_bin_idx++) {
        NamedFunc nb_bin_cuts = "";
        std::string nb_filename = "";
        switch (nb_bin_idx) {
	  case 3:
            nb_bin_cuts = "nbm==3&&nbl==3";
	    nb_filename = "nb3";
	    break;
	  default:
            nb_bin_cuts = "nbm>=3&&nbl>=4";
	    nb_filename = "nb4";
	}
	for (int njet_bin_idx = 4; njet_bin_idx <= 5; njet_bin_idx++) {
          NamedFunc njet_bin_cuts = "";
          std::string njet_filename = "";
          switch (njet_bin_idx) {
	    case 4:
              njet_bin_cuts = "njet==4";
	      njet_filename = "_nj4";
	      continue;
	      break;
	    default:
              njet_bin_cuts = "njet>=4";
	  }
	  //code here
          pm.Push<Hist1D>(Axis(20, 0, 1000, mt2, "M_{T2} [GeV]", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&njet==4" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_lin).Weight(mixed_model_weight).Tag(("FixName:varstudies_mt2_lin_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
          pm.Push<Hist1D>(Axis(20, 0, 1000, mt2, "M_{T2} [GeV]", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140&&njet==4" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_log).Weight(mixed_model_weight).Tag(("FixName:varstudies_mt2_log_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
          pm.Push<Hist1D>(Axis(20, 0, 500, mtmin, "m_{Tmin} [GeV]", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_lin).Weight(mixed_model_weight).Tag(("FixName:varstudies_mtmin_lin_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
          pm.Push<Hist1D>(Axis(20, 0, 500, mtmin, "m_{Tmin} [GeV]", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_log).Weight(mixed_model_weight).Tag(("FixName:varstudies_mtmin_log_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
          pm.Push<Hist1D>(Axis(20, 0., 1., max_top_tag, "max top tag score", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_lin).Weight(mixed_model_weight).Tag(("FixName:varstudies_maxtopscore_lin_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
          pm.Push<Hist1D>(Axis(20, 0., 1., max_top_tag, "max top tag score", {}),
              cms_baseline && "100<=hig_cand_am[0]&&hig_cand_am[0]<140" && met_bin_cuts && drmax_bin_cuts && nb_bin_cuts && njet_bin_cuts,
              search_procs, plt_log).Weight(mixed_model_weight).Tag(("FixName:varstudies_maxtopscore_log_"+met_filename+drmax_filename+nb_filename+njet_filename).c_str()).LuminosityTag(total_luminosity_string);
	} // /jet loop
      } // /bjet loop
    } // /drmax loop
  } // /met loop


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
