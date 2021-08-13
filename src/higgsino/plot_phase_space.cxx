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
#include "Math/Vector4D.h"

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
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

Double_t DeltaR(const ROOT::Math::PtEtaPhiMVector & v1, const ROOT::Math::PtEtaPhiMVector & v2)
{
   Double_t deta = v1.Eta()-v2.Eta();
   Double_t dphi = TVector2::Phi_mpi_pi(v1.Phi()-v2.Phi());
   return TMath::Sqrt( deta*deta+dphi*dphi );
}

const NamedFunc mc_true_avg_dr("mc_true_avg_dr", [](const Baby &b) ->NamedFunc::ScalarType{
  // Find b pair relaed with higgs
  // index_higgs_to_b[higgs_index] = [b_index]
  map<int, vector<int> > index_higgs_to_b;
  for (unsigned iParticle = 0; iParticle < (*b.mc_id()).size(); ++iParticle) {
    if (abs((*b.mc_id())[iParticle]) != 5) continue;
    if ((*b.mc_mom())[iParticle] != 25) continue;
    index_higgs_to_b[(*b.mc_momidx())[iParticle]].push_back(iParticle);
  }
  //for(auto it = index_higgs_to_b.begin(); it != index_higgs_to_b.end(); ++it) {
  //  if (it->second.size()!=2) continue;
  //  cout<<(*b.mc_id())[it->first]<<" "<<(*b.mc_id())[it->second[0]]<<" "<<(*b.mc_id())[it->second[1]]<<endl;
  //}
  vector<double> dr_bb;
  for(auto it = index_higgs_to_b.begin(); it != index_higgs_to_b.end(); ++it) {
    if (it->second.size()!=2) continue;
    //int higgs_index = it->first;
    int b0_index = it->second[0];
    int b1_index = it->second[1];
    ROOT::Math::PtEtaPhiMVector higgs_b0((*b.mc_pt())[b0_index], (*b.mc_eta())[b0_index], (*b.mc_phi())[b0_index], (*b.mc_mass())[b0_index]);
    ROOT::Math::PtEtaPhiMVector higgs_b1((*b.mc_pt())[b1_index], (*b.mc_eta())[b1_index], (*b.mc_phi())[b1_index], (*b.mc_mass())[b1_index]);
    dr_bb.push_back(DeltaR(higgs_b0,higgs_b1));
  };
  if (dr_bb.size()!=2) return -1;
  return (dr_bb[0]+dr_bb[1])/2;
});

const NamedFunc mc_true_leadH_dr("mc_true_leadH_dr", [](const Baby &b) ->NamedFunc::ScalarType{
  // Find b pair relaed with higgs
  // index_higgs_to_b[higgs_index] = [b_index]
  map<int, vector<int> > index_higgs_to_b;
  for (unsigned iParticle = 0; iParticle < (*b.mc_id()).size(); ++iParticle) {
    if (abs((*b.mc_id())[iParticle]) != 5) continue;
    if ((*b.mc_mom())[iParticle] != 25) continue;
    index_higgs_to_b[(*b.mc_momidx())[iParticle]].push_back(iParticle);
  }
  vector<double> dr_bb;
  vector<double> higgs_pt;
  for(auto it = index_higgs_to_b.begin(); it != index_higgs_to_b.end(); ++it) {
    if (it->second.size()!=2) continue;
    int higgs_index = it->first;
    int b0_index = it->second[0];
    int b1_index = it->second[1];
    ROOT::Math::PtEtaPhiMVector higgs_b0((*b.mc_pt())[b0_index], (*b.mc_eta())[b0_index], (*b.mc_phi())[b0_index], (*b.mc_mass())[b0_index]);
    ROOT::Math::PtEtaPhiMVector higgs_b1((*b.mc_pt())[b1_index], (*b.mc_eta())[b1_index], (*b.mc_phi())[b1_index], (*b.mc_mass())[b1_index]);
    dr_bb.push_back(DeltaR(higgs_b0,higgs_b1));
    higgs_pt.push_back((*b.mc_pt())[higgs_index]);
  };
  if (dr_bb.size()!=2) return -1;
  if (higgs_pt[0]>higgs_pt[1]) return dr_bb[0];
  else return dr_bb[1];
});
const NamedFunc mc_true_subleadH_dr("mc_true_subleadH_dr", [](const Baby &b) ->NamedFunc::ScalarType{
  // Find b pair relaed with higgs
  // index_higgs_to_b[higgs_index] = [b_index]
  map<int, vector<int> > index_higgs_to_b;
  for (unsigned iParticle = 0; iParticle < (*b.mc_id()).size(); ++iParticle) {
    if (abs((*b.mc_id())[iParticle]) != 5) continue;
    if ((*b.mc_mom())[iParticle] != 25) continue;
    index_higgs_to_b[(*b.mc_momidx())[iParticle]].push_back(iParticle);
  }
  vector<double> dr_bb;
  vector<double> higgs_pt;
  for(auto it = index_higgs_to_b.begin(); it != index_higgs_to_b.end(); ++it) {
    if (it->second.size()!=2) continue;
    int higgs_index = it->first;
    int b0_index = it->second[0];
    int b1_index = it->second[1];
    ROOT::Math::PtEtaPhiMVector higgs_b0((*b.mc_pt())[b0_index], (*b.mc_eta())[b0_index], (*b.mc_phi())[b0_index], (*b.mc_mass())[b0_index]);
    ROOT::Math::PtEtaPhiMVector higgs_b1((*b.mc_pt())[b1_index], (*b.mc_eta())[b1_index], (*b.mc_phi())[b1_index], (*b.mc_mass())[b1_index]);
    dr_bb.push_back(DeltaR(higgs_b0,higgs_b1));
    higgs_pt.push_back((*b.mc_pt())[higgs_index]);
  };
  if (dr_bb.size()!=2) return -1;
  if (higgs_pt[0]<=higgs_pt[1]) return dr_bb[0];
  else return dr_bb[1];
});
const NamedFunc mc_true_leadH_pt("mc_true_leadH_pt", [](const Baby &b) ->NamedFunc::ScalarType{
  // Find b pair relaed with higgs
  // index_higgs_to_b[higgs_index] = [b_index]
  map<int, vector<int> > index_higgs_to_b;
  for (unsigned iParticle = 0; iParticle < (*b.mc_id()).size(); ++iParticle) {
    if (abs((*b.mc_id())[iParticle]) != 5) continue;
    if ((*b.mc_mom())[iParticle] != 25) continue;
    index_higgs_to_b[(*b.mc_momidx())[iParticle]].push_back(iParticle);
  }
  vector<double> dr_bb;
  vector<double> higgs_pt;
  for(auto it = index_higgs_to_b.begin(); it != index_higgs_to_b.end(); ++it) {
    if (it->second.size()!=2) continue;
    int higgs_index = it->first;
    int b0_index = it->second[0];
    int b1_index = it->second[1];
    ROOT::Math::PtEtaPhiMVector higgs_b0((*b.mc_pt())[b0_index], (*b.mc_eta())[b0_index], (*b.mc_phi())[b0_index], (*b.mc_mass())[b0_index]);
    ROOT::Math::PtEtaPhiMVector higgs_b1((*b.mc_pt())[b1_index], (*b.mc_eta())[b1_index], (*b.mc_phi())[b1_index], (*b.mc_mass())[b1_index]);
    dr_bb.push_back(DeltaR(higgs_b0,higgs_b1));
    higgs_pt.push_back((*b.mc_pt())[higgs_index]);
  };
  if (dr_bb.size()!=2) return -1;
  if (higgs_pt[0]>higgs_pt[1]) return higgs_pt[0];
  else return higgs_pt[1];
});

namespace{
  bool single_thread = false;
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
  //PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::lumi_shapes).Bottom(BottomType::ratio).RatioMaximum(0.5).RatioMinimum(0.0);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::info)};

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};

  // Set options
  string mc_base_folder = "/net/cms25/cms25r0/jbkim/pico/NanoAODv5/higgsino_humboldt/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";

  //string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";

  set<int> years;
  years = {2016, 2017, 2018};
  //years = {2016};

  NamedFunc weight = Higfuncs::final_weight;

  vector<shared_ptr<Process> > procs;

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
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
  //// Set mc processes
  //procs.push_back(Process::MakeShared<Baby_pico>("tt+x", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  //procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  //procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  //procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  //procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  // Set signal mc
  //vector<string> signal_mass = {"200", "450", "700", "950"}; 
  vector<string> signal_mass = {"950"}; 
  //vector<string> signal_mass = {}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need signal_mass.size() >= sig_colors.size()
  // 1D signal
  for (unsigned isig(0); isig<signal_mass.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+signal_mass[isig]+",1)", Process::Type::background, 
      sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+signal_mass[isig]+"_mLSP-0*.root"}), "1"));
  }
  //// 2D signal
  ////procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(500,150)", Process::Type::background, 
  ////  sig_colors[0], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-500_mLSP-150*.root"}), "1"));
  //procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(300,150)", Process::Type::background, 
  //  //sig_colors[0], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-950_mLSP-600*.root"}), "1"));
  //  sig_colors[0], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-*.root"}), "1"));

  NamedFunc base_filters = Higfuncs::final_pass_filters;

  // resolved cuts
  NamedFunc base_resolved = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)||(nbt==2&&nbm==2))";
  NamedFunc signal_resolved = "(nbt>=2&&nbm>=3&&nbl>=4)&&hig_cand_am[0]<140&&hig_cand_am[0]>100";

  // boosted cuts
  NamedFunc base_boosted = "met>300 && ht>300 && !low_dphi_met && nvlep==0 && ntk==0 && ht>600 && nfjet>=2 &&"
                        "fjet_pt[0]>300 && fjet_pt[1]>300 && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 &&"
                        "fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  NamedFunc signal_boosted = "fjet_deep_md_hbb_btv[0]>0.7 && fjet_deep_md_hbb_btv[1]>0.7 &&"
                          "fjet_msoftdrop[0]>95 && fjet_msoftdrop[0]<=145&&"
                          "fjet_msoftdrop[1]>95 && fjet_msoftdrop[1]<=145";

  PlotMaker pm;

  //// Plot average mc_true deltaR_bb
  ////pm.Push<Hist1D>(Axis(16, 0, 3.2, mc_true_avg_dr, "true avg. \\Delta R_{bb}", {}),
  ////  base_filters&&base_resolved&&signal_resolved, procs, plt_lin).Weight(weight).Tag("ShortName:ResolvedSignalRegion");
  //pm.Push<Hist2D>(Axis(16, 0, 3.2, mc_true_leadH_dr, "Leading H \\Delta R_{bb}", {}),
  //  Axis(20, 0, 500, mc_true_leadH_pt, "Leading H pT [GeV]", {}),
  //  base_filters, procs, plt_2D).Weight(weight);
  pm.Push<Hist2D>(Axis(16, 0, 3.2, mc_true_leadH_dr, "Leading H \\DeltaR_{bb}", {}),
    Axis(16, 0, 3.2, mc_true_subleadH_dr, "Sub-leading H \\DeltaR_{bb}", {}),
    base_filters, procs, plt_2D).Weight(weight).Tag("FixName:leadHDeltaR_vs_subleadHDeltaR");
  //pm.Push<Hist2D>(Axis(11, 0, 2.2, mc_true_leadH_dr, "Leading H \\Delta R_{bb}", {}),
  //  Axis(11, 0, 2.2, mc_true_subleadH_dr, "Subleading H \\Delta R_{bb}", {}),
  //  base_filters&&base_resolved&&signal_resolved, procs, plt_2D).Weight(weight).Tag("ShortName:ResolvedSignalRegion");
  //pm.Push<Hist2D>(Axis(11, 0, 2.2, mc_true_leadH_dr, "Leading H \\Delta R_{bb}", {}),
  //  Axis(11, 0, 2.2, mc_true_subleadH_dr, "Subleading H \\Delta R_{bb}", {}),
  //  base_filters&&base_boosted&&signal_boosted, procs, plt_2D).Weight(weight).Tag("ShortName:BoostedSignalRegion");

  //// Plot drmax vs met
  //pm.Push<Hist2D>(
  //  Axis(20, 0, 1000, "met", "p_{T}^{miss} [GeV]", {150, 200, 300, 400}),
  //  Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
  //  base_filters&&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<3.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)||(nbt==2&&nbm==2))"&&
  //  signal_resolved, procs, plt_2D).Weight(weight).Tag("FixName:drmax_vs_met");

  ////string signalName = "*TChiHH_mChi-450_mLSP-0*.root";
  ////string signalName = "*TChiHH_mChi-900_mLSP-0*.root";
  //string signalName = "*TChiHH_mChi-1400_mLSP-0*.root";

  //vector<shared_ptr<Process> > procs_analysis;
  //procs_analysis.push_back(Process::MakeShared<Baby_pico>("all", Process::Type::background,sig_colors[0],
  //                attach_folder(sig_base_folder, years, sig_skim_folder, {signalName}),
  //                base_filters&&"met>150&&ntk==0&&!low_dphi_met&&nvlep==0"));
  //procs_analysis.push_back(Process::MakeShared<Baby_pico>("resolved", Process::Type::background,sig_colors[1],
  //                attach_folder(sig_base_folder, years, sig_skim_folder, {signalName}),
  //                base_filters&&base_resolved&&signal_resolved));
  //procs_analysis.push_back(Process::MakeShared<Baby_pico>("boosted", Process::Type::background,sig_colors[2],
  //                attach_folder(sig_base_folder, years, sig_skim_folder, {signalName}),
  //                base_filters&&base_boosted&&signal_boosted));

  //pm.Push<Hist1D>(Axis(16, 0, 3.2, mc_true_avg_dr, "true avg. \\Delta R_{bb}", {}),
  //  base_filters, procs_analysis, plt_shapes_info).Weight(weight);
  // Make process for boosted
  // Make process for resolved
  // Plot mc_true deltaR_bb
  
  //// Resolved critical variables: mass, nb, met, maxdr, 
  //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
  //  base_resolved&&signal_resolved, procs, plt_lin).Weight(weight).Tag("ShortName:ResolvedSignalRegion");
  //pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::higd_bcat, "Number b tags", {}),
  //  base_resolved&&signal_resolved, procs, plt_lin).Weight(weight).Tag("ShortName:ResolvedSignalRegion");
  //pm.Push<Hist1D>(Axis(20,0,1000,"met", "p_{T}^{miss} [GeV]", {150, 200, 300, 400}),
  //  base_resolved&&signal_resolved, procs, plt_lin).Weight(weight).Tag("ShortName:ResolvedSignalRegion");
  //pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1,2.2}),
  //  base_resolved&&signal_resolved, procs, plt_lin).Weight(weight).Tag("ShortName:ResolvedSignalRegion");

  //// Boosted critical variablees: mass, nb, met
  //pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[0]", "Leading AK8 softdrop mass [GeV]", {95, 145}),
  //  base_boosted&&signal_boosted, procs, plt_lin).Weight(weight).Tag("ShortName:BoostedSignalRegion");
  //pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "Subleading AK8 softdrop mass [GeV]", {95, 145}),
  //  base_boosted&&signal_boosted, procs, plt_lin).Weight(weight).Tag("ShortName:BoostedSignalRegion");
  //pm.Push<Hist1D>(Axis(10, 0, 1, "fjet_deep_md_hbb_btv[0]", "Leading AK8 b discriminator", {0.7}),
  //  base_boosted&&signal_boosted, procs, plt_lin).Weight(weight).Tag("ShortName:BoostedSignalRegion");
  //pm.Push<Hist1D>(Axis(10, 0, 1, "fjet_deep_md_hbb_btv[1]", "Subleading AK8 b discriminator", {0.7}),
  //  base_boosted&&signal_boosted, procs, plt_lin).Weight(weight).Tag("ShortName:BoostedSignalRegion");
  //pm.Push<Hist1D>(Axis(20,0,1000,"met", "p_{T}^{miss} [GeV]", {300, 500, 700}),
  //  base_boosted&&signal_boosted, procs, plt_lin).Weight(weight).Tag("ShortName:BoostedSignalRegion");

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  // Make custom pretty plot
  float leftMargin = 0.13;
  float rightMargin = 0.175;
  float bottomMargin = 0.14;
  float topMargin = 0.075;
  TCanvas canvas("canvas", "canvas", 850, 600);
  canvas.SetTicks(1,1);
  canvas.SetFillStyle(4000);
  canvas.SetMargin(leftMargin, rightMargin, bottomMargin, topMargin);
  // Draw histogram
  Hist2D * drPhaseSpace = static_cast<Hist2D*>(pm.GetFigure("leadHDeltaR_vs_subleadHDeltaR").get());
  TH2D hist_drPhaseSpace = drPhaseSpace->GetBkgHist(1);
  hist_drPhaseSpace.Draw("colz");
  // Draw labels
  TLatex cmsLabel(leftMargin+0.03-0.013, 1-topMargin-0.03+0.005+0.08,"#font[62]{CMS}");
  cmsLabel.SetNDC();
  cmsLabel.SetTextAlign(13);
  cmsLabel.Draw();
  TLatex supplementaryLabel(leftMargin+0.03+0.07, 1-topMargin-0.03-0.005+0.08, "#scale[0.69]{#font[52]{Simulation Supplementary}}");
  supplementaryLabel.SetNDC();
  supplementaryLabel.SetTextAlign(13);
  supplementaryLabel.Draw();
  TLatex luminosityLabel(1-rightMargin-0.03+0.02, 1-topMargin-0.03+0.005+0.08, "137 fb^{-1} (13 TeV)");
  luminosityLabel.SetNDC();
  luminosityLabel.SetTextAlign(33);
  luminosityLabel.SetTextFont(42);
  luminosityLabel.SetTextSize(0.043);
  luminosityLabel.Draw();
  TLatex physicsLabel(1-rightMargin-0.03, 1-topMargin-0.03+0.005, "TChiHH-G(950,1)");
  physicsLabel.SetNDC();
  physicsLabel.SetTextAlign(33);
  physicsLabel.SetTextFont(42);
  physicsLabel.SetTextSize(0.043);
  physicsLabel.Draw();
  // Draw boxes
  TBox resolvedBox(0.4,0.4,2.2,2.2);
  resolvedBox.SetLineColor(kBlue);
  resolvedBox.SetFillStyle(0);
  resolvedBox.SetLineWidth(5);
  resolvedBox.Draw();
  TLatex resolvedLabel(0.8, 1.4, "Resolved");
  resolvedLabel.SetTextColor(kBlue);
  resolvedLabel.SetTextFont(42);
  resolvedLabel.SetTextSize(0.06);
  resolvedLabel.Draw();
  TBox boostedBox(0,0,0.8,0.8);
  boostedBox.SetLineColor(kRed);
  boostedBox.SetFillStyle(0);
  boostedBox.SetLineWidth(5);
  boostedBox.Draw();
  TLatex boostedLabel(0.07, 0.15, "Boosted");
  boostedLabel.SetTextColor(kRed);
  boostedLabel.SetTextFont(42);
  boostedLabel.SetTextSize(0.06);
  boostedLabel.Draw();
  TLine background_left1(0,2.3,0,3.2);
  background_left1.SetLineWidth(5);
  background_left1.Draw();
  TLine background_left2(2.3,0,2.3,2.3);
  background_left2.SetLineWidth(5);
  background_left2.Draw();
  TLine background_top(0,3.18,3.18,3.18);
  background_top.SetLineWidth(5);
  background_top.Draw();
  TLine background_right(3.18,3.18,3.19,0);
  background_right.SetLineWidth(5);
  background_right.Draw();
  TLine background_bottom1(2.3,0,3.19,0);
  background_bottom1.SetLineWidth(5);
  background_bottom1.Draw();
  TLine background_bottom2(0,2.3,2.3,2.3);
  background_bottom2.SetLineWidth(5);
  background_bottom2.Draw();
  TLatex backgroundLabel(1.3, 2.6, "#splitline{Background}{Dominated}");
  backgroundLabel.SetTextFont(42);
  backgroundLabel.SetTextSize(0.06);
  backgroundLabel.Draw();
  // Save canvas
  canvas.Print("plots/leadHDeltaR_vs_subleadHDeltaR.pdf");

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
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
      if(false){
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
