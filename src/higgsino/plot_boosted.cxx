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
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = true;
}

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double wgt = 1;
  if (b.type()==106000) {
    return 137.;
  }
  if (b.SampleType()==2016){
    return wgt*35.9;
  } else if (b.SampleType()==2017){
    return wgt*41.5;
  } else {
    return wgt*59.6;
  }
});

const NamedFunc fjet1_dphi("fjet1_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  size_t ifjet = 0;
  if ((*b.fjet_phi()).size()<=ifjet) return 4;
  else return fabs(TVector2::Phi_mpi_pi((*b.fjet_phi())[ifjet]-b.met_phi()));
});
const NamedFunc fjet2_dphi("fjet2_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  size_t ifjet = 1;
  if ((*b.fjet_phi()).size()<=ifjet) return 4;
  else return fabs(TVector2::Phi_mpi_pi((*b.fjet_phi())[ifjet]-b.met_phi()));
});
const NamedFunc min_fjet_dphi("min_fjet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ifjet = 0; ifjet < (*b.fjet_phi()).size(); ++ifjet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.fjet_phi())[ifjet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});
const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  //vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    //if ((*b.jet_deepcsv())[ijet]<btag_wpts_2016[2]) continue;
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

const NamedFunc boosted_mass_region("boosted_mass_region", [](const Baby &b) -> NamedFunc::ScalarType{
  if ((*b.fjet_msoftdrop())[0]>95 && (*b.fjet_msoftdrop())[0]<=145 && (*b.fjet_msoftdrop())[1]>95 && (*b.fjet_msoftdrop())[1]<=145) return 1;
  else if ((*b.fjet_msoftdrop())[0]>50 && (*b.fjet_msoftdrop())[0]<=250 && (*b.fjet_msoftdrop())[1]>50 && (*b.fjet_msoftdrop())[1]<=250) return 0;
  else return -1;
});

const NamedFunc boosted_average_mass_region("boosted_average_mass_region", [](const Baby &b) -> NamedFunc::ScalarType{
  if (((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2>95 && ((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2<=145) return 1;
  else if (((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2>50 && ((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2<=250) return 0;
  else return -1;
});

const NamedFunc boosted_h1_mass_region("boosted_h1_mass_region", [](const Baby &b) -> NamedFunc::ScalarType{
  if ((*b.fjet_msoftdrop())[1]>95 && (*b.fjet_msoftdrop())[1]<=145) return 1;
  else if ((*b.fjet_msoftdrop())[1]>50 && (*b.fjet_msoftdrop())[1]<=250) return 0;
  else return -1;
});

const NamedFunc dmsoftdrop("dmsoftdrop", [](const Baby &b) -> NamedFunc::ScalarType{
  return fabs((*b.fjet_msoftdrop())[0]-(*b.fjet_msoftdrop())[1]);
});

const NamedFunc boo_0htag("boo_0htag", [](const Baby &b) -> NamedFunc::ScalarType{
  return (*b.fjet_deep_md_hbb_btv())[0]<=0.7&&(*b.fjet_deep_md_hbb_btv())[1]<=0.7;
});
const NamedFunc boo_1htag("boo_1htag", [](const Baby &b) -> NamedFunc::ScalarType{
  return ((*b.fjet_deep_md_hbb_btv())[0]>0.7&&(*b.fjet_deep_md_hbb_btv())[1]<=0.7)||((*b.fjet_deep_md_hbb_btv())[1]>0.7&&(*b.fjet_deep_md_hbb_btv())[0]<=0.7);
});
const NamedFunc boo_2htag("boo_2htag", [](const Baby &b) -> NamedFunc::ScalarType{
  return ((*b.fjet_deep_md_hbb_btv())[0]>0.7&&(*b.fjet_deep_md_hbb_btv())[1]>0.7);
});

//const NamedFunc boo_0htag("boo_0htag", [](const Baby &b) -> NamedFunc::ScalarType{
//  return (*b.fjet_mva_hbb_btv())[0]<=0.3&&(*b.fjet_mva_hbb_btv())[1]<=0.3;
//});
//const NamedFunc boo_1htag("boo_1htag", [](const Baby &b) -> NamedFunc::ScalarType{
//  return ((*b.fjet_mva_hbb_btv())[0]>0.3&&(*b.fjet_mva_hbb_btv())[1]<=0.3)||((*b.fjet_mva_hbb_btv())[1]>0.3&&(*b.fjet_mva_hbb_btv())[0]<=0.3);
//});
//const NamedFunc boo_2htag("boo_2htag", [](const Baby &b) -> NamedFunc::ScalarType{
//  return ((*b.fjet_mva_hbb_btv())[0]>0.3&&(*b.fjet_mva_hbb_btv())[1]>0.3);
//});

const NamedFunc boo_mass("boo_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  return (*b.fjet_msoftdrop())[0]>50&&(*b.fjet_msoftdrop())[0]<=250&&(*b.fjet_msoftdrop())[1]>50&&(*b.fjet_msoftdrop())[1]<=250;
  //return (*b.fjet_msoftdrop())[0]>1;
  //return ((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2>50&&((*b.fjet_msoftdrop())[0]+(*b.fjet_msoftdrop())[1])/2<=250;
});

// Combine the 0H, 1H, 2H plots
void add_fjet_2mass_plots(PlotMaker & pm, NamedFunc common_cuts, NamedFunc wgt, string sample_name, string base_folder, set<int> years, string mc_folder, set<string> mctag, vector<PlotOpt> & plt_shapes) {
  vector<shared_ptr<Process> > procs_htag;
  vector<int> colors = {kGreen+1, kRed, kBlue, kOrange};
  procs_htag.push_back(Process::MakeShared<Baby_pico>(sample_name+" 0h", Process::Type::background,colors[0],
                  attach_folder(base_folder, years, mc_folder, mctag),
                  common_cuts&&boo_0htag));
  procs_htag.push_back(Process::MakeShared<Baby_pico>(sample_name+" 1h", Process::Type::background,colors[1],
                  attach_folder(base_folder, years, mc_folder, mctag),
                  common_cuts&&boo_1htag));
  procs_htag.push_back(Process::MakeShared<Baby_pico>(sample_name+" 2h", Process::Type::background,colors[2],
                  attach_folder(base_folder, years, mc_folder, mctag),
                  common_cuts&&boo_2htag));

  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    1, procs_htag, plt_shapes).Weight(wgt);
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    1, procs_htag, plt_shapes).Weight(wgt);
  pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    1, procs_htag, plt_shapes).Weight(wgt);
}

void add_fjet_2d_plots(PlotMaker & pm, NamedFunc common_cuts, NamedFunc wgt, vector<shared_ptr<Process> > & procs, vector<PlotOpt> & plt_shapes_info, vector<PlotOpt> & plt_2D) {
  // lead mass 0H, 1H, 2H
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
  // sublead mass 0H, 1H, 2H
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
  // Sum for ABCD 0H, 1H, 2H
  pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    common_cuts&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    common_cuts&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
  pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    common_cuts&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
  // 2D lead-sublead mass 0H, 1H, 2H
  pm.Push<Hist2D>(Axis(15, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_2htag, procs, plt_2D).Weight(wgt);
  pm.Push<Hist2D>(Axis(15, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_1htag, procs, plt_2D).Weight(wgt);
  pm.Push<Hist2D>(Axis(15, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    common_cuts&&boo_0htag, procs, plt_2D).Weight(wgt);
}

const NamedFunc hcut_0mid("hcut_0mid", [](const Baby &b) -> NamedFunc::ScalarType{
  //return (*b.fjet_deep_md_hbb_btv())[0]>0.7 &&(*b.fjet_deep_md_hbb_btv())[0]<0.86 && (*b.fjet_deep_md_hbb_btv())[1]>0.7 && (*b.fjet_deep_md_hbb_btv())[1]<0.86;
  //return ((*b.fjet_deep_md_hbb_btv())[0]>0.7 &&(*b.fjet_deep_md_hbb_btv())[0]<0.86 && (*b.fjet_deep_md_hbb_btv())[1]<0.86) || ((*b.fjet_deep_md_hbb_btv())[1]>0.7 &&(*b.fjet_deep_md_hbb_btv())[1]<0.86 && (*b.fjet_deep_md_hbb_btv())[0]<0.86);
  return (*b.fjet_deep_md_hbb_btv())[0]>0.5 &&(*b.fjet_deep_md_hbb_btv())[0]<0.86 && (*b.fjet_deep_md_hbb_btv())[1]>0.5 && (*b.fjet_deep_md_hbb_btv())[1]<0.86;
});
const NamedFunc hcut_1mid("hcut_1mid", [](const Baby &b) -> NamedFunc::ScalarType{
  return (((*b.fjet_deep_md_hbb_btv())[0]>0.86 && (*b.fjet_deep_md_hbb_btv())[1]>0.7 && (*b.fjet_deep_md_hbb_btv())[1]<0.86) || ((*b.fjet_deep_md_hbb_btv())[1]>0.86 && (*b.fjet_deep_md_hbb_btv())[0]>0.7 && (*b.fjet_deep_md_hbb_btv())[0]<0.86));
});
const NamedFunc hcut_2mid("hcut_2mid", [](const Baby &b) -> NamedFunc::ScalarType{
  return ((*b.fjet_deep_md_hbb_btv())[0]>0.86  && (*b.fjet_deep_md_hbb_btv())[1]>0.86);
});

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  string hostname = execute("echo $HOSTNAME");
  string bfolder("");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-") || Contains(hostname, "physics.ucsb"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  set<int> years;
  years = {2016, 2017, 2018};
  //years = {2016};
  //years = {2017};
  //years = {2018};

  string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string sig_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/";
  string sig_skim_folder = "TChiHH/merged_higmc_unskimmed/";

  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                  "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"
                                 });
  mctags["ttonly"]     = set<string>({"*TTJets_*Lept*"});
  //mctags["topx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
  //                                   "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mctags["all"] = set<string>({
                               "*TTJets_*Lept*",
                               "*TTJets_DiLept_Tune*", "*TTJets_HT*", "*TTJets_SingleLeptFromTbar_Tune*", "*TJets_SingleLeptFromT_Tune*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  //mctags["all"] = mctags["wjets"];
  //mctags["all"] = set<string>({
  //                             "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
  //                             "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
  //                             "*_QCD_HT2000toInf_*",
  //                             "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
  //                             "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  //});

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("tt", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),base_filters&&"stitch")); 


  //vector<string> sigm = {}; 
  vector<string> sigm = {"450", "700", "950"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  // extended list
  // vector<string> sigm = {"127","300","400","500","600","700","850","1000"}; 
  // vector<int> sig_colors = {kAzure, kGreen+2, kRed, kViolet-6, kYellow, kMagenta+1, kCyan+1, kOrange-3};
  for (unsigned isig(0); isig<sigm.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
      sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_*.root"}), base_filters));
  }

  vector<shared_ptr<Process> > procs_onlytt;
  procs_onlytt.push_back(Process::MakeShared<Baby_pico>("ttonly", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),base_filters&&"stitch"));
  //sigm = {"900"}; 
  //for (unsigned isig(0); isig<sigm.size(); isig++){
  //  procs_onlytt.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
  //    sig_colors[isig], {foldersig+"*TChiHH_mChi-"+sigm[isig]+"_*.root"}, base_filters));
  //}

  vector<shared_ptr<Process> > procs_onlyz;
  procs_onlyz.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,  mctags["zjets"]),base_filters&&"stitch"));

  vector<shared_ptr<Process> > procs_onlyw;
  procs_onlyw.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["wjets"]),base_filters&&"stitch"));

  vector<shared_ptr<Process> > procs_tzw;
  procs_tzw.push_back(Process::MakeShared<Baby_pico>("tt", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),base_filters&&"stitch"));
  procs_tzw.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),base_filters&&"stitch"));
  procs_tzw.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),base_filters&&"stitch"));


  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    //.Stack(StackType::shapes).LegendColumns(3);
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
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

  //double luminosity = 137;

  NamedFunc wgt = "w_lumi*w_isr"*w_years*Higfuncs::eff_higtrig;

  //string common = "!lowDphiFix && nvlep==0 && ntk==0 && met>300" ;
  //string common = "!lowDphiFix && nvlep==0 && ntk==0";
  string common = "!lowDphiFix && nvlep==0 && ntk == 0 && met>300" ;

  //NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600 && nbm>=2";
  //NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600 "&&dmsoftdrop<40;
  //NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600 && nbm>=2"&&dmsoftdrop<40;
  //NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600 ";
  //string common_notrk = "!lowDphiFix && nvlep==0";
  //NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600 && nbm>=2"&&dmsoftdrop<40;
  NamedFunc boo_base = "nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600"&&dmsoftdrop<40;

  //string boo_mass = "fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  //string boo_2htag = "fjet_deep_md_hbb_btv[0]>0.7 && fjet_deep_md_hbb_btv[1]>0.7";
  //string boo_1htag = "(fjet_deep_md_hbb_btv[0]>0.7 && fjet_deep_md_hbb_btv[1]<=0.7)||(fjet_deep_md_hbb_btv[0]<=0.7 && fjet_deep_md_hbb_btv[1]>0.7)";
  //string boo_0htag = "fjet_deep_md_hbb_btv[0]<=0.7 && fjet_deep_md_hbb_btv[1]<=0.7";
  //// boostedRegionIdx: 0 = D, 1 = B1, 2 = B2, 3 = C, 4 = A1, 5 = A2
  //string boo_SR = "boostedRegionIdx==5";

  bool plot_fjet_vs_htag = false;
  bool plot_hmass_vs_htag = false;
  bool plot_h2mass_vs_htag = false;
  bool plot_htag_relations = true;

  PlotMaker pm;

  //pm.Push<Hist1D>(Axis(20, 150, 1000, "met", "met [GeV]", {}),
  //  common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);
  //pm.Push<Hist1D>(Axis(20, 0, 3.1416,fjet1_dphi, "fjet[0] dphi", {}),
  //  common&&boo_base&&boo_mass&&fjet1_dphi!=4, procs, plt_lin).Weight(wgt);
  //pm.Push<Hist1D>(Axis(20, 0, 3.1416,fjet2_dphi, "fjet[1] dphi", {}),
  //  common&&boo_base&&boo_mass&&fjet2_dphi!=4, procs, plt_lin).Weight(wgt);
  //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_fjet_dphi, "min fjet dphi", {}),
  //  common&&boo_base&&boo_mass&&min_fjet_dphi!=4, procs, plt_lin).Weight(wgt);
  //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
  //  common&&boo_base&&boo_mass&&min_jet_dphi!=4, procs, plt_lin).Weight(wgt);


  if (plot_fjet_vs_htag) {
    ////// Against met
    //pm.Push<Hist1D>(Axis({300, 500, 700, 1000},"met", "p_{T}^{miss} [GeV]", {}),
    //  common&&boo_base&&boo_mass&&"boostedRegionIdx==4"&&"met>300", procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis({300, 500, 700, 1000},"met", "p_{T}^{miss} [GeV]", {}),
    //  common&&boo_base&&boo_mass&&"boostedRegionIdx==5"&&"met>300", procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis({300, 500, 700, 1000},"met", "p_{T}^{miss} [GeV]", {}),
      common&&boo_base&&boo_mass&&"met>300", procs, plt_lin).Weight(wgt);

    add_fjet_2d_plots(pm, common&&boo_base&&boo_mass, wgt, procs, plt_lin, plt_2D);
    //add_fjet_2d_plots(pm, common&&boo_base&&boo_mass, wgt, procs_onlytt, plt_shapes_info, plt_2D);
    //add_fjet_2d_plots(pm, common&&boo_base&&boo_mass, wgt, procs_onlyz, plt_shapes_info, plt_2D);
    //add_fjet_2d_plots(pm, common&&boo_base&&boo_mass, wgt, procs_onlyw, plt_shapes_info, plt_2D);

     add_fjet_2mass_plots(pm, base_filters&&"stitch"&&common&&boo_base&&boo_mass, wgt, "mc", mc_base_folder, years, mc_skim_folder, mctags["all"], plt_shapes);
     //add_fjet_2mass_plots(pm, base_filters&&"stitch"&&common&&boo_base&&boo_mass, wgt, "mc", foldermc, mctags["all"], plt_shapes);
     //add_fjet_2mass_plots(pm, base_filters&&"stitch"&&common&&boo_base&&boo_mass, wgt, "ttonly", foldermc, mctags["ttonly"], plt_shapes);

    //// All with met>300
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 50, 250,"fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist1D>(Axis({0,1,2}, boosted_mass_region, "Boosted mass region", {}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_shapes_info).Weight(wgt);
    //pm.Push<Hist2D>(Axis(10, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_2D).Weight(wgt);
    //pm.Push<Hist2D>(Axis(10, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_2D).Weight(wgt);
    //pm.Push<Hist2D>(Axis(10, 50, 250, "fjet_msoftdrop[0]", "fjet[0] softdrop [GeV]", {95, 145}),
    //  Axis(10, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95, 145}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_2D).Weight(wgt);

  }

  if (plot_hmass_vs_htag) {

    // Show where to do dmsoftdrop mass cut
    pm.Push<Hist1D>(Axis(8, 0, 80, dmsoftdrop, "Delta fjet mass [GeV]", {40}),
      common&&boo_base&&boo_mass&&boo_2htag, procs, plt_lin).Weight(wgt);

    // Show full distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);

    // 0H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_0htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_lin).Weight(wgt);

    // 1H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_1htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_lin).Weight(wgt);

    // 2H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_2htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_lin).Weight(wgt);

    // 0H, 1H, 2H comparison
    vector<shared_ptr<Process> > procs_htag;
    procs_htag.push_back(Process::MakeShared<Baby_pico>("0h", Process::Type::background,kGreen+1,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_0htag));
    procs_htag.push_back(Process::MakeShared<Baby_pico>("1h", Process::Type::background,kRed,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_1htag));
    procs_htag.push_back(Process::MakeShared<Baby_pico>("2h", Process::Type::background,kBlue,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_2htag));
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95,145}),
      1, procs_htag, plt_shapes).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {40}),
    //  1, procs_htag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis({0,1,2}, boosted_average_mass_region, "Boosted average mass region", {}),
      1, procs_htag, plt_shapes).Weight(wgt);

    //vector<shared_ptr<Process> > procs_htag_tt;
    //procs_htag_tt.push_back(Process::MakeShared<Baby_pico>("tt 0h", Process::Type::background,kGreen+1,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_0htag));
    //procs_htag_tt.push_back(Process::MakeShared<Baby_pico>("tt 1h", Process::Type::background,kRed,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_1htag));
    //procs_htag_tt.push_back(Process::MakeShared<Baby_pico>("tt 2h", Process::Type::background,kBlue,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_2htag));
    //pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {85,135}),
    //  1, procs_htag_tt, plt_shapes).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {40}),
    //  1, procs_htag_tt, plt_shapes).Weight(wgt);

    //vector<shared_ptr<Process> > procs_htag_wjets;
    //procs_htag_wjets.push_back(Process::MakeShared<Baby_pico>("W+Jets 0h", Process::Type::background,kGreen+1,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["wjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_0htag));
    //procs_htag_wjets.push_back(Process::MakeShared<Baby_pico>("W+Jets 1h", Process::Type::background,kRed,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["wjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_1htag));
    //procs_htag_wjets.push_back(Process::MakeShared<Baby_pico>("W+Jets 2h", Process::Type::background,kBlue,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["wjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_2htag));
    //pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {85,135}),
    //  1, procs_htag_wjets, plt_shapes).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {40}),
    //  1, procs_htag_wjets, plt_shapes).Weight(wgt);

    //vector<shared_ptr<Process> > procs_htag_zjets;
    //procs_htag_zjets.push_back(Process::MakeShared<Baby_pico>("Z+Jets 0h", Process::Type::background,kGreen+1,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["zjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_0htag));
    //procs_htag_zjets.push_back(Process::MakeShared<Baby_pico>("Z+Jets 1h", Process::Type::background,kRed,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["zjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_1htag));
    //procs_htag_zjets.push_back(Process::MakeShared<Baby_pico>("Z+Jets 2h", Process::Type::background,kBlue,
    //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["zjets"]),
    //                base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_2htag));
    //pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {85,135}),
    //  1, procs_htag_zjets, plt_shapes).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {40}),
    //  1, procs_htag_zjets, plt_shapes).Weight(wgt);


  }

  if (plot_h2mass_vs_htag) {

    // Show full distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "fjet1 mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);

    // 0H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "fjet1 mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_0htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_0htag, procs, plt_lin).Weight(wgt);

    // 1H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "fjet1 mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_1htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_1htag, procs, plt_lin).Weight(wgt);

    // 2H distribution
    pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "fjet1 mass [GeV]", {95, 145}),
      common&&boo_base&&boo_mass&&boo_2htag, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(8, 0, 40, dmsoftdrop, "Delta fjet mass [GeV]", {}),
    //  common&&boo_base&&boo_mass&&boo_2htag, procs, plt_lin).Weight(wgt);

    // 0H, 1H, 2H comparison
    vector<shared_ptr<Process> > procs_htag;
    procs_htag.push_back(Process::MakeShared<Baby_pico>("0h", Process::Type::background,kGreen+1,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_0htag));
    procs_htag.push_back(Process::MakeShared<Baby_pico>("1h", Process::Type::background,kRed,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_1htag));
    procs_htag.push_back(Process::MakeShared<Baby_pico>("2h", Process::Type::background,kBlue,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&boo_2htag));
    pm.Push<Hist1D>(Axis(15, 50, 250, "fjet_msoftdrop[1]", "fjet[1] softdrop [GeV]", {95,145}),
      1, procs_htag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis({0,1,2}, boosted_h1_mass_region, "Boosted fjet1 mass region", {}),
      1, procs_htag, plt_shapes).Weight(wgt);
  }

  if (plot_htag_relations) {
    // H tag relation
    pm.Push<Hist2D>(Axis(20, 0, 1, "fjet_deep_md_hbb_btv[0]", "leading fjet Htag", {0.7, 0.86, 0.91}),
      Axis(20, 0, 1, "fjet_deep_md_hbb_btv[1]", "subleading fjet Htag", {0.7, 0.86, 0.91}),
      common&&boo_base&&boo_mass, procs, plt_lin).Weight(wgt);

   //// 1 loose H, 0 mid H
   //string hcut_0mid = "fjet_deep_md_hbb_btv[0]>0.7 &&fjet_deep_md_hbb_btv[0]<0.86 && fjet_deep_md_hbb_btv[1]>0.7 &&  fjet_deep_md_hbb_btv[1]>0.86";
   //// 2 loose H, 1  mid H
   //string hcut_1mid = "((fjet_deep_md_hbb_btv[0]>0.86 && fjet_deep_md_hbb_btv[1]>0.7 &&  fjet_deep_md_hbb_btv[1]>0.86) || (fjet_deep_md_hbb_btv[1]>0.86 && fjet_deep_md_hbb_btv[0]>0.7 &&  fjet_deep_md_hbb_btv[0]>0.86))";
   //// 2 mid H
   //string hcut_2mid = "(fjet_deep_md_hbb_btv[0]>0.86  && fjet_deep_md_hbb_btv[1]>0.86)";

   pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
     common&&boo_base&&boo_mass&&hcut_0mid, procs, plt_lin).Weight(wgt);
   pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
     common&&boo_base&&boo_mass&&hcut_1mid, procs, plt_lin).Weight(wgt);
   pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95, 145}),
     common&&boo_base&&boo_mass&&hcut_2mid, procs, plt_lin).Weight(wgt);

    // 0mid, 1mid, 2mid comparison
    vector<shared_ptr<Process> > procs_hcut;
    procs_hcut.push_back(Process::MakeShared<Baby_pico>("0 mid h", Process::Type::background,kGreen+1,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&hcut_0mid));
    procs_hcut.push_back(Process::MakeShared<Baby_pico>("1 mid h", Process::Type::background,kRed,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&hcut_1mid));
    procs_hcut.push_back(Process::MakeShared<Baby_pico>("2 mid h", Process::Type::background,kBlue,
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&common&&boo_base&&boo_mass&&hcut_2mid));
    pm.Push<Hist1D>(Axis(15, 50, 250, "(fjet_msoftdrop[0]+fjet_msoftdrop[1])/2", "Average fjet mass [GeV]", {95,145}),
      1, procs_hcut, plt_shapes).Weight(wgt);
  }


  pm.multithreaded_ = true;
  pm.min_print_ = true;
  //pm.MakePlots(luminosity);
  pm.MakePlots(1);

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
