#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  // signal points to include and their colors
  vector<string> sigm = {"225","400"}; 
  vector<int> sig_colors = {kCyan, kGreen}; // need sigm.size() >= sig_colors.size()
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::log).Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  vector<PlotOpt> linplot = {lin_lumi};
  vector<PlotOpt> logplot = {log_lumi};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  set<int> years; 
  years = {2016, 2017, 2018};

  string base_dir(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/");
  string mc_skim_dir("mc/merged_higmc_met150/");
  string sig_skim_dir("SMS-TChiHH_2D/skim_met150/");

  vector<string> mctags = {
    "*QCD_HT200to300*",  "*QCD_HT300to500*", 
    "*QCD_HT500to700*",  
    "*QCD_HT700to1000*", 
    "*QCD_HT1000to1500*", "*QCD_HT1500to2000*", "*QCD_HT2000toInf*"
  };
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baseline cuts ///////////////////////////////////////
  string baseline_s = "pass && stitch && nvlep==0 && ntk==0 && met>150 && njet>=4 && njet<=5 && nbt>=2";
  // baseline_s += " && hig_cand_dm[0]<=40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<=2.2";
  NamedFunc baseline = baseline_s;
  NamedFunc wgt = "w_lumi"*HigUtilities::w_years;

  vector<int> colors = {kRed+3, kRed+1, kOrange-3, kOrange, kSpring+10, kGreen+1, kAzure+1};
  vector<shared_ptr<Process> > procs;
  for (unsigned i(0); i<mctags.size(); i++)
    procs.push_back(Process::MakeShared<Baby_pico>("QCD "+mctags[i].substr(mctags[i].find("HT")+2, mctags[i].length()-8),        
    Process::Type::background, colors[i],    
    attach_folder(base_dir, years, mc_skim_dir, {mctags[i]}),  baseline));

  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1) x10", Process::Type::signal, 
      sig_colors[isig], attach_folder(base_dir, years, sig_skim_dir, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), baseline));

  NamedFunc v_dphi("v_dphi",[](const Baby &b) -> NamedFunc::VectorType{ 
    vector<double> v_dphi_;
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (fabs(b.jet_eta()->at(i))<=2.4 && !b.jet_islep()->at(i)) 
        v_dphi_.push_back(b.jet_met_dphi()->at(i));
    }
    return v_dphi_;
  });

  PlotMaker pm;


  pm.Push<Hist1D>(Axis(30,0.,3.,v_dphi[0.], "#Delta#phi(j_{0}, MET)",{0.5}), baseline_s, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(30,0.,3.,v_dphi[1], "#Delta#phi(j_{1}, MET)",{0.5}), baseline_s && v_dphi[0.]>0.5, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(30,0.,3.,v_dphi[2], "#Delta#phi(j_{2}, MET)",{0.3}), baseline_s && v_dphi[0.]>0.5 && v_dphi[1]>0.5, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(30,0.,3.,v_dphi[3], "#Delta#phi(j_{3}, MET)",{0.3}), baseline_s && v_dphi[0.]>0.5 && v_dphi[1]>0.5 && v_dphi[2]>0.3, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(30,0.,3.,v_dphi[4], "#Delta#phi(j_{4}, MET)",{}),baseline_s && v_dphi[0.]>0.5 && v_dphi[1]>0.5 && v_dphi[2]>0.3 && v_dphi[3]>0.3 && "njet>4", procs, linplot).Weight(wgt);
  
  pm.Push<Hist1D>(Axis(40,-2.,2.,"(met-mht)/met", "(MET-MHT)/MET",{-0.5,0.5}), "1", procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(40,-2.,2.,"(met-mht)/met", "(MET-MHT)/MET",{-0.5,0.5}), "hig_cand_dm[0]<=40 && hig_cand_am[0]<200 && hig_cand_drmax[0]<=2.2", procs, linplot).Weight(wgt);

  pm.Push<Hist1D>(Axis(60,0.,6.,"met/mht", "p^{miss}_{T}/H^{miss}_{T}",{2.}), baseline_s, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(60,0.,6.,"met/mht", "p^{miss}_{T}/H^{miss}_{T}",{2.}), baseline_s+"&& !low_dphi_met", procs, linplot).Weight(wgt);

  pm.Push<Hist1D>(Axis(60,0.,6.,"met/met_calo", "p^{miss}_{T}/p^{miss}_{T,calo}",{2., 5}), baseline_s, procs, linplot).Weight(wgt);
  pm.Push<Hist1D>(Axis(60,0.,6.,"met/met_calo", "p^{miss}_{T}/p^{miss}_{T,calo}",{2., 5}), baseline_s+"&& !low_dphi_met", procs, linplot).Weight(wgt);

  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
