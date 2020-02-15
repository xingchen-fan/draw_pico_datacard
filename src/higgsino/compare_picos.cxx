#include <iostream>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby_pico.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
  double lumi = 137;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  NamedFunc filters = "1";//HigUtilities::pass_2016 && "stitch"; 

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("2016", Process::Type::background, kBlack, 
      {bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_higloose/*TTJets_*SingleLeptFromT_T*"}, filters));
  procs.push_back(Process::MakeShared<Baby_pico>("2017", Process::Type::background, kAzure+1, 
      {bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/2017/mc/merged_higmc_higloose/*TTJets_*SingleLeptFromT_T*"}, filters));
  procs.push_back(Process::MakeShared<Baby_pico>("2018", Process::Type::background, kRed+1, 
      {bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/2018/mc/merged_higmc_higloose/*TTJets_*SingleLeptFromT_T*"}, filters));

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info).YAxis(YAxisType::linear).Stack(StackType::lumi_shapes);
  //.Bottom(BottomType::ratio);//.YAxis(YAxisType::linear).Stack(StackType::lumi_shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> all_plot_types = {log_shapes};

  PlotMaker pm;

  string baseline = "1";
  vector<NamedFunc> weights = {"w_lumi"};
  for(unsigned iw=0; iw<weights.size(); iw++){
    pm.Push<Hist1D>(Axis(20, 0, 2000., "ht", "H_{T} [GeV]", {500.}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(14, 150, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njet", "N_{jets}", {6.5, 8.5}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nvlep", "N_{leps}", {5.5, 8.5}),
    //     "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(5, -0.5, 5.5, "nbm", "N_{b}", {1.5, 2.5,3.5}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
        "njet>=4", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
        "njet>=4", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20,0,200,"fjet_msoftdrop[0]", "M_{SD}(leading-p_{T} jet) [GeV]", {100., 140.}),
        "nfjet>=1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(41,-1,1,"fjet_mva_hbb_btv[0]", "DoubleB tagger", {}),
        "nfjet>=1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20,0,1,"fjet_deep_md_hbb_btv[0]", "Deep DoubleB tagger", {}),
        "nfjet>=1", procs, all_plot_types).Weight(weights[iw]);

    pm.Push<Hist1D>(Axis(40, 0., 2., "w_btag", "b-tag weight", {}),
       "nbm==0", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(40, 0., 2., "w_btag", "b-tag weight", {}),
       "nbm==1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_lep", "lepton weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_isr", "ISR weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_pu", "Pileup weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
  }
   
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
