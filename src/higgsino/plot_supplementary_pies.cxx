//This script generates the pie charts that may be used in the supplementary plots
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options (default indicates it is an AN plot): 

#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TLatex.h"
#include "TPad.h"
#include "TPie.h"
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
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  Palette colors("txt/colors.txt","default");
  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "plot_variables,plot_nminus1,make_cutflow,make_pies");

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
  PlotOpt lin_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt lin_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt log_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt log_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_info_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_paper_data = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::supplementary);

  vector<PlotOpt> plt_lin = {lin_norm_paper};
  vector<PlotOpt> plt_log = {log_norm_paper};

  set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  int lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string, lumi_precision);

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
  mctags["other_and_single_t"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root","*_ST_*.root"});
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

  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";

  std::vector<std::shared_ptr<Process>> search_procs;
  std::shared_ptr<Process> proc_ttbar = Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch");
  std::shared_ptr<Process> proc_zjets = Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch");
  std::shared_ptr<Process> proc_wjets = Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch");
  std::shared_ptr<Process> proc_qcd = Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch"); 
  std::shared_ptr<Process> proc_other = Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other_and_single_t"]),"stitch");

  search_procs.push_back(proc_ttbar);
  search_procs.push_back(proc_zjets);
  search_procs.push_back(proc_wjets);
  search_procs.push_back(proc_qcd);
  search_procs.push_back(proc_other);

  std::vector<std::string> proc_names = {"t#bar{t}+X","Z+jets","W+jets","QCD","Other"};
  std::vector<Int_t> proc_colors = {TColor::GetColor(102,102,255),kOrange+1,kGreen+1,TColor::GetColor(255,200,0),kGray+2};

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc filters = Higfuncs::final_pass_filters;
  NamedFunc weight = Higfuncs::final_weight;
  NamedFunc sr_baseline = filters&&script_utilities::search_resolved;
  NamedFunc metnbbin("metnbbin",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbm()<3) {
      if (b.met() < 300) return 1;
      else return 2;
    }
    else if (b.nbm()==3 && b.nbl()<4) {
      if (b.met() < 300) return 3;
      else return 4;
    }
    if (b.met() < 300) return 5;
    return 6;
  });
  NamedFunc drmaxnbbin("drmaxnbbin",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbm()<3) {
      if (b.hig_cand_drmax()->at(0) > 1.1) return 1;
      else return 2;
    }
    else if (b.nbm()==3 && b.nbl()<4) {
      if (b.hig_cand_drmax()->at(0) > 1.1) return 3;
      else return 4;
    }
    if (b.hig_cand_drmax()->at(0) > 1.1) return 5;
    return 6;
  });

  //------------------------------------------------------------------------------------
  //                                     make plots and pie charts
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(6, 0.5, 6.5, metnbbin, "MET-Nb category", {}),
      sr_baseline,
      search_procs, plt_lin).Weight(weight)
      .Tag("FixName:supppies__metnbcat")
      .LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(6, 0.5, 6.5, drmaxnbbin, "Drmax-Nb category", {}),
      sr_baseline,
      search_procs, plt_lin).Weight(weight)
      .Tag("FixName:supppies__drmaxnbcat")
      .LuminosityTag(total_luminosity_string);
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  TCanvas pie_canvas("pie_canvas","pie_canvas",700,1000);
  Hist1D* metnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("supppies__metnbcat").get());
  std::vector<TH1D*> metnb_proc_hist;
  Hist1D* drmaxnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("supppies__drmaxnbcat").get());
  std::vector<TH1D*> drmaxnb_proc_hist;
  for (auto proc : search_procs) {
    metnb_proc_hist.push_back(new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(proc.get())))->scaled_hist_));
    drmaxnb_proc_hist.push_back(new TH1D((static_cast<Hist1D::SingleHist1D*>(drmaxnb_hist1d->GetComponent(proc.get())))->scaled_hist_));
  }

  pie_canvas.cd();
  TPad cms_title_pad("cms_title_pad","cms_title_pad",0,0.9,0.571,1.0);
  cms_title_pad.Draw();
  cms_title_pad.cd();
  TLatex cms_label;
  cms_label.SetTextSize(0.3);
  cms_label.SetNDC(kTRUE);
  cms_label.SetTextAlign(22);
  cms_label.DrawLatex(0.5, 0.5,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}");

  pie_canvas.cd();
  TPad lowmet_2b_pad("lowmet_2b_pad","lowmet_2b_pad",0,0.7,0.286,0.9);
  lowmet_2b_pad.Draw();
  lowmet_2b_pad.cd();
  std::vector<double> lowmet_2b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowmet_2b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(1));
  }
  TPie lowmet_2b_pie("lowmet_2b_pie","2b Low p_{T}^{miss}",5,&lowmet_2b_proc_comp[0],&proc_colors[0]);
  lowmet_2b_pie.SetLabelFormat("%perc");
  lowmet_2b_pie.SetCircle(0.5, 0.48, 0.35);
  lowmet_2b_pie.Draw();

  pie_canvas.cd();
  TPad lowmet_3b_pad("lowmet_3b_pad","lowmet_3b_pad",0.286,0.7,0.571,0.9);
  lowmet_3b_pad.Draw();
  lowmet_3b_pad.cd();
  std::vector<double> lowmet_3b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowmet_3b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(3));
  }
  TPie lowmet_3b_pie("lowmet_3b_pie","3b Low p_{T}^{miss}",5,&lowmet_3b_proc_comp[0],&proc_colors[0]);
  lowmet_3b_pie.SetLabelFormat("%perc");
  lowmet_3b_pie.SetCircle(0.5, 0.48, 0.35);
  lowmet_3b_pie.Draw();

  pie_canvas.cd();
  TPad lowmet_4b_pad("lowmet_4b_pad","lowmet_4b_pad",0.571,0.7,0.857,0.9);
  lowmet_4b_pad.Draw();
  lowmet_4b_pad.cd();
  std::vector<double> lowmet_4b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowmet_4b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(5));
  }
  TPie lowmet_4b_pie("lowmet_4b_pie","4b Low p_{T}^{miss}",5,&lowmet_4b_proc_comp[0],&proc_colors[0]);
  lowmet_4b_pie.SetLabelFormat("%perc");
  lowmet_4b_pie.SetCircle(0.5, 0.48, 0.35);
  lowmet_4b_pie.Draw();

  pie_canvas.cd();
  TPad legend_pad("legend_pad","legend_pad",0.857,0.7,1.0,0.9);
  legend_pad.Draw();
  legend_pad.cd();
  TLegend leg(0.0, 0.0, 1.0, 1.0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    leg.AddEntry(metnb_proc_hist[proc_idx], proc_names[proc_idx].c_str(), "f");
  }
  leg.Draw();

  pie_canvas.cd();
  TPad highmet_2b_pad("highmet_2b_pad","highmet_2b_pad",0,0.5,0.286,0.7);
  highmet_2b_pad.Draw();
  highmet_2b_pad.cd();
  std::vector<double> highmet_2b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highmet_2b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(2));
  }
  TPie highmet_2b_pie("highmet_2b_pie","2b High p_{T}^{miss}",5,&highmet_2b_proc_comp[0],&proc_colors[0]);
  highmet_2b_pie.SetLabelFormat("%perc");
  highmet_2b_pie.SetCircle(0.5, 0.48, 0.35);
  highmet_2b_pie.Draw();

  pie_canvas.cd();
  TPad highmet_3b_pad("highmet_3b_pad","highmet_3b_pad",0.286,0.5,0.571,0.7);
  highmet_3b_pad.Draw();
  highmet_3b_pad.cd();
  std::vector<double> highmet_3b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highmet_3b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(4));
  }
  TPie highmet_3b_pie("highmet_3b_pie","3b High p_{T}^{miss}",5,&highmet_3b_proc_comp[0],&proc_colors[0]);
  highmet_3b_pie.SetLabelFormat("%perc");
  highmet_3b_pie.SetCircle(0.5, 0.48, 0.35);
  highmet_3b_pie.Draw();

  pie_canvas.cd();
  TPad highmet_4b_pad("highmet_4b_pad","highmet_4b_pad",0.571,0.5,0.857,0.7);
  highmet_4b_pad.Draw();
  highmet_4b_pad.cd();
  std::vector<double> highmet_4b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highmet_4b_proc_comp.push_back(metnb_proc_hist[proc_idx]->GetBinContent(6));
  }
  TPie highmet_4b_pie("highmet_4b_pie","4b High p_{T}^{miss}",5,&highmet_4b_proc_comp[0],&proc_colors[0]);
  highmet_4b_pie.SetLabelFormat("%perc");
  highmet_4b_pie.SetCircle(0.5, 0.48, 0.35);
  highmet_4b_pie.Draw();

  pie_canvas.cd();
  TPad highdrmax_2b_pad("highdrmax_2b_pad","highdrmax_2b_pad",0,0.2,0.286,0.4);
  highdrmax_2b_pad.Draw();
  highdrmax_2b_pad.cd();
  std::vector<double> highdrmax_2b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highdrmax_2b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(1));
  }
  TPie highdrmax_2b_pie("highdrmax_2b_pie","2b High #DeltaR_{max}",5,&highdrmax_2b_proc_comp[0],&proc_colors[0]);
  highdrmax_2b_pie.SetLabelFormat("%perc");
  highdrmax_2b_pie.SetCircle(0.5, 0.48, 0.35);
  highdrmax_2b_pie.Draw();

  pie_canvas.cd();
  TPad highdrmax_3b_pad("highdrmax_3b_pad","highdrmax_3b_pad",0.286,0.2,0.571,0.4);
  highdrmax_3b_pad.Draw();
  highdrmax_3b_pad.cd();
  std::vector<double> highdrmax_3b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highdrmax_3b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(3));
  }
  TPie highdrmax_3b_pie("highdrmax_3b_pie","3b High #DeltaR_{max}",5,&highdrmax_3b_proc_comp[0],&proc_colors[0]);
  highdrmax_3b_pie.SetLabelFormat("%perc");
  highdrmax_3b_pie.SetCircle(0.5, 0.48, 0.35);
  highdrmax_3b_pie.Draw();

  pie_canvas.cd();
  TPad highdrmax_4b_pad("highdrmax_4b_pad","highdrmax_4b_pad",0.571,0.2,0.857,0.4);
  highdrmax_4b_pad.Draw();
  highdrmax_4b_pad.cd();
  std::vector<double> highdrmax_4b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    highdrmax_4b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(5));
  }
  TPie highdrmax_4b_pie("highdrmax_4b_pie","4b High #DeltaR_{max}",5,&highdrmax_4b_proc_comp[0],&proc_colors[0]);
  highdrmax_4b_pie.SetLabelFormat("%perc");
  highdrmax_4b_pie.SetCircle(0.5, 0.48, 0.35);
  highdrmax_4b_pie.Draw();

  pie_canvas.cd();
  TPad lowdrmax_2b_pad("lowdrmax_2b_pad","lowdrmax_2b_pad",0,0.0,0.286,0.2);
  lowdrmax_2b_pad.Draw();
  lowdrmax_2b_pad.cd();
  std::vector<double> lowdrmax_2b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowdrmax_2b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(2));
  }
  TPie lowdrmax_2b_pie("lowdrmax_2b_pie","2b Low #DeltaR_{max}",5,&lowdrmax_2b_proc_comp[0],&proc_colors[0]);
  lowdrmax_2b_pie.SetLabelFormat("%perc");
  lowdrmax_2b_pie.SetCircle(0.5, 0.48, 0.35);
  lowdrmax_2b_pie.Draw();

  pie_canvas.cd();
  TPad lowdrmax_3b_pad("lowdrmax_3b_pad","lowdrmax_3b_pad",0.286,0.0,0.571,0.2);
  lowdrmax_3b_pad.Draw();
  lowdrmax_3b_pad.cd();
  std::vector<double> lowdrmax_3b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowdrmax_3b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(4));
  }
  TPie lowdrmax_3b_pie("lowdrmax_3b_pie","3b Low #DeltaR_{max}",5,&lowdrmax_3b_proc_comp[0],&proc_colors[0]);
  lowdrmax_3b_pie.SetLabelFormat("%perc");
  lowdrmax_3b_pie.SetCircle(0.5, 0.48, 0.35);
  lowdrmax_3b_pie.Draw();

  pie_canvas.cd();
  TPad lowdrmax_4b_pad("lowdrmax_4b_pad","lowdrmax_4b_pad",0.571,0.0,0.857,0.2);
  lowdrmax_4b_pad.Draw();
  lowdrmax_4b_pad.cd();
  std::vector<double> lowdrmax_4b_proc_comp;
  for (unsigned int proc_idx = 0; proc_idx < search_procs.size(); proc_idx++) {
    lowdrmax_4b_proc_comp.push_back(drmaxnb_proc_hist[proc_idx]->GetBinContent(6));
  }
  TPie lowdrmax_4b_pie("lowdrmax_4b_pie","4b Low #DeltaR_{max}",5,&lowdrmax_4b_proc_comp[0],&proc_colors[0]);
  lowdrmax_4b_pie.SetLabelFormat("%perc");
  lowdrmax_4b_pie.SetCircle(0.5, 0.48, 0.35);
  lowdrmax_4b_pie.Draw();

  pie_canvas.Print("plots/supp_res_pies.pdf");
  std::cout << "Opening plots/supp_res_pies.pdf\n";

  //for (TH1D* hist : metnb_proc_hist) {
  //  delete hist;
  //}
  //for (TH1D* hist : drmaxnb_proc_hist) {
  //  delete hist;
  //}

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

