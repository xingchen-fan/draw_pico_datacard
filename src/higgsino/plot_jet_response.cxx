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
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TGraphErrors.h"

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
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  // Options: 
  string string_options = "";
}

Double_t DeltaR(const ROOT::Math::PtEtaPhiMVector & v1, const ROOT::Math::PtEtaPhiMVector & v2)
{
   Double_t deta = v1.Eta()-v2.Eta();
   Double_t dphi = TVector2::Phi_mpi_pi(v1.Phi()-v2.Phi());
   return TMath::Sqrt( deta*deta+dphi*dphi );
}

// Get matched rec. jet pt
const NamedFunc all_matched_jet_pt("all_matched_jet_pt", [](const Baby &b) -> NamedFunc::VectorType{
  vector<int> genjet_indices(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    genjet_indices[iJet] = (*b.jet_genjet_idx())[iJet];
  }

  vector<double> jet_pt;
  jet_pt.reserve(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_pt.push_back((*b.jet_pt())[iJet]);
  }
  return jet_pt;
});
// Get matched gen. jet pt
const NamedFunc all_matched_genjet_pt("all_matched_genjet_pt", [](const Baby &b) -> NamedFunc::VectorType{
  vector<int> genjet_indices(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    genjet_indices[iJet] = (*b.jet_genjet_idx())[iJet];
  }

  vector<double> jet_pt;
  jet_pt.reserve(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_pt.push_back((*b.genjet_pt())[genjet_indices.at(iJet)]);
  }
  return jet_pt;
});
// Get all matched (rec jet pt) / (gen. jets pt)
const NamedFunc all_matched_response("all_matched_response", [](const Baby &b) -> NamedFunc::VectorType{
  vector<int> genjet_indices(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    genjet_indices[iJet] = (*b.jet_genjet_idx())[iJet];
  }

  vector<double> jet_response;
  jet_response.reserve(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_response.push_back((*b.jet_pt())[iJet] / (*b.genjet_pt())[genjet_indices.at(iJet)] );
  }
  return jet_response;
});
// Get all matched deltaR(rec jet, gen. jet)
const NamedFunc all_matched_deltaR("all_matched_deltaR", [](const Baby &b) -> NamedFunc::VectorType{
  vector<int> genjet_indices(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    genjet_indices[iJet] = (*b.jet_genjet_idx())[iJet];
  }

  vector<double> deltaR;
  deltaR.reserve(b.njet());
  for (unsigned iJet = 0; iJet < unsigned(b.njet()); ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    ROOT::Math::PtEtaPhiMVector recjet ((*b.jet_pt())[iJet], (*b.jet_eta())[iJet], (*b.jet_phi())[iJet], (*b.jet_m())[iJet]);
    ROOT::Math::PtEtaPhiMVector genjet ((*b.genjet_pt())[genjet_indices.at(iJet)], (*b.genjet_eta())[genjet_indices.at(iJet)], (*b.genjet_phi())[genjet_indices.at(iJet)], (*b.genjet_m())[genjet_indices.at(iJet)]);
    deltaR.push_back(DeltaR(recjet, genjet));
  }
  return deltaR;
});

// Get matched rec. jet pt
const NamedFunc matched_jet_pt("matched_jet_pt", [](const Baby &b) -> NamedFunc::VectorType{
  vector<unsigned> bbjet_indices = Higfuncs::get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  vector<int> genjet_indices(4);
  genjet_indices[0] = (*b.jet_genjet_idx())[bbjet_indices.at(0)];
  genjet_indices[1] = (*b.jet_genjet_idx())[bbjet_indices.at(1)];
  genjet_indices[2] = (*b.jet_genjet_idx())[bbjet_indices.at(2)];
  genjet_indices[3] = (*b.jet_genjet_idx())[bbjet_indices.at(3)];

  vector<double> jet_pt;
  jet_pt.reserve(4);
  for (unsigned iJet = 0; iJet < 4; ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_pt.push_back((*b.jet_pt())[bbjet_indices.at(iJet)]);
  }
  return jet_pt;
});
// Get matched gen. jet pt
const NamedFunc matched_genjet_pt("matched_genjet_pt", [](const Baby &b) -> NamedFunc::VectorType{
  vector<unsigned> bbjet_indices = Higfuncs::get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  vector<int> genjet_indices(4);
  genjet_indices[0] = (*b.jet_genjet_idx())[bbjet_indices.at(0)];
  genjet_indices[1] = (*b.jet_genjet_idx())[bbjet_indices.at(1)];
  genjet_indices[2] = (*b.jet_genjet_idx())[bbjet_indices.at(2)];
  genjet_indices[3] = (*b.jet_genjet_idx())[bbjet_indices.at(3)];

  vector<double> jet_pt;
  jet_pt.reserve(4);
  for (unsigned iJet = 0; iJet < 4; ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_pt.push_back((*b.genjet_pt())[genjet_indices.at(iJet)]);
  }
  return jet_pt;
});
// Get matched (rec jet pt) / (gen. jets pt)
const NamedFunc matched_response("matched_response", [](const Baby &b) -> NamedFunc::VectorType{
  vector<unsigned> bbjet_indices = Higfuncs::get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  vector<int> genjet_indices(4);
  genjet_indices[0] = (*b.jet_genjet_idx())[bbjet_indices.at(0)];
  genjet_indices[1] = (*b.jet_genjet_idx())[bbjet_indices.at(1)];
  genjet_indices[2] = (*b.jet_genjet_idx())[bbjet_indices.at(2)];
  genjet_indices[3] = (*b.jet_genjet_idx())[bbjet_indices.at(3)];

  vector<double> jet_response;
  jet_response.reserve(4);
  for (unsigned iJet = 0; iJet < 4; ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    jet_response.push_back((*b.jet_pt())[bbjet_indices.at(iJet)] / (*b.genjet_pt())[genjet_indices.at(iJet)] );
  }
  return jet_response;
});

// Get matched deltaR(rec jet, gen. jet)
const NamedFunc matched_deltaR("matched_deltaR", [](const Baby &b) -> NamedFunc::VectorType{
  vector<unsigned> bbjet_indices = Higfuncs::get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi(), *b.jet_isgood());
  vector<int> genjet_indices(4);
  genjet_indices[0] = (*b.jet_genjet_idx())[bbjet_indices.at(0)];
  genjet_indices[1] = (*b.jet_genjet_idx())[bbjet_indices.at(1)];
  genjet_indices[2] = (*b.jet_genjet_idx())[bbjet_indices.at(2)];
  genjet_indices[3] = (*b.jet_genjet_idx())[bbjet_indices.at(3)];

  vector<double> deltaR;
  deltaR.reserve(4);
  for (unsigned iJet = 0; iJet < 4; ++iJet) {
    if (genjet_indices[iJet] == -1) continue;
    ROOT::Math::PtEtaPhiMVector recjet ((*b.jet_pt())[bbjet_indices.at(iJet)], (*b.jet_eta())[bbjet_indices.at(iJet)], (*b.jet_phi())[bbjet_indices.at(iJet)], (*b.jet_m())[bbjet_indices.at(iJet)]);
    ROOT::Math::PtEtaPhiMVector genjet ((*b.genjet_pt())[genjet_indices.at(iJet)], (*b.genjet_eta())[genjet_indices.at(iJet)], (*b.genjet_phi())[genjet_indices.at(iJet)], (*b.genjet_m())[genjet_indices.at(iJet)]);
    deltaR.push_back(DeltaR(recjet, genjet));
  }
  return deltaR;
});

// Collect histograms per process from plotIndex
void get_histograms_2d(PlotMaker & pm, int plotIndex, map<string, TH2D *> & histograms_2d) {
  Hist2D * hist2d = static_cast<Hist2D*> (pm.Figures()[plotIndex].get());
  for (auto & process : hist2d->GetProcesses()) {
    //cout<<process->name_<<endl;
    Hist2D::SingleHist2D * singleHist2D = static_cast<Hist2D::SingleHist2D*>(hist2d->GetComponent(process));
    TH2D * th2d;
    if (process->type_ == Process::Type::data) {
      th2d = static_cast<TH2D*>(singleHist2D->clusterizer_.GetHistogram().Clone());
    } else {
      th2d = static_cast<TH2D*>(singleHist2D->clusterizer_.GetHistogram(hist2d->GetLuminosity()).Clone());
    }
    histograms_2d[process->name_] = th2d;
  }
}

// Plot weighted average across y for each x
void plot_weighted_average(TH2D * response_2d, string const & name, string const & title) {
  int nbinx = response_2d->GetNbinsX();
  int nbiny = response_2d->GetNbinsY();
  vector<Double_t> xArray(nbinx);
  vector<Double_t> yArray(nbinx);
  vector<Double_t> xErrorArray(nbinx);
  vector<Double_t> yErrorArray(nbinx);
  for (int iX = 1; iX <= nbinx; ++iX) {
    // Calculate mean
    double total_entry = 0;
    double total_weight = 0;
    for (int iY = 1; iY <= nbiny; ++iY) {
      total_entry += response_2d->GetBinContent(iX,iY) * response_2d->GetYaxis()->GetBinCenter(iY);
      total_weight += response_2d->GetBinContent(iX,iY);
      //cout<<"x: "<<iX<<" y: "<<iY<<endl;
      //cout<<"  "<<response_2d->GetXaxis()->GetBinLowEdge(iX)<<" "<<response_2d->GetXaxis()->GetBinCenter(iX)<<" "<<response_2d->GetXaxis()->GetBinUpEdge(iX)<<endl;
      //cout<<"  "<<response_2d->GetYaxis()->GetBinLowEdge(iY)<<" "<<response_2d->GetYaxis()->GetBinCenter(iY)<<" "<<response_2d->GetYaxis()->GetBinUpEdge(iY)<<endl;
      //cout<<"  "<<response_2d->GetBinContent(iX,iY)<<endl;
    }
    double average = total_entry/total_weight;
    double deltaX = 0;
    double total_deltaXSquare = 0;
    for (int iY = 1; iY <= nbiny; ++iY) {
      deltaX = response_2d->GetYaxis()->GetBinCenter(iY) - average;
      total_deltaXSquare += response_2d->GetBinContent(iX,iY) * deltaX * deltaX;
    }
    double variance = total_deltaXSquare/total_weight;
    double std_dev = sqrt(variance);
    //cout<<"x bin: "<<response_2d->GetXaxis()->GetBinLowEdge(iX)<<" "<<response_2d->GetXaxis()->GetBinCenter(iX)<<" "<<response_2d->GetXaxis()->GetBinUpEdge(iX)<<endl;
    //cout<<"  total_entry: "<<total_entry<<" total_weight: "<<total_weight<<endl; 
    //cout<<"  average: "<<average<<" std_dev: "<<std_dev<<endl;
    xArray[iX-1] = response_2d->GetXaxis()->GetBinCenter(iX);
    yArray[iX-1] = average;
    xErrorArray[iX-1] = response_2d->GetXaxis()->GetBinCenter(iX) - response_2d->GetXaxis()->GetBinLowEdge(iX);
    yErrorArray[iX-1] = std_dev/sqrt(total_weight);
    yErrorArray[iX-1] = 0;
  }
  // Plot
  TCanvas * c1 = new TCanvas("c1","c1",500, 500);
  TGraphErrors * graph = new TGraphErrors(nbinx, &xArray[0], &yArray[0], &xErrorArray[0], &yErrorArray[0]);
  graph->SetTitle((title+";gen. AK4 p_{T}[GeV];rec. AK4 p_{T}/gen. AK4 p_{T}").c_str());
  graph->SetMarkerStyle(21);
  graph->SetMaximum(1.10);
  graph->SetMinimum(0.90);
  graph->Draw("AP goff");
  c1->SetLeftMargin(0.13);
  c1->SaveAs(string("plots/"+name+".pdf").c_str());
  cout<<"open plots/"<<name<<".pdf"<<endl;
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
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt lin_lumi = lin_norm_info.Title(TitleType::data).Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::data_norm).RatioMaximum(1.86);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_lin_lumi = {lin_lumi};
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

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

  // Set folders
  string mc_production_folder = nanoAodFolder+"/"+production;
  string data_production_folder = nanoAodFolder+"/"+production;
  string signal_production_folder = nanoAodFolder+"/"+production;

  signal_production_folder = "/homes/jbkim/analysis/nano2pico.klamath/out";
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
  string mc_skim_folder = "mc/merged_higmc_preselect/";
  string data_skim_folder = "data/merged_higdata_preselect/";
  string signal_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_preselect/";

  signal_skim_folder = "unskimmed/";


  // Data cuts
  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  NamedFunc triggers_data = "1";
  //if (sample_name == "zll") triggers_data = lepton_triggers;
  //else if (sample_name == "ttbar") triggers_data = lepton_triggers || met_triggers;
  //else if (sample_name == "qcd") triggers_data = met_triggers;
  //else  triggers_data = met_triggers;
  triggers_data = met_triggers;

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

  // Weights
  NamedFunc weight = Higfuncs::final_weight; //"weight"*eff_higtrig_run2*w_years*Functions::w_pileup;

  // Filters for each sample
  NamedFunc search_filters = Higfuncs::final_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2&&weight<1.5"
  NamedFunc ttbar_filters = Higfuncs::final_ttbar_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc zll_filters = Higfuncs::final_zll_pass_filters; //pass_filters&& "met/met_calo<5&&weight<1.5"
  NamedFunc qcd_filters = Higfuncs::final_qcd_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2"


  // Set procs
  map<string, vector<shared_ptr<Process> > > procsDict;
  //procsDict["mc_ttbar"].push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_production_folder, years, mc_skim_folder, 
  //                {"*TTJets_*Lept*","*_TTZ*.root", "*_TTW*.root","*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"}
  //                ),"stitch"));
  //procsDict["signal"].push_back(Process::MakeShared<Baby_pico>("TChiHH(500,0)", Process::Type::signal, kGreen+1, 
  //                attach_folder(signal_production_folder, years, signal_skim_folder, 
  //                {"*TChiHH_mChi-500_mLSP-0_*.root"} ), "1"));

  //procsDict["signal"].push_back(Process::MakeShared<Baby_pico>("TChiHH(500,0)", Process::Type::background, kGreen+1, 
  //                attach_folder(signal_production_folder, years, signal_skim_folder, 
  //                {"*TChiHH_mChi-500_mLSP-0_*.root"} ), "1"));
  procsDict["signal"].push_back(Process::MakeShared<Baby_pico>("TChiHH(950,0)", Process::Type::background, kGreen+1, 
                  attach_folder(signal_production_folder, years, signal_skim_folder, 
                  {"*TChiHH_mChi-950_mLSP-0_*.root"} ), "1"));


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

  {
    pm.Push<Hist2D>(
      //Axis(25, 30, 780, matched_genjet_pt, "Matched gen. AK4 p_{T} (from rec. Higgs) [GeV]", {}),
      Axis(16, 0, 800, matched_genjet_pt, "Matched gen. AK4 p_{T} (from rec. Higgs) [GeV]", {}),
      Axis(100, 0.6, 1.6, matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T} (from rec. Higgs)", {}), 
      search_filters&&search_resolved_cuts, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:rec_higgs_response_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(
      Axis(25, 30, 780, matched_genjet_pt, "Matched gen. AK4 p_{T} (from rec. Higgs) [GeV]", {}),
      Axis(25, 30, 780, matched_jet_pt, "Matched higg rec. AK4 p_{T} (from rec. Higgs) [GeV]", {}), 
      search_filters&&search_resolved_cuts, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:rec_higgs_recPt_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //// 1D Axis({,}, , ) doesn't scale bin entries correctly
    //pm.Push<Hist1D>(Axis(25, 30, 780, matched_jet_pt, "Matched higg rec. AK4 p_{T} (from rec. Higgs) [GeV]", {}),
    //  search_filters&&search_resolved_cuts,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:rec_higgs_matched_recjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 30, 780, matched_genjet_pt, "Matched gen. AK4 p_{T} (from rec. Higgs) [GeV]", {}),
    //  search_filters&&search_resolved_cuts,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:rec_higgs_matched_genjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 2, matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T} [from rec. Higgs]", {}),
    //  search_filters&&search_resolved_cuts,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:rec_higgs_matched_response__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 0.6, matched_deltaR, "Delta R(rec. AK4 jet, gen. AK4 jet) [from rec. Higgs]", {}),
    //  search_filters&&search_resolved_cuts,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:rec_higgs_matched_deltaR__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);


    //// With baseline cuts
    ////pm.Push<Hist2D>(
    ////  Axis(25, 30, 780, all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
    ////  Axis(100, 0.6, 1.6, all_matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T}", {}), 
    ////  search_filters&&search_resolved_cuts, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:response_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    ////pm.Push<Hist2D>(
    ////  Axis(25, 30, 780, all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
    ////  Axis(25, 30, 780, all_matched_jet_pt, "Matched higg rec. AK4 p_{T} [GeV]", {}), 
    ////  search_filters&&search_resolved_cuts, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:recPt_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    ////pm.Push<Hist1D>(Axis(25, 30, 780., all_matched_jet_pt, "Matched higg rec. AK4 p_{T} [GeV]", {}),
    ////  search_filters&&search_resolved_cuts,
    ////  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_recjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    ////pm.Push<Hist1D>(Axis(25, 30, 780., all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
    ////  search_filters&&search_resolved_cuts,
    ////  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_genjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    ////pm.Push<Hist1D>(Axis(20, 0, 2, all_matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T}", {}),
    ////  search_filters&&search_resolved_cuts,
    ////  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_response__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    ////pm.Push<Hist1D>(Axis(20, 0, 0.6, all_matched_deltaR, "Delta R(rec. AK4 jet, gen. AK4 jet)", {}),
    ////  search_filters&&search_resolved_cuts,
    ////  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_deltaR__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

    // No baseline cuts
    pm.Push<Hist2D>(
      //Axis(25, 30, 780., all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
      Axis(16, 0, 800., all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
      Axis(100, 0.6, 1.6, all_matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T}", {}), 
      search_filters, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:response_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    pm.Push<Hist2D>(
      Axis(25, 30, 780, all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
      Axis(25, 30, 780, all_matched_jet_pt, "Matched higg rec. AK4 p_{T} [GeV]", {}), 
      search_filters, procsDict["signal"], plt_2D).Weight(weight).Tag("FixName:recPt_vs_genPt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 30, 780., all_matched_jet_pt, "Matched higg rec. AK4 p_{T} [GeV]", {}),
    //  search_filters,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_recjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(25, 30, 780., all_matched_genjet_pt, "Matched gen. AK4 p_{T} [GeV]", {}),
    //  search_filters,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_genjet_pt__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 2, all_matched_response, "rec. AK4 p_{T} / gen. AK4 p_{T}", {}),
    //  search_filters,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_response__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 0.6, all_matched_deltaR, "Delta R(rec. AK4 jet, gen. AK4 jet)", {}),
    //  search_filters,
    //  procsDict["signal"], plt_lin).Weight(weight).Tag("FixName:matched_deltaR__search_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

    pm.multithreaded_ = !single_thread;
    pm.min_print_ = true;
    pm.MakePlots(1.);
  }

  {
    map<string, TH2D *> histograms_2d;
    get_histograms_2d(pm, /*plotIndex*/ 0, histograms_2d);
    TH2D * response_2d = histograms_2d.begin()->second;
    plot_weighted_average(response_2d, "resolved_cuts", "Resolved cuts");
  }
  {
    map<string, TH2D *> histograms_2d;
    get_histograms_2d(pm, /*plotIndex*/ 2, histograms_2d);
    TH2D * response_2d = histograms_2d.begin()->second;
    plot_weighted_average(response_2d, "no_cuts", "No cuts");
  }


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
