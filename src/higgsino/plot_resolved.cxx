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

bool greater_btag(pair<float, unsigned> bjet_1, pair<float, unsigned> bjet_2) {
  return bjet_1.first > bjet_2.first;
}

bool smaller_dm(tuple<double, double, vector<unsigned> > higgs_1, tuple<double, double, vector<unsigned> > higgs_2) {
  return get<0>(higgs_1) < get<0>(higgs_2);
}

vector<unsigned> get_higgs_bbjet_indices(vector<float> const & jet_m, vector<float> const & jet_deepcsv, vector<float> const & jet_pt, vector<float> const & jet_eta, vector<float> const & jet_phi) {
  //cout<<"event start"<<endl;
  // Reconstruct resolved
  // Sort by btag
  vector<pair<float, unsigned> > bjet;
  for (unsigned ijet = 0; ijet < jet_m.size(); ijet++) {
    //cout<<"ijet ["<<ijet<<"] btag: "<<jet_deepcsv[ijet]<<" pt:"<<jet_pt[ijet]<<endl;
    bjet.push_back({jet_deepcsv[ijet],ijet});
  }
  sort(bjet.begin(), bjet.end(), greater_btag);
  for (unsigned ibjet = 0; ibjet < bjet.size(); ibjet++) {
    //cout<<"ibjet ["<<ibjet<<"] btag: "<<bjet[ibjet].first<<" jet index: "<<bjet[ibjet].second<<endl;
  }
  // Make three higgs candidates
  vector<vector<unsigned> > higgs_bjet_indices = {{0,1,2,3}, {0,2,1,3}, {0,3,1,2}};
  // higgs_info = [higgs_dm, higgs_am, higgs_jet_indices]
  vector<tuple<double, double, vector<unsigned> > > higgs_info;
  for (unsigned ihiggs = 0; ihiggs < higgs_bjet_indices.size(); ihiggs++) {
    //cout<<"ihiggs ["<<ihiggs<<"] bjet index: "<<higgs_bjet_indices[ihiggs][0]<<" "<<higgs_bjet_indices[ihiggs][1]<<" "<<higgs_bjet_indices[ihiggs][2]<<" "<<higgs_bjet_indices[ihiggs][3]<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] jet index: "<<bjet[higgs_bjet_indices[ihiggs][0]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][1]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][2]].second<<" "<<bjet[higgs_bjet_indices[ihiggs][3]].second<<endl;
    ROOT::Math::PtEtaPhiMVector higgs1_b1 (jet_pt[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][0]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][0]].second]);
    ROOT::Math::PtEtaPhiMVector higgs1_b2 (jet_pt[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][1]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][1]].second]);
    ROOT::Math::PtEtaPhiMVector higgs2_b1 (jet_pt[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][2]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][2]].second]);
    ROOT::Math::PtEtaPhiMVector higgs2_b2 (jet_pt[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_eta[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_phi[bjet[higgs_bjet_indices[ihiggs][3]].second], jet_m[bjet[higgs_bjet_indices[ihiggs][3]].second]);
    //cout<<"ihiggs ["<<ihiggs<<"] higg1_b1 e: "<<higgs1_b1.E()<<" pt: "<<higgs1_b1.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg1_b2 e: "<<higgs1_b2.E()<<" pt: "<<higgs1_b2.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2_b1 e: "<<higgs2_b1.E()<<" pt: "<<higgs2_b1.Pt()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2_b2 e: "<<higgs2_b2.E()<<" pt: "<<higgs2_b2.Pt()<<endl;

    ROOT::Math::PtEtaPhiMVector higgs1 = higgs1_b1 + higgs1_b2;
    ROOT::Math::PtEtaPhiMVector higgs2 = higgs2_b1 + higgs2_b2;
    //cout<<"ihiggs ["<<ihiggs<<"] higg1: "<<higgs1.E()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] higg2: "<<higgs2.E()<<endl;
    //cout<<"ihiggs ["<<ihiggs<<"] am: "<<(higgs1.M() + higgs2.M())/2<<" dm: "<<fabs(higgs1.M()-higgs2.M())<<endl;

    vector<unsigned> higgs_indices;
    // Sort by higgs pt
    if (higgs1.Pt()>higgs2.Pt()) {
      higgs_indices = {bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second, bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second};
    } else {
      higgs_indices = {bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second, bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second};
    }
    //// Sort by dr
    //if (DeltaR(higgs1_b1, higgs1_b2) < DeltaR(higgs2_b1, higgs2_b2)) {
    //  higgs_indices = {bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second, bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second};
    //} else {
    //  higgs_indices = {bjet[higgs_bjet_indices[ihiggs][2]].second, bjet[higgs_bjet_indices[ihiggs][3]].second, bjet[higgs_bjet_indices[ihiggs][0]].second, bjet[higgs_bjet_indices[ihiggs][1]].second};
    //}

    higgs_info.push_back({fabs(higgs1.M()-higgs2.M()), (higgs1.M() + higgs2.M())/2, higgs_indices});
  }
  // Sort by dm
  sort(higgs_info.begin(), higgs_info.end(), smaller_dm);
  //cout<<"higgs [0] am: "<<get<1>(higgs_info[0])<<" dm: "<<get<0>(higgs_info[0])<<endl;
  //cout<<"higgs [1] am: "<<get<1>(higgs_info[1])<<" dm: "<<get<0>(higgs_info[1])<<endl;
  //cout<<"higgs [2] am: "<<get<1>(higgs_info[2])<<" dm: "<<get<0>(higgs_info[2])<<endl;
  //cout<<"event end"<<endl;

  return get<2>(higgs_info[0]);
}

const NamedFunc boost_low_dphi("boost_low_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  bool low_dphi = false;
  if ((*b.jet_met_dphi()).size()>3) {
    low_dphi = (*b.jet_met_dphi())[0]<0.5 || (*b.jet_met_dphi())[1]<0.5 ||
                          (*b.jet_met_dphi())[2]<0.3 || (*b.jet_met_dphi())[3]<0.3;
  }
  if ((*b.jet_met_dphi()).size()==3) {
    low_dphi = (*b.jet_met_dphi())[0]<0.5 || (*b.jet_met_dphi())[1]<0.5 ||
                          (*b.jet_met_dphi())[2]<0.3;
  }
  return low_dphi;
});
extern const NamedFunc higgs_h1_dr("higgs_h1_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  return DeltaR(h1b1, h1b2);
});
extern const NamedFunc higgs_h2_dr("higgs_h2_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  return DeltaR(h2b1, h2b2);
});
extern const NamedFunc higgs_min_dr("higgs_min_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  float h1_dr = DeltaR(h1b1, h1b2);
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  float h2_dr = DeltaR(h2b1, h2b2);
  return min(h1_dr, h2_dr);
});
extern const NamedFunc higgs_h1_mass("higgs_h1_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  ROOT::Math::PtEtaPhiMVector h1b1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector h1b2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector higgs1 = h1b1 + h1b2;
  return higgs1.M();
});
extern const NamedFunc higgs_h2_mass("higgs_h2_mass", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  ROOT::Math::PtEtaPhiMVector h2b1 ((*b.jet_pt())[bbjet_indices.at(2)], (*b.jet_eta())[bbjet_indices.at(2)], (*b.jet_phi())[bbjet_indices.at(2)], (*b.jet_m())[bbjet_indices.at(2)]);
  ROOT::Math::PtEtaPhiMVector h2b2 ((*b.jet_pt())[bbjet_indices.at(3)], (*b.jet_eta())[bbjet_indices.at(3)], (*b.jet_phi())[bbjet_indices.at(3)], (*b.jet_m())[bbjet_indices.at(3)]);
  ROOT::Math::PtEtaPhiMVector higgs2 = h2b1 + h2b2;
  return higgs2.M();
});
extern const NamedFunc higgs_h1b1_bcsv("higgs_h1b1_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  return (*b.jet_deepcsv())[bbjet_indices.at(0)];
});
extern const NamedFunc higgs_h1b2_bcsv("higgs_h1b2_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  return (*b.jet_deepcsv())[bbjet_indices.at(1)];
});
extern const NamedFunc higgs_h1_bcsv_diff("higgs_h1_bcsv_diff", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  float h1b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(0)];
  float h1b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(1)];
  return fabs(h1b1_bcsv - h1b2_bcsv);
});
extern const NamedFunc higgs_h2b1_bcsv("higgs_h2b1_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  return (*b.jet_deepcsv())[bbjet_indices.at(2)];
});
extern const NamedFunc higgs_h2b2_bcsv("higgs_h2b2_bcsv", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  return (*b.jet_deepcsv())[bbjet_indices.at(3)];
});
extern const NamedFunc higgs_h2_bcsv_diff("higgs_h2_bcsv_diff", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  float h2b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(2)];
  float h2b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(3)];
  return fabs(h2b1_bcsv - h2b2_bcsv);
});
extern const NamedFunc higgs_btag_good("higgs_btag_good", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  float h1b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(0)];
  float h1b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(1)];
  float h2b1_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(2)];
  float h2b2_bcsv = (*b.jet_deepcsv())[bbjet_indices.at(3)];
  return h1b1_bcsv>0.8953&&h1b2_bcsv>0.2217&&h2b1_bcsv>0.8953&&h2b2_bcsv>0.2217;
});

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

extern const NamedFunc jet0_p("jet0_p", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.jet_phi()->size() < 1) return -9999;
  else return ROOT::Math::PtEtaPhiMVector ((*b.jet_pt())[0], (*b.jet_eta())[0], (*b.jet_phi())[0], (*b.jet_m())[0]).P();
});

extern const NamedFunc average_btag("average_btag", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<unsigned> bbjet_indices = get_higgs_bbjet_indices(*b.jet_m(), *b.jet_deepcsv(), *b.jet_pt(), *b.jet_eta(), *b.jet_phi());
  return ((*b.jet_deepcsv())[bbjet_indices.at(0)]+(*b.jet_deepcsv())[bbjet_indices.at(1)]+(*b.jet_deepcsv())[bbjet_indices.at(2)]+(*b.jet_deepcsv())[bbjet_indices.at(3)])/4;
});

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double wgt = 1;
  if (b.type()==106000) {
    return 35.9;
  }
  if (b.SampleType()==2016){
    return wgt*35.9;
  } else if (b.SampleType()==2017){
    return wgt*41.5;
  } else {
    return wgt*59.6;
  }
});

  
namespace{
  bool single_thread = true;
  //bool do_twiki = true;
  double luminosity = 1.;
  //double luminosity = 137;
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

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  //vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

  set<int> years;
  //years = {2016, 2017, 2018};
  years = {2016};
  //years = {2017};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if((Contains(hostname, "cms") || Contains(hostname, "compute-")) && !Contains(hostname, "cms37"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  //string mc_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string sig_base_folder = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";

  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...
  //NamedFunc wgt = "w_lumi*w_isr"*w_years*Higfuncs::eff_higtrig;
  //NamedFunc wgt = "w_lumi*w_isr*137."*Higfuncs::eff_higtrig;
  NamedFunc wgt = "w_lumi*w_isr"*Higfuncs::eff_higtrig;
  if (years.size()==1 && *years.begin()==2016) wgt *= "137.";
  else wgt *= w_years;

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                  "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"
                                 });
  mctags["ttonly"]     = set<string>({"*TTJets_*Lept*"});
  mctags["topx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root"});
  //mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
  //                                 "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
  //                                 "*_QCD_HT2000toInf_*"});
  mctags["qcd"]     = set<string>({
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"});
  mctags["all"] = set<string>({"*TTJets_*Lept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               //"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               //"*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("tt", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),base_filters&&"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),base_filters&&"stitch")); 
  //procs.push_back(Process::MakeShared<Baby_pico>("tt+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["topx"]),base_filters&&"stitch"));
  //procs.push_back(Process::MakeShared<Baby_pico>("ttonly", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["ttonly"]),base_filters&&"stitch"));

  bool plot_btag_comparison = false;
  bool plot_btag_relation = false;
  bool plot_btag_category = false;
  bool plot_ABCD = false;
  bool plot_dphi = true;

  vector<string> sigm = {"450", "700", "950"}; 
  //vector<string> sigm = {"950"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  if (plot_ABCD || plot_dphi) {
    for (unsigned isig(0); isig<sigm.size(); isig++){
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), base_filters));
    }
  }


  string baseline = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150";
  string baseline_fourjet = "njet>=4&&njet<=5&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&nbm>=2"; 
  string drcut = "hig_cand_drmax[0]>=0 && hig_cand_drmax[0]<2.2";
  //string drcut = "hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<2.2";
  //string drcut = "hig_cand_drmax[0]>=0 && hig_cand_drmax[0]<1.1";
  //string drcut = "hig_cand_drmax[0]>=0.8 && hig_cand_drmax[0]<1.2";
  string bcut_2b = "(nbt==2&&nbm==2)";
  string bcut_3b = "(nbt>=2&&nbm==3&&nbl==3)";
  string bcut_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
  string bcuts = bcut_2b+"||"+bcut_3b+"||"+bcut_4b;
  //baseline_fourjet += "&&("+bcut_3b+"||"+bcut_4b+")";
  //baseline_fourjet += "&&"+bcut_4b;
  //baseline_fourjet += "&&"+bcut_2b;

  PlotMaker pm;

  if (plot_btag_comparison) {
    // 2b, 3b, 4b comparison
    vector<shared_ptr<Process> > procs_btag;
    procs_btag.push_back(Process::MakeShared<Baby_pico>("2b", Process::Type::background,sig_colors[0],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("3b", Process::Type::background,sig_colors[1],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("4b", Process::Type::background,sig_colors[2],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b));
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      1, procs_btag, plt_shapes).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1,2.2}),
      1, procs_btag, plt_shapes).Weight(wgt);

    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, higgs_h2_mass, "subleading mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_2b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_3b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "drmax", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&bcut_4b, procs, plt_lin).Weight(wgt);

    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(22, 0, 2.2, higgs_h1_dr, "dr", {1.1, 2.2}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);

    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, higgs_h2_dr, "dr", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);

    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    //pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {}),
    //  base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);
  }

  if (plot_btag_relation) {
    // b tag relation
    pm.Push<Hist2D>(Axis(20, 0, 1, higgs_h1b1_bcsv, "higgs1 b1 btag", {0.2217, 0.6321, 0.8953}),
      Axis(15, 0, 1, higgs_h1b2_bcsv, "higgs1 b2 btag", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1_bcsv_diff, "Btag diff of leading Higg", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist2D>(Axis(20, 0, 1, higgs_h2b1_bcsv, "higgs2 b1 btag", {0.2217, 0.6321, 0.8953}),
      Axis(15, 0, 1, higgs_h2b2_bcsv, "higgs2 b2 btag", {0.2217, 0.6321, 0.8953}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2_bcsv_diff, "Btag diff of sub-leading Higg", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    // b tags
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b1_bcsv, "Leading pt Higg, leading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h1b2_bcsv, "Leading pt Higg, subleading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b1_bcsv, "SubLeading pt Higg, leading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, higgs_h2b2_bcsv, "SubLeading pt Higg, subleading Btag jet, Btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
  }

  if (plot_btag_category) {
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&higgs_btag_good, procs, plt_lin).Weight(wgt);
  }

  if (plot_ABCD) {
    // 2b, 3b, 4b comparison
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);
    vector<shared_ptr<Process> > procs_btag;
    procs_btag.push_back(Process::MakeShared<Baby_pico>("2b", Process::Type::background,sig_colors[0],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("3b", Process::Type::background,sig_colors[1],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b));
    procs_btag.push_back(Process::MakeShared<Baby_pico>("4b", Process::Type::background,sig_colors[2],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b));
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      1, procs_btag, plt_shapes).Weight(wgt).Tag("nbtagVsAm");


    // Average btag vs average higgs mass
    // Correlated because of physics?
    pm.Push<Hist2D>(Axis(20, 0, 1, average_btag, "average_btag", {0.2217, 0.6321, 0.8953}),
      Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, average_btag, "average btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    // See how mass shape changes depending on average btag
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.2217&&average_btag<=0.6321, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.6321&&average_btag<=0.8953, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100, 140}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.8953, procs, plt_lin).Weight(wgt);
    // See how mass shape changes depending on average btag
    vector<shared_ptr<Process> > procs_average_btag;
    procs_average_btag.push_back(Process::MakeShared<Baby_pico>("low b", Process::Type::background,sig_colors[0],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.2217&&average_btag<=0.6321));
    procs_average_btag.push_back(Process::MakeShared<Baby_pico>("mid b", Process::Type::background,sig_colors[1],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.6321&&average_btag<=0.8953));
    procs_average_btag.push_back(Process::MakeShared<Baby_pico>("high b", Process::Type::background,sig_colors[2],
                    attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),
                    base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.8953));
    pm.Push<Hist1D>(Axis(15, 50, 250, "hig_cand_am[0]", "Average mass [GeV]", {100,140}),
      1, procs_average_btag, plt_shapes).Weight(wgt).Tag("btagVsAm");

    // Average b tag vs dr max
    pm.Push<Hist2D>(Axis(20, 0, 1, average_btag, "average_btag", {0.2217, 0.6321, 0.8953}),
      Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max dela R_{bb}", {1.1,2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_2D).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 1, average_btag, "average btag", {}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    // See how drmax changes depending on average btag
    pm.Push<Hist1D>(Axis(16, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.2217&&average_btag<=0.6321, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.6321&&average_btag<=0.8953, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&average_btag>0.8953, procs, plt_lin).Weight(wgt);
    // See how mass shape changes depending on average btag
    pm.Push<Hist1D>(Axis(15, 0, 3.2, "hig_cand_drmax[0]", "Max \\Delta R_{bb}", {1.1, 2.2}),
      1, procs_average_btag, plt_shapes).Weight(wgt).Tag("btagVsDrmax");
  }

  if (plot_dphi) {
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_2b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_3b, procs, plt_lin).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, 0, 3.1416,min_jet_dphi, "min jet dphi", {0.5}),
      base_filters&&"stitch"&&baseline&&baseline_fourjet&&drcut&&bcut_4b, procs, plt_lin).Weight(wgt);
  }

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(luminosity);

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
