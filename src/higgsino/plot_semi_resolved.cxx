///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "Math/Vector4D.h"
#include "TVector2.h"
#include "TMath.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  //float lumi = 35.9;
  float lumi = 137;
  bool doSignal = true;
  bool showData = false;
  string mass_points_string = "";
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


Double_t DeltaR(const ROOT::Math::PtEtaPhiMVector & v1, const ROOT::Math::PtEtaPhiMVector & v2)
{
   Double_t deta = v1.Eta()-v2.Eta();
   Double_t dphi = TVector2::Phi_mpi_pi(v1.Phi()-v2.Phi());
   return TMath::Sqrt( deta*deta+dphi*dphi );
}

int get_semi_resolved_fjet_index(vector<float> const & fjet_msoftdrop, vector<float> const & fjet_mva_hbb_btv)
{
  int best_fjet_index = -1;
  // Find event that has one fjet with mass [85,150] and htag>0.6
  vector<unsigned> pass_softdrop_indices;
  for (unsigned ifjet = 0; ifjet < fjet_msoftdrop.size(); ifjet++)
  {
    if (fjet_msoftdrop[ifjet]>85 && fjet_msoftdrop[ifjet]<150) pass_softdrop_indices.push_back(ifjet);
  }
  vector<unsigned> pass_htagger_indices;
  for (unsigned index = 0; index < pass_softdrop_indices.size(); index++)
  {
    unsigned ifjet = pass_softdrop_indices[index];
    if (fjet_mva_hbb_btv[ifjet]>0.6) pass_htagger_indices.push_back(ifjet);
  }
  if (pass_htagger_indices.size() == 1) best_fjet_index = pass_htagger_indices[0];
  return best_fjet_index;
}

extern const NamedFunc semi_fjet_m("semi_fjet_m", [](const Baby &b) -> NamedFunc::ScalarType{
  int best_fjet_index = get_semi_resolved_fjet_index(*b.fjet_msoftdrop(), *b.fjet_mva_hbb_btv());
  if (best_fjet_index==-1) return 0;
  //return ROOT::Math::PtEtaPhiMVector((*b.fjet_pt())[best_fjet_index], (*b.fjet_eta())[best_fjet_index], (*b.fjet_phi())[best_fjet_index], (*b.fjet_m())[best_fjet_index]).M();
  return (*b.fjet_msoftdrop())[best_fjet_index];
});

vector<int> get_semi_bbjet_indices(vector<float> const & fjet_msoftdrop, vector<float> const & fjet_mva_hbb_btv, vector<int> const & jet_fjet_idx, vector<float> const & jet_deepcsv)
{
  vector<int> bbjet_index {-1, -1};
  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
  //cout<<"event start"<<endl;
  int best_fjet_index = get_semi_resolved_fjet_index(fjet_msoftdrop, fjet_mva_hbb_btv);
  if (best_fjet_index==-1) return bbjet_index;
  // Reconstruct resolved
  // Find best btag jets
  pair<float, float> best_btag = {-100,-100};
  pair<int, int> best_btag_index = {-1, -1};
  for (unsigned ijet = 0; ijet < jet_fjet_idx.size(); ijet++)
  {
    if (best_fjet_index == jet_fjet_idx[ijet]) continue;
    float btag = jet_deepcsv[ijet];
    //cout<<"["<<ijet<<"] btag: "<<btag<<endl;
    if (btag > best_btag.first) 
    { 
      best_btag.first = btag;
      best_btag_index.first = ijet;
    }
    else if (btag > best_btag.second) 
    {
      best_btag.second = btag;
      best_btag_index.second = ijet;
    }
    //cout<<"Best ["<<best_btag.first<<", "<<best_btag_index.first<<"]"<<"["<<best_btag.second<<", "<<best_btag_index.second<<"]"<<endl;
  }
  //cout<<"event end"<<endl;
  if (best_btag_index.first == -1 || best_btag_index.second == -1) return bbjet_index;
  if (best_btag.second < btag_wpts_2016[2]) return bbjet_index;
  bbjet_index.at(0) = best_btag_index.first;
  bbjet_index.at(1) = best_btag_index.second;
  return bbjet_index;
}

ROOT::Math::PtEtaPhiMVector get_semi_bbjet(vector<float> const & fjet_msoftdrop, vector<float> const & fjet_mva_hbb_btv, vector<int> const & jet_fjet_idx, vector<float> const & jet_deepcsv, vector<float> const & jet_pt, vector<float> const & jet_eta, vector<float> const & jet_phi, vector<float> const & jet_m)
{
  vector<int> bbjet_indices = get_semi_bbjet_indices(fjet_msoftdrop, fjet_mva_hbb_btv, jet_fjet_idx, jet_deepcsv);
  if (bbjet_indices.at(0) == -1 || bbjet_indices.at(1) == -1) return ROOT::Math::PtEtaPhiMVector();
  ROOT::Math::PtEtaPhiMVector bjet_1 (jet_pt[bbjet_indices.at(0)], jet_eta[bbjet_indices.at(0)], jet_phi[bbjet_indices.at(0)], jet_m[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector bjet_2 (jet_pt[bbjet_indices.at(1)], jet_eta[bbjet_indices.at(1)], jet_phi[bbjet_indices.at(1)], jet_m[bbjet_indices.at(1)]);
  ROOT::Math::PtEtaPhiMVector higgs = bjet_1 + bjet_2;
  return higgs;
}

extern const NamedFunc semi_bbjet_m("semi_bbjet_m", [](const Baby &b) -> NamedFunc::ScalarType{
  return get_semi_bbjet((*b.fjet_msoftdrop()), (*b.fjet_mva_hbb_btv()), (*b.jet_fjet_idx()), (*b.jet_deepcsv()), (*b.jet_pt()), (*b.jet_eta()), (*b.jet_phi()), (*b.jet_m())).M();
});

extern const NamedFunc semi_bbjet_dr("semi_bbjet_dr", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> bbjet_indices = get_semi_bbjet_indices((*b.fjet_msoftdrop()), (*b.fjet_mva_hbb_btv()), (*b.jet_fjet_idx()), (*b.jet_deepcsv()));
  if (bbjet_indices.at(0) == -1 || bbjet_indices.at(1) == -1) return -999;
  ROOT::Math::PtEtaPhiMVector bjet_1 ((*b.jet_pt())[bbjet_indices.at(0)], (*b.jet_eta())[bbjet_indices.at(0)], (*b.jet_phi())[bbjet_indices.at(0)], (*b.jet_m())[bbjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector bjet_2 ((*b.jet_pt())[bbjet_indices.at(1)], (*b.jet_eta())[bbjet_indices.at(1)], (*b.jet_phi())[bbjet_indices.at(1)], (*b.jet_m())[bbjet_indices.at(1)]);
  return DeltaR(bjet_1, bjet_2);
});

//extern const NamedFunc semi_fjet_m("semi_fjet_m", [](const Baby &b) -> NamedFunc::ScalarType{
//  // Find event that has one fjet with mass [85,135] and htag>0.3
//  vector<unsigned> pass_softdrop_indices;
//  for (unsigned ifjet = 0; ifjet < (*b.fjet_msoftdrop()).size(); ifjet++)
//  {
//    if ((*b.fjet_msoftdrop())[ifjet]>85 && (*b.fjet_msoftdrop())[ifjet]<135) pass_softdrop_indices.push_back(ifjet);
//  }
//  vector<unsigned> pass_htagger_indices;
//  for (unsigned index = 0; index < pass_softdrop_indices.size(); index++)
//  {
//    unsigned ifjet = pass_softdrop_indices[index];
//    if ((*b.fjet_mva_hbb_btv())[ifjet]>0.6) pass_htagger_indices.push_back(ifjet);
//  }
//  if (pass_htagger_indices.size() != 1) return 0;
//  unsigned best_fjet_index = pass_htagger_indices[0];
//  return ROOT::Math::PtEtaPhiMVector((*b.fjet_pt())[best_fjet_index], (*b.fjet_eta())[best_fjet_index], (*b.fjet_phi())[best_fjet_index], (*b.fjet_m())[best_fjet_index]).M();
//});

//extern const NamedFunc semi_bbjet_m("semi_bbjet_m", [](const Baby &b) -> NamedFunc::ScalarType{
//  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
//  // Find event that has one fjet with mass [85,135] and htag>0.3
//  vector<unsigned> pass_softdrop_indices;
//  for (unsigned ifjet = 0; ifjet < (*b.fjet_msoftdrop()).size(); ifjet++)
//  {
//    if ((*b.fjet_msoftdrop())[ifjet]>85 && (*b.fjet_msoftdrop())[ifjet]<135) pass_softdrop_indices.push_back(ifjet);
//  }
//  vector<unsigned> pass_htagger_indices;
//  for (unsigned index = 0; index < pass_softdrop_indices.size(); index++)
//  {
//    unsigned ifjet = pass_softdrop_indices[index];
//    if ((*b.fjet_mva_hbb_btv())[ifjet]>0.6) pass_htagger_indices.push_back(ifjet);
//  }
//  if (pass_htagger_indices.size() != 1) return 0;
//
//  // Reconstruct resolved
//  // Find best btag
//  pair<float, float> best_btag = {-100,-100};
//  pair<int, int> best_btag_index = {-1, -1};
//  for (unsigned ijet = 0; ijet < (*b.jet_fjet_idx()).size(); ijet++)
//  {
//    if (pass_htagger_indices[0] == static_cast<unsigned>((*b.jet_fjet_idx())[ijet])) continue;
//    float btag = (*b.jet_deepcsv())[ijet];
//    if (btag > best_btag.first) 
//    { 
//      best_btag.first = btag;
//      best_btag_index.first = ijet;
//    }
//    else if (btag > best_btag.second) 
//    {
//      best_btag.second = btag;
//      best_btag_index.second = ijet;
//    }
//  }
//  if (best_btag_index.first == -1 || best_btag_index.second == -1) return 0;
//  if (best_btag_index.second < btag_wpts_2016[2]) return 0;
//  // Combine best btag
//  ROOT::Math::PtEtaPhiMVector bjet_1 ((*b.jet_pt())[best_btag_index.first], (*b.jet_eta())[best_btag_index.first], (*b.jet_phi())[best_btag_index.first], (*b.jet_m())[best_btag_index.first]);
//  ROOT::Math::PtEtaPhiMVector bjet_2 ((*b.jet_pt())[best_btag_index.second], (*b.jet_eta())[best_btag_index.second], (*b.jet_phi())[best_btag_index.second], (*b.jet_m())[best_btag_index.second]);
//  ROOT::Math::PtEtaPhiMVector h = bjet_1 + bjet_2;
//  return h.M();
//});

//extern const NamedFunc semi_bbjet_dr("semi_bbjet_dr", [](const Baby &b) -> NamedFunc::ScalarType{
//  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
//  // Find event that has one fjet with mass [85,135] and htag>0.3
//  vector<unsigned> pass_softdrop_indices;
//  for (unsigned ifjet = 0; ifjet < (*b.fjet_msoftdrop()).size(); ifjet++)
//  {
//    if ((*b.fjet_msoftdrop())[ifjet]>85 && (*b.fjet_msoftdrop())[ifjet]<135) pass_softdrop_indices.push_back(ifjet);
//  }
//  vector<unsigned> pass_htagger_indices;
//  for (unsigned index = 0; index < pass_softdrop_indices.size(); index++)
//  {
//    unsigned ifjet = pass_softdrop_indices[index];
//    if ((*b.fjet_mva_hbb_btv())[ifjet]>0.6) pass_htagger_indices.push_back(ifjet);
//  }
//  if (pass_htagger_indices.size() != 1) return 0;
//
//  // Reconstruct resolved
//  // Find best btag
//  pair<float, float> best_btag = {-100,-100};
//  pair<int, int> best_btag_index = {-1, -1};
//  for (unsigned ijet = 0; ijet < (*b.jet_fjet_idx()).size(); ijet++)
//  {
//    if (pass_htagger_indices[0] == static_cast<unsigned>((*b.jet_fjet_idx())[ijet])) continue;
//    float btag = (*b.jet_deepcsv())[ijet];
//    if (btag > best_btag.first) 
//    { 
//      best_btag.first = btag;
//      best_btag_index.first = ijet;
//    }
//    else if (btag > best_btag.second) 
//    {
//      best_btag.second = btag;
//      best_btag_index.second = ijet;
//    }
//  }
//  if (best_btag_index.first == -1 || best_btag_index.second == -1) return 0;
//  if (best_btag_index.second < btag_wpts_2016[2]) return 0;
//  // Combine best btag
//  ROOT::Math::PtEtaPhiMVector bjet_1 ((*b.jet_pt())[best_btag_index.first], (*b.jet_eta())[best_btag_index.first], (*b.jet_phi())[best_btag_index.first], (*b.jet_m())[best_btag_index.first]);
//  ROOT::Math::PtEtaPhiMVector bjet_2 ((*b.jet_pt())[best_btag_index.second], (*b.jet_eta())[best_btag_index.second], (*b.jet_phi())[best_btag_index.second], (*b.jet_m())[best_btag_index.second]);
//
//  return DeltaR(bjet_1, bjet_2);
//});


//extern const NamedFunc pass_semi_cand("pass_semi_cand", [](const Baby &b) -> NamedFunc::ScalarType{
//  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
//  // Find event that has one fjet with mass [85,135] and htag>0.3
//  vector<unsigned> pass_softdrop_indices;
//  for (unsigned ifjet = 0; ifjet < (*b.fjet_msoftdrop()).size(); ifjet++)
//  {
//    if ((*b.fjet_msoftdrop())[ifjet]>85 && (*b.fjet_msoftdrop())[ifjet]<135) pass_softdrop_indices.push_back(ifjet);
//  }
//  vector<unsigned> pass_htagger_indices;
//  for (unsigned index = 0; index < pass_softdrop_indices.size(); index++)
//  {
//    unsigned ifjet = pass_softdrop_indices[index];
//    if ((*b.fjet_mva_hbb_btv())[ifjet]>0.6) pass_htagger_indices.push_back(ifjet);
//  }
//  if (pass_htagger_indices.size() != 1) return 0;
//
//  // Reconstruct resolved
//  // Find best btag
//  pair<float, float> best_btag = {-100,-100};
//  pair<int, int> best_btag_index = {-1, -1};
//  for (unsigned ijet = 0; ijet < (*b.jet_fjet_idx()).size(); ijet++)
//  {
//    if (pass_htagger_indices[0] == static_cast<unsigned>((*b.jet_fjet_idx())[ijet])) continue;
//    float btag = (*b.jet_deepcsv())[ijet];
//    if (btag > best_btag.first) 
//    { 
//      best_btag.first = btag;
//      best_btag_index.first = ijet;
//    }
//    else if (btag > best_btag.second) 
//    {
//      best_btag.second = btag;
//      best_btag_index.second = ijet;
//    }
//  }
//  if (best_btag_index.first == -1 || best_btag_index.second == -1) return 0;
//  if (best_btag_index.second < btag_wpts_2016[2]) return 0;
//  // Combine best btag
//  ROOT::Math::PtEtaPhiMVector bjet_1 ((*b.jet_pt())[best_btag_index.first], (*b.jet_eta())[best_btag_index.first], (*b.jet_phi())[best_btag_index.first], (*b.jet_m())[best_btag_index.first]);
//  ROOT::Math::PtEtaPhiMVector bjet_2 ((*b.jet_pt())[best_btag_index.second], (*b.jet_eta())[best_btag_index.second], (*b.jet_phi())[best_btag_index.second], (*b.jet_m())[best_btag_index.second]);
//
//  if (DeltaR(bjet_1, bjet_2)>2.2) return 0;
//  ROOT::Math::PtEtaPhiMVector h = bjet_1 + bjet_2;
//  if (h.M()>100 && h.M()<140) return 1;
//
//  else return 0;
//});

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

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


  string baseFolder = HigUtilities::getBaseFolder();
  set<int> years = {2016};
  map<string, string> samplePaths;
  samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higloose/";
  //samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higloose/";
  //samplePaths["mc_2016"] = "/net/cms1/cms1r0/jbkim/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higloose/";
  samplePaths["signal_2016"] = "/net/cms1/cms1r0/jbkim/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higloose/";

  NamedFunc metFilters = HigUtilities::pass_2016;
  NamedFunc weight = "w_lumi";

  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, metFilters&&"stitch", sampleProcesses);

  vector<pair<string, string> > massPoints = { {"225", "1"}, {"400", "1"}, {"700", "1"}, {"900", "1"}};
  vector<int> signal_colors = {kGreen, kRed, kCyan, kBlack};
  map<string, set<string> > signalPaths;
  for (auto massPoint : massPoints)
  {
   for (auto year : years)
   {
     signalPaths[HigUtilities::setProcessName("TChiHH",atoi(massPoint.first.c_str()),atoi(massPoint.second.c_str()))].insert(samplePaths["signal_"+to_string(year)]+"*mChi-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root");
     cout<<"Adding "<<samplePaths["signal_"+to_string(year)]+"*mChi-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root"<<endl;
   }
  }
  sampleProcesses["signal"] = vector<shared_ptr<Process> > ();
  int iColor = 0;
  for (auto signalPath : signalPaths)
  {
    sampleProcesses["signal"].push_back(Process::MakeShared<Baby_pico>(signalPath.first, Process::Type::signal, signal_colors[iColor++], signalPath.second, metFilters));
    for (auto path : signalPath.second) cout<<path<<endl;
  }

  vector<shared_ptr<Process> > mcSignalProcesses = sampleProcesses["mc"];
  for (auto proc : sampleProcesses["signal"]) mcSignalProcesses.push_back(proc);

  NamedFunc semi_baseline = !boost_low_dphi&&"njet>=3&&njet<=5 && nvlep==0 && ntk==0"&&semi_fjet_m!=0.&&semi_bbjet_m!=0.&&semi_bbjet_dr<2.2;

  PlotMaker pm;

  pm.Push<Hist1D>(Axis(20,85,150,semi_fjet_m, "semi_fjet_m [GeV]", {}),
    semi_baseline,mcSignalProcesses, plt_lin).Weight(weight);
  pm.Push<Hist1D>(Axis(40,0,500,semi_bbjet_m, "semi_bbjet_m [GeV]", {60, 160}),
    semi_baseline,mcSignalProcesses, plt_lin).Weight(weight);
  pm.Push<Hist1D>(Axis(20,0,2.2,semi_bbjet_dr, "semi_bbjet_dr", {}),
    semi_baseline,mcSignalProcesses, plt_lin).Weight(weight);

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;

}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"no_signal", no_argument, 0, 'n'},    
      {"luminosity", required_argument, 0, 'l'},    
      {"showData", no_argument, 0, 'd'},    
      {"mass_points", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "nl:dp:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      doSignal = false;
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 'd':
      showData = true;
      break;
    case 'p': 
      mass_points_string = optarg; 
      break;
    case 0:
      // optname = long_options[option_index].name;
      // if(optname == "note"){
      //   note = true;
      // }else{
      //   printf("Bad option! Found option name %s\n", optname.c_str());
      //   exit(1);
      // }
      // break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
