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
  string mass_points_string = "225_1,400_1,700_1,1000_1";
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
  if ((*b.jet_met_dphi()).size()==2) {
    low_dphi = (*b.jet_met_dphi())[0]<0.5 || (*b.jet_met_dphi())[1]<0.5;
  }
  return low_dphi;
});


Double_t DeltaR(const ROOT::Math::PtEtaPhiMVector & v1, const ROOT::Math::PtEtaPhiMVector & v2)
{
   Double_t deta = v1.Eta()-v2.Eta();
   Double_t dphi = TVector2::Phi_mpi_pi(v1.Phi()-v2.Phi());
   return TMath::Sqrt( deta*deta+dphi*dphi );
}

vector<int> get_twoj_indices(vector<float> const & jet_m, vector<float> const & jet_deepcsv, vector<float> const & jet_pt)
{
  // Find event that has two jet with mass [85,150]. Select highest b tagged.
  vector<int> hjet_indices {-1, -1};
  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
  //cout<<"event start"<<endl;
  // Find best btag jets
  pair<float, float> best_btag = {-100,-100};
  pair<int, int> best_btag_index = {-1, -1};
  for (unsigned ijet = 0; ijet < jet_m.size(); ijet++)
  {
    if (jet_m[ijet]<85 || jet_m[ijet]>150) continue;
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
  if (best_btag_index.first == -1 || best_btag_index.second == -1) return hjet_indices;
  if (best_btag.second < btag_wpts_2016[1]) return hjet_indices;
  // Sort by pt
  if (jet_pt[best_btag_index.first] > jet_pt[best_btag_index.second])
  {
    hjet_indices.at(0) = best_btag_index.first;
    hjet_indices.at(1) = best_btag_index.second;
  } else {
    hjet_indices.at(1) = best_btag_index.first;
    hjet_indices.at(0) = best_btag_index.second;
  }
  return hjet_indices;
}

extern const NamedFunc twoj_leadingH_m("twoj_leadingH_m", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> hjet_indices = get_twoj_indices((*b.jet_m()), (*b.jet_deepcsv()), (*b.jet_pt()));
  if (hjet_indices.at(0) == -1 || hjet_indices.at(1) == -1) return 0;
  ROOT::Math::PtEtaPhiMVector hjet_0 ((*b.jet_pt())[hjet_indices.at(0)], (*b.jet_eta())[hjet_indices.at(0)], (*b.jet_phi())[hjet_indices.at(0)], (*b.jet_m())[hjet_indices.at(0)]);
  return hjet_0.M();
});

extern const NamedFunc twoj_subleadingH_m("twoj_subleadingH_m", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> hjet_indices = get_twoj_indices((*b.jet_m()), (*b.jet_deepcsv()), (*b.jet_pt()));
  if (hjet_indices.at(0) == -1 || hjet_indices.at(1) == -1) return 0;
  ROOT::Math::PtEtaPhiMVector hjet_1 ((*b.jet_pt())[hjet_indices.at(1)], (*b.jet_eta())[hjet_indices.at(1)], (*b.jet_phi())[hjet_indices.at(1)], (*b.jet_m())[hjet_indices.at(1)]);
  return hjet_1.M();
});

extern const NamedFunc twoj_dm("twoj_dm", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<int> hjet_indices = get_twoj_indices((*b.jet_m()), (*b.jet_deepcsv()), (*b.jet_pt()));
  if (hjet_indices.at(0) == -1 || hjet_indices.at(1) == -1) return -1;
  ROOT::Math::PtEtaPhiMVector hjet_0 ((*b.jet_pt())[hjet_indices.at(0)], (*b.jet_eta())[hjet_indices.at(0)], (*b.jet_phi())[hjet_indices.at(0)], (*b.jet_m())[hjet_indices.at(0)]);
  ROOT::Math::PtEtaPhiMVector hjet_1 ((*b.jet_pt())[hjet_indices.at(1)], (*b.jet_eta())[hjet_indices.at(1)], (*b.jet_phi())[hjet_indices.at(1)], (*b.jet_m())[hjet_indices.at(1)]);
  return fabs(hjet_0.M() - hjet_1.M());
});

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
  samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higloose/";

  NamedFunc metFilters = HigUtilities::pass_2016;
  NamedFunc weight = "w_lumi";

  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, metFilters&&"stitch", sampleProcesses);

  vector<pair<string, string> > massPoints = { {"225","1"}, {"400", "1"}, {"700", "1"}, {"900", "1"} };
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

  //vector<shared_ptr<Process> > mcSignalProcesses = sampleProcesses["signal"];

  NamedFunc twoj_baseline = !boost_low_dphi&&"njet>=2&&njet<=3 && nvlep==0 && ntk==0"&&twoj_leadingH_m!=0.&&twoj_subleadingH_m!=0.;

  PlotMaker pm;

  //pm.Push<Hist1D>(Axis(20,0,30,"njet", "njets [GeV]", {}),
  //  "1",mcSignalProcesses, plt_lin).Weight(weight);
  pm.Push<Hist1D>(Axis(20,85,150,twoj_leadingH_m, "twoj_leadingH_m [GeV]", {}),
    twoj_baseline,mcSignalProcesses, plt_lin).Weight(weight);
  pm.Push<Hist1D>(Axis(40,85,150,twoj_subleadingH_m, "twoj_subleadingH_m [GeV]", {80,160}),
    twoj_baseline,mcSignalProcesses, plt_lin).Weight(weight);
  pm.Push<Hist1D>(Axis(20,0,40,twoj_dm, "twoj_dm [GeV]", {}),
    twoj_baseline,mcSignalProcesses, plt_lin).Weight(weight);

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
