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
  bool csv = false;
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
  // Find event that has one fjet with mass [85,150] and htag>0.3
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
  if (best_btag_index.second < btag_wpts_2016[2]) return bbjet_index;
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


//extern const NamedFunc pass_semi_cand("pass_semi_cand", [](const Baby &b) -> NamedFunc::ScalarType{
//  vector<float> btag_wpts_2016 = {0.2217, 0.6321, 0.8953};
//  // Find event that has one fjet with mass [85,135] and htag>0.3
//  vector<unsigned> pass_softdrop_indices;
//  for (unsigned ifjet = 0; ifjet < (*b.fjet_msoftdrop()).size(); ifjet++)
//  {
//    if ((*b.fjet_msoftdrop())[ifjet]>85 && (*b.fjet_msoftdrop())[ifjet]<150) pass_softdrop_indices.push_back(ifjet);
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

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higloose/");
  //string folderdata(bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/");
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higloose/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  string c_ps = "stitch";
  NamedFunc base_filters = HigUtilities::pass_2016;

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGreen+1,
							    attach_folder(foldermc, mctags["other"]),c_ps&&base_filters));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
							    attach_folder(foldermc, mctags["qcd"]),c_ps&&base_filters)); 
  procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
							    attach_folder(foldermc,mctags["vjets"]),c_ps&&base_filters));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc, mctags["ttx"]),c_ps&&base_filters));

  ////auto data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, 1,
  ////                {folderdata+"*root"}, Higfuncs::trig_hig>0. && "pass"&&filters);
  //auto data = Process::MakeShared<Baby_pico>("Data", Process::Type::data, 1,
  //                {folderdata+"*root"}, "pass"&&base_filters);
  //if (showData) procs.push_back(data);

  if (doSignal) {
    vector<pair<string, string> > massPoints;
    HigUtilities::parseMassPoints(mass_points_string, massPoints);
    vector<string> sigm(massPoints.size());
    for(unsigned iMass=0; iMass<massPoints.size(); ++iMass) sigm[iMass] = massPoints[iMass].first;
    //vector<string> sigm({"225","400", "700"});
    for (unsigned isig(0); isig<sigm.size(); isig++)
      // TODO There is no pass_fsmet 
      //procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", 
      //  Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+sigm[isig]+"*.root"}, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet"&&base_filters));
      procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", 
        Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+sigm[isig]+"*.root"}, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso"&&base_filters));
  }
 
  string filters = "pass_ra2_badmu && met/met_calo<5";
  NamedFunc baseline = HigUtilities::pass_nhig_cand&&(filters+"&& njet>=4 && njet<=5 && nvlep==0 && ntk==0 && !low_dphi");
  string c_2b = "nbt==2&&nbm==2";
  string c_3b = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b = "nbt>=2&&nbm>=3&&nbl>=4";
  string hig = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && (hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  string sbd = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && !(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  //NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;//*(1+Higfuncs::wgt_syst_ttx);// * Higfuncs::wgt_comp;
  NamedFunc wgt = "w_lumi";

  string met_1 = "met>150 && met<=200";
  string met_2 = "met>200 && met<=300";
  string met_3 = "met>300 && met<=450";
  string met_4 = "met>450";

  NamedFunc boosted_baseline = !boost_low_dphi && "met>150 && ht>300 && nvlep==0 && ntk==0 && nfjet>1";
  string boosted_baseline_ak8 = "fjet_pt[0]>300 && fjet_pt[1]>300 && fjet_msoftdrop[0]>85 && fjet_msoftdrop[0]<135 && fjet_msoftdrop[1]>85 && fjet_msoftdrop[1]<135";
  //string boosted_bcuts = "(fjet_mva_hbb_btv[0]>0.3||fjet_mva_hbb_btv[1]>0.3)";
  string boosted_doubleh = "((met<300&&fjet_mva_hbb_btv[0]>0.6&&fjet_mva_hbb_btv[1]>0.6)||(met>300&&fjet_mva_hbb_btv[0]>0.3&&fjet_mva_hbb_btv[1]>0.3))";
  string boosted_singleh = "((fjet_mva_hbb_btv[0]>0.3&&fjet_mva_hbb_btv[1]<=0.3)||(fjet_mva_hbb_btv[0]<=0.3&&fjet_mva_hbb_btv[1]>0.3))";
  string boosted_bcuts = "("+boosted_doubleh + "||" + boosted_singleh+")";

  NamedFunc semi_baseline = !boost_low_dphi && "njet>=3&&njet<=5 && nvlep==0 && ntk==0"&&semi_fjet_m!=0.&&semi_bbjet_m>60&&semi_bbjet_m<160&&semi_bbjet_dr<1.1;
  
  //        Cutflow table
  //-------------------------------- 
  PlotMaker pm;
  string tabname = "regions_combined";
  if (csv) tabname +="_csv";
  pm.Push<Table>(tabname, vector<TableRow>{
  TableRow("Resolved, 3b"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", baseline && c_3b+"&&"+hig+"&&"+met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", baseline && c_3b+"&&"+hig+"&&"+met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", baseline && c_3b+"&&"+hig+"&&"+met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",              baseline && c_3b+"&&"+hig+"&&"+met_4,0,1, wgt),

  TableRow("Resolved, 4b"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", baseline && c_4b+"&&"+hig+"&&"+met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", baseline && c_4b+"&&"+hig+"&&"+met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", baseline && c_4b+"&&"+hig+"&&"+met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",              baseline && c_4b+"&&"+hig+"&&"+met_4,0,1, wgt),

  TableRow("Boosted, 1H"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_singleh+"&&"+met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_singleh+"&&"+met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_singleh+"&&"+met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",              boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_singleh+"&&"+met_4,0,1, wgt),

  TableRow("Boosted, 2H"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_doubleh+"&&"+met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_doubleh+"&&"+met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_doubleh+"&&"+met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",              boosted_baseline&&boosted_baseline_ak8+"&&"+boosted_doubleh+"&&"+met_4,0,1, wgt),

  TableRow("Semi-resolved: resolved-overlap"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", semi_baseline&&met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", semi_baseline&&met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", semi_baseline&&met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",             semi_baseline&&met_4,0,1, wgt),

  TableRow("Semi-resolved: no-overlap"),
  TableRow("$150 < p_{T}^{\\rm{miss}} \\leq 200$", semi_baseline&&!(baseline &&(c_3b || c_4b)&&hig)&&met_1,0,1, wgt),
  TableRow("$200 < p_{T}^{\\rm{miss}} \\leq 300$", semi_baseline&&!(baseline &&(c_3b || c_4b)&&hig)&&met_2,0,1, wgt),
  TableRow("$300 < p_{T}^{\\rm{miss}} \\leq 450$", semi_baseline&&!(baseline &&(c_3b || c_4b)&&hig)&&met_3,0,1, wgt),
  TableRow("$p_{T}^{\\rm{miss}}>450$",             semi_baseline&&!(baseline &&(c_3b || c_4b)&&hig)&&met_4,0,1, wgt),

  },procs,0);

  //string listname = "eventlist";
  //if (csv) listname +="_csv";
  //pm.Push<EventScan>(listname, baseline && "met>150 && nbt>=2 && nbm>=3" && hig, 
  //                   vector<NamedFunc>{"run", "lumiblock", "event"},
  //                   vector<shared_ptr<Process> >{data});

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
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
