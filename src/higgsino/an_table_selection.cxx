///// verify_signal - cutflow maker for verifying signal productions

#include <set>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

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
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  // vector<string> sigm = {}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  float lumi = 35.9;
  string bfolder = "/net/cms25";
  //If you modify mchi_sigm, modify mchi_sigm_int in the NamedFunc below
  vector<string> mchi_sigm = {"175","500","900"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  vector<string> mchi_sigm2d = {"250","350","450"}; 
  vector<string> mlsp_sigm2d = {"50","200","100"}; 

  string mc_production = "higgsino_humboldt"; // higgsino_eldorado
  string year = "2016"; // 2016, 2017, 2018

  string foldermc(bfolder+"/cms25r5/pico/NanoAODv5/"+mc_production+"/"+year+"/mc/skim_met150/");
  string foldersig(bfolder+"/cms25r5/pico/NanoAODv5/"+mc_production+"/"+year+"/SMS-TChiHH_2D/unskimmed/");
  if (year == "2017") {
	  lumi = 41.5;
  }
  if (year == "2018") {
	  lumi = 60.0;
  }

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lep*","*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
  string c_ps = "stitch";
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, 1,
              attach_folder(foldermc, mctags["other"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background, 1,
              attach_folder(foldermc,mctags["singlet"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, 1,
              attach_folder(foldermc, mctags["qcd"]),c_ps)); 
  procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, 1,
              attach_folder(foldermc,mctags["vjets"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,1,
              attach_folder(foldermc, mctags["ttx"]),c_ps));
  for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
    procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "1"));
      std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
  }
  for (unsigned isig(0); isig<mchi_sigm2d.size(); isig++) {
    procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "1"));
      std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
  }

  //give the right weight for 2D Higgsino scans, 1D scans, and other stuff
  //also applies trigger efficiencies
  const NamedFunc mixed_model_weight("mixed_model_weight",[](const Baby &b) -> NamedFunc::ScalarType{
    set<int> mchi_sigm_int = {175, 500, 900}; 
    double trig_eff = 1.0;
    if (b.SampleType()==2016)     	trig_eff = Higfuncs::get_0l_trigeff2016.GetScalar(b);
    else if (b.SampleType()==2017) 	trig_eff = Higfuncs::get_0l_trigeff2017.GetScalar(b);
    else if (b.SampleType()==2018) 	trig_eff = Higfuncs::get_0l_trigeff2018.GetScalar(b);
    if (b.mprod() == -999) {
      //not signal
      return b.weight()*trig_eff;
    }
    if (mchi_sigm_int.count(b.mprod()) > 0) {
      //1d signal
      return b.weight()*trig_eff;
    }
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return b.weight()/xsec1d*xsec2d*trig_eff;
  });

  string c_2bt  = "nbt>=2";
  string c_2b   = "nbt==2&&nbm==2";
  string c_3b   = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b   = "nbt>=2&&nbm>=3&&nbl>=4";
  string c_ge3b = "nbt>=2&&nbm>=3";

  string c_hig_dm = "hig_cand_dm[0]<=40";
  string c_drmax  = "hig_cand_drmax[0]<=2.2";
  string hig = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && (hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  string sbd = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && !(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  string baseline = "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2";

  NamedFunc wgt = "w_lumi*w_isr";//Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  string sigonly = "type>100e3";

  string ncols = to_string(procs.size()+2);  
  PlotMaker pm;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("No selection", 
    "met>0",0,0, mixed_model_weight),
  TableRow("$p_\\text{t}^\\text{miss}>150 \\text{ GeV}$", 
    "met>150",0,0,mixed_model_weight),
  TableRow("$N_\\text{vl}=0$, $4\\leq N_\\text{j}\\leq 5$", 
    "met>150 && nvlep==0 && njet>=4 && njet<=5",0,0,mixed_model_weight),
  TableRow("$N_\\text{b}\\geq 2$", 
    "met>150 && nvlep==0 && njet>=4 && njet<=5 && nbt>=2",0,0,mixed_model_weight),
  TableRow("$N_\\text{tk}=0$", 
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0",0,0,mixed_model_weight),
  TableRow("Fake MET Cuts",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2",0,0,mixed_model_weight),
  TableRow("$\\Delta m<40 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40",0,0,mixed_model_weight),
  TableRow("$\\Delta R_\\text{max}<2.2$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2",0,0,mixed_model_weight),
  TableRow("$100<\\langle m\\rangle <140 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140",0,0,mixed_model_weight),
  TableRow("$N_\\text{b}\\geq 3$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3",0,0,mixed_model_weight),
  TableRow("$N_\\text{b}=4$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
  TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>200 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
  TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>300 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
  TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>400 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
  //TableRow("$\\Delta R_\\text{max}>1.1$",        
  //  "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>450 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
  },procs,false,true,false,true,false,true);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
