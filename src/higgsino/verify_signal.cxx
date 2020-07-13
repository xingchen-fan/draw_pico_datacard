///// verify_signal - cutflow maker for verifying signal productions

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
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  // vector<string> sigm = {}; 
  vector<string> mchi_sigm = {"225","400","400","700","700","700"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "150","0",  "150","450"}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  float lumi = 35.9;
  string bfolder = "/net/cms25";

  string mc_production = "higgsino_humboldt"; // higgsino_eldorado
  string year = "2016"; // 2017, 2018

  string foldermc(bfolder+"/cms25r5/pico/NanoAODv5/"+mc_production+"/"+year+"/mc/unskimmed/");
  string foldersig(bfolder+"/cms25r5/pico/NanoAODv5/"+mc_production+"/"+year+"/SMS-TChiHH_2D/unskimmed/");
  if (year == "2017") {
	  lumi = 41.5;
  }
  if (year == "2018") {
	  lumi = 60.0;
  }

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lep*"});
  mctags["ttx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
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
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}", Process::Type::background,1,
              attach_folder(foldermc, mctags["tt"]),c_ps));
  for (unsigned isig(0); isig<mchi_sigm.size(); isig++)
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+mchi_sigm[isig]+","+mlsp_sigm[isig]+")", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "1"));


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
    "met>0",0,0, "weight"),
  TableRow("$N_\\text{vl}=0$, $4\\leq N_\\text{j}\\geq 5$", 
    "nvlep==0 && njet>=4 && njet<=5",0,0,"weight"),
  TableRow("$N_\\text{b}\\geq 2$", 
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2",0,0,"weight"),
  TableRow("$E_\\text{t}^\\text{miss}>150 \\text{ GeV}$", 
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150",0,0,"weight"),
  TableRow("$N_\\text{tk}=0$", 
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0",0,0,"weight"),
  TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met",0,0,"weight"),
  TableRow("$\\Delta m<40 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && hig_cand_dm[0]<=40",0,0,"weight"),
  TableRow("$\\Delta R_\\text{max}<2.2$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2",0,0,"weight"),
  TableRow("$100<\\langle m\\rangle <140 \\text{ GeV}$",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140",0,0,"weight"),
  TableRow("N_\\text{b}=4",        
    "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>150 && ntk==0 && !low_dphi_met && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,"weight"),

  //TableRow("\\multicolumn{"+ncols+"}{c}{HIG: $100<\\left< m \\right>\\leq140$}\\\\%", 
  //  "met>1e6",0,1, "weight"),

  //TableRow("$100<\\left< m \\right>\\leq140$ GeV", 
  //  baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !low_dphi &&"+hig,0,0, wgt),
  //TableRow("3b + 4b", 
  //  baseline +"  && ntk==0 && met>150 &&"+c_ge3b+"&& !low_dphi &&"+hig,0,0,wgt),
  //TableRow("4b", 
  //  baseline + " && ntk==0 && met>150 &&"+c_4b+"  && !low_dphi &&"+hig,0,0, wgt),
  //TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", 
  //  baseline + " && ntk==0 && met>200 &&"+c_4b+"  && !low_dphi &&"+hig,0,0,wgt),
  //TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", 
  //  baseline + " && ntk==0 && met>300 &&"+c_4b+"  && !low_dphi &&"+hig,0,0,wgt),
  //TableRow("$p_{\\rm T}^{\\rm miss}>450$ GeV", 
  //  baseline + " && ntk==0 && met>450 &&"+c_4b+"  && !low_dphi &&"+hig,0,1,wgt),

  // TableRow("\\multicolumn{"+ncols+"}{c}{SBD: $\\left< m\\right> <100$ or $140<\\left< m\\right>\\leq200$}\\\\%", 
  //   "met>1e6",0,1, wgt),

  // TableRow("CHECK:: SBD 2b, MET 150-200", 
  //   baseline + " && ntk==0 && met>150 && met<=200 &&"+c_2b+"   && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("$\\left< m\\right> <100$ or $140<\\left< m\\right>\\leq200$ GeV", 
  //   baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("3b + 4b", 
  //   baseline +"  && ntk==0 && met>150 &&"+c_ge3b+"&& !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("4b", 
  //   baseline + " && ntk==0 && met>150 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", 
  //   baseline + " && ntk==0 && met>200 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", 
  //   baseline + " && ntk==0 && met>300 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>450$ GeV", 
  //   baseline + " && ntk==0 && met>450 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),

  },procs,0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
