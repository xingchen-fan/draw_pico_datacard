///// table_preds: Makes piecharts

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
  float lumi = 35.9;
  // vector<string> sigm = {}; 
  vector<string> sigm = {"225","400","700"}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  //string mc_production = "higgsino_angeles"; // higgsino_eldorado
  //string mc_production = "higgsino_eldorado";
  string mc_production = "higgsino_humboldt";
  string year = "2016"; // 2017, 2018

  string foldermc(bfolder+"/cms29r0/pico/NanoAODv5/"+mc_production+"/"+year+"/mc1/merged_higmc_higloose/");
  //string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_unskimmed/");
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/"+mc_production+"/2016/SMS-TChiHH_2D/unskimmed/");

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
  for (unsigned isig(0); isig<sigm.size(); isig++) {
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}, "1"));
    cout<<foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"<<endl;
  }


  string c_2bt  = "nbt>=2";
  string c_2b   = "nbt==2&&nbm==2";
  string c_3b   = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b   = "nbt>=2&&nbm>=3&&nbl>=4";
  string c_ge3b = "nbt>=2&&nbm>=3";

  string c_hig_dm = "hig_cand_dm[0]<=40";
  string c_drmax  = "hig_cand_drmax[0]<=2.2";
  string hig = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && (hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  string sbd = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && !(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";

  NamedFunc wgt = "w_lumi*w_isr";//Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  string baseline = "nvlep==0 && njet>=4 && njet<=5"; //pass_muon_jet && met/met_calo<5
  string sigonly = "type>100e3";

  string ncols = to_string(procs.size()+2);  
  PlotMaker pm;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("1", 
    "1",0,0, wgt),
  TableRow("$0\\ell$, $\\text{4-5 jets}$", 
    baseline,0,0, wgt),
  TableRow("nbt>=2", 
    baseline+"&&"+c_2bt,0,0, wgt),
  TableRow("MET>150", 
    baseline+"&&"+c_2bt+"&&met>150",0,0, wgt),
  TableRow("$0\\ell$, $\\text{4-5 jets}$, $N_{\\text{b,T}}\\geq 2$, $p_{\\rm T}^{\\rm miss}>150$ GeV", 
    baseline+"&& met>150 &&" + c_2bt,0,0, wgt),
  TableRow("Track veto", 
    baseline + " && ntk==0 && met>150 &&" + c_2bt,0,0, wgt),
  TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$",        
    baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix",0,0, wgt),
  TableRow("$|\\Delta m| < 40$ GeV",     
    baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix && "+c_hig_dm,0,0, wgt),
  TableRow("$\\Delta R_{\\text{max}} < 2.2$",                    
    baseline +"  && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix && "+c_hig_dm+" && "+c_drmax,0,1,wgt),

  TableRow("\\multicolumn{"+ncols+"}{c}{HIG: $100<\\left< m \\right>\\leq140$}\\\\%", 
    "met>1e6",0,1, wgt),

  TableRow("$100<\\left< m \\right>\\leq140$ GeV", 
    baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix &&"+hig,0,0, wgt),
  TableRow("3b + 4b", 
    baseline +"  && ntk==0 && met>150 &&"+c_ge3b+"&& !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("4b", 
    baseline + " && ntk==0 && met>150 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0, wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", 
    baseline + " && ntk==0 && met>200 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", 
    baseline + " && ntk==0 && met>300 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>450$ GeV", 
    baseline + " && ntk==0 && met>450 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,1,wgt),

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
