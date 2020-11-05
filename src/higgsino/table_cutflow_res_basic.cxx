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
#include "higgsino/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

// Modified to use no stitch
// Modified to use only singleLept no extension
// Modified lumi to be 1
// Removes weights

namespace{
  //float lumi = 35.9;
  float lumi = 1;
  // vector<string> sigm = {}; 
  vector<string> sigm = {"225","400","700"}; 
}

NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
  current_cut = current_cut && additional_cut;
  return current_cut;
}
string add_string(string & current_string, string additional_string) {
  if (current_string == "") current_string = additional_string;
  else current_string += ", "+additional_string;
  return current_string;
}

const NamedFunc Flag_goodVertices("Flag_goodVertices", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_goodv();
});
const NamedFunc Flag_globalSuperTightHalo2016Filter("Flag_globalSuperTightHalo2016Filter", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()>=100e3) return 1;// FastSim
  else return b.pass_cschalo_tight();
});
const NamedFunc Flag_HBHENoiseFilter("Flag_HBHENoiseFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_hbhe();
});
const NamedFunc Flag_HBHENoiseIsoFilter("Flag_HBHENoiseIsoFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_hbheiso();
});
const NamedFunc Flag_EcalDeadCellTriggerPrimitiveFilter("Flag_EcalDeadCellTriggerPrimitiveFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_ecaldeadcell();
});
const NamedFunc Flag_BadPFMuonFilter("Flag_BadPFMuonFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_badpfmu();
});
const NamedFunc Flag_eeBadScFilter("Flag_eeBadScFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()<1000 && b.type()>0) return b.pass_eebadsc();// Data
  else return 1;
});
const NamedFunc Flag_ecalBadCalibFilterV2("Flag_ecalBadCalibFilterV2", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()==2016) return 1;
  else return b.pass_badcalib();
});
const NamedFunc pass_susyCleanFastsimJets("pass_susyCleanFastsimJets", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()>=100e3) return b.pass_jets();// FastSim
  else return 1;
});
//const NamedFunc pass_ra4St("pass_ra4St", [](const Baby &b) -> NamedFunc::ScalarType{
//  if (b.st()<10000) return 0;
//  else return 1;
//});
const NamedFunc pass_ra4MuonJet("pass_ra4MuonJet", [](const Baby &b) -> NamedFunc::ScalarType{
  return (b.pass_ra2_badmu() || b.pass_muon_jet());
});
const NamedFunc pass_ra4MetMetCalo("pass_ra4MetMetCalo", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.met()/b.met_calo()<5;
});
const NamedFunc pass_ra2bLowNeutralJetFilter("pass_ra2bLowNeutralJetFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_low_neutral_jet();
});
const NamedFunc pass_ra2bHTRatioDPhiTightFilter("pass_ra2bHTRatioDPhiTightFilter", [](const Baby &b) -> NamedFunc::ScalarType{
  return b.pass_htratio_dphi_tight();
});

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-") || Contains(hostname, "physics.ucsb.edu"))
    // TODO
    bfolder = "/net/cms29/cms29r0/"; // In laptops, you can't create a /net folder
    //bfolder = "/net/cms25/cms25r5/jbkim/"; // In laptops, you can't create a /net folder

  // TODO
  //string production = "higgsino_angeles"; // higgsino_eldorado
  string production = "higgsino_eldorado";
  //string production = "higgsino_humboldt";
  string year = "2016"; // 2017, 2018
  //string skim = "unskimmed";
  string skim = "skim_met150";
  //string skim = "merged_higmc_preselect";

  string foldermc(bfolder+"/pico/NanoAODv5/"+production+"/"+year+"/mc/"+skim+"/");
  string foldersig(bfolder+"/pico/NanoAODv5/"+production+"/"+year+"/TChiHH/"+skim+"/");

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_SingleLept*",
                                  "*TTJets_DiLept*",
                                });
  mctags["tt_sl"] = set<string>({"*TTJets_*SingleLeptFromT_Tune*__100*"});
  mctags["ttx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  // TODO none
  //string c_ps = "1";
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
  // TODO none
  //procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t} {\\rm{ SingleLeptFromT\\_Tune}}", Process::Type::background,1,
  //            attach_folder(foldermc, mctags["tt_sl"]),c_ps));

  string lsp = (production == "higgsino_angeles" ? "1" : "0");
  // TODO
  for (unsigned isig(0); isig<sigm.size(); isig++)
  {
    //procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", 
    //  Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-"+lsp+"_*.root"}, "1"));
    //cout<<foldersig+"*TChiHH_mChi-"+sigm[isig]+"_mLSP-"+lsp+"_*.root"<<endl;
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

  // TODO none
  NamedFunc wgt = "w_isr" * Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  //NamedFunc wgt = "1";//w_isr * Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  //string baseline = "nvlep==0 && njet>=4 && njet<=5"; //pass_muon_jet && met/met_calo<5

  string ncols = to_string(procs.size()+1);  
  PlotMaker pm;
  NamedFunc current_cut = "1";
  pm.Push<Table>("cutflow", vector<TableRow>{

  //TableRow("(No cut) (weight = 1)", current_cut,0,0, wgt),
  //TableRow("$p_{\\rm T}^{\\rm miss}>150$ (weight = cross-section)", current_cut,0,0, wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>150$ GeV", add_cut(current_cut, "met>150"),0,0, wgt),

  TableRow("goodVerticies", add_cut(current_cut, Flag_goodVertices),0,0, wgt),
  TableRow("globalSuperTightHalo2016Filter", add_cut(current_cut, Flag_globalSuperTightHalo2016Filter),0,0, wgt),
  TableRow("HBHENoiseFilter", add_cut(current_cut, Flag_HBHENoiseFilter),0,0, wgt),
  TableRow("HBHENoiseIsoFilter", add_cut(current_cut, Flag_HBHENoiseIsoFilter),0,0, wgt),
  TableRow("EcalDeadCellTriggerPrimitiveFilter", add_cut(current_cut, Flag_EcalDeadCellTriggerPrimitiveFilter),0,0, wgt),
  TableRow("BadPFMuonFilter", add_cut(current_cut, Flag_BadPFMuonFilter),0,0, wgt),
  //TableRow("eeBadScFilter", add_cut(current_cut, Flag_eeBadScFilter),0,0, wgt),
  TableRow("ecalBadCalibFilterV2", add_cut(current_cut, Flag_ecalBadCalibFilterV2),0,0, wgt),
  TableRow("susyCleanFastsimJets", add_cut(current_cut, pass_susyCleanFastsimJets),0,0, wgt),
  //TableRow("pass_ra4St", add_cut(current_cut, pass_ra4St),0,0, wgt), // There is currently no st in nano
  TableRow("pass\\_ra4MuonJet", add_cut(current_cut, pass_ra4MuonJet),0,0, wgt),
  TableRow("pass\\_ra4MetMetCalo", add_cut(current_cut, pass_ra4MetMetCalo),0,0, wgt),
  //TableRow("pass\\_ra2bLowNeutralJetFilter", add_cut(current_cut, pass_ra2bLowNeutralJetFilter),0,0, wgt),
  //TableRow("pass\\_ra2bHTRatioDPhiTightFilter", add_cut(current_cut, pass_ra2bHTRatioDPhiTightFilter),0,0, wgt),

  //TableRow("$0\\ell$", add_cut(current_cut, "nvlep==0"),0,0, wgt),
  TableRow("0 muon", add_cut(current_cut, "nvmu==0"),0,0, wgt),
  TableRow("0 electrons", add_cut(current_cut, "nvel==0"),0,0, wgt),
  TableRow("Track veto", add_cut(current_cut, "ntk==0"),0,0, wgt),

  TableRow("jets $>$ 1", add_cut(current_cut, "njet>1"),0,0, wgt),
  TableRow("$\\text{4-5 jets}$", add_cut(current_cut, "njet>=4 && njet<=5"),0,0, wgt),

  TableRow("$N_{\\text{b,T}}\\geq 2$", add_cut(current_cut, c_2bt),0,0, wgt),
  TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$", add_cut(current_cut, "!lowDphiFix"),0,0, wgt),
  TableRow("$|\\Delta m| < 40$ GeV", add_cut(current_cut, c_hig_dm),0,0, wgt),
  TableRow("$\\Delta R_{\\text{max}} < 2.2$", add_cut(current_cut, c_drmax),0,1,wgt),

  TableRow("\\multicolumn{"+ncols+"}{c}{HIG: $100<\\left< m \\right>\\leq140$}\\\\%", "met>1e6",0,1, wgt),

  TableRow("$100<\\left< m \\right>\\leq140$ GeV", add_cut(current_cut, hig),0,0, wgt),
  TableRow("3b + 4b", add_cut(current_cut, c_ge3b),0,0,wgt),
  TableRow("4b", add_cut(current_cut, c_4b),0,0, wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", add_cut(current_cut, "met>200"),0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", add_cut(current_cut, "met>300"),0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>450$ GeV", add_cut(current_cut, "met>450"),0,1,wgt),

  },procs,0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
