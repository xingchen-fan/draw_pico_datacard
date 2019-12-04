///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 35.9;
  string tag = "";
}

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  string foldermc(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higbase/");
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higtight/");

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*"});
  mctags["topx"]     = set<string>({"*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"});

  //    Define processes, including intersections
  //--------------------------------------------------
  NamedFunc base_filters = HigUtilities::pass_2016 && "met/met_calo<5"; //since pass_fsjets is not quite usable...

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(foldermc, mctags["other"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z$+$jets", Process::Type::background, kOrange+1,
                  attach_folder(foldermc,mctags["zjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W$+$jets", Process::Type::background, kGreen+1,
                  attach_folder(foldermc,mctags["wjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("$t\\bar{t}$", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc, mctags["tt"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("$t/t\\bar{t}+X$", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc, mctags["topx"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(foldermc, mctags["qcd"]),"stitch")); 

  vector<int> sigm = {700};
  for (unsigned isig(0); isig<sigm.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+to_string(sigm[isig])+",1)", Process::Type::signal, 
      1, {foldersig+"*TChiHH_mChi-"+to_string(sigm[isig])+"_*.root"}));
    cout<<"Adding: "<<foldersig+"*TChiHH_mChi-"+to_string(sigm[isig])+"_*.root"<<endl;
  }

  NamedFunc wgt = "w_lumi*w_isr";
  vector<string> cuts;
  // cuts.push_back("met>150 && ht>300 ");
  // cuts.push_back("!lowDphiFix");
  // cuts.push_back("nvlep==0");
  // cuts.push_back("ntk==0");
  cuts.push_back("met>300 && ht>300 && !lowDphiFix && nvlep==0 && ntk==0");
  cuts.push_back("ht>600");
  cuts.push_back("nfjet>=2");
  cuts.push_back("fjet_pt[0]>300 && fjet_pt[1]>300");
  cuts.push_back("fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250");
  cuts.push_back("fjet_deep_md_hbb_btv[0]>0.7 && fjet_deep_md_hbb_btv[1]>0.7");
  cuts.push_back("fjet_msoftdrop[0]>85 && fjet_msoftdrop[0]<=135 && fjet_msoftdrop[1]>85 && fjet_msoftdrop[1]<=135");

  vector<TableRow> tablerows;
  string tmp ="";
  for (unsigned i(0); i<cuts.size(); i++) {
    if (i>0) tmp += "&&"; 
    tmp += cuts[i];
    tablerows.push_back(TableRow("$"+CodeToLatex(cuts[i])+"$", base_filters && tmp, 0,0, wgt));
  }

  PlotMaker pm;
  TString tabname = "table";
  pm.Push<Table>(tabname.Data(), tablerows, procs, 0);
  pm.min_print_ = true;
  pm.MakePlots(lumi);



  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
