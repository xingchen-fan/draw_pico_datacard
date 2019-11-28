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

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 137;
  bool doSignal = false;
  bool res_view = true; 
  bool quick = false;
  string tag = "";
  string mass_points_string = "900_1";//"225_1,400_1,900_1";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higtight/");
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_unskimmed/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  set<string> all_bkg;
  if (!quick) {
    all_bkg.insert(mctags["ttx"].begin(),   mctags["ttx"].end());
    all_bkg.insert(mctags["vjets"].begin(), mctags["vjets"].end());
    all_bkg.insert(mctags["qcd"].begin(),   mctags["qcd"].end());
  }
  all_bkg.insert(mctags["other"].begin(), mctags["other"].end());



  //    Baseline definitions
  //-----------------------------------------
  NamedFunc wgt = "w_lumi*w_isr";//"weight"; 
  string common = "!lowDphiFix && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  string res_base = common+"&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  string boo_base = common+"&& nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300 && ht>600";
  boo_base += " && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";


  //    Useful boosted cut definitions
  //------------------------------------------
  string boo_SR = boo_base+"&& boostedRegionIdx>=4";

  //    Useful resolved cut definitions
  //------------------------------------------
  string c_2b = "nbm==2";
  string c_3b = "nbm==3&&nbl==3";
  string c_4b = "nbm>=3&&nbl>=4";
  string hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
  string sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  vector<string> res_abcd = {c_2b+"&&"+sbd, c_3b+"&&"+sbd, c_4b+"&&"+sbd, 
                             c_2b+"&&"+hig, c_3b+"&&"+hig, c_4b+"&&"+hig};

  string res_SR = res_base+"&& nbm>=3 &&"+hig;

  //    Define processes, including intersections
  //--------------------------------------------------
  string anaA(res_view?"Res":"Boo"), anaB(res_view?"Boo":"Res");
  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  auto proc_bkg_anaA = Process::MakeShared<Baby_pico>("MC bkg", 
    Process::Type::background, kGreen+1,attach_folder(foldermc, all_bkg),
    base_filters && "stitch" && (res_view ? res_base : boo_base));
  auto proc_bkg_anaA_anaB_crsr = Process::MakeShared<Baby_pico>("MC bkg $\\cap$ "+anaB+" (CR+SR)", 
    Process::Type::background, kGreen+1, attach_folder(foldermc, all_bkg),
    base_filters && "stitch &&"+res_base+"&&"+boo_base);
  auto proc_bkg_anaA_anaB_sr = Process::MakeShared<Baby_pico>("MC bkg $\\cap$ "+anaB+" SR only", 
    Process::Type::background, kGreen+1, attach_folder(foldermc, all_bkg),
    base_filters && "stitch" && (res_view ? (res_base+"&&"+boo_SR) : (boo_base+"&&"+res_SR)));
  auto proc_sig_anaA = Process::MakeShared<Baby_pico>("TChiHH(900,1)", 
    Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-900*.root"}, 
    base_filters && (res_view ? res_base : boo_base));
  auto proc_sig_anaA_anaB_sr = Process::MakeShared<Baby_pico>("TChiHH(900,1) $\\cap$ "+anaB+" SR only", 
    Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-900*.root"}, 
    base_filters && "stitch" && (res_view ? (res_base+"&&"+boo_SR) : (boo_base+"&&"+res_SR)));

  vector<shared_ptr<Process> > procs = {proc_bkg_anaA, proc_bkg_anaA_anaB_crsr, proc_bkg_anaA_anaB_sr, 
                                        proc_sig_anaA, proc_sig_anaA_anaB_sr};


  //    Useful binning definitions
  //------------------------------------------
  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");

  // assume common MET binning
  vector<TString> vc_met;
  if (res_view) {
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=450");
    vc_met.push_back("met>450");
  } else {
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=500");
    vc_met.push_back("met>500 && met<=700");
    vc_met.push_back("met>700");
  }   

  //        Define and fill tables
  //------------------------------------
  PlotMaker pm;
  vector<TableRow> tablerows;
  unsigned nbins(0), nheads(0);
  if (res_view) {
    nheads = vc_met.size()*vc_drmax.size();
    for (auto &imet: vc_met) {
      for (auto &idrmax: vc_drmax) {
        tablerows.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+"$"));
        for (auto &iabcd: res_abcd) {
          NamedFunc res_reg_ = imet+"&&"+idrmax+"&&"+iabcd;
          tablerows.push_back(TableRow("$"+CodeToLatex(iabcd)+"$", res_reg_, 0,0, wgt));
          nbins++;
        }
      }
    }
  } else { //boosted priority
    nheads = vc_met.size();
    for (auto &imet: vc_met) {
      tablerows.push_back(TableRow("$"+CodeToLatex(imet.Data())+"$"));
      for (unsigned ireg(0); ireg<6; ireg++) {
        NamedFunc boo_reg_ = imet && "boostedRegionIdx=="+to_string(ireg);
        tablerows.push_back(TableRow(to_string(ireg), boo_reg_, 0,0, wgt));
        nbins++;
      }
    }
  }
  TString tabname = "table_"; tabname += (res_view ? "resolved" : "boosted");
  pm.Push<Table>(tabname.Data(), tablerows, procs, 0);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //     Retrieve the relevant yields
  //-----------------------------------------
  vector<vector<GammaParams>> yields; 
  Table* itable = static_cast<Table*>(pm.Figures()[0].get());
  yields.push_back(itable->Yield(proc_bkg_anaA.get(), lumi)); 
  yields.push_back(itable->Yield(proc_bkg_anaA_anaB_crsr.get(), lumi)); 
  yields.push_back(itable->Yield(proc_bkg_anaA_anaB_sr.get(), lumi)); 
  yields.push_back(itable->Yield(proc_sig_anaA.get(), lumi)); 
  yields.push_back(itable->Yield(proc_sig_anaA_anaB_sr.get(), lumi)); 

  //      Make the plot
  //-------------------------------------
  PlotOpt opts("txt/plot_styles.txt", "Bins");
  if (!res_view) opts.BottomMargin(0.13);
  setPlotStyle(opts);

  vector<int> colors = {EColor::kBlack, EColor::kGray+2, EColor::kAzure+1, EColor::kRed, EColor::kMagenta};
  vector<int> fills = {3003, 3004, 3005, 0, 0};

  TCanvas can("can","");
  can.SetFillStyle(4000);

  TLegend legA(0.2, 0.75, .5, .88,"","NDC");
  legA.SetFillStyle(0); legA.SetBorderSize(0); legA.SetTextSize(opts.LegendEntryHeight()); legA.SetTextFont(42);
  TLegend legB(0.5, 0.75, .9, .88,"","NDC");
  legB.SetFillStyle(0); legB.SetBorderSize(0); legB.SetTextSize(opts.LegendEntryHeight()); legB.SetTextFont(42);

  vector<TH1D> histo;
  TString htitle = res_view? "Resolved regions":"Boosted regions";
  for (unsigned ihist(0); ihist<yields.size(); ihist++) {
    histo.push_back(TH1D("", htitle+"; ; Events", nbins, 0.5, nbins+0.5));
  }
  histo[0].SetTitleOffset(0.6,"y");
  histo[0].SetLabelOffset(10,"x");
  histo[0].SetTitleSize(0.06,"y");
  histo[0].SetMinimum(0.1);
  histo[0].GetYaxis()->CenterTitle(true);
  for (unsigned ihist(0); ihist<yields.size(); ihist++) {
    TString label = procs[ihist].get()->name_;
    label.ReplaceAll("\\cap", "#cap").ReplaceAll("$", "");
    for (unsigned ibin(0); ibin<nbins; ibin++){
      int iyield = ibin+ibin/(nbins/nheads)+1; //skip header rows
      // cout<<iyield<<"\t"<<yields[ihist][iyield]<<endl;
      histo[ihist].SetBinContent(ibin+1, yields[ihist][iyield].Yield());
      histo[ihist].SetBinError(ibin+1, yields[ihist][iyield].Uncertainty());   
    }
    histo[ihist].SetLineColor(colors[ihist]);
    histo[ihist].SetFillColor(colors[ihist]);
    if (label.Contains("TChi")) {
      if (label.Contains(anaB)) histo[ihist].SetLineStyle(2);
      histo[ihist].SetLineWidth(2);
      histo[ihist].SetFillStyle(0);
    }
    if (ihist==0) {
      histo[0].SetMaximum(histo[0].GetMaximum()*15);
      histo[ihist].SetFillStyle(3003);
    }
    TString draw_opt = "hist";
    histo[ihist].Draw(draw_opt+(ihist==0?"":"same"));
    if (label.Contains(anaB))
      legB.AddEntry(&histo[ihist], label, label.Contains("TChi") ? "L":"LF");
    else
      legA.AddEntry(&histo[ihist], label, label.Contains("TChi") ? "L":"LF");
  }
  legA.Draw("same");
  legB.Draw("same");

  //       Adding annotations
  //-------------------------------------
  unsigned nmet = res_view ? nbins/12 : nbins/6;
  unsigned nbin_per_met = res_view ? 12 : 6;
  unsigned nbin_per_drmax = 6;
  float space = (1-opts.LeftMargin() - opts.RightMargin())/nmet;
  TLatex latex; latex.SetNDC();
  latex.SetTextAlign(21);
  TLine line; line.SetNDC();
  line.SetLineStyle(3); line.SetLineWidth(2);
  for (unsigned i(0); i<nmet; i++) { 
    TString metlabel = vc_met[i];
    // metlatex.ReplaceAll("&&","-").ReplaceAll("met","").ReplaceAll("<=","").ReplaceAll(">","");
    latex.SetTextSize(0.045);
    float met_pos_ = opts.LeftMargin()+(i+0.5)*space;
    latex.SetTextColor(kBlack);
    latex.DrawLatex(met_pos_, 0.05, CodeToRootTex(metlabel.Data()).c_str());
    if (i>0) {
      line.SetLineColor(kBlack);
      line.DrawLine(i*nbin_per_met+0.5, 0, i*nbin_per_met+0.5, histo[0].GetMaximum()/15);
    }
    if (res_view){
      for (unsigned j(0); j<vc_drmax.size();j++){
        float pos_ = met_pos_ + (j-0.5)*space/2;
        latex.SetTextColor(j ? kOrange:kOrange-3);
        latex.DrawLatex(pos_, 0.12, CodeToRootTex(vc_drmax[j].Data()).c_str());
        if (j>0) {
          line.SetLineColor(j ? kOrange:kOrange-3);
          line.DrawLine(i*nbin_per_met+j*nbin_per_drmax+0.5, 0,
                        i*nbin_per_met+j*nbin_per_drmax+0.5, histo[0].GetMaximum()/30);
        }
      }
    }
  }
  latex.SetTextSize(0.045);
  // label.DrawLatex(nbins/2, 0.0005, "p_{T}^{miss} [GeV]");

  TString recotag = (res_view ? "resolved" : "boosted");
  TString pdfname = tag+"_"+recotag+".pdf";
  if (tag=="") pdfname = recotag+".pdf";
  can.SaveAs(pdfname);
  can.SetLogy(); 
  can.SaveAs(pdfname.ReplaceAll(".pdf","_log.pdf"));
  cout<<" open "<<pdfname<<endl;

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"no_signal", no_argument, 0, 'n'},    
      {"luminosity", required_argument, 0, 'l'},    
      {"mass_points", required_argument, 0, 'p'},
      {"tag", required_argument, 0, 't'},
      {"boo", no_argument, 0, 0},
      {"quick", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "nl:p:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      doSignal = false;
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 'p': 
      mass_points_string = optarg; 
      break;
    case 't': 
      mass_points_string = optarg; 
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "boo"){
        res_view = false;
      }else if(optname == "quick"){
        quick = true;
      } else {
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
