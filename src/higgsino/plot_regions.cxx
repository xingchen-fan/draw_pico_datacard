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
  string mass_points_string = "700_1";//"225_1,400_1,700_1";
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
  string foldersig(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higtight/");
  string folderdata(bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/data/merged_higmc_higtight/");

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

  NamedFunc base_filters = HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  auto proc_bkg = Process::MakeShared<Baby_pico>("Background", Process::Type::background, kGreen+1,
                  attach_folder(foldermc, all_bkg),base_filters && "stitch");
  // auto proc_bkg = Process::MakeShared<Baby_pico>("Background", Process::Type::background, kGreen+1,
  //                 {foldermc+"*TTJets_*SingleLeptFromT_genM*root"},pass);
  // procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGreen+1,
  //                 attach_folder(foldermc, mctags["other"]),pass));
  // procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                 attach_folder(foldermc, mctags["qcd"]),pass)); 
  // procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                 attach_folder(foldermc,mctags["vjets"]),pass));
  // procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                 attach_folder(foldermc, mctags["ttx"]),pass));

  // if (doSignal) {
  //   vector<pair<string, string> > sigm;
  //   HigUtilities::parseMassPoints(mass_points_string, sigm);
  //   for (unsigned isig(0); isig<sigm.size(); isig++)
  //     procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig].first+",1)", 
  //       Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-"+sigm[isig].first+"*.root"}, pass_sig));
  // }
  auto proc_sig700 = Process::MakeShared<Baby_pico>("TChiHH(700,1)", 
        Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-700*.root"}, base_filters);
  auto proc_sig1000 = Process::MakeShared<Baby_pico>("TChiHH(1000,1)", 
        Process::Type::signal, 1, {foldersig+"*TChiHH_mChi-1000*.root"}, base_filters);

  vector<shared_ptr<Process> > procs = {proc_bkg, proc_sig700, proc_sig1000};

  NamedFunc lowDphi("lowDphi",[](const Baby &b) -> NamedFunc::ScalarType{
    bool lowDphi(false);
    if (b.jet_pt()->size()<4) return false;
    for (unsigned i(0); i<b.jet_pt()->size(); i++){
      if (i>3) break;
      if (i<=1 && b.jet_met_dphi()->at(i)<0.5) lowDphi = true;
      else if (i>=2 && b.jet_met_dphi()->at(i)<0.3) lowDphi = true;
    }
    return lowDphi;
  });

  NamedFunc nGoodFatJets("nGoodFatJets",[](const Baby &b) -> NamedFunc::ScalarType{ // AR = SR+CR
    int pass_pt_m_loose(0);
    for (unsigned i(0); i<b.fjet_pt()->size(); i++){
      if (i>1) break;
      if (b.fjet_pt()->at(i)<=300) continue;
      if (b.fjet_msoftdrop()->at(i)<=50 || b.fjet_msoftdrop()->at(i)>250) continue;
      pass_pt_m_loose++;
    } 
    return pass_pt_m_loose;
  });

  NamedFunc boostedRegionIdx("boostedRegionIdx",[](const Baby &b) -> NamedFunc::ScalarType{
    float bb_wp = 0.3;
    if (b.met()<300) bb_wp = 0.6;
    // 0 = D, 1 = B1, 2 = B2, 3 = C, 4 = A1, 5 = A2
    if (b.fjet_pt()->size()<2) return 0;
    bool j1mass = b.fjet_msoftdrop()->at(0)>85 && b.fjet_msoftdrop()->at(0)<=135;
    bool j2mass = b.fjet_msoftdrop()->at(1)>85 && b.fjet_msoftdrop()->at(1)<=135;
    bool j1bb = b.fjet_mva_hbb_btv()->at(0)>bb_wp;
    bool j2bb = b.fjet_mva_hbb_btv()->at(1)>bb_wp;

    if (j1mass && j2mass) {
      if (j1bb && j2bb) return 5; 
      else if (j1bb || j2bb) return 4; 
      else return 3;
    } else { 
      if (j1bb && j2bb) return 2;
      else if (j1bb || j2bb) return 1;
      else return 0;
    }
  });

  //    Baseline definitions
  //-----------------------------------------
  NamedFunc wgt = "w_lumi";//"weight"; 
  NamedFunc common = !lowDphi && "nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40";
  NamedFunc res_base = common && "nbt>=2 && njet>=4 && njet<=5 &&"+higtrim;
  NamedFunc boo_base = common && nGoodFatJets>=2 && "ht>300";

  //    Useful boosted cut definitions
  //------------------------------------------
  NamedFunc boo_SR = boo_base && boostedRegionIdx>=4;
  NamedFunc boo_CR = boo_base && boostedRegionIdx>=0. && boostedRegionIdx<=3;

  //    Useful resolved cut definitions
  //------------------------------------------
  string c_2b = "nbm==2";
  string c_3b = "nbm==3&&nbl==3";
  string c_4b = "nbm>=3&&nbl>=4";
  string hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
  string sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  vector<string> res_abcd = {c_2b+"&&"+sbd, c_3b+"&&"+sbd, c_4b+"&&"+sbd, 
                             c_2b+"&&"+hig, c_3b+"&&"+hig, c_4b+"&&"+hig};

  NamedFunc res_SR = res_base && ("nbm>=3 &&"+hig);
  NamedFunc res_CR = res_base && ("(nbm==2 &&"+hig+") || ("+sbd+")");

  //    Useful binning definitions
  //------------------------------------------
  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");

  // assume common MET binning
  vector<TString> vc_met;
  vc_met.push_back("met>150 && met<=200");
  vc_met.push_back("met>200 && met<=300");
  vc_met.push_back("met>300 && met<=450");
  if (res_view) {
    vc_met.push_back("met>450");
  } else {
    vc_met.push_back("met>450 && met<=700");
    vc_met.push_back("met>700");
  }   

  //        Define and fill tables
  //------------------------------------
  TString anaA(res_view?"resolved":"boosted"), anaB(res_view?"boosted":"resolved");
  PlotMaker pm;
  vector<vector<TableRow>> tablerows;
  tablerows.resize(3, vector<TableRow>());
  // Table 0: Analysis A
  // Table 1: Analysis B SR yields in the Analysis A regions
  // Table 2: Analysis B CR yields in the Analysis A regions
  unsigned nbins(0), nheads(0);
  if (res_view) {
    nheads = vc_met.size()*vc_drmax.size();
    for (auto &imet: vc_met) {
      for (auto &idrmax: vc_drmax) {
        for (auto &itab: tablerows)
          itab.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+"$"));
        for (auto &iabcd: res_abcd) {
          NamedFunc res_reg_ = res_base && (imet+"&&"+idrmax+"&&"+iabcd);
          tablerows[0].push_back(TableRow("$"+CodeToLatex(iabcd)+"$", res_reg_, 0,0, wgt));
          tablerows[1].push_back(TableRow("$"+CodeToLatex(iabcd)+"\\cap$ Boo (CR+SR)", res_reg_ && (boo_CR || boo_SR), 0,0, wgt));
          tablerows[2].push_back(TableRow("$"+CodeToLatex(iabcd)+"\\cap$ Boo SR", res_reg_ && boo_SR, 0,0, wgt));
          nbins++;
        }
      }
    }
  } else { //boosted priority
    nheads = vc_met.size();
    for (auto &imet: vc_met) {
      for (auto &itab: tablerows)
        itab.push_back(TableRow("$"+CodeToLatex(imet.Data())+"$"));
      for (unsigned ireg(0); ireg<6; ireg++) {
        NamedFunc boo_reg_ = boo_base && imet && boostedRegionIdx==ireg;
        tablerows[0].push_back(TableRow(to_string(ireg), boo_reg_, 0,0, wgt));
        tablerows[1].push_back(TableRow(to_string(ireg)+"$\\cap$ Res (CR+SR)", boo_reg_ && (res_CR || res_SR), 0,0, wgt));
        tablerows[2].push_back(TableRow(to_string(ireg)+"$\\cap$ Res SR", boo_reg_ && res_SR, 0,0, wgt));
        nbins++;
      }
    }
  }

  for (unsigned itab(0); itab<tablerows.size(); itab++) 
    pm.Push<Table>("table_"+to_string(itab)+"_"+anaA.Data(), tablerows[itab], procs, 0);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //     Retrieve the relevant yields
  //-----------------------------------------
  vector<vector<GammaParams>> yields; 
  vector<TString> labels;
  // 0: Analysis A background -> tab 0 bkg 
  labels.push_back("MC bkg.");
  Table* itable = static_cast<Table*>(pm.Figures()[0].get());
  yields.push_back(itable->BackgroundYield(lumi)); 
  // 1: Analysis B - CR bkg -> tab 2 bkg
  labels.push_back("#cap "+anaB+" (CR+SR), MC bkg.");
  itable = static_cast<Table*>(pm.Figures()[1].get());
  yields.push_back(itable->BackgroundYield(lumi)); 
  // 2: Analysis B - SR bkg -> tab 1 bkg 
  labels.push_back("#cap "+anaB+" SR only, MC bkg.");
  itable = static_cast<Table*>(pm.Figures()[2].get());
  yields.push_back(itable->BackgroundYield(lumi)); 
  // 3: Analysis A signal -> tab 0 sig
  labels.push_back(proc_sig700.get()->name_);
  itable = static_cast<Table*>(pm.Figures()[0].get()); 
  yields.push_back(itable->Yield(proc_sig700.get(), lumi)); 
  // 4: Analysis B - SR sig -> tab 1 sig
  labels.push_back("#cap "+anaB+" SR only, "+proc_sig700.get()->name_);
  itable = static_cast<Table*>(pm.Figures()[2].get());
  yields.push_back(itable->Yield(proc_sig700.get(), lumi)); 

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
  // histo[0].SetMinimum(miny);
  // histo[0].SetMaximum(maxy);
  // histo[0].GetXaxis()->SetLabelOffset(0.008);
  // histo[0].GetXaxis()->SetLabelSize(0.055);
  histo[0].SetTitleOffset(0.6,"y");
  histo[0].SetLabelOffset(10,"x");
  histo[0].SetTitleSize(0.06,"y");
  histo[0].SetMinimum(0.1);
  histo[0].GetYaxis()->CenterTitle(true);
  // histo[0].GetYaxis()->SetLabelOffset(0.4);
  // histo[0].GetYaxis()->SetLabelSize(0.06);
  // histo[0].GetXaxis()->SetTickLength(0);
  for (unsigned ihist(0); ihist<yields.size(); ihist++) {
    for (unsigned ibin(0); ibin<nbins; ibin++){
      int iyield = ibin+ibin/(nbins/nheads)+1; //skip header rows
      // cout<<iyield<<"\t"<<yields[ihist][iyield]<<endl;
      histo[ihist].SetBinContent(ibin+1, yields[ihist][iyield].Yield());
      histo[ihist].SetBinError(ibin+1, yields[ihist][iyield].Uncertainty());   
    }
    if (labels[ihist].Contains("TChi")) {
      if (labels[ihist].Contains(anaB)) histo[ihist].SetLineStyle(2);
      histo[ihist].SetLineWidth(2);
    }
    histo[ihist].SetLineColor(colors[ihist]);
    histo[ihist].SetFillColor(colors[ihist]);
    if (ihist==0) histo[ihist].SetFillStyle(3003);
    if (labels[ihist].Contains("TChi")) histo[ihist].SetFillStyle(0);
  }
  TString draw_opt = "hist";
  histo[0].SetMaximum(histo[0].GetMaximum()*15);
  for (unsigned ihist(0); ihist<yields.size(); ihist++){
    histo[ihist].Draw(draw_opt+(ihist==0?"":"same"));
    if (labels[ihist].Contains(anaB))
      legB.AddEntry(&histo[ihist], labels[ihist], labels[ihist].Contains("TChi") ? "L":"LF");
    else
      legA.AddEntry(&histo[ihist], labels[ihist], labels[ihist].Contains("TChi") ? "L":"LF");
  }
  legA.Draw("same");
  legB.Draw("same");

  //       Adding annotations
  //-------------------------------------
  unsigned nmet = res_view ? nbins/12 : nbins/6;
  unsigned nbin_per_met = res_view ? 12 : 6;
  unsigned nbin_per_drmax = 6;
  float space = (1-opts.LeftMargin() - opts.RightMargin())/nmet;
  TLatex label; label.SetNDC();
  label.SetTextAlign(21);
  TLine line; line.SetNDC();
  line.SetLineStyle(3); line.SetLineWidth(2);
  for (unsigned i(0); i<nmet; i++) { 
    TString metlabel = vc_met[i];
    // metlabel.ReplaceAll("&&","-").ReplaceAll("met","").ReplaceAll("<=","").ReplaceAll(">","");
    label.SetTextSize(0.045);
    float met_pos_ = opts.LeftMargin()+(i+0.5)*space;
    label.SetTextColor(kBlack);
    label.DrawLatex(met_pos_, 0.05, CodeToRootTex(metlabel.Data()).c_str());
    if (i>0) {
      line.SetLineColor(kBlack);
      line.DrawLine(i*nbin_per_met+0.5, 0, i*nbin_per_met+0.5, histo[0].GetMaximum()/15);
    }
    if (res_view){
      for (unsigned j(0); j<vc_drmax.size();j++){
        float pos_ = met_pos_ + (j-0.5)*space/2;
        label.SetTextColor(j ? kOrange:kOrange-3);
        label.DrawLatex(pos_, 0.12, CodeToRootTex(vc_drmax[j].Data()).c_str());
        if (j>0) {
          line.SetLineColor(j ? kOrange:kOrange-3);
          line.DrawLine(i*nbin_per_met+j*nbin_per_drmax+0.5, 0,
                        i*nbin_per_met+j*nbin_per_drmax+0.5, histo[0].GetMaximum()/30);
        }
      }
    }
  }
  label.SetTextSize(0.045);
  // label.DrawLatex(nbins/2, 0.0005, "p_{T}^{miss} [GeV]");

  TString pdfname = tag+"_"+anaA+".pdf";
  if (tag=="") pdfname = anaA+".pdf";
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
