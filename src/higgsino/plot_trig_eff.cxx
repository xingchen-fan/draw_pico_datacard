#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TBox.h"
#include "TFile.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "higgsino/hig_functions.hpp"

using namespace std;

double trig_eff(double ht, double met);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baby //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  PlotOpt opts("txt/plot_styles.txt", "Eff2D");

  setPlotStyle(opts);
  gStyle->SetPaintTextFormat("5.2f");

  TCanvas can("can","");


  TH2D hist = TH2D("","",30,150,300,5,0,1000);

  for (int iht(0); iht<hist.GetNbinsY(); iht++) {
    for (int imet(0); imet<hist.GetNbinsX(); imet++) {
      hist.SetBinContent(imet+1, iht+1, trig_eff(iht*200+1, imet*5+151));
      cout<<"MET = "<<imet*5+1<<",HT = "<<iht*200+1<<" = "<<trig_eff(iht*200+1, imet*5+151)<<endl;
    }
  }
  ////// Setting additional histogram style
  hist.GetXaxis()->SetLabelOffset(0.012);
  hist.GetXaxis()->SetTitle("MET [GeV]");

  hist.GetYaxis()->SetTitle("H_{T} [GeV]");

  hist.GetZaxis()->SetLabelSize(0.04);
  hist.GetZaxis()->SetRangeUser(0.,1.);
  hist.GetZaxis()->SetTitleSize(0.045);
  hist.GetZaxis()->CenterTitle(true);
  hist.GetZaxis()->SetTitle("Efficiency");

  ////// Dividing and plot histogram
  hist.Draw("COLZTEXT");

  TString hname = "trig_eff.pdf";

  can.SaveAs(hname);
  cout<<endl<<" open "<<hname<<endl<<endl;

  hname.ReplaceAll(".pdf", ".root");
  hname.ReplaceAll("plots/","");
  TFile file(hname, "recreate");
  file.cd();
  hname.ReplaceAll(".root","");
  hist.Write(hname);
  file.Close();

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

double trig_eff(double ht, double met){
  if(ht>   0 && ht<= 200 && met> 150 && met<= 155) return 0.532;
  else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) return 0.612;
  else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) return 0.589;
  else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) return 0.515;
  else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) return 0.588;
  else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) return 0.591;
  else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) return 0.684;
  else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) return 0.678;
  else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) return 0.537;
  else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) return 0.511;
  else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) return 0.619;
  else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) return 0.727;
  else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) return 0.699;
  else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) return 0.690;
  else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) return 0.568;
  else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) return 0.678;
  else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) return 0.769;
  else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) return 0.732;
  else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) return 0.609;
  else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) return 0.685;
  else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) return 0.670;
  else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) return 0.811;
  else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) return 0.779;
  else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) return 0.736;
  else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) return 0.663;
  else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) return 0.730;
  else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) return 0.838;
  else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) return 0.820;
  else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) return 0.819;
  else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) return 0.736;
  else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) return 0.745;
  else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) return 0.874;
  else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) return 0.848;
  else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) return 0.869;
  else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) return 0.759;
  else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) return 0.777;
  else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) return 0.903;
  else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) return 0.850;
  else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) return 0.839;
  else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) return 0.847;
  else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) return 0.792;
  else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) return 0.907;
  else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) return 0.884;
  else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) return 0.870;
  else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) return 0.781;
  else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) return 0.757;
  else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) return 0.924;
  else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) return 0.921;
  else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) return 0.936;
  else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) return 0.803;
  else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) return 0.841;
  else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) return 0.949;
  else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) return 0.927;
  else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) return 0.894;
  else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) return 0.839;
  else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) return 0.850;
  else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) return 0.966;
  else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) return 0.952;
  else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) return 0.919;
  else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) return 0.959;
  else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) return 0.896;
  else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) return 0.973;
  else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) return 0.979;
  else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) return 0.956;
  else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) return 0.971;
  else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) return 0.844;
  else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) return 0.983;
  else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) return 0.976;
  else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) return 0.983;
  else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) return 0.942;
  else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) return 0.880;
  else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) return 0.985;
  else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) return 0.992;
  else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) return 0.989;
  else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) return 0.931;
  else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) return 0.915;
  else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) return 0.989;
  else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) return 0.992;
  else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) return 0.984;
  else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) return 0.965;
  else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) return 0.862;
  else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) return 0.992;
  else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) return 0.989;
  else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) return 0.963;
  else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) return 0.991;
  else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) return 0.744;
  else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) return 0.994;
  else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) return 0.996;
  else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) return 1.000;
  else return 0.987;
}
