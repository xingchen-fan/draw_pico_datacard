#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>
#include <dirent.h>

#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/palette.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/write_kappa_datacards.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
namespace 
{
  string outFolder = getenv("PWD");
  string mass_points_string = ""; // run over all the points, otherwise specify, e.g. "175_1,200_1"
  string years_string = "2016,2017,2018";
  float luminosity = 1.;
  //float luminosity = 35.9;
  string dimensionFilePath = "";
  bool unblind = false;
  string tag = "resolved";
  bool do_met_average = true;
  string higgsino_model = "N1N2";
}

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  //if (b.type()==106000) {
  //  return 35.9;
  //}
  if (b.SampleType()==2016){
    return weight*35.9;
  } else if (b.SampleType()==2017){
    return weight*41.5;
  } else {
    return weight*59.6;
  }
});

const NamedFunc min_jet_dphi("min_jet_dphi", [](const Baby &b) -> NamedFunc::ScalarType{
  float min_dphi = 4;
  for (size_t ijet = 0; ijet < (*b.jet_phi()).size(); ++ijet) {
    float dphi = fabs(TVector2::Phi_mpi_pi((*b.jet_phi())[ijet]-b.met_phi()));
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
});

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  HigWriteDataCards::GetOptions(argc, argv);

  string baseFolder = ""; string hostName = execute("echo $HOSTNAME");
  if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = "/net/cms29";
  //baseFolder = "";

  cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_sys_abcd!!"<<endl;
  gSystem->mkdir(outFolder.c_str(), kTRUE);

  set<int> years;
  HigUtilities::parseYears(years_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  map<string, string> samplePaths;
  //samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  samplePaths["mc_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/mc/merged_higmc_preselect/";
  samplePaths["signal_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  samplePaths["mc_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/mc/merged_higmc_preselect/";
  samplePaths["signal_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  samplePaths["mc_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/mc/merged_higmc_preselect/";
  samplePaths["signal_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  // massPoints = { {"1000","1"} }
  vector<pair<string, string> > massPoints;
  if (mass_points_string == "") HigUtilities::findMassPoints(samplePaths["signal_2016"], massPoints);
  else HigUtilities::parseMassPoints(mass_points_string, massPoints);

  //NamedFunc filters = HigUtilities::pass_2016;
  NamedFunc filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, filters && "stitch", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  //HigUtilities::setDataProcesses(years, samplePaths, filters&&Higfuncs::trig_hig>0., sampleProcesses);
  HigUtilities::setDataProcesses(years, samplePaths, filters, sampleProcesses);
  HigUtilities::setSignalProcesses(massPoints, years, samplePaths, filters, sampleProcesses);
  
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years*Functions::w_pileup;
  if (higgsino_model=="N1N2") weight *= HigUtilities::w_CNToN1N2;
  string baseline = "!low_dphi_met && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  if (tag=="resolved") baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  else if (tag=="boosted") {
    baseline += " && ht>600 && nfjet>1 && fjet_pt[0]>300 && fjet_pt[1]>300";
    baseline += " && fjet_msoftdrop[0]>50 && fjet_msoftdrop[0]<=250 && fjet_msoftdrop[1]>50 && fjet_msoftdrop[1]<=250";
  }

  // Set ABCDbins
  // "bkg" and "sig" should be kept for setDataCardBackground()
  map<string, string> xBins;
  if (tag=="resolved") {
    xBins["bkg"] = "nbm==2";
    xBins["sig0"] = "nbm==3&&nbl==3";
    xBins["sig1"] = "nbm>=3&&nbl>=4";
  } else if (tag=="boosted") {
    string j1bb = "((met<=300 && fjet_mva_hbb_btv[0]>0.6)||(met>300 && fjet_mva_hbb_btv[0]>0.3))";
    string j2bb = "((met<=300 && fjet_mva_hbb_btv[1]>0.6)||(met>300 && fjet_mva_hbb_btv[1]>0.3))";
    xBins["bkg"] = "(!"+j1bb+"&& !"+j2bb+")";
    xBins["sig0"] = "((!"+j1bb+"&&"+j2bb+")||("+j1bb+"&& !"+j2bb+"))";
    xBins["sig1"] = j1bb+"&&"+j2bb;
  }
  
  map<string, string> yBins; // Shares sideband
  if (tag=="resolved") {
    yBins["sig"] = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  } else if (tag=="boosted") {
    yBins["sig"] = "fjet_msoftdrop[0]>95 && fjet_msoftdrop[0]<=145 && fjet_msoftdrop[1]>95 && fjet_msoftdrop[1]<=145";
    yBins["bkg"] = "!("+yBins["sig"]+")";
  }
  
  map<string, vector<pair<string, string> > > dimensionBins;
  if (dimensionFilePath==""){
    if (tag=="resolved") {
      //// Original
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300 && met<=400"});
      //dimensionBins["met"].push_back({"met3", "met>400"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      // Slightly lower met binning
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=250"});
      dimensionBins["met"].push_back({"met2", "met>250 && met<=325"});
      dimensionBins["met"].push_back({"met3", "met>325"});
      dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// 3 met bin
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      //dimensionBins["met"].push_back({"met2", "met>300"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// mix binning
      //dimensionBins["mix"].push_back({"met0_drmax1", "met>150 && met<=200 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax1", "met>200 && met<=300 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax1", "met>300 && met<=400 &&hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met3_drmax1", "met>400 && hig_cand_drmax[0]>1.1"});
      //dimensionBins["mix"].push_back({"met0_drmax0", "met>150 && met<=200 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met1_drmax0", "met>200 && met<=300 &&hig_cand_drmax[0]<=1.1"});
      //dimensionBins["mix"].push_back({"met2_drmax0", "met>300 &&hig_cand_drmax[0]<=1.1"});

    } else {
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=500"});
      dimensionBins["met"].push_back({"met3", "met>500 && met<=700"});
      dimensionBins["met"].push_back({"met4", "met>700"});
    }
  } else {
    HigWriteDataCards::readDimensionFile(dimensionFilePath, dimensionBins);
    for (auto & dimension : dimensionBins){
      cout<<"dimension bins"<<endl;
      cout<<dimension.first<<endl;
      for (auto & entry : dimension.second){
        cout<<" "<<entry.first<<" "<<entry.second<<endl;
      }
    }
  }

  // sampleBins = { {label, cut} }
  vector<pair<string, string> > sampleBins;
  HigUtilities::setABCDBins(xBins, yBins, dimensionBins, sampleBins);
  //sampleBins = {{"test","1"}};

  // Set systematics somehow..
  // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = systematic value
  map<string, map<string, float> > controlSystematics;
  HigWriteDataCards::setControlSystematics(controlSystematics);
  // Consider weight for nonHH somehow..

  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  //HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);

  PlotMaker pm;
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  bool verbose = false;
  HigUtilities::makePlots(cutTable, sampleProcesses, luminosity, pm, verbose);

  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  //HigUtilities::fillDataYields(pm, cutTable["data"], mYields);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields);

  // Calculate kappas from mc
  // kappas[xsigX_ysigY_metA_drmaxB]
  map<string, double> kappas;
  // kappa_uncertainties[metX_drmaxY]
  map<string, double> kappa_uncertainties;
  HigWriteDataCards::calculateKappas(mYields, dimensionBins, kappas, kappa_uncertainties);
  //// For testing
  //for (auto & item : kappas) item.second = 2;
  //for (auto & item : kappa_uncertainties) item.second = 0.001;

  //cout<<"cut name: "<<cutTable["signal"].tableRows[5].cut_.Name()<<endl;

  for (auto & process : sampleProcesses["signal"]){
    cout<<"Process name: "<<process->name_<<endl;
    string model;
    int mGluino=0, mLSP=0;
    HigUtilities::getInfoFromProcessName(process->name_, model, mGluino, mLSP);
    TString outPath = outFolder+"/datacard-"+HigUtilities::setProcessNameLong(model, mGluino, mLSP)+"_"+years_string+"_"+tag+".txt";
    cout<<"open "<<outPath<<endl;
    ofstream cardFile(outPath);
    HigWriteDataCards::writeDataCardHeader(sampleBins,cardFile);

    vector<vector<string> > tableValues;
    if (unblind) HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "mc", tableValues);
    if (do_met_average) {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
    } else {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signal", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signal", mYields, sampleBins, tableValues);      
    }
    //HigWriteDataCards::setDataCardControlSystematics(controlSystematics, sampleBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile);
    tableValues.clear();
    HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile, true);
    tableValues.clear();
    HigWriteDataCards::setDataCardKappa(kappas, kappa_uncertainties, dimensionBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues, cardFile);
    tableValues.clear();
    //writeDataCardClosure()
    //writeDataCardUncertainties()
    //writeDataCardParam()
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

namespace HigWriteDataCards{
  void calculateKappas(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins, std::map<std::string, double> & kappas, std::map<std::string, double> & kappa_uncertainties) {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        GammaParams xsig_ysig = mYields.at("mc_xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xsig_ybkg = mYields.at("mc_xsig"+to_string(iXbin)+"_ybkg_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ysig = mYields.at("mc_xbkg_ysig_"+(dimension->second)[iBin].first).first;
        GammaParams xbkg_ybkg = mYields.at("mc_xbkg_ybkg_"+(dimension->second)[iBin].first).first;
        vector<vector<float> > entries = {{static_cast<float>(xbkg_ybkg.NEffective())}, {static_cast<float>(xsig_ybkg.NEffective())}, {static_cast<float>(xbkg_ysig.NEffective())}, {static_cast<float>(xsig_ysig.NEffective())}};
        vector<vector<float> > entry_weights = {{static_cast<float>(xbkg_ybkg.Weight())}, {static_cast<float>(xsig_ybkg.Weight())}, {static_cast<float>(xbkg_ysig.Weight())}, {static_cast<float>(xsig_ysig.Weight())}};
        vector<float> powers = {1, -1, -1, 1};
        float kappa_up = 0, kappa_down = 0;
        double kappa = calcKappa(entries, entry_weights, powers, kappa_down, kappa_up);
        double kappa_uncertainty = kappa_up>kappa_down ? kappa_up:kappa_down;
        //cout<<xsig_ysig.Yield()<<" "<<xsig_ybkg.Yield()<<" "<<xbkg_ysig.Yield()<<" "<<xbkg_ybkg.Yield()<<endl;
        //cout<<"  kappa: "<<kappa<<" stat. uncertainty: "<<kappa_uncertainty<<endl;
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        kappas[kappaName] = kappa;
        kappa_uncertainties[kappaName] = kappa_uncertainty;
      }
    }

    //for (unsigned iMet = 0; iMet < dimensionBins["met"].size(); iMet++) {
    //  for (unsigned iDrmax = 0 ; iDrmax < dimensionBins["drmax"].size(); iDrmax++) {
    //    for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
    //      GammaParams xsig_ysig = mYields.at("mc_xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax)).first;
    //      GammaParams xsig_ybkg = mYields.at("mc_xsig"+to_string(iXbin)+"_ybkg_met"+to_string(iMet)+"_drmax"+to_string(iDrmax)).first;
    //      GammaParams xbkg_ysig = mYields.at("mc_xbkg_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax)).first;
    //      GammaParams xbkg_ybkg = mYields.at("mc_xbkg_ybkg_met"+to_string(iMet)+"_drmax"+to_string(iDrmax)).first;
    //      vector<vector<float> > entries = {{static_cast<float>(xbkg_ybkg.NEffective())}, {static_cast<float>(xsig_ybkg.NEffective())}, {static_cast<float>(xbkg_ysig.NEffective())}, {static_cast<float>(xsig_ysig.NEffective())}};
    //      vector<vector<float> > entry_weights = {{static_cast<float>(xbkg_ybkg.Weight())}, {static_cast<float>(xsig_ybkg.Weight())}, {static_cast<float>(xbkg_ysig.Weight())}, {static_cast<float>(xsig_ysig.Weight())}};
    //      vector<float> powers = {1, -1, -1, 1};
    //      float kappa_up = 0, kappa_down = 0;
    //      double kappa = calcKappa(entries, entry_weights, powers, kappa_down, kappa_up);
    //      double kappa_uncertainty = kappa_up>kappa_down ? kappa_up:kappa_down;
    //      //cout<<xsig_ysig.Yield()<<" "<<xsig_ybkg.Yield()<<" "<<xbkg_ysig.Yield()<<" "<<xbkg_ybkg.Yield()<<endl;
    //      //cout<<"  kappa: "<<kappa<<" stat. uncertainty: "<<kappa_uncertainty<<endl;
    //      string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
    //      kappas[kappaName] = kappa;
    //      kappa_uncertainties[kappaName] = kappa_uncertainty;
    //    }
    //  }
    //}
  }

  void setControlSystematics(map<string, map<string, float> > & controlSystematics) {
    controlSystematics["xsig0_ybkg_met0_drmax0"]["ttbar"] = 1.09; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["ttbar"] = 1.19;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["ttbar"] = 1.05;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["ttbar"] = 1.15;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["ttbar"] = 1.08;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["ttbar"] = 1.20;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["ttbar"] = 1.04;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["ttbar"] = 1.15;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["ttbar"] = 1.05;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["ttbar"] = 1.04;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["ttbar"] = 1.12;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["ttbar"] = 1.10;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["ttbar"] = 1.03;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["ttbar"] = 1.09;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["vjets"] = 1.02; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["vjets"] = 1.05;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["vjets"] = 1.04;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["vjets"] = 1.14;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["vjets"] = 1.09;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["vjets"] = 1.10;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["vjets"] = 1.12;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["vjets"] = 1.01;
  }

  void readDimensionFile(std::string const & dimensionFilePath, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins){
    map<string, int> dimensionIndex;
    ifstream dimensionFile(dimensionFilePath.c_str());
    if (dimensionFile.is_open()){
      string dimension, dimensionCut;
      while (dimensionFile >> dimension >> dimensionCut){
        dimensionBins[dimension].push_back(std::make_pair(dimension+to_string(dimensionIndex[dimension]), dimensionCut));
        dimensionIndex[dimension]++;
      }
      dimensionFile.close();
    }
  }

  void writeDataCardHeader(vector<pair<string, string> > sampleBins, ofstream & cardFile)
  {
    cardFile<<"imax "<<sampleBins.size()<<"  number of channels\n";
    cardFile<<"jmax 1  number of backgrounds\n";
    cardFile<<"kmax *  number of nuisance parameters\n";
    cardFile<<"shapes * * FAKE\n\n";
  }

  void setDataCardObserved(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & dataTag, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "Observation";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())); // Integerize. Round down
      //row[iBin*2+3] = to_string(int(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield())+1); // Integerize. Round up
      //row[iBin*2+3] = RoundNumber(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield(), 0, 1).Data(); // Integerize. Round
      row[iBin*2+3] = to_string(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield());
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }

  void setDataCardSignalBackground(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues){
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = sampleBins[iBin].first;
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "sig";
      row[iBin*2+3] = "bkg";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      row[iBin*2+2] = "0";
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "rate";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin){
      string label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[iBin].first;
      row[iBin*2+2] = to_string(mYields.at(label).first.Yield());
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }
  
  void setDataCardSignalStatistics(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    string label;
  
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      setRow(row,"-");
      row[0] = "stat_"+sampleBins[rBin].first;
      row[1] = "lnN";
      label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[rBin].first;
      string uncertainty;
      if (mYields.at(label).first.Yield() == 0) uncertainty = "2.00";
      else uncertainty = to_string(1+mYields.at(label).first.Uncertainty()/mYields.at(label).first.Yield());
      row[2+2*rBin] = uncertainty;
      tableValues.push_back(row);
    }
  
    setRow(row,"");
    tableValues.push_back(row);
  }

  void setDataCardControlSystematics(map<string, map<string, float> > controlSystematics, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues) 
  {
    // controlSystematics[xsig0_ybkg_met0_drmax0]["ttbar"] = 1.03

    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    
    // controlSystematicNames[0] = ttbar
    vector<string> controlSystematicNames;
    for (auto it = controlSystematics.begin()->second.begin(); it != controlSystematics.begin()->second.end(); ++it) controlSystematicNames.push_back(it->first);

    // Set control sample systematics
    for (unsigned iControlSystematic = 0; iControlSystematic < controlSystematicNames.size(); ++iControlSystematic) {
      string systematicName = controlSystematicNames[iControlSystematic]; // ex: ttbar
      setRow(row,"-");
      row[0] = "cr_"+systematicName;
      row[1] = "lnN";
      for (unsigned rBin=0;rBin<sampleBins.size();++rBin) {
        string regionName = sampleBins[rBin].first;
        if (controlSystematics.find(regionName) != controlSystematics.end()) {
          row[2+2*rBin+1] = to_string(controlSystematics[regionName][systematicName]);
        }
      }
      tableValues.push_back(row);
    }

    setRow(row, "");
    tableValues.push_back(row);
  }

  
  // returns count of non-overlapping occurrences of 'sub' in 'str'
  int countSubstring(const std::string& str, const std::string& sub)
  {
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos; offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
  }

  void setDataCardBackground(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & mcTag, vector<vector<string> > & tableValues)
  {
    // title + type + tag + tagType + (yield,equation), equation arguments
    vector<string> row(6);
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      row[0] = "rp_"+sampleBins[rBin].first;
      row[1] = "rateParam";
      row[2] = sampleBins[rBin].first;
      row[3] = "bkg";
      if (countSubstring(sampleBins[rBin].first,"sig") != 2) 
      {
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())); //Integerize. Round down
        //row[4] = to_string(int(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield())+1); //Integerize. Round up
        //row[4] = RoundNumber(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield(), 0, 1).Data(); //Integerize. Round
        row[4] = to_string(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield());
        // Infinity defined in https://root.cern.ch/doc/master/RooNumber_8cxx_source.html
        row[5] = "[0,1.0e30]";
      }
      else 
      {
        row[4] = "(@0*@1/@2)*@3";
        vector<string> xydimension;
        HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        vector<string> aNameVector = xydimension;
        aNameVector[0] = "xbkg";
        aNameVector[1] = "ybkg";
        string aName;
        HigUtilities::vectorStringToString(aNameVector, aName, "_");
        vector<string> bNameVector = xydimension;
        bNameVector[1] = "ybkg";
        string bName;
        HigUtilities::vectorStringToString(bNameVector, bName, "_");
        vector<string> cNameVector = xydimension;
        cNameVector[0] = "xbkg";
        string cName;
        HigUtilities::vectorStringToString(cNameVector, cName, "_");
        string kappaName;
        vector<string> kappaNameVector = xydimension;
        HigUtilities::vectorStringToString(kappaNameVector, kappaName, "_");
        //cout<<aName<<" "<<bName<<" "<<cName<<endl;
        row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName+",kappa_"+kappaName;
      }
      tableValues.push_back(row);
      setRow(row,"");
    }
  
  }

  void setDataCardKappa(map<string, double > & kappas, map<string, double > & kappa_uncertainties, map<string, vector<pair<string, string> > > & dimensionBins, vector<vector<string> > & tableValues)
  {
    map<string, vector<pair<string, string> > > t_dimensionBins = dimensionBins;
    // Combine dimensions to one list of dimensions.
    HigUtilities::combineDimensionBins(t_dimensionBins);

    auto dimension = t_dimensionBins.begin();
    vector<string> row(4);
    for (unsigned iBin =0; iBin < dimension->second.size(); ++iBin) {
      for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
        //string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
        string kappaName = "xsig"+to_string(iXbin)+"_ysig_"+(dimension->second)[iBin].first;
        row[0] = "kappa_"+kappaName;
        row[1] = "param";
        row[2] = to_string(kappas[kappaName]);
        row[3] = to_string(kappa_uncertainties[kappaName]);
        tableValues.push_back(row);
        setRow(row,"");
      }
    }
    cout<<"End"<<endl;

    //for (unsigned iMet = 0; iMet < dimensionBins["met"].size(); iMet++) {
    //  for (unsigned iDrmax = 0 ; iDrmax < dimensionBins["drmax"].size(); iDrmax++) {
    //    for (unsigned iXbin = 0; iXbin < 2; iXbin++) {
    //      string kappaName = "xsig"+to_string(iXbin)+"_ysig_met"+to_string(iMet)+"_drmax"+to_string(iDrmax);
    //      row[0] = "kappa_"+kappaName;
    //      row[1] = "param";
    //      row[2] = to_string(kappas[kappaName]);
    //      row[3] = to_string(kappa_uncertainties[kappaName]);
    //      tableValues.push_back(row);
    //      setRow(row,"");
    //    }
    //  }
    //}
  }
  
  void writeTableValues(vector<vector<string> > & tableValues, ofstream & cardFile, bool alignLeft)
  {
    // Set space values
    vector<unsigned> xSpace(tableValues.back().size());
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (xSpace[xIndex] < tableValues[yIndex][xIndex].size()) xSpace[xIndex] = tableValues[yIndex][xIndex].size();
      }
    }
  
    // Write values
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (alignLeft) cardFile << setw(xSpace[xIndex]) << std::left << tableValues[yIndex][xIndex] << " ";
        else cardFile << setw(xSpace[xIndex]) << tableValues[yIndex][xIndex] << " ";
      }
      cardFile<<endl;
    }
  }
  
  void setRow(vector<string> & row, string const & value)
  {
    for(auto & item : row) item = value;
  }

  void GetOptions(int argc, char *argv[])
  {
    string blah;
    while(true){
      static struct option long_options[] = {
        {"output_folder", required_argument, 0, 'o'},
        {"mass_points", required_argument, 0, 'p'},
        {"years", required_argument, 0, 'y'},
        {"luminosity", required_argument, 0, 'l'},
        {"dimension", required_argument, 0, 'd'},
        {"tag", required_argument, 0, 't'},
        {"unblind", no_argument, 0, 'u'},
        {"recomet", no_argument, 0, 0},
        {"higgsino_model", required_argument, 0, 'm'},
        {0, 0, 0, 0}
      };
  
      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "o:p:y:l:d:t:m:nu", long_options, &option_index);
      if( opt == -1) break;
  
      string optname;
      switch(opt){
        case 'o': outFolder = optarg; break;
        case 'p': mass_points_string = optarg; break;
        case 'y': years_string = optarg; break;
        case 'l': luminosity = atof(optarg); break;
        case 't': tag = optarg; break;
        case 'm': higgsino_model = optarg; break;
        case 'd': 
          dimensionFilePath = optarg; 
          if (!FileExists(dimensionFilePath)) 
          {
            cout<<"[Error] No file caled "<<dimensionFilePath<<endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 'u':
          unblind = true;
          break;
        case 0:
          optname = long_options[option_index].name;
          if(optname == "recomet"){
            do_met_average = false;
          }else{
            printf("Bad option! Found option name %s\n", optname.c_str());
          }
          break;
        default: 
          printf("Bad option! getopt_long returned character code 0%o\n", opt); 
          break;
      }
    }
  }

}
