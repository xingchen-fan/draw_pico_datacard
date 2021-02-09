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
#include <utility>

#include "TH1.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/write_kappa_datacards.hpp"

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
  ////samplePaths["mc_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/mc/merged_higmc_preselect/";
  ////samplePaths["signal_2016"] = baseFolder + "/cms29r0/pico/NanoAODv5/higgsino_eldorado/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  ////samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  //samplePaths["mc_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/mc/merged_higmc_preselect/";
  //samplePaths["signal_2016"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/mc/merged_higmc_preselect/";
  //samplePaths["signal_2017"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldtv3/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  //samplePaths["mc_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/mc/merged_higmc_preselect/";
  //samplePaths["signal_2018"] = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  samplePaths["mc_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2016/mc/merged_higmc_preselect/";
  samplePaths["signal_2016"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2016/SMS-TChiHH_2D/merged_higmc_preselect/";
  samplePaths["mc_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2017/mc/merged_higmc_preselect/";
  samplePaths["signal_2017"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2017/SMS-TChiHH_2D/merged_higmc_preselect/";
  samplePaths["mc_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2018/mc/merged_higmc_preselect/";
  samplePaths["signal_2018"] = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2018/SMS-TChiHH_2D/merged_higmc_preselect/";

  // massPoints = { {"1000","1"} }
  vector<pair<string, string> > massPoints;
  if (higgsino_model == "N1N2") {
    string signal_folder = samplePaths["signal_2016"];
    set<string> signal_files = Glob(signal_folder+"/*.root");
    string mLSP, mChi;
    for (string signal_file : signal_files) {
      HigUtilities::filenameToMassPoint(signal_file, mChi, mLSP);
      if (stoi(mChi)>800) continue; // Ingore points above 800 GeV
      massPoints.push_back({mChi, mLSP});
    }
  }
  else if (mass_points_string == "") HigUtilities::findMassPoints(samplePaths["signal_2016"], massPoints);
  else HigUtilities::parseMassPoints(mass_points_string, massPoints);

  //NamedFunc filters = HigUtilities::pass_2016;
  //NamedFunc filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";
  NamedFunc filters = Higfuncs::final_pass_filters;

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, filters && "stitch", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  //HigUtilities::setDataProcesses(years, samplePaths, filters&&Higfuncs::trig_hig>0., sampleProcesses);
  HigUtilities::setDataProcesses(years, samplePaths, filters, sampleProcesses);
  HigUtilities::setSignalProcesses(massPoints, years, samplePaths, filters, sampleProcesses);
  
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  //NamedFunc weight = Higfuncs::final_weight;
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;

  if (higgsino_model=="N1N2") weight *= HigUtilities::w_CNToN1N2;
  string baseline = "!low_dphi_met && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200 && met/mht<2 && met/met_calo<2&& weight<1.5";
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
      // Original
      dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      dimensionBins["met"].push_back({"met1", "met>200 && met<=300"});
      dimensionBins["met"].push_back({"met2", "met>300 && met<=400"});
      dimensionBins["met"].push_back({"met3", "met>400"});
      dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

      //// Slightly lower met binning
      //dimensionBins["met"].push_back({"met0", "met>150 && met<=200"});
      //dimensionBins["met"].push_back({"met1", "met>200 && met<=250"});
      //dimensionBins["met"].push_back({"met2", "met>250 && met<=325"});
      //dimensionBins["met"].push_back({"met3", "met>325"});
      //dimensionBins["drmax"].push_back({"drmax0", "hig_cand_drmax[0]<=1.1"});
      //dimensionBins["drmax"].push_back({"drmax1", "hig_cand_drmax[0]>1.1"});

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
  
  //vector holding systematic variations as a name and weights. May need to be extended to a class
  vector<pair<string, vector<NamedFunc>>> systematics_vector;
  systematics_vector.push_back(make_pair(string("lumi"),
      vector<NamedFunc>({weight*Higfuncs::wgt_syst_lumi_up, weight*Higfuncs::wgt_syst_lumi_down})));
  systematics_vector.push_back(make_pair(string("trigeff"),
      vector<NamedFunc>({Higfuncs::final_weight_notrgeff*Higfuncs::eff_higtrig_run2_syst_up, 
      Higfuncs::final_weight_notrgeff*Higfuncs::eff_higtrig_run2_syst_down})));
  systematics_vector.push_back(make_pair(string("fsjetid"),
      vector<NamedFunc>({weight*1.01, weight*0.99})));
  systematics_vector.push_back(make_pair(string("bctag"),
      vector<NamedFunc>({weight*"sys_bchig[0]/w_bhig", weight*"sys_bchig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("fs_bctag"),
      vector<NamedFunc>({weight*"sys_fs_bchig[0]/w_bhig", weight*"sys_fs_bchig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("udsgtag"),
      vector<NamedFunc>({weight*"sys_udsghig[0]/w_bhig", weight*"sys_udsghig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("fs_udsgtag"),
      vector<NamedFunc>({weight*"sys_fs_udsghig[0]/w_bhig", weight*"sys_fs_udsghig[1]/w_bhig"})));
  systematics_vector.push_back(make_pair(string("prefire"),
      vector<NamedFunc>({weight*"sys_prefire[0]/w_prefire", weight*"sys_prefire[1]/w_prefire"})));
  systematics_vector.push_back(make_pair(string("isr"),
      vector<NamedFunc>({weight*"sys_isr[0]/w_isr", weight*"sys_isr[1]/w_isr"})));
  systematics_vector.push_back(make_pair(string("fs_met"),
      vector<NamedFunc>({})));
  systematics_vector.push_back(make_pair(string("PUsig"),
      vector<NamedFunc>({weight*"(npv<=20)", weight*"(npv>=21)"})));
  //TODO: add the following systematics? (from RA4, Old HH+MET, and HH+MET AN)
  //PU weights
  //fs_btag
  //renorm & factorization scales
  //JECS and JER

  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  map<string, HigUtilities::HistInformation > histInfo;
  //HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);
  for (pair<string, vector<NamedFunc>> sys : systematics_vector) {
    for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
      string sys_name = sys.first+to_string(wgt_idx);
      HigUtilities::addBinCuts(sampleBins, baseline, sys.second[wgt_idx], 
                               "signal_"+sys_name, cutTable["signal"]);
      HigUtilities::addBinCuts(sampleBins, baseline, sys.second[wgt_idx], 
                               "signalGenMet_"+sys_name, HigUtilities::nom2genmet,
                               cutTable["signal"]);
    }
  }
  HigWriteDataCards::make_npv_plots("signal", histInfo);
  if (unblind) HigWriteDataCards::make_npv_plots("data", histInfo);
  else HigWriteDataCards::make_npv_plots("mc", histInfo);

  PlotMaker pm;
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  bool verbose = false;
  HigUtilities::makePlots(cutTable, histInfo, sampleProcesses, luminosity, pm, verbose);

  // fill mYields
  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  //HigUtilities::fillDataYields(pm, cutTable["data"], mYields);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields, true);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields);
  for (pair<string, vector<NamedFunc>> sys : systematics_vector) {
    for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
      string sys_name = sys.first+to_string(wgt_idx);
      HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, 
          "signal_"+sys_name, "signalGenMet_"+sys_name, 
          "signalAverageGenMet_"+sys_name, mYields);
    }
  }
  
  map<string, TH1D> h_sig_npv;
  TH1D* h_data_npv;
  HigWriteDataCards::fill_npv_hists(pm, histInfo, sampleProcesses, h_sig_npv, &h_data_npv);

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
      HigWriteDataCards::setDataCardSignalSystematics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues, systematics_vector, h_sig_npv, h_data_npv);
    } else {
      HigWriteDataCards::setDataCardSignalBackground(process->name_, "signal", mYields, sampleBins, tableValues);
      HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signal", mYields, sampleBins, tableValues);      
      HigWriteDataCards::setDataCardSignalSystematics(process->name_, "signal", mYields, sampleBins, tableValues, systematics_vector, h_sig_npv, h_data_npv);
    }
    HigWriteDataCards::setDataCardControlSystematics(controlSystematics, sampleBins, tableValues);
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

  delete h_data_npv;
  HigWriteDataCards::delete_hist_info(histInfo);

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
  }

  void setControlSystematics(map<string, map<string, float> > & controlSystematics) {
    controlSystematics["xsig0_ybkg_met0_drmax0"]["ttbar"] = 1.12; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["ttbar"] = 1.15;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["ttbar"] = 1.06;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["ttbar"] = 1.10;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["ttbar"] = 1.14;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["ttbar"] = 1.02;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["ttbar"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["ttbar"] = 1.11;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["ttbar"] = 1.05;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["ttbar"] = 1.06;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["ttbar"] = 1.09;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["ttbar"] = 1.01;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["ttbar"] = 1.05;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["vjets"] = 1.02; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["vjets"] = 1.04;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["vjets"] = 1.05;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["vjets"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["vjets"] = 1.09;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["vjets"] = 1.06;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["vjets"] = 1.01;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["vjets"] = 1.07;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["vjets"] = 1.08;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["vjets"] = 1.02;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["vjets"] = 1.00;

    controlSystematics["xsig0_ybkg_met0_drmax0"]["qcd"] = 1.00; 
    controlSystematics["xsig1_ybkg_met0_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met0_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met1_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax0"]["qcd"] = 1.02;
    controlSystematics["xsig0_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met1_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met2_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax0"]["qcd"] = 1.00;
    controlSystematics["xsig0_ybkg_met3_drmax1"]["qcd"] = 1.00;
    controlSystematics["xsig1_ybkg_met3_drmax1"]["qcd"] = 1.01;
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

  void make_npv_plots(string tag, map<string, HigUtilities::HistInformation> & hist_info) {
    //add npv plot with no cuts past higloose+filters
    PlotOpt plot_opt("txt/plot_styles.txt", "CMSPaper");
    plot_opt = plot_opt.Title(PlotOptTypes::TitleType::info)
      .Bottom(PlotOptTypes::BottomType::off)
      .Overflow(PlotOptTypes::OverflowType::none)
      .YAxis(PlotOptTypes::YAxisType::linear)
      .Stack(PlotOptTypes::StackType::lumi_shapes).LegendColumns(3);
    //if (tag=="mc") plot_opt = plot_opt.Stack(StackType::signal_overlay);
    HigUtilities::HistInformation this_hist_info;
    this_hist_info.axis_ = new Axis(100, -0.5, 99.5, "npv", "N_{pv}", {});
    this_hist_info.cut_ = new NamedFunc("1");
    this_hist_info.plot_opt_ = new PlotOpt(plot_opt);
    this_hist_info.figure_index = -1;
    hist_info[tag] = this_hist_info;
  }

  void fill_npv_hists(PlotMaker& pm, map<string, HigUtilities::HistInformation> & hist_info, map<string, vector<shared_ptr<Process>>> & sample_processes, map<string, TH1D> &h_mc_npv, TH1D **h_data_npv) {
    Hist1D * sig_h1d = static_cast<Hist1D*>(pm.Figures()[hist_info["signal"].figure_index].get());
    for (auto & process : sample_processes["signal"]) {
      Hist1D::SingleHist1D* sig_sh1d = static_cast<Hist1D::SingleHist1D*>(sig_h1d->GetComponent(process.get()));
      h_mc_npv[process->name_] = sig_sh1d->scaled_hist_;
      //h_mc_npv[process->name_].Scale(1.0/h_mc_npv[process->name_].Integral());
    }

    TH1D* h_ret_npv = new TH1D("h_ret_npv","(Pseudo)data N_{pv}", 100, -0.5, 99.5);

    string signal_name = "mc";
    if (unblind) {
      signal_name = "data";
    }
    Hist1D * data_h1d = static_cast<Hist1D*>(pm.Figures()[hist_info[signal_name].figure_index].get());
    for (auto & process : sample_processes[signal_name]) {
      Hist1D::SingleHist1D* data_sh1d = static_cast<Hist1D::SingleHist1D*>(data_h1d->GetComponent(process.get()));
      h_ret_npv->Add(&(data_sh1d->scaled_hist_));
    }
    //h_ret_npv->Scale(1.0/h_ret_npv->Integral());
    (*h_data_npv) = h_ret_npv;
  }

  void delete_hist_info(map<string, HigUtilities::HistInformation> & hist_info) {
    for (auto single_hist_info : hist_info) {
      string const & type = single_hist_info.first;
      delete hist_info[type].axis_;
      delete hist_info[type].cut_;
      delete hist_info[type].plot_opt_;
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

  void setDataCardSignalSystematics(string const & processName, string const & signalAverageGenMetTag, 
      map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, 
      vector<vector<string> > & tableValues, vector<pair<string, vector<NamedFunc>>> systematics_vector,
      map<string, TH1D> h_sig_npv, TH1D* h_data_npv)
  {

    vector<string> row(2+2*sampleBins.size());

    //get information from pileup histograms for pileup systematic
    TH1D this_h_sig_npv = h_sig_npv[processName];
    float avr_npv_highpu = 0.0;
    float total_entries_highpu = 0.0;
    float avr_npv_lowpu = 0.0;
    float total_entries_lowpu = 0.0;
    for (int npv_bin = 0; npv_bin < 20; npv_bin++) {
      avr_npv_lowpu += this_h_sig_npv.GetBinContent(npv_bin+1)*static_cast<float>(npv_bin);
      total_entries_lowpu += this_h_sig_npv.GetBinContent(npv_bin+1);
    }
    for (int npv_bin = 0; npv_bin < this_h_sig_npv.GetNbinsX(); npv_bin++) {
      avr_npv_highpu += this_h_sig_npv.GetBinContent(npv_bin+1)*static_cast<float>(npv_bin);
      total_entries_highpu += this_h_sig_npv.GetBinContent(npv_bin+1);
    }
    avr_npv_lowpu = avr_npv_lowpu/total_entries_lowpu;
    avr_npv_highpu = avr_npv_highpu/total_entries_highpu;
    float total_entries_sig = this_h_sig_npv.Integral();
    float total_entries_dat = h_data_npv->Integral();

    for (pair<string, vector<NamedFunc>> sys : systematics_vector) {
      setRow(row,"-");
      row[0] = sys.first;
      row[1] = "lnN";
      for (unsigned int bin_idx = 0; bin_idx < (sampleBins.size()); bin_idx++) {
        if (sys.first == "fs_met") {
          //reco MET vs gen MET systematic
          string reclabel = processName + "_signal_" + sampleBins[bin_idx].first;
          string genlabel = processName + "_signalGenMet_" + sampleBins[bin_idx].first;
          string avrlabel = processName + "_signalAverageGenMet_" + sampleBins[bin_idx].first;
          float recmet_value = mYields.at(reclabel).first.Yield();
          float avrmet_value = mYields.at(avrlabel).first.Yield();
          float genmet_value = mYields.at(genlabel).first.Yield();
          //arbitrary convention for sign
          row[2+2*bin_idx] = to_string((recmet_value-genmet_value)/2.0/avrmet_value+1.0);
        }
        else if (sys.first == "PUsig") {
          //this systematic only uses reco met
          string label_prefix = processName + "_signal_";
          string label_suffix = "_" + sampleBins[bin_idx].first;
          float lowpu_cut_yield = mYields.at(label_prefix+sys.first+"0"+label_suffix).first.Yield();
          float highpu_cut_yield = mYields.at(label_prefix+sys.first+"1"+label_suffix).first.Yield();
          float lowpu_eff = lowpu_cut_yield/total_entries_lowpu;
          float highpu_eff = highpu_cut_yield/total_entries_highpu;
          float m = (highpu_eff-lowpu_eff)/(avr_npv_highpu-avr_npv_lowpu);
          float b = (lowpu_eff*avr_npv_highpu-highpu_eff*avr_npv_lowpu)/(avr_npv_highpu-avr_npv_lowpu);
          float eff_data = 0, eff_mc = 0;
          for (int npv_bin = 0; npv_bin < h_data_npv->GetNbinsX(); npv_bin++) {
            //compare linear fit
            float fx = m*static_cast<float>(npv_bin)+b;
            eff_data += fx*h_data_npv->GetBinContent(npv_bin+1);
            eff_mc += fx*this_h_sig_npv.GetBinContent(npv_bin+1);
          }
          eff_data = eff_data/total_entries_dat;
          eff_mc = eff_mc/total_entries_sig;
          row[2+2*bin_idx] = to_string((eff_data-eff_mc)/eff_mc+1.0);
        }
        else {
          //normal up/down systematic
          string label = processName + "_" + signalAverageGenMetTag + "_" + sampleBins[bin_idx].first;
          float nominal_value = mYields.at(label).first.Yield();
          if (nominal_value == 0) {
            //no signal, poor stats
            row[2+2*bin_idx] = "2.00";
          }
          else {
            float max_variation = 0.0;
            float variation_sign = 1.0;
            bool up_variation = true;
            for (unsigned wgt_idx = 0; wgt_idx < sys.second.size(); wgt_idx++) {
              label = processName + "_" + signalAverageGenMetTag + "_" + sys.first + to_string(wgt_idx) + "_" + sampleBins[bin_idx].first;
              float variation_value = mYields.at(label).first.Yield();
              //float variation = fabs(variation_value-nominal_value)/nominal_value;
              float variation = variation_value/nominal_value-1.0;
              //take reference sign from the first variation
              if (up_variation == true) {
                if (variation < 0) variation_sign = -1.0;
              }
              up_variation = false;
              if (fabs(variation) > max_variation)
                max_variation = fabs(variation);
            }
            row[2+2*bin_idx] = to_string(variation_sign*max_variation+1.0);
          }
        }
      }
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
