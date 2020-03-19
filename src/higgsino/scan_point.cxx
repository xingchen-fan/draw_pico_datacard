#include "higgsino/scan_point.hpp"

#include <cstdlib>

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>

#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "core/utilities.hpp"
#include "core/cross_sections.hpp"

using namespace std;

namespace{
  string in_dir = "$PWD";
  string file_name = "";
  string model = "N1";
  bool do_signif = false;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  if(file_name == "") ERROR("Must supply an input file name");

  //// Parsing the gluino and LSP masses
  int mchi, mlsp;
  parseMasses(file_name, mchi, mlsp);
  double xsec, xsec_unc;
  if (model=="N1N2") xsec::higgsino2DCrossSection(mchi, xsec, xsec_unc);
  else xsec::higgsinoCrossSection(mchi, xsec, xsec_unc);
  string glu_lsp("mChi-"+to_string(mchi)+"_mLSP-"+to_string(mlsp));

  string workdir = in_dir+"/scan_point_"+model+"_"+glu_lsp+"/";
  gSystem->mkdir(workdir.c_str(), kTRUE);
 
  ostringstream command;
  string done = " < /dev/null &> /dev/null; ";
  done = "; ";
  command 
    << "ln -s $(readlink -f " << in_dir <<"/"<< file_name << ") " << workdir << done
    << "cd " << workdir << done
    << "combine -M AsymptoticLimits " << file_name << done;
  if(do_signif){
    command
      << "combine -M Significance --significance --expectSignal=1 --verbose=999999 --rMin=-10. --uncapped=1 " << file_name
      << " < /dev/null &> signif_obs.log; "
      << "combine -M Significance --significance --expectSignal=1 -t -1 --verbose=999999 --rMin=-10. --uncapped=1 --toysFreq " << file_name
      << " < /dev/null &> signif_exp.log; ";
  }
  command << flush;
  execute(command.str());
  
  string limits_file_name = workdir+"/higgsCombineTest.AsymptoticLimits.mH120.root";
  TFile limits_file(limits_file_name.c_str(), "read");
  if(!limits_file.IsOpen()) ERROR("Could not open limits file "+limits_file_name);
  TTree *tree = static_cast<TTree*>(limits_file.Get("limit"));
  if(tree == nullptr) ERROR("Could not get limits tree");
  double limit;
  tree->SetBranchAddress("limit", &limit);
  int num_entries = tree->GetEntries();
  if(num_entries != 6) ERROR("Expected 6 tree entries. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double exp_2down = limit;
  tree->GetEntry(1);
  double exp_down = limit;
  tree->GetEntry(2);
  double exp = limit;
  tree->GetEntry(3);
  double exp_up = limit;
  tree->GetEntry(4);
  double exp_2up = limit;
  tree->GetEntry(5);
  double obs = limit;
  limits_file.Close();

  double sig_obs, sig_exp;
  if(do_signif){
    sig_obs = GetSignif(workdir+"/signif_obs.log");
    sig_exp = GetSignif(workdir+"/signif_exp.log");
  }

  cout
    << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << mchi
    << ' ' << mlsp
    << ' ' << xsec
    << ' ' << xsec_unc
    << ' ' << obs
    << ' ' << exp
    << ' ' << exp_up
    << ' ' << exp_down
    << ' ' << exp_2up
    << ' ' << exp_2down;
  if(do_signif){
    cout
      << ' ' << sig_obs
      << ' ' << sig_exp;
  }
  cout << endl;

  string txtname(workdir+"/limits_"+model+"_"+glu_lsp+".txt");
  ofstream txtfile(txtname);
  txtfile
    << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << mchi
    << ' ' << mlsp
    << ' ' << xsec
    << ' ' << xsec_unc
    << ' ' << obs
    << ' ' << exp
    << ' ' << exp_up
    << ' ' << exp_down
    << ' ' << exp_2up
    << ' ' << exp_2down;
  if(do_signif){
    txtfile
      << ' ' << sig_obs
      << ' ' << sig_exp;
  }
  txtfile << endl;
}

double GetSignif(const string &filename){
  double signif = 0.;
  ifstream file(filename);
  string line;
  while(getline(file, line)){
    auto pos = line.find("Significance: ");
    if(pos != 0) continue;
    string val = line.substr(14);
    signif = stod(val);
  }
  return signif;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"in_dir", required_argument, 0, 'i'},
      {"file_name", required_argument, 0, 'f'},
      {"signif", required_argument, 0, 's'},
      {"model", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "i:f:m:s", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      file_name = optarg;
      break;
    case 'i':
      in_dir = optarg;
      break;
    case 's':
      do_signif = true;
      break;
    case 'm':
      model = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      cerr << "Bad option! getopt_long returned character code " << static_cast<int>(opt) << endl;
      break;
    }
  }
}
