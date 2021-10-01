#ifndef H_SCRIPT_UTILITIES
#define H_SCRIPT_UTILITIES

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <getopt.h>
#include <unistd.h>

#include "TColor.h"
#include "TError.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/cross_sections.hpp"
#include "core/event_scan.hpp"
#include "core/functions.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

namespace script_utilities { 

  //struct for holding options passed via command line arguments
  struct ArgStruct {
    bool single_thread;
    std::string year_string;
    std::string string_options;
    std::string tag;
    bool unblind;
  };

  //functions
  
  //function to parse command line arguments and return an ArgStruct
  //generic blurb to add to beginning of scripts
  //Arguments
  // --single_thread (-s) to run single thread for debugging
  // --unblind (-u) to unblind
  // --year (-y) yearname to select data year (2016, 2017, 2018, run2)
  // --tag (-t) to add a tag to produced plots
  // --string_options (-o) to specify what to plot
  ArgStruct get_options(int argc, char* argv[], std::string default_string_options="");

  //returns a vector containing a single plotopt corresponding to a data title linear plot
  //(with a ratio plot if unblind==true)
  std::vector<PlotOpt> plot_lin(bool unblind=false);

  //returns a vector containing a single plotopt corresponding to a data title log plot
  //(with a ratio plot if unblind==true)
  std::vector<PlotOpt> plot_log(bool unblind=false);

  //returns a vector containing a single plotopt corresponding to a data title linear shapes plot
  std::vector<PlotOpt> plot_shapes();

  //returns a vector containing a single plotopt corresponding to a data title log shapes plot
  std::vector<PlotOpt> plot_log_shapes();
}

#endif
