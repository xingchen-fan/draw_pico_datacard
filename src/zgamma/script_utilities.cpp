//collection of things commonly used in plotting scripts
//maybe should be merged with higutilities.cpp

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
#include "zgamma/script_utilities.hpp"

namespace script_utilities { 

  ArgStruct get_options(int argc, char* argv[], std::string default_string_options) {
    //
    //initialize to defaults
    ArgStruct options;
    options.single_thread = false;
    options.year_string = "run2";
    options.string_options = default_string_options;
    options.tag = "";
    options.unblind = false;

    //parse arguments
    while(true){
      static struct option long_options[] = {
        {"single_thread", no_argument, 0, 's'},
        {"year", required_argument, 0, 'y'},
        {"unblind", no_argument, 0, 'u'},
        {"tag", required_argument, 0, 't'},
        {"string_options", required_argument, 0, 'o'},
        {0, 0, 0, 0}
      };

      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "sy:ut:o:", long_options, &option_index);
      if (opt == -1) break;

      std::string optname;
      switch(opt){
      case 's':
        options.single_thread = true;
        break;
      case 'u':
        options.unblind = true;
        break;
      case 'y':
        options.year_string = optarg;
        break;
      case 't':
        options.tag = optarg;
        break;
      case 'o':
        options.string_options = optarg;
        break;
      case 0:
        //handle cases with no short argument form
        std::cout << "Bad option! Found option name " << optname << std::endl;
        break;
      default:
        //std::cout << "Bad option! getopt_long returned character code 0" << opt << std::endl;
        printf("Bad option! getopt_long returned character code 0%o\n", opt);
        break;
      }
    }
    return options;
  }

  std::vector<PlotOpt> plot_lin(bool unblind) {
    PlotOpt lin_norm("txt/plot_styles.txt", "CMSPaper");
    lin_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .LegendColumns(3);
    if (unblind) {
      lin_norm.Stack(PlotOptTypes::StackType::data_norm);
      lin_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {lin_norm};
  }

  std::vector<PlotOpt> plot_log(bool unblind) {
    PlotOpt log_norm("txt/plot_styles.txt", "CMSPaper");
    log_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::signal_overlay)
        .LegendColumns(3);
    if (unblind) {
      log_norm.Stack(PlotOptTypes::StackType::data_norm);
      log_norm.Bottom(PlotOptTypes::BottomType::ratio);
    }
    return {log_norm};
  }

  std::vector<PlotOpt> plot_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::linear)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }

  std::vector<PlotOpt> plot_log_shapes() {
    PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
    shapes_norm.Title(PlotOptTypes::TitleType::info)   
        .Bottom(PlotOptTypes::BottomType::off)
        .YAxis(PlotOptTypes::YAxisType::log)
        .LogMinimum(.2)
        .Stack(PlotOptTypes::StackType::shapes)
        .LegendColumns(3);
    return {shapes_norm};
  }

}
