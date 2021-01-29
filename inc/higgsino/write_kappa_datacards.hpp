#ifndef H_HIG_WRITE_DATACARDS_LUMINOSITY
#define H_HIG_WRITE_DATACARDS_LUMINOSITY

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <map>

#include "core/gamma_params.hpp"
#include "core/named_func.hpp"
#include "core/table_row.hpp"

namespace HigWriteDataCards {
  void calculateKappas(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins, std::map<std::string, double> & kappas, std::map<std::string, double> & kappa_uncertainties);
  void setControlSystematics(std::map<std::string, std::map< std::string, float> > & controlSystematics);
  void writeDataCardHeader(std::vector<std::pair<std::string, std::string> > sampleBins, std::ofstream & cardFile);
  void setDataCardObserved(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::string const & dataTag, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardSignalBackground(std::string const & processName, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardSignalStatistics(std::string const & processName, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardSignalSystematics(std::string const & processName, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues, std::vector<std::pair<std::string, std::vector<NamedFunc>>> systematics_vector);
  void setDataCardControlSystematics(std::map<std::string, std::map<std::string, float> > controlSystematics ,std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues);
  // Counts number of "sig" to find tags for equations. Uses "bkg" tag. Relies on "xname_yname_dimensionname" tag 
  int countSubstring(const std::string& str, const std::string& sub);
  void setDataCardBackground(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::string const & mcTag, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardKappa(std::map<std::string, double > & kappas, std::map<std::string, double > & kappa_uncertainties, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins, std::vector<std::vector<std::string> > & tableValues);
  void writeTableValues(std::vector<std::vector<std::string> > & tableValues, std::ofstream & cardFile, bool alignLeft=false);
  void setRow(std::vector<std::string> & row, std::string const & value);
  // dimensionFile should be in "dimensioName dimensionCut" format
  void readDimensionFile(std::string const & dimensionFilePath, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins);
  void GetOptions(int argc, char *argv[]);
}

#endif
