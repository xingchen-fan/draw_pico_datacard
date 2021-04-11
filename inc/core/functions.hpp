#ifndef H_FUNCTIONS
#define H_FUNCTIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
 
namespace Functions{
  
  extern const NamedFunc hem_veto;
  extern const NamedFunc lowDphiFix;
  extern const NamedFunc boostedRegionIdx;
  extern const NamedFunc amBoostedRegionIdx;
  extern const NamedFunc SampleType;
  extern const NamedFunc ntrub;
  extern const NamedFunc w_pileup;
  extern const NamedFunc w_syst_pileup_up;
  extern const NamedFunc w_syst_pileup_down;

  bool IsGoodJet(const Baby &b, std::size_t ijet);

  // these are used for doing systematic variations of weights
  enum class Variation{central, up, down};
  NamedFunc MismeasurementCorrection(const std::string &file_path,
                                     const std::string &mismeas_scenario,
                                     Variation variation = Variation::central);
  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario);

  extern const NamedFunc boostSignalRegion;
  extern const NamedFunc boostControlRegion;
}

#endif
