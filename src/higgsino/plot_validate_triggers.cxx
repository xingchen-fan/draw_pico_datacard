#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

const NamedFunc get_0l_trigeff2016("get_0l_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 250 && met> 150 && met<= 155) {eff = 0.493894; errup = 0.0366069; errdown = 0.0365979;}
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.616091; errup = 0.040626; errdown = 0.0406625;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.630731; errup = 0.0417236; errdown = 0.0417726;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.654605; errup = 0.0436736; errdown = 0.0437577;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.658147; errup = 0.0459552; errdown = 0.0461807;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.62381; errup = 0.0528005; errdown = 0.0536457;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.640449; errup = 0.0678307; errdown = 0.070652;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.55665; errup = 0.0409086; errdown = 0.0410074;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.672944; errup = 0.0442055; errdown = 0.0442663;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.693822; errup = 0.0455462; errdown = 0.0456194;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.691021; errup = 0.0458807; errdown = 0.0459876;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.621951; errup = 0.0445279; errdown = 0.0447372;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.711628; errup = 0.0552757; errdown = 0.0565272;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.579545; errup = 0.0677244; errdown = 0.0694046;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.62419; errup = 0.032886; errdown = 0.0332858;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.717157; errup = 0.029367; errdown = 0.0295024;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.733951; errup = 0.0298227; errdown = 0.0299625;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.740271; errup = 0.0306371; errdown = 0.0308365;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.726916; errup = 0.0337854; errdown = 0.0343537;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.702247; errup = 0.0442877; errdown = 0.0462883;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.703704; errup = 0.0714117; errdown = 0.0794391;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.674352; errup = 0.0361548; errdown = 0.0369257;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.768562; errup = 0.0310229; errdown = 0.0311954;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.771899; errup = 0.0313409; errdown = 0.0315393;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.771282; errup = 0.0317878; errdown = 0.0320403;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.774947; errup = 0.0349064; errdown = 0.0356142;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.710983; errup = 0.0446227; errdown = 0.0467699;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.655738; errup = 0.0699558; errdown = 0.0753537;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.784946; errup = 0.0386148; errdown = 0.0400792;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.80334; errup = 0.0323132; errdown = 0.0325395;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.818095; errup = 0.0327434; errdown = 0.0329806;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.827869; errup = 0.0334729; errdown = 0.0337948;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.808786; errup = 0.0363654; errdown = 0.0373441;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.798658; errup = 0.0450919; errdown = 0.0484874;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.740741; errup = 0.0689347; errdown = 0.0783647;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.765625; errup = 0.0394371; errdown = 0.0410056;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.829294; errup = 0.0332461; errdown = 0.0335221;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.859719; errup = 0.0338653; errdown = 0.0341199;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.868649; errup = 0.0344868; errdown = 0.0348334;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.813299; errup = 0.0363376; errdown = 0.0373078;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.858156; errup = 0.0438374; errdown = 0.0478987;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.703704; errup = 0.0714117; errdown = 0.0794391;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.784884; errup = 0.0374544; errdown = 0.0406106;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.867628; errup = 0.0242716; errdown = 0.0248058;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.880851; errup = 0.0237118; errdown = 0.0241066;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.868759; errup = 0.0245513; errdown = 0.0251527;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.860724; errup = 0.0278056; errdown = 0.0292921;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.837209; errup = 0.0389778; errdown = 0.0440665;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.76; errup = 0.0663161; errdown = 0.0778814;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.759124; errup = 0.0421401; errdown = 0.0459079;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.89931; errup = 0.0243848; errdown = 0.024947;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.897311; errup = 0.0240833; errdown = 0.0245572;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.904271; errup = 0.024543; errdown = 0.0251564;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.880117; errup = 0.0276601; errdown = 0.0292663;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.857143; errup = 0.0387544; errdown = 0.0446173;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.826923; errup = 0.0578011; errdown = 0.0717659;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.873874; errup = 0.0383798; errdown = 0.0449809;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.890302; errup = 0.0248328; errdown = 0.0255271;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.914483; errup = 0.0243483; errdown = 0.0248965;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.932874; errup = 0.0247376; errdown = 0.0254601;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.886435; errup = 0.0279373; errdown = 0.0297191;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.85567; errup = 0.0419407; errdown = 0.0493746;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.84; errup = 0.0570871; errdown = 0.0722736;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.777778; errup = 0.0456141; errdown = 0.0508714;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.934489; errup = 0.0249785; errdown = 0.0258267;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.949025; errup = 0.0243608; errdown = 0.0249119;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.926621; errup = 0.0247748; errdown = 0.0255011;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.903226; errup = 0.0281424; errdown = 0.0302625;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.844444; errup = 0.0441926; errdown = 0.0520767;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.880952; errup = 0.0547064; errdown = 0.0755438;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.844961; errup = 0.0492065; errdown = 0.0533571;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.942308; errup = 0.0417775; errdown = 0.0419929;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.94065; errup = 0.0415114; errdown = 0.0416419;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.952336; errup = 0.0419764; errdown = 0.0421269;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.952813; errup = 0.0424617; errdown = 0.0428832;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.9; errup = 0.0449583; errdown = 0.0474417;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.836066; errup = 0.0609471; errdown = 0.071605;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.932584; errup = 0.0484419; errdown = 0.0556785;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.958926; errup = 0.0410742; errdown = 0.0428166;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.967811; errup = 0.0321888; errdown = 0.0426961;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.975419; errup = 0.024581; errdown = 0.0429401;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.960784; errup = 0.0392157; errdown = 0.0433312;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.945652; errup = 0.0444227; errdown = 0.0468189;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.940299; errup = 0.0498009; errdown = 0.0606091;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.955556; errup = 0.0444444; errdown = 0.0694684;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.970588; errup = 0.0294118; errdown = 0.0434143;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.984869; errup = 0.0151307; errdown = 0.0433006;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.983247; errup = 0.0167526; errdown = 0.043227;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.968668; errup = 0.0313316; errdown = 0.0437547;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.966216; errup = 0.0337838; errdown = 0.0475711;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.970149; errup = 0.0298507; errdown = 0.0568301;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.928571; errup = 0.0471241; errdown = 0.0871465;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.973134; errup = 0.0138291; errdown = 0.0161141;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.984686; errup = 0.0118413; errdown = 0.0126288;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.98609; errup = 0.0117738; errdown = 0.0125537;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.982143; errup = 0.0129215; errdown = 0.0150878;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.964286; errup = 0.0200601; errdown = 0.0293301;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.0384858;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0978447;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.985915; errup = 0.0132899; errdown = 0.0173362;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.992537; errup = 0.00746269; errdown = 0.0124079;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.983986; errup = 0.0120326; errdown = 0.0130306;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.983221; errup = 0.0130238; errdown = 0.0155758;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.019218;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.969697; errup = 0.0272627; errdown = 0.0671161;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0976169;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.994065; errup = 0.00593472; errdown = 0.0116793;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.990919; errup = 0.00908059; errdown = 0.009761;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.00892534;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.994774; errup = 0.00522648; errdown = 0.0100829;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.995708; errup = 0.00429185; errdown = 0.0131247;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.98913; errup = 0.0108696; errdown = 0.0260365;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.1852;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.992188; errup = 0.0078125; errdown = 0.0197536;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.998155; errup = 0.00184502; errdown = 0.00972075;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.997222; errup = 0.00277778; errdown = 0.00947604;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.00960046;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.993464; errup = 0.00653595; errdown = 0.0172321;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.986486; errup = 0.0135135; errdown = 0.0315945;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.458651;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0280804;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.997773; errup = 0.00222717; errdown = 0.00588578;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00371801;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00425727;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.987124; errup = 0.00757714; errdown = 0.0127016;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.990291; errup = 0.00854428; errdown = 0.022162;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2016("get_0l_fakemet_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff = 0.404605; errup = 0.0281502; errdown = 0.0281502;}
  else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff = 0.571429; errup = 0.0374088; errdown = 0.0374088;}
  else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff = 0.613636; errup = 0.0423806; errdown = 0.0423806;}
  else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff = 0.716418; errup = 0.0550662; errdown = 0.0550662;}
  else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff = 0.666667; errup = 0.0727393; errdown = 0.0727393;}
  else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff = 0.656716; errup = 0.0580067; errdown = 0.0580067;}
  else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff = 0.846154; errup = 0.0707589; errdown = 0.0707589;}
  else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff = 0.73913; errup = 0.0915605; errdown = 0.0915605;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.457086; errup = 0.0222559; errdown = 0.0222559;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.509589; errup = 0.0261664; errdown = 0.0261664;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.52795; errup = 0.0278203; errdown = 0.0278203;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.60084; errup = 0.0317442; errdown = 0.0317442;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.636364; errup = 0.0332746; errdown = 0.0332746;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.517647; errup = 0.0383244; errdown = 0.0383244;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.661417; errup = 0.0419922; errdown = 0.0419922;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.694444; errup = 0.0443253; errdown = 0.0443253;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.621053; errup = 0.0497728; errdown = 0.0497728;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.659091; errup = 0.0505302; errdown = 0.0505302;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.778626; errup = 0.0362737; errdown = 0.0362737;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.738636; errup = 0.0468378; errdown = 0.0468378;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.746269; errup = 0.0531615; errdown = 0.0531615;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.758621; errup = 0.0561886; errdown = 0.0561886;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.744186; errup = 0.0665378; errdown = 0.0665378;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff = 0.754098; errup = 0.0389865; errdown = 0.0389865;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.755556; errup = 0.0640644; errdown = 0.0640644;}
  else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff = 0.432544; errup = 0.0158991; errdown = 0.0158991;}
  else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff = 0.467866; errup = 0.0178888; errdown = 0.0178888;}
  else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff = 0.510542; errup = 0.0193994; errdown = 0.0193994;}
  else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff = 0.525773; errup = 0.0206981; errdown = 0.0206981;}
  else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff = 0.569476; errup = 0.0236322; errdown = 0.0236322;}
  else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff = 0.604369; errup = 0.0240906; errdown = 0.0240906;}
  else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff = 0.629747; errup = 0.0271637; errdown = 0.0271637;}
  else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff = 0.620567; errup = 0.028896; errdown = 0.028896;}
  else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff = 0.626609; errup = 0.0316885; errdown = 0.0316885;}
  else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff = 0.633663; errup = 0.0338995; errdown = 0.0338995;}
  else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff = 0.680645; errup = 0.0264799; errdown = 0.0264799;}
  else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff = 0.744589; errup = 0.0286928; errdown = 0.0286928;}
  else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff = 0.658163; errup = 0.0338804; errdown = 0.0338804;}
  else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff = 0.706667; errup = 0.0371743; errdown = 0.0371743;}
  else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff = 0.700855; errup = 0.0423314; errdown = 0.0423314;}
  else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff = 0.642173; errup = 0.0270951; errdown = 0.0270951;}
  else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff = 0.609626; errup = 0.035674; errdown = 0.035674;}
  else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff = 0.926471; errup = 0.0316513; errdown = 0.0316513;}
  else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff = 0.404868; errup = 0.0109406; errdown = 0.0109406;}
  else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff = 0.437424; errup = 0.0122272; errdown = 0.0122272;}
  else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff = 0.488841; errup = 0.0134125; errdown = 0.0134125;}
  else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff = 0.510135; errup = 0.014528; errdown = 0.014528;}
  else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff = 0.534884; errup = 0.0158603; errdown = 0.0158603;}
  else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff = 0.554398; errup = 0.0169094; errdown = 0.0169094;}
  else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff = 0.565807; errup = 0.0182575; errdown = 0.0182575;}
  else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff = 0.625; errup = 0.0189018; errdown = 0.0189018;}
  else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff = 0.590909; errup = 0.0209647; errdown = 0.0209647;}
  else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff = 0.655031; errup = 0.0215405; errdown = 0.0215405;}
  else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff = 0.616162; errup = 0.0172806; errdown = 0.0172806;}
  else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff = 0.635556; errup = 0.0185242; errdown = 0.0185242;}
  else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff = 0.642147; errup = 0.021374; errdown = 0.021374;}
  else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff = 0.663529; errup = 0.0229197; errdown = 0.0229197;}
  else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff = 0.617391; errup = 0.0261666; errdown = 0.0261666;}
  else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff = 0.606112; errup = 0.0201328; errdown = 0.0201328;}
  else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff = 0.602532; errup = 0.0246231; errdown = 0.0246231;}
  else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff = 0.704819; errup = 0.017701; errdown = 0.017701;}
  else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff = 0.955; errup = 0.00732931; errdown = 0.00732931;}
  else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff = 0.403256; errup = 0.00518037; errdown = 0.00518037;}
  else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff = 0.433503; errup = 0.00573216; errdown = 0.00573216;}
  else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff = 0.476213; errup = 0.00629978; errdown = 0.00629978;}
  else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff = 0.504982; errup = 0.00679124; errdown = 0.00679124;}
  else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff = 0.529961; errup = 0.00735406; errdown = 0.00735406;}
  else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff = 0.56639; errup = 0.00786624; errdown = 0.00786624;}
  else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff = 0.577905; errup = 0.00836627; errdown = 0.00836627;}
  else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff = 0.600321; errup = 0.00877644; errdown = 0.00877644;}
  else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff = 0.619677; errup = 0.00930155; errdown = 0.00930155;}
  else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff = 0.623563; errup = 0.0098163; errdown = 0.0098163;}
  else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff = 0.66071; errup = 0.0074079; errdown = 0.0074079;}
  else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff = 0.683814; errup = 0.00823533; errdown = 0.00823533;}
  else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff = 0.702176; errup = 0.00893583; errdown = 0.00893583;}
  else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff = 0.724342; errup = 0.00952245; errdown = 0.00952245;}
  else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff = 0.745802; errup = 0.0104774; errdown = 0.0104774;}
  else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff = 0.775077; errup = 0.00732399; errdown = 0.00732399;}
  else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff = 0.808443; errup = 0.00838429; errdown = 0.00838429;}
  else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff = 0.844642; errup = 0.00697529; errdown = 0.00697529;}
  else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff = 0.904545; errup = 0.00748778; errdown = 0.00748778;}
  else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff = 0.953261; errup = 0.00695909; errdown = 0.00695909;}
  else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff = 0.969849; errup = 0.00699865; errdown = 0.00699865;}
  else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.990217; errup = 0.00324488; errdown = 0.00324488;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.393245; errup = 0.00244599; errdown = 0.00244599;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.43438; errup = 0.00271824; errdown = 0.00271824;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.480271; errup = 0.00299093; errdown = 0.00299093;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.51571; errup = 0.00327196; errdown = 0.00327196;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.553577; errup = 0.00351606; errdown = 0.00351606;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.592697; errup = 0.00377356; errdown = 0.00377356;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.619873; errup = 0.00400325; errdown = 0.00400325;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.651545; errup = 0.00424702; errdown = 0.00424702;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.682976; errup = 0.00444573; errdown = 0.00444573;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.704728; errup = 0.00462451; errdown = 0.00462451;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.734168; errup = 0.00347758; errdown = 0.00347758;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.777321; errup = 0.00366278; errdown = 0.00366278;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.803227; errup = 0.00394339; errdown = 0.00394339;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.827035; errup = 0.00416252; errdown = 0.00416252;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.84477; errup = 0.00445204; errdown = 0.00445204;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.859084; errup = 0.00316684; errdown = 0.00316684;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.890933; errup = 0.00359947; errdown = 0.00359947;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.905381; errup = 0.00318308; errdown = 0.00318308;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.927529; errup = 0.00414893; errdown = 0.00414893;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.947605; errup = 0.00497749; errdown = 0.00497749;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.967402; errup = 0.00549868; errdown = 0.00549868;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.985313; errup = 0.0030399; errdown = 0.0030399;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.294527; errup = 0.00165956; errdown = 0.00165956;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.323816; errup = 0.00185338; errdown = 0.00185338;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.352023; errup = 0.00204674; errdown = 0.00204674;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.378594; errup = 0.00225679; errdown = 0.00225679;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.402409; errup = 0.00245637; errdown = 0.00245637;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.433744; errup = 0.00267275; errdown = 0.00267275;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.463049; errup = 0.00289395; errdown = 0.00289395;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.48537; errup = 0.00310318; errdown = 0.00310318;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.513993; errup = 0.00332597; errdown = 0.00332597;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.533194; errup = 0.00356519; errdown = 0.00356519;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.560839; errup = 0.00275238; errdown = 0.00275238;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.599586; errup = 0.00306054; errdown = 0.00306054;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.629718; errup = 0.00337431; errdown = 0.00337431;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.669613; errup = 0.00367228; errdown = 0.00367228;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.69982; errup = 0.00397398; errdown = 0.00397398;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.733677; errup = 0.00285148; errdown = 0.00285148;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.775371; errup = 0.00333625; errdown = 0.00333625;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.805159; errup = 0.00296283; errdown = 0.00296283;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.847872; errup = 0.00376673; errdown = 0.00376673;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.875209; errup = 0.00478206; errdown = 0.00478206;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.888848; errup = 0.00603014; errdown = 0.00603014;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.92654; errup = 0.00382917; errdown = 0.00382917;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1el_trigeff2016("get_1el_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.10209; errup = 0.00591574; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.508521; errup = 0.0133826; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.540469; errup = 0.00506853; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.197796; errup = 0.00954788; errdown = 0.00922175;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.522956; errup = 0.0172376; errdown = 0.0172907;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.592268; errup = 0.00609976; errdown = 0.00612797;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.343518; errup = 0.00991985; errdown = 0.0097879;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.632466; errup = 0.0136146; errdown = 0.0138232;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.676785; errup = 0.00477725; errdown = 0.00481403;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.51545; errup = 0.00513898; errdown = 0.00514221;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.734544; errup = 0.0056898; errdown = 0.00576763;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.756135; errup = 0.00212987; errdown = 0.00214247;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.76007; errup = 0.0128345; errdown = 0.0133055;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.85582; errup = 0.0129388; errdown = 0.0139281;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.860034; errup = 0.00541667; errdown = 0.00559399;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.851656; errup = 0.00672062; errdown = 0.00697494;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.884747; errup = 0.00665368; errdown = 0.00699398;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.90492; errup = 0.00282916; errdown = 0.00290529;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.207113; errup = 0.0202127; errdown = 0.0189198;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.564516; errup = 0.0268312; errdown = 0.0271955;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.564735; errup = 0.00789192; errdown = 0.00792435;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.371345; errup = 0.0279443; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.664093; errup = 0.0305449; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.651934; errup = 0.00869838; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.523629; errup = 0.0225912; errdown = 0.0226843;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.722334; errup = 0.0206066; errdown = 0.021545;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.728115; errup = 0.00644023; errdown = 0.0065357;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.656417; errup = 0.0101746; errdown = 0.0103171;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.792651; errup = 0.00856531; errdown = 0.00882788;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.802124; errup = 0.00282306; errdown = 0.00285346;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.829493; errup = 0.0261432; errdown = 0.0294446;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.852459; errup = 0.0231571; errdown = 0.0262989;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.879075; errup = 0.00753224; errdown = 0.00794474;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.866873; errup = 0.013536; errdown = 0.0147365;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.90665; errup = 0.0104923; errdown = 0.0115955;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.91716; errup = 0.00409646; errdown = 0.0042842;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.283465; errup = 0.0456264; errdown = 0.0418996;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.75; errup = 0.0330874; errdown = 0.036026;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.664792; errup = 0.00995901; errdown = 0.0101047;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.578512; errup = 0.0481009; errdown = 0.0495105;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.741722; errup = 0.0370708; errdown = 0.0405552;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.735936; errup = 0.0107074; errdown = 0.0109856;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.703883; errup = 0.033095; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.814545; errup = 0.0239949; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.78274; errup = 0.00762671; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.784091; errup = 0.0140963; errdown = 0.014768;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.871389; errup = 0.0103228; errdown = 0.0110463;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.849313; errup = 0.00324759; errdown = 0.00330547;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.870968; errup = 0.0354651; errdown = 0.0446107;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.966942; errup = 0.015755; errdown = 0.0253664;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.91099; errup = 0.00864424; errdown = 0.00942993;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.893333; errup = 0.0208625; errdown = 0.0247257;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.926509; errup = 0.013463; errdown = 0.0159179;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.944102; errup = 0.00452597; errdown = 0.00488294;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.5; errup = 0.0780973; errdown = 0.0780973;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.739726; errup = 0.0541361; errdown = 0.0615126;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.746437; errup = 0.0120965; errdown = 0.0124777;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.54; errup = 0.0786763; errdown = 0.080479;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.794872; errup = 0.0475795; errdown = 0.0561273;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.819343; errup = 0.0117716; errdown = 0.0123769;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.843478; errup = 0.0347781; errdown = 0.0414801;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.908397; errup = 0.0254458; errdown = 0.0325889;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.842222; errup = 0.00867254; errdown = 0.00906496;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.864662; errup = 0.0173909; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.926448; errup = 0.0103967; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.89891; errup = 0.00351003; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.966667; errup = 0.0275914; errdown = 0.072517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.942529; errup = 0.0246014; errdown = 0.0370083;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.9377; errup = 0.00971033; errdown = 0.011223;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.966292; errup = 0.0182896; errdown = 0.0317032;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.950495; errup = 0.0152082; errdown = 0.0203749;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.965256; errup = 0.00460922; errdown = 0.0052372;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.851852; errup = 0.0695101; errdown = 0.101731;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.936508; errup = 0.0301398; errdown = 0.0473598;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.853598; errup = 0.0126086; errdown = 0.013529;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.909091; errup = 0.0490632; errdown = 0.0805904;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.871069; errup = 0.0134541; errdown = 0.0146888;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.883333; errup = 0.0420083; errdown = 0.0571758;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.929825; errup = 0.0332829; errdown = 0.0520157;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.876818; errup = 0.00970247; errdown = 0.0103755;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.920354; errup = 0.0181519; errdown = 0.0223173;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.948276; errup = 0.0103093; errdown = 0.0124414;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.935926; errup = 0.00348682; errdown = 0.00366757;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.952381; errup = 0.0394264; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.981481; errup = 0.0153245; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.971074; errup = 0.00757748; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.981481; errup = 0.0153245; errdown = 0.0412969;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.951613; errup = 0.0190274; errdown = 0.0277831;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.97007; errup = 0.00505453; errdown = 0.00595382;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.857143; errup = 0.0917089; errdown = 0.158271;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0595223;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.925703; errup = 0.0118357; errdown = 0.0136893;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.932367; errup = 0.0124155; errdown = 0.0147054;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.810811; errup = 0.0670146; errdown = 0.0869581;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.980769; errup = 0.0159141; errdown = 0.0428344;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.950556; errup = 0.00764631; errdown = 0.00884972;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.972603; errup = 0.0130667; errdown = 0.0211322;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.977564; errup = 0.00823774; errdown = 0.0118755;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.960628; errup = 0.00329173; errdown = 0.00356596;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.115502;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0409781;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.975385; errup = 0.00847391; errdown = 0.0119156;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.990654; errup = 0.00773258; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988032; errup = 0.00390372; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.123222;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.95572; errup = 0.00884752; errdown = 0.0107028;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.115502;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.96875; errup = 0.025866; errdown = 0.0682225;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.974895; errup = 0.00709912; errdown = 0.00937376;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.102638;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0376284;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.982256; errup = 0.00466259; errdown = 0.00604897;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.993333; errup = 0.00551564; errdown = 0.0151622;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.990291; errup = 0.00527929; errdown = 0.00935364;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.981423; errup = 0.00219908; errdown = 0.00247001;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0419109;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.994505; errup = 0.00354816; errdown = 0.00720076;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0879414;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0180628;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.992408; errup = 0.00279574; errdown = 0.00406535;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.168149;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.108691;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.991718; errup = 0.00395933; errdown = 0.00649962;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.264229;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0839348;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.986143; errup = 0.00548298; errdown = 0.00818471;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.115502;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0559083;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.99409; errup = 0.00255052; errdown = 0.00397839;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0173807;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.996441; errup = 0.00294413; errdown = 0.00813542;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.994857; errup = 0.0011116; errdown = 0.00138053;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.458642;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0498539;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 0.995455; errup = 0.00293541; errdown = 0.00596358;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0149771;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998031; errup = 0.00127137; errdown = 0.0025904;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1mu_trigeff2016("get_1mu_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.309278; errup = 0.012177; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.896958; errup = 0.00960512; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.902498; errup = 0.00336926; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.433962; errup = 0.0171325; errdown = 0.016981;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.867868; errup = 0.0132858; errdown = 0.0144527;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.908151; errup = 0.00408863; errdown = 0.00425486;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.631429; errup = 0.0117375; errdown = 0.0118913;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.929134; errup = 0.00660721; errdown = 0.00719829;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.935745; errup = 0.0022462; errdown = 0.00232037;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.938377; errup = 0.00596572; errdown = 0.00652796;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.962165; errup = 0.00418374; errdown = 0.00465172;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.961258; errup = 0.00142423; errdown = 0.00147532;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.379603; errup = 0.0275657; errdown = 0.0268474;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.91954; errup = 0.0147031; errdown = 0.0173473;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.898394; errup = 0.00492584; errdown = 0.00514117;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.595238; errup = 0.0324355; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.91791; errup = 0.0169135; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.928719; errup = 0.00479442; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.734667; errup = 0.0164532; errdown = 0.0171043;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.927536; errup = 0.00906659; errdown = 0.0101678;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.94032; errup = 0.0026961; errdown = 0.00281232;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.948035; errup = 0.00793019; errdown = 0.00915617;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.965458; errup = 0.00537137; errdown = 0.00623961;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.963843; errup = 0.00173444; errdown = 0.00181629;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.553333; errup = 0.0433853; errdown = 0.0441515;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.930233; errup = 0.0194967; errdown = 0.0252262;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.915444; errup = 0.00580904; errdown = 0.00618079;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.675214; errup = 0.0457257; errdown = 0.0489983;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.951515; errup = 0.016605; errdown = 0.0230425;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.94313; errup = 0.00521512; errdown = 0.00568264;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.848397; errup = 0.019717; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.973843; errup = 0.00711092; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.950236; errup = 0.00288025; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.946067; errup = 0.0107415; errdown = 0.0129539;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.980488; errup = 0.00479992; errdown = 0.006119;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.967188; errup = 0.00193449; errdown = 0.00204772;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.671233; errup = 0.0587767; errdown = 0.0640049;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.914634; errup = 0.0309441; errdown = 0.0429513;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.933731; errup = 0.00643465; errdown = 0.00703905;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.775862; errup = 0.0574412; errdown = 0.068324;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.948454; errup = 0.0220862; errdown = 0.0333658;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.958838; errup = 0.00565631; errdown = 0.00644963;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.918478; errup = 0.0203311; errdown = 0.0254746;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.947368; errup = 0.0132382; errdown = 0.0167863;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.966378; errup = 0.00292606; errdown = 0.00318206;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.983264; errup = 0.00799298; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.985102; errup = 0.00513915; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.97674; errup = 0.00198722; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.828571; errup = 0.0658118; errdown = 0.0883336;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.955556; errup = 0.0286549; errdown = 0.0556317;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.970874; errup = 0.00551699; errdown = 0.00663155;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.973684; errup = 0.02178; errdown = 0.0579268;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.964286; errup = 0.0230346; errdown = 0.0451714;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.976055; errup = 0.00514512; errdown = 0.00635002;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.94382; errup = 0.0240536; errdown = 0.0362177;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.987952; errup = 0.00777825; errdown = 0.0156694;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.983293; errup = 0.00258303; errdown = 0.0030079;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.987261; errup = 0.0082239; errdown = 0.0165543;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.988604; errup = 0.00544618; errdown = 0.00891882;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.982685; errup = 0.00206541; errdown = 0.00232236;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0923495;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0683597;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.987118; errup = 0.0044459; errdown = 0.00629175;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0659133;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.978873; errup = 0.00598041; errdown = 0.00791061;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.952381; errup = 0.0258048; errdown = 0.044158;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0149771;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.987744; errup = 0.00258316; errdown = 0.00318481;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988764; errup = 0.00929678; errdown = 0.0253616;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.986425; errup = 0.0073789; errdown = 0.013028;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.994496; errup = 0.00135936; errdown = 0.00174222;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.123222;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.973684; errup = 0.02178; errdown = 0.0579268;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.982704; errup = 0.00511276; errdown = 0.00685994;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0512411;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.998138; errup = 0.00154055; errdown = 0.00426903;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0361508;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.991379; errup = 0.00713254; errdown = 0.019543;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.997409; errup = 0.00111857; errdown = 0.00174877;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.985294; errup = 0.0121686; errdown = 0.0330032;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.00811302;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.996846; errup = 0.000979628; errdown = 0.00134259;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 0.996764; errup = 0.00209004; errdown = 0.00425238;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.23126;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0559083;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.994709; errup = 0.0028782; errdown = 0.00511987;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0498539;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0172182;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.997996; errup = 0.000958848; errdown = 0.0015817;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0335184;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00783675;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.999413; errup = 0.000378822; errdown = 0.00077304;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2el_trigeff2016("get_2el_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.938144; errup = 0.0141214; errdown = 0.0141214;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.949468; errup = 0.0112961; errdown = 0.0112961;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.966443; errup = 0.00851778; errdown = 0.00851778;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.962617; errup = 0.0082014; errdown = 0.0082014;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.972632; errup = 0.00748604; errdown = 0.00748604;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.973626; errup = 0.00751234; errdown = 0.00751234;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.960688; errup = 0.00963289; errdown = 0.00963289;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.983213; errup = 0.00629125; errdown = 0.00629125;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.987531; errup = 0.00554136; errdown = 0.00554136;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.973988; errup = 0.00855701; errdown = 0.00855701;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.985876; errup = 0.00627181; errdown = 0.00627181;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.974522; errup = 0.00889224; errdown = 0.00889224;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.977199; errup = 0.00851926; errdown = 0.00851926;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.985294; errup = 0.00729868; errdown = 0.00729868;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.995914; errup = 0.000770677; errdown = 0.000770677;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2mu_trigeff2016("get_2mu_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.978495; errup = 0.0075211; errdown = 0.0075211;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.993151; errup = 0.0030526; errdown = 0.0030526;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.99209; errup = 0.000276502; errdown = 0.000276502;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_trigeff2017("get_0l_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 250 && met> 150 && met<= 155) {eff = 0.178698; errup = 0.0225437; errdown = 0.0220535;}
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.247156; errup = 0.0263155; errdown = 0.0262127;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.313043; errup = 0.0331531; errdown = 0.0330574;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.359584; errup = 0.0378066; errdown = 0.0377303;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.350832; errup = 0.0388109; errdown = 0.0386298;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.412698; errup = 0.0524354; errdown = 0.0519716;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.571429; errup = 0.0793748; errdown = 0.0806037;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.197836; errup = 0.0257522; errdown = 0.0251198;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.291738; errup = 0.0307667; errdown = 0.0306785;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.372247; errup = 0.0389705; errdown = 0.0389039;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.412712; errup = 0.0431359; errdown = 0.0430838;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.457831; errup = 0.0492426; errdown = 0.0491888;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.424; errup = 0.0534294; errdown = 0.0530271;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.42268; errup = 0.0695934; errdown = 0.0681992;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.208633; errup = 0.0216225; errdown = 0.0206788;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.350923; errup = 0.0223027; errdown = 0.0221885;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.427288; errup = 0.0266864; errdown = 0.0266198;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.451982; errup = 0.0281213; errdown = 0.0280742;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.508527; errup = 0.0335266; errdown = 0.0335433;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.423645; errup = 0.0435078; errdown = 0.0428184;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.5; errup = 0.0635302; errdown = 0.0635302;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.291572; errup = 0.0278471; errdown = 0.0270174;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.376812; errup = 0.0243114; errdown = 0.0241923;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.488679; errup = 0.0300519; errdown = 0.0300403;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.516588; errup = 0.0313049; errdown = 0.0313213;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.515094; errup = 0.0351539; errdown = 0.0351921;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.593909; errup = 0.0482373; errdown = 0.0490218;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.551282; errup = 0.0679688; errdown = 0.0692868;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.29316; errup = 0.0321667; errdown = 0.0309177;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.449911; errup = 0.028093; errdown = 0.0280427;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.562955; errup = 0.0337113; errdown = 0.0337775;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.583333; errup = 0.0345929; errdown = 0.0346779;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.553942; errup = 0.037316; errdown = 0.0374639;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.670103; errup = 0.049731; errdown = 0.0510836;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.662162; errup = 0.0682863; errdown = 0.072533;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.33617; errup = 0.0379333; errdown = 0.0366224;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.493603; errup = 0.0308202; errdown = 0.0308125;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.575588; errup = 0.0345602; errdown = 0.0346466;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.616114; errup = 0.0365028; errdown = 0.0366377;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.691837; errup = 0.0420474; errdown = 0.0424717;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.664706; errup = 0.051484; errdown = 0.0530429;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.671875; errup = 0.0720834; errdown = 0.0773989;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.424893; errup = 0.0489984; errdown = 0.0485139;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.52524; errup = 0.0462387; errdown = 0.046263;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.598726; errup = 0.051834; errdown = 0.051925;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.672794; errup = 0.0571534; errdown = 0.0572841;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.659955; errup = 0.0583602; errdown = 0.0586616;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.654545; errup = 0.0658697; errdown = 0.0670822;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.764706; errup = 0.0822803; errdown = 0.088313;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.384181; errup = 0.0507023; errdown = 0.0496024;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.568249; errup = 0.0501677; errdown = 0.0502506;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.672603; errup = 0.0574368; errdown = 0.0575908;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.682065; errup = 0.0581009; errdown = 0.0582584;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.75; errup = 0.0648599; errdown = 0.0653311;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.683544; errup = 0.0677075; errdown = 0.069177;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.71875; errup = 0.0834579; errdown = 0.0891295;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.397163; errup = 0.0556659; errdown = 0.0543881;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.601351; errup = 0.0530601; errdown = 0.0532006;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.702663; errup = 0.0598218; errdown = 0.0600126;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.741214; errup = 0.0627835; errdown = 0.0630176;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.744318; errup = 0.06499; errdown = 0.0655452;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.738255; errup = 0.0707083; errdown = 0.0726148;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.72; errup = 0.0894686; errdown = 0.0972978;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.552632; errup = 0.0673381; errdown = 0.0680904;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.629699; errup = 0.0555027; errdown = 0.0557027;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.7552; errup = 0.0637816; errdown = 0.0640219;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.790774; errup = 0.0663718; errdown = 0.0666336;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.784741; errup = 0.0673963; errdown = 0.0679512;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.761905; errup = 0.0733868; errdown = 0.0759424;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.857143; errup = 0.0844708; errdown = 0.0937921;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.535484; errup = 0.046646; errdown = 0.0470994;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.702381; errup = 0.0289533; errdown = 0.0292343;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.798263; errup = 0.0301237; errdown = 0.0303795;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.823843; errup = 0.0304719; errdown = 0.0307028;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.82148; errup = 0.0317993; errdown = 0.032267;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.887446; errup = 0.0370015; errdown = 0.0392256;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.789474; errup = 0.0512002; errdown = 0.0570999;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.607143; errup = 0.0611004; errdown = 0.0637535;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.769231; errup = 0.0322748; errdown = 0.032899;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.845422; errup = 0.0315886; errdown = 0.0319421;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.891209; errup = 0.0322527; errdown = 0.0325563;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.864407; errup = 0.0332178; errdown = 0.0339039;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.921951; errup = 0.0367854; errdown = 0.0393666;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.9; errup = 0.0474975; errdown = 0.0584863;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.673913; errup = 0.0782297; errdown = 0.0865062;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.823666; errup = 0.0338487; errdown = 0.0347553;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.883436; errup = 0.0328181; errdown = 0.0333208;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.902622; errup = 0.0326708; errdown = 0.0330305;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.929314; errup = 0.0339327; errdown = 0.0346537;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.914634; errup = 0.0382833; errdown = 0.0418475;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.940476; errup = 0.0410651; errdown = 0.0500144;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.73913; errup = 0.103; errdown = 0.127009;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.865517; errup = 0.0416831; errdown = 0.0430847;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.933712; errup = 0.040711; errdown = 0.0412243;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.92515; errup = 0.040197; errdown = 0.0405667;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.942356; errup = 0.0412866; errdown = 0.0420479;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.942197; errup = 0.0433697; errdown = 0.0461073;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.981481; errup = 0.0185185; errdown = 0.0583584;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 0.631579; errup = 0.127314; errdown = 0.14351;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.89083; errup = 0.0428662; errdown = 0.0448225;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.946602; errup = 0.0412909; errdown = 0.0420009;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.948998; errup = 0.0409667; errdown = 0.0414155;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.971722; errup = 0.0282776; errdown = 0.0423138;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.972222; errup = 0.0277778; errdown = 0.0461203;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.969697; errup = 0.030303; errdown = 0.0561038;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 0.8125; errup = 0.10127; errdown = 0.15011;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.934783; errup = 0.0224715; errdown = 0.0249545;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.959212; errup = 0.018765; errdown = 0.0193679;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 0.970817; errup = 0.0182131; errdown = 0.018526;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.976934; errup = 0.0183933; errdown = 0.0188786;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.987755; errup = 0.0122449; errdown = 0.0212926;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.967742; errup = 0.0232094; errdown = 0.0302634;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.308548;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.95614; errup = 0.0254774; errdown = 0.0333458;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.985366; errup = 0.0146341; errdown = 0.0196978;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.989101; errup = 0.0108992; errdown = 0.0185522;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.994286; errup = 0.00571429; errdown = 0.0186984;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.994764; errup = 0.0052356; errdown = 0.0214907;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.980198; errup = 0.019802; errdown = 0.0310067;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.460004;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 0.962963; errup = 0.037037; errdown = 0.0578629;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.990132; errup = 0.00986842; errdown = 0.0362957;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.993781; errup = 0.00621891; errdown = 0.0354061;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.995399; errup = 0.00460123; errdown = 0.0354958;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0361871;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.99115; errup = 0.00884956; errdown = 0.0403945;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2017("get_0l_fakemet_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff = 0.0974304; errup = 0.00970318; errdown = 0.00970318;}
  else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff = 0.14405; errup = 0.016044; errdown = 0.016044;}
  else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff = 0.184874; errup = 0.025163; errdown = 0.025163;}
  else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff = 0.321101; errup = 0.0447209; errdown = 0.0447209;}
  else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff = 0.348837; errup = 0.0513934; errdown = 0.0513934;}
  else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff = 0.545455; errup = 0.0567443; errdown = 0.0567443;}
  else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff = 0.62069; errup = 0.0901022; errdown = 0.0901022;}
  else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff = 0.818182; errup = 0.116291; errdown = 0.116291;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.148239; errup = 0.0101691; errdown = 0.0101691;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.170517; errup = 0.012474; errdown = 0.012474;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.214286; errup = 0.0150635; errdown = 0.0150635;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.241993; errup = 0.0180663; errdown = 0.0180663;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.277533; errup = 0.0210154; errdown = 0.0210154;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.369444; errup = 0.0254381; errdown = 0.0254381;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.395833; errup = 0.0315667; errdown = 0.0315667;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.424731; errup = 0.036244; errdown = 0.036244;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.486188; errup = 0.0371505; errdown = 0.0371505;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.51049; errup = 0.0418029; errdown = 0.0418029;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.622642; errup = 0.0332911; errdown = 0.0332911;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.669173; errup = 0.0407985; errdown = 0.0407985;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.728261; errup = 0.0463795; errdown = 0.0463795;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.709677; errup = 0.0576468; errdown = 0.0576468;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.757576; errup = 0.0746009; errdown = 0.0746009;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff = 0.916031; errup = 0.0242315; errdown = 0.0242315;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.954545; errup = 0.0314022; errdown = 0.0314022;}
  else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff = 0.189972; errup = 0.00841328; errdown = 0.00841328;}
  else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff = 0.228241; errup = 0.0103543; errdown = 0.0103543;}
  else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff = 0.243159; errup = 0.0119953; errdown = 0.0119953;}
  else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff = 0.323308; errup = 0.0153295; errdown = 0.0153295;}
  else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff = 0.377577; errup = 0.0174026; errdown = 0.0174026;}
  else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff = 0.407874; errup = 0.0195022; errdown = 0.0195022;}
  else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff = 0.418468; errup = 0.0218655; errdown = 0.0218655;}
  else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff = 0.486339; errup = 0.0261257; errdown = 0.0261257;}
  else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff = 0.454286; errup = 0.0266142; errdown = 0.0266142;}
  else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff = 0.582143; errup = 0.0294747; errdown = 0.0294747;}
  else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff = 0.590571; errup = 0.0244947; errdown = 0.0244947;}
  else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff = 0.671141; errup = 0.0272147; errdown = 0.0272147;}
  else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff = 0.795349; errup = 0.0275148; errdown = 0.0275148;}
  else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff = 0.823944; errup = 0.0319617; errdown = 0.0319617;}
  else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff = 0.796875; errup = 0.0355608; errdown = 0.0355608;}
  else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff = 0.911972; errup = 0.0168129; errdown = 0.0168129;}
  else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff = 0.975309; errup = 0.0121923; errdown = 0.0121923;}
  else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff = 0.966667; errup = 0.0231741; errdown = 0.0231741;}
  else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff = 0.245086; errup = 0.00879666; errdown = 0.00879666;}
  else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff = 0.279679; errup = 0.0103794; errdown = 0.0103794;}
  else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff = 0.319121; errup = 0.0118475; errdown = 0.0118475;}
  else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff = 0.363946; errup = 0.0140301; errdown = 0.0140301;}
  else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff = 0.376166; errup = 0.0155941; errdown = 0.0155941;}
  else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff = 0.445813; errup = 0.0174432; errdown = 0.0174432;}
  else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff = 0.507184; errup = 0.0189505; errdown = 0.0189505;}
  else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff = 0.532957; errup = 0.021651; errdown = 0.021651;}
  else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff = 0.599109; errup = 0.0231283; errdown = 0.0231283;}
  else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff = 0.631043; errup = 0.02434; errdown = 0.02434;}
  else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff = 0.685155; errup = 0.0187591; errdown = 0.0187591;}
  else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff = 0.694382; errup = 0.0218378; errdown = 0.0218378;}
  else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff = 0.740413; errup = 0.0238111; errdown = 0.0238111;}
  else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff = 0.8; errup = 0.0245718; errdown = 0.0245718;}
  else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff = 0.856436; errup = 0.0246715; errdown = 0.0246715;}
  else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff = 0.880342; errup = 0.0173238; errdown = 0.0173238;}
  else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff = 0.903382; errup = 0.0205343; errdown = 0.0205343;}
  else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff = 0.933333; errup = 0.0117589; errdown = 0.0117589;}
  else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff = 0.978102; errup = 0.00510454; errdown = 0.00510454;}
  else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff = 0.322707; errup = 0.00710557; errdown = 0.00710557;}
  else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff = 0.366185; errup = 0.00810742; errdown = 0.00810742;}
  else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff = 0.406324; errup = 0.00905657; errdown = 0.00905657;}
  else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff = 0.439822; errup = 0.00997531; errdown = 0.00997531;}
  else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff = 0.463355; errup = 0.0110594; errdown = 0.0110594;}
  else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff = 0.515652; errup = 0.0119226; errdown = 0.0119226;}
  else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff = 0.554072; errup = 0.0128428; errdown = 0.0128428;}
  else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff = 0.598348; errup = 0.0134323; errdown = 0.0134323;}
  else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff = 0.642009; errup = 0.0144877; errdown = 0.0144877;}
  else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff = 0.696538; errup = 0.0146713; errdown = 0.0146713;}
  else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff = 0.715486; errup = 0.0110539; errdown = 0.0110539;}
  else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff = 0.781983; errup = 0.0113604; errdown = 0.0113604;}
  else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff = 0.827618; errup = 0.0114987; errdown = 0.0114987;}
  else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff = 0.839763; errup = 0.0115368; errdown = 0.0115368;}
  else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff = 0.884467; errup = 0.0114531; errdown = 0.0114531;}
  else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff = 0.898712; errup = 0.00730038; errdown = 0.00730038;}
  else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff = 0.930855; errup = 0.00691768; errdown = 0.00691768;}
  else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff = 0.941468; errup = 0.00522821; errdown = 0.00522821;}
  else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff = 0.959596; errup = 0.00528902; errdown = 0.00528902;}
  else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff = 0.976963; errup = 0.00485453; errdown = 0.00485453;}
  else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff = 0.991124; errup = 0.00360739; errdown = 0.00360739;}
  else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.973333; errup = 0.00515956; errdown = 0.00515956;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.489496; errup = 0.00369257; errdown = 0.00369257;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.522675; errup = 0.00395583; errdown = 0.00395583;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.568861; errup = 0.00421634; errdown = 0.00421634;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.592914; errup = 0.00450707; errdown = 0.00450707;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.635472; errup = 0.00469609; errdown = 0.00469609;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.674262; errup = 0.00491711; errdown = 0.00491711;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.686402; errup = 0.00517974; errdown = 0.00517974;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.726592; errup = 0.00534945; errdown = 0.00534945;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.737255; errup = 0.00551234; errdown = 0.00551234;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.775837; errup = 0.00559582; errdown = 0.00559582;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.802235; errup = 0.00410939; errdown = 0.00410939;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.83659; errup = 0.00423424; errdown = 0.00423424;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.849854; errup = 0.00442152; errdown = 0.00442152;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.880527; errup = 0.00454664; errdown = 0.00454664;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.893832; errup = 0.00466479; errdown = 0.00466479;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.916677; errup = 0.00306379; errdown = 0.00306379;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.936391; errup = 0.00323515; errdown = 0.00323515;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.941998; errup = 0.0028653; errdown = 0.0028653;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.952878; errup = 0.00376834; errdown = 0.00376834;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.963603; errup = 0.00465149; errdown = 0.00465149;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.954495; errup = 0.00694311; errdown = 0.00694311;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.983751; errup = 0.00328979; errdown = 0.00328979;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.387888; errup = 0.001472; errdown = 0.001472;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.417927; errup = 0.00160609; errdown = 0.00160609;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.451373; errup = 0.0017459; errdown = 0.0017459;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.475537; errup = 0.00189276; errdown = 0.00189276;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.512738; errup = 0.00204099; errdown = 0.00204099;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.53824; errup = 0.00218953; errdown = 0.00218953;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.567544; errup = 0.0023476; errdown = 0.0023476;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.598143; errup = 0.00249669; errdown = 0.00249669;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.624385; errup = 0.00264351; errdown = 0.00264351;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.649303; errup = 0.00279531; errdown = 0.00279531;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.681851; errup = 0.002131; errdown = 0.002131;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.729521; errup = 0.00231462; errdown = 0.00231462;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.762535; errup = 0.00251095; errdown = 0.00251095;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.798438; errup = 0.00267273; errdown = 0.00267273;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.821241; errup = 0.00285038; errdown = 0.00285038;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.856744; errup = 0.00196232; errdown = 0.00196232;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.886954; errup = 0.00223653; errdown = 0.00223653;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.916497; errup = 0.00185908; errdown = 0.00185908;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.936176; errup = 0.00234582; errdown = 0.00234582;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.947549; errup = 0.00292349; errdown = 0.00292349;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.954474; errup = 0.0037065; errdown = 0.0037065;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.967807; errup = 0.00238727; errdown = 0.00238727;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1el_trigeff2017("get_1el_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.00606061; errup = 0.00297522; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.252727; errup = 0.0137849; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.47488; errup = 0.0101789; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.0580311; errup = 0.00853854; errdown = 0.00756;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.273902; errup = 0.0169699; errdown = 0.0163626;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.512061; errup = 0.0119679; errdown = 0.0119814;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.239788; errup = 0.0123203; errdown = 0.0119141;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.432836; errup = 0.0161873; errdown = 0.0160495;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.610241; errup = 0.00967965; errdown = 0.00976557;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.368855; errup = 0.0057901; errdown = 0.00575302;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.58404; errup = 0.00788592; errdown = 0.0079285;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.718103; errup = 0.00425673; errdown = 0.00429573;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.358079; errup = 0.0165383; errdown = 0.0162179;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.574176; errup = 0.0270377; errdown = 0.0274661;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.721596; errup = 0.01343; errdown = 0.0138268;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.727555; errup = 0.00887959; errdown = 0.00906034;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.776404; errup = 0.014199; errdown = 0.0148447;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.830705; errup = 0.00707675; errdown = 0.00731432;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.0128205; errup = 0.00857983; errdown = 0.00552679;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.251185; errup = 0.0228961; errdown = 0.0216526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.515088; errup = 0.0135744; errdown = 0.0135961;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.07; errup = 0.0180052; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.359756; errup = 0.0284325; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.58209; errup = 0.0154441; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.313953; errup = 0.0239721; errdown = 0.0230631;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.564612; errup = 0.0229432; errdown = 0.0232113;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.668443; errup = 0.0116455; errdown = 0.0118501;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.478528; errup = 0.0103138; errdown = 0.0102959;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.690422; errup = 0.0112097; errdown = 0.0114324;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.76283; errup = 0.00509721; errdown = 0.00517276;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.404682; errup = 0.0303426; errdown = 0.0296689;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.722973; errup = 0.0383999; errdown = 0.0416756;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.766862; errup = 0.0165056; errdown = 0.0173229;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.73375; errup = 0.0159396; errdown = 0.0165468;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.814136; errup = 0.0203184; errdown = 0.0220684;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.843819; errup = 0.00860726; errdown = 0.0089989;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.0526316; errup = 0.0271984; errdown = 0.019208;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.295082; errup = 0.0375578; errdown = 0.0351452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.577465; errup = 0.0167048; errdown = 0.0168781;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.196429; errup = 0.0448753; errdown = 0.0388449;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.467105; errup = 0.0438353; errdown = 0.0433692;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.620134; errup = 0.0182735; errdown = 0.018608;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.495098; errup = 0.037384; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.59375; errup = 0.03452; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.738619; errup = 0.0125052; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.584366; errup = 0.0156025; errdown = 0.0157683;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.737811; errup = 0.0147447; errdown = 0.0152789;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.796101; errup = 0.00566351; errdown = 0.00578094;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.607843; errup = 0.0518085; errdown = 0.0541227;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.77551; errup = 0.0439117; errdown = 0.050173;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.801205; errup = 0.0182366; errdown = 0.0195126;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.767516; errup = 0.0244905; errdown = 0.0263047;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.791878; errup = 0.0297937; errdown = 0.0330066;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.863636; errup = 0.00987432; errdown = 0.0104889;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.0338983; errup = 0.0429671; errdown = 0.021865;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.357895; errup = 0.0557905; errdown = 0.0525014;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.63606; errup = 0.0202397; errdown = 0.0207128;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.230769; errup = 0.0733859; errdown = 0.0614831;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.532468; errup = 0.0624335; errdown = 0.0633651;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.701613; errup = 0.0211067; errdown = 0.0219585;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.645161; errup = 0.0455298; errdown = 0.0480775;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.700787; errup = 0.0426733; errdown = 0.0461241;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.807339; errup = 0.0135565; errdown = 0.0142912;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.716075; errup = 0.0211522; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.797227; errup = 0.0170527; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.84023; errup = 0.0061067; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.785714; errup = 0.0664487; errdown = 0.0823485;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.891892; errup = 0.0510106; errdown = 0.077259;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.827922; errup = 0.0219648; errdown = 0.0242503;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.841667; errup = 0.0342095; errdown = 0.0405813;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.893443; errup = 0.0283253; errdown = 0.0356815;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.905502; errup = 0.0102026; errdown = 0.0112294;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.2; errup = 0.130119; errdown = 0.0931225;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.604167; errup = 0.0776099; errdown = 0.0825214;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.801508; errup = 0.0204268; errdown = 0.0220347;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.222222; errup = 0.109068; errdown = 0.0843892;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.575; errup = 0.0872414; errdown = 0.0915181;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.788406; errup = 0.0225238; errdown = 0.0243027;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.745098; errup = 0.0646476; errdown = 0.0756564;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.813333; errup = 0.0466366; errdown = 0.0561452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.85; errup = 0.0150446; errdown = 0.0163204;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.815686; errup = 0.0248726; errdown = 0.027542;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.876254; errup = 0.0193393; errdown = 0.0220715;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.895428; errup = 0.00633502; errdown = 0.00668102;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.892857; errup = 0.0577336; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.885965; errup = 0.0213678; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.982456; errup = 0.0145177; errdown = 0.0391869;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.970588; errup = 0.0189746; errdown = 0.0374789;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.93956; errup = 0.0102474; errdown = 0.0120003;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.461538; errup = 0.172004; errdown = 0.164847;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.612903; errup = 0.0974957; errdown = 0.105934;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.454545; errup = 0.189662; errdown = 0.179582;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.883178; errup = 0.0223037; errdown = 0.0262543;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.75; errup = 0.0943497; errdown = 0.119341;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.906977; errup = 0.0439838; errdown = 0.0674608;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.914286; errup = 0.0132412; errdown = 0.0152122;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.873418; errup = 0.026946; errdown = 0.0322169;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929648; errup = 0.0182155; errdown = 0.0231163;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929748; errup = 0.00601652; errdown = 0.00650973;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.809524; errup = 0.0888157; errdown = 0.125184;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.907801; errup = 0.0246102; errdown = 0.0312058;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.942857; errup = 0.0368224; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.975; errup = 0.0206905; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.951351; errup = 0.011195; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.8; errup = 0.106751; errdown = 0.157061;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0536391; errdown = 0.0995072;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.955645; errup = 0.0130229; errdown = 0.017253;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.769231; errup = 0.122762; errdown = 0.174724;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.869565; errup = 0.0701228; errdown = 0.110814;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.920833; errup = 0.0175659; errdown = 0.0214819;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.944444; errup = 0.046004; errdown = 0.116415;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.97619; errup = 0.0197048; errdown = 0.0526298;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.961039; errup = 0.00899222; errdown = 0.0112252;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.950413; errup = 0.0194949; errdown = 0.0284436;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.978022; errup = 0.0104893; errdown = 0.0170363;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.969503; errup = 0.00381731; errdown = 0.00430889;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0690403; errdown = 0.16652;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.964706; errup = 0.013912; errdown = 0.0204859;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.966667; errup = 0.0275914; errdown = 0.072517;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0384134;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.980815; errup = 0.00661202; errdown = 0.00932515;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.913043; errup = 0.0559625; errdown = 0.103371;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.981651; errup = 0.00725455; errdown = 0.0107985;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0769247;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.984674; errup = 0.00732055; errdown = 0.0119517;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0659133;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.042887;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.997938; errup = 0.00170573; errdown = 0.00472518;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.990099; errup = 0.00819202; errdown = 0.0223979;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00964279;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.992439; errup = 0.00171528; errdown = 0.00215231;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.115502;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00804214;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.108691;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0283562;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998397; errup = 0.00132575; errdown = 0.0036754;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1mu_trigeff2017("get_1mu_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.0708592; errup = 0.00848378; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.386643; errup = 0.0149712; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.79144; errup = 0.00809149; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.271298; errup = 0.0170541; errdown = 0.0164307;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.550075; errup = 0.0198882; errdown = 0.0200436;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.843672; errup = 0.0091337; errdown = 0.00957456;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.431312; errup = 0.012669; errdown = 0.0125819;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.704082; errup = 0.0120989; errdown = 0.0123844;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.886586; errup = 0.00491944; errdown = 0.00510812;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.887264; errup = 0.00786673; errdown = 0.00835697;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.92204; errup = 0.00652681; errdown = 0.00704377;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.936299; errup = 0.00313807; errdown = 0.00328508;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.0614251; errup = 0.0143316; errdown = 0.0119574;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.414583; errup = 0.0236888; errdown = 0.0233169;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.814067; errup = 0.00972747; errdown = 0.0101237;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.325843; errup = 0.0311932; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.614198; errup = 0.0281699; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.857381; errup = 0.0102039; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.543956; errup = 0.0190702; errdown = 0.0191953;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.780458; errup = 0.0146207; errdown = 0.0153252;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.893736; errup = 0.0056122; errdown = 0.00587802;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.921429; errup = 0.0093482; errdown = 0.0104152;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.952077; errup = 0.00699194; errdown = 0.00802783;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.940981; errup = 0.0035208; errdown = 0.00372265;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.15; errup = 0.0337209; errdown = 0.0288809;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.45098; errup = 0.0332565; errdown = 0.0328488;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.842715; errup = 0.0105919; errdown = 0.0111815;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.520325; errup = 0.0487155; errdown = 0.0490739;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.647668; errup = 0.0360589; errdown = 0.037699;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.861354; errup = 0.0115486; errdown = 0.0123748;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.666667; errup = 0.025072; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.806922; errup = 0.0171572; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.903462; errup = 0.00599439; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.941414; errup = 0.0106018; errdown = 0.0125537;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.954397; errup = 0.00843586; errdown = 0.0100577;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.95625; errup = 0.00337964; errdown = 0.00363757;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.14; errup = 0.0672437; errdown = 0.0501519;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.556452; errup = 0.0479412; errdown = 0.0489288;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.868726; errup = 0.0122566; errdown = 0.0132553;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.619048; errup = 0.0664109; errdown = 0.0706357;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.816514; errup = 0.0382778; errdown = 0.044765;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.919672; errup = 0.0110898; errdown = 0.0125669;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.73262; errup = 0.0336118; errdown = 0.0362961;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.855596; errup = 0.0215151; errdown = 0.0242958;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.931906; errup = 0.00593101; errdown = 0.0064273;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.956364; errup = 0.0122799; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.972569; errup = 0.00808855; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.968242; errup = 0.00341417; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.423077; errup = 0.116981; errdown = 0.110076;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.0654592; errdown = 0.0717057;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.912046; errup = 0.0124958; errdown = 0.0141913;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.724138; errup = 0.089299; errdown = 0.107524;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.850746; errup = 0.044684; errdown = 0.0568246;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.943262; errup = 0.0112896; errdown = 0.0136025;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.94186; errup = 0.0248848; errdown = 0.0374166;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.928571; errup = 0.0191788; errdown = 0.024547;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.956884; errup = 0.00560047; errdown = 0.0063382;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.980392; errup = 0.0106526; errdown = 0.0187052;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.973913; errup = 0.0102998; errdown = 0.0152564;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.981845; errup = 0.00295218; errdown = 0.0034642;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.533333; errup = 0.152935; errdown = 0.15827;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.882353; errup = 0.0554381; errdown = 0.0832913;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.929412; errup = 0.01398; errdown = 0.0167685;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.769231; errup = 0.122762; errdown = 0.174724;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.904762; errup = 0.0450174; errdown = 0.0689193;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.981308; errup = 0.0073897; errdown = 0.0109973;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0283562;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.94898; errup = 0.0218627; errdown = 0.0330405;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.973348; errup = 0.00524904; errdown = 0.00636007;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988636; errup = 0.00940245; errdown = 0.0256444;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.985469; errup = 0.00306041; errdown = 0.00377044;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.888889; errup = 0.0598486; errdown = 0.0963981;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.98374; errup = 0.00643122; errdown = 0.00958564;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.947368; errup = 0.0435805; errdown = 0.110836;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.988473; errup = 0.00550887; errdown = 0.00902056;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0392319;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0182418;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.990955; errup = 0.00295229; errdown = 0.00410363;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0249041;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.993711; errup = 0.0052034; errdown = 0.0143129;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.995131; errup = 0.00168324; errdown = 0.00239245;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.205568;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0709947;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 0.992788; errup = 0.00392226; errdown = 0.00696502;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0972223;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.99505; errup = 0.00319693; errdown = 0.00649192;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0542609;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0224723;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.994958; errup = 0.00199818; errdown = 0.00299934;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0461088;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.012972;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.997174; errup = 0.00112047; errdown = 0.0016842;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2el_trigeff2017("get_2el_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.949309; errup = 0.0148916; errdown = 0.0148916;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.939502; errup = 0.0142222; errdown = 0.0142222;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.930233; errup = 0.0158603; errdown = 0.0158603;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.919014; errup = 0.0161885; errdown = 0.0161885;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.931818; errup = 0.0155131; errdown = 0.0155131;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.94382; errup = 0.0140922; errdown = 0.0140922;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.973913; errup = 0.0105101; errdown = 0.0105101;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.941704; errup = 0.01569; errdown = 0.01569;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.960591; errup = 0.0136558; errdown = 0.0136558;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.919075; errup = 0.0207345; errdown = 0.0207345;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.941176; errup = 0.0180462; errdown = 0.0180462;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.972222; errup = 0.0122488; errdown = 0.0122488;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.934211; errup = 0.0201085; errdown = 0.0201085;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.983325; errup = 0.00197631; errdown = 0.00197631;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2mu_trigeff2017("get_2mu_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.948413; errup = 0.0139338; errdown = 0.0139338;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.963134; errup = 0.0090451; errdown = 0.0090451;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.983369; errup = 0.000433117; errdown = 0.000433117;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_trigeff2018("get_0l_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 250 && met> 150 && met<= 155) {eff = 0.0967153; errup = 0.0150137; errdown = 0.0145224;}
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.167497; errup = 0.0211652; errdown = 0.0210723;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.241638; errup = 0.0301657; errdown = 0.0300795;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.296256; errup = 0.0366391; errdown = 0.03657;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.329746; errup = 0.0417888; errdown = 0.0416636;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.324786; errup = 0.0467917; errdown = 0.0461939;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.355705; errup = 0.0603952; errdown = 0.0589062;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.113527; errup = 0.0180531; errdown = 0.0174079;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.206503; errup = 0.0258924; errdown = 0.0258034;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.303987; errup = 0.0375402; errdown = 0.0374743;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.343447; errup = 0.0422467; errdown = 0.0421913;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.386041; errup = 0.0484809; errdown = 0.0483928;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.408163; errup = 0.0557802; errdown = 0.0554969;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.557692; errup = 0.0842001; errdown = 0.0849625;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.150068; errup = 0.0212125; errdown = 0.0205461;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.238877; errup = 0.0268564; errdown = 0.0267496;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.320253; errup = 0.035451; errdown = 0.0353715;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.378999; errup = 0.0415648; errdown = 0.0415101;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.450299; errup = 0.050112; errdown = 0.0500683;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.47076; errup = 0.0566384; errdown = 0.056548;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.442478; errup = 0.0690247; errdown = 0.0682026;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.167954; errup = 0.0251208; errdown = 0.0241645;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.305413; errup = 0.0337355; errdown = 0.0336594;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.400266; errup = 0.0436087; errdown = 0.0435682;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.453343; errup = 0.0491396; errdown = 0.0491197;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.50641; errup = 0.0558154; errdown = 0.055821;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.52459; errup = 0.0622959; errdown = 0.0623782;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.514019; errup = 0.074971; errdown = 0.0751716;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.233546; errup = 0.0321916; errdown = 0.0314222;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.354263; errup = 0.0389686; errdown = 0.0389045;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.447485; errup = 0.048566; errdown = 0.0485431;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.488976; errup = 0.0528465; errdown = 0.0528416;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.556686; errup = 0.0610976; errdown = 0.0611526;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.580645; errup = 0.0678513; errdown = 0.0681328;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.574074; errup = 0.0786405; errdown = 0.0796228;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.213115; errup = 0.0323576; errdown = 0.0311651;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.381985; errup = 0.0422952; errdown = 0.042226;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.4975; errup = 0.0538146; errdown = 0.0538134;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.56453; errup = 0.060499; errdown = 0.0605259;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.599411; errup = 0.0652644; errdown = 0.0653553;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.626506; errup = 0.0725999; errdown = 0.0730831;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.577982; errup = 0.078725; errdown = 0.0797419;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.242537; errup = 0.0365664; errdown = 0.0349584;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.468608; errup = 0.0456795; errdown = 0.0456595;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.54512; errup = 0.0523274; errdown = 0.0523529;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.628895; errup = 0.0596188; errdown = 0.0596832;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.678516; errup = 0.0649901; errdown = 0.0651589;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.708333; errup = 0.0724515; errdown = 0.0733919;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.72619; errup = 0.0840402; errdown = 0.0878291;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.256637; errup = 0.0400091; errdown = 0.038164;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.476334; errup = 0.0466341; errdown = 0.0466169;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.604627; errup = 0.0577631; errdown = 0.0578274;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.686888; errup = 0.0646849; errdown = 0.0647724;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.722523; errup = 0.0690384; errdown = 0.0692797;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.776699; errup = 0.077237; errdown = 0.0784073;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.8; errup = 0.0868232; errdown = 0.0916411;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.318182; errup = 0.0486661; errdown = 0.0469038;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.507812; errup = 0.0496426; errdown = 0.0496489;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.664678; errup = 0.0631778; errdown = 0.0632866;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.710407; errup = 0.0669595; errdown = 0.0670759;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.76234; errup = 0.0723177; errdown = 0.072583;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.693548; errup = 0.0727106; errdown = 0.0738229;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.785714; errup = 0.088376; errdown = 0.094047;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.392593; errup = 0.0587338; errdown = 0.0573817;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.575663; errup = 0.0564932; errdown = 0.0565812;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.717647; errup = 0.0678563; errdown = 0.0680036;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.768065; errup = 0.0719207; errdown = 0.072056;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.791075; errup = 0.0748891; errdown = 0.0752093;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.784211; errup = 0.0782049; errdown = 0.0795381;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.814815; errup = 0.0870468; errdown = 0.0918832;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.447115; errup = 0.0460682; errdown = 0.0456331;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.669557; errup = 0.0436146; errdown = 0.043731;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.73572; errup = 0.0467041; errdown = 0.0468034;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.796218; errup = 0.0499228; errdown = 0.0500172;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.844789; errup = 0.053137; errdown = 0.0533263;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.797251; errup = 0.0544585; errdown = 0.0554635;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.849315; errup = 0.0601831; errdown = 0.0630276;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.483051; errup = 0.0582314; errdown = 0.0579628;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.707895; errup = 0.0464989; errdown = 0.046708;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.823691; errup = 0.0517625; errdown = 0.0519056;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.853846; errup = 0.0531989; errdown = 0.053307;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.86014; errup = 0.0542692; errdown = 0.0545357;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.924528; errup = 0.0589093; errdown = 0.0599973;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.881818; errup = 0.0624045; errdown = 0.0667303;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.5875; errup = 0.0696658; errdown = 0.0717458;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.762238; errup = 0.0500787; errdown = 0.0504354;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.890011; errup = 0.055502; errdown = 0.0556851;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.905991; errup = 0.0561762; errdown = 0.0563066;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.897516; errup = 0.0562536; errdown = 0.0565501;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.877049; errup = 0.0577854; errdown = 0.0591199;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.921569; errup = 0.0624077; errdown = 0.0671418;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.759259; errup = 0.0821659; errdown = 0.0906687;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.809045; errup = 0.0616139; errdown = 0.062172;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.9; errup = 0.0657743; errdown = 0.0659915;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.915713; errup = 0.0665311; errdown = 0.0666614;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.965074; errup = 0.0349265; errdown = 0.0701522;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.948357; errup = 0.0516432; errdown = 0.07113;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.936364; errup = 0.0636364; errdown = 0.0748697;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 0.757576; errup = 0.0961604; errdown = 0.111793;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.888158; errup = 0.066507; errdown = 0.0673262;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.945899; errup = 0.0541012; errdown = 0.0690012;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.949398; errup = 0.0506024; errdown = 0.0689052;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.953722; errup = 0.0462777; errdown = 0.0696014;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.956098; errup = 0.0439024; errdown = 0.0715069;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.952381; errup = 0.047619; errdown = 0.0774568;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 0.75; errup = 0.0952577; errdown = 0.12006;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.923469; errup = 0.0210724; errdown = 0.0226616;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.968597; errup = 0.0178036; errdown = 0.0181361;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 0.975101; errup = 0.0175303; errdown = 0.0177027;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.977588; errup = 0.0177678; errdown = 0.0181074;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.983092; errup = 0.0169082; errdown = 0.0194056;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.970588; errup = 0.0211717; errdown = 0.0257862;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 0.571429; errup = 0.222712; errdown = 0.248043;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.960265; errup = 0.0229593; errdown = 0.0284651;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.988679; errup = 0.0113208; errdown = 0.0185487;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.985045; errup = 0.0149551; errdown = 0.0179142;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.989262; errup = 0.0107383; errdown = 0.018086;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.993174; errup = 0.00682594; errdown = 0.0195361;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.981481; errup = 0.0185185; errdown = 0.0294305;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.369032;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 0.97561; errup = 0.0188677; errdown = 0.0329513;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.994859; errup = 0.00514139; errdown = 0.0125718;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.996457; errup = 0.00354296; errdown = 0.0109897;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.990426; errup = 0.00957447; errdown = 0.0114221;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.997326; errup = 0.0026738; errdown = 0.0122737;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0144134;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2018("get_0l_fakemet_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 350 && met> 150 && met<= 160) {eff = 0.078; errup = 0.011993; errdown = 0.011993;}
  else if (ht> 0 && ht<= 350 && met> 160 && met<= 170) {eff = 0.102837; errup = 0.0180878; errdown = 0.0180878;}
  else if (ht> 0 && ht<= 350 && met> 170 && met<= 180) {eff = 0.186207; errup = 0.0323274; errdown = 0.0323274;}
  else if (ht> 0 && ht<= 350 && met> 180 && met<= 190) {eff = 0.253731; errup = 0.0531615; errdown = 0.0531615;}
  else if (ht> 0 && ht<= 350 && met> 190 && met<= 200) {eff = 0.363636; errup = 0.0725204; errdown = 0.0725204;}
  else if (ht> 0 && ht<= 350 && met> 200 && met<= 225) {eff = 0.428571; errup = 0.0763604; errdown = 0.0763604;}
  else if (ht> 0 && ht<= 350 && met> 225 && met<= 250) {eff = 0.529412; errup = 0.121058; errdown = 0.121058;}
  else if (ht> 0 && ht<= 350 && met> 250 && met<= 9999) {eff = 0.823529; errup = 0.0924594; errdown = 0.0924594;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.153846; errup = 0.0114786; errdown = 0.0114786;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.164521; errup = 0.014023; errdown = 0.014023;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.227273; errup = 0.0178692; errdown = 0.0178692;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.239374; errup = 0.0201823; errdown = 0.0201823;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.285294; errup = 0.024489; errdown = 0.024489;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.309353; errup = 0.0277225; errdown = 0.0277225;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.402655; errup = 0.0326231; errdown = 0.0326231;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.411111; errup = 0.0366741; errdown = 0.0366741;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.503876; errup = 0.0440212; errdown = 0.0440212;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.473684; errup = 0.0512278; errdown = 0.0512278;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.649351; errup = 0.0384517; errdown = 0.0384517;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.666667; errup = 0.0488824; errdown = 0.0488824;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.860759; errup = 0.0389502; errdown = 0.0389502;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.722222; errup = 0.0609519; errdown = 0.0609519;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.837209; errup = 0.0562986; errdown = 0.0562986;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 300) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.973684; errup = 0.0259672; errdown = 0.0259672;}
  else if (ht> 450 && ht<= 550 && met> 150 && met<= 155) {eff = 0.185145; errup = 0.0090772; errdown = 0.0090772;}
  else if (ht> 450 && ht<= 550 && met> 155 && met<= 160) {eff = 0.224373; errup = 0.0111693; errdown = 0.0111693;}
  else if (ht> 450 && ht<= 550 && met> 160 && met<= 165) {eff = 0.253235; errup = 0.0132203; errdown = 0.0132203;}
  else if (ht> 450 && ht<= 550 && met> 165 && met<= 170) {eff = 0.296253; errup = 0.0156247; errdown = 0.0156247;}
  else if (ht> 450 && ht<= 550 && met> 170 && met<= 175) {eff = 0.351145; errup = 0.0186508; errdown = 0.0186508;}
  else if (ht> 450 && ht<= 550 && met> 175 && met<= 180) {eff = 0.397363; errup = 0.0212361; errdown = 0.0212361;}
  else if (ht> 450 && ht<= 550 && met> 180 && met<= 185) {eff = 0.430233; errup = 0.0238763; errdown = 0.0238763;}
  else if (ht> 450 && ht<= 550 && met> 185 && met<= 190) {eff = 0.512894; errup = 0.0267555; errdown = 0.0267555;}
  else if (ht> 450 && ht<= 550 && met> 190 && met<= 195) {eff = 0.553191; errup = 0.0296056; errdown = 0.0296056;}
  else if (ht> 450 && ht<= 550 && met> 195 && met<= 200) {eff = 0.571429; errup = 0.033065; errdown = 0.033065;}
  else if (ht> 450 && ht<= 550 && met> 200 && met<= 210) {eff = 0.626263; errup = 0.0243116; errdown = 0.0243116;}
  else if (ht> 450 && ht<= 550 && met> 210 && met<= 220) {eff = 0.731707; errup = 0.0282492; errdown = 0.0282492;}
  else if (ht> 450 && ht<= 550 && met> 220 && met<= 230) {eff = 0.797753; errup = 0.0301069; errdown = 0.0301069;}
  else if (ht> 450 && ht<= 550 && met> 230 && met<= 240) {eff = 0.758865; errup = 0.0360249; errdown = 0.0360249;}
  else if (ht> 450 && ht<= 550 && met> 240 && met<= 250) {eff = 0.824074; errup = 0.0366384; errdown = 0.0366384;}
  else if (ht> 450 && ht<= 550 && met> 250 && met<= 300) {eff = 0.90458; errup = 0.0181507; errdown = 0.0181507;}
  else if (ht> 450 && ht<= 550 && met> 300 && met<= 400) {eff = 0.948387; errup = 0.0177708; errdown = 0.0177708;}
  else if (ht> 450 && ht<= 550 && met> 400 && met<= 9999) {eff = 0.959459; errup = 0.0229267; errdown = 0.0229267;}
  else if (ht> 550 && ht<= 650 && met> 150 && met<= 155) {eff = 0.2286; errup = 0.00942298; errdown = 0.00942298;}
  else if (ht> 550 && ht<= 650 && met> 155 && met<= 160) {eff = 0.262344; errup = 0.0107295; errdown = 0.0107295;}
  else if (ht> 550 && ht<= 650 && met> 160 && met<= 165) {eff = 0.332286; errup = 0.0132019; errdown = 0.0132019;}
  else if (ht> 550 && ht<= 650 && met> 165 && met<= 170) {eff = 0.356405; errup = 0.0153936; errdown = 0.0153936;}
  else if (ht> 550 && ht<= 650 && met> 170 && met<= 175) {eff = 0.418685; errup = 0.0167548; errdown = 0.0167548;}
  else if (ht> 550 && ht<= 650 && met> 175 && met<= 180) {eff = 0.471642; errup = 0.0192856; errdown = 0.0192856;}
  else if (ht> 550 && ht<= 650 && met> 180 && met<= 185) {eff = 0.472222; errup = 0.0222374; errdown = 0.0222374;}
  else if (ht> 550 && ht<= 650 && met> 185 && met<= 190) {eff = 0.543033; errup = 0.02255; errdown = 0.02255;}
  else if (ht> 550 && ht<= 650 && met> 190 && met<= 195) {eff = 0.584699; errup = 0.0257577; errdown = 0.0257577;}
  else if (ht> 550 && ht<= 650 && met> 195 && met<= 200) {eff = 0.635328; errup = 0.0256919; errdown = 0.0256919;}
  else if (ht> 550 && ht<= 650 && met> 200 && met<= 210) {eff = 0.669922; errup = 0.0207819; errdown = 0.0207819;}
  else if (ht> 550 && ht<= 650 && met> 210 && met<= 220) {eff = 0.737113; errup = 0.0223478; errdown = 0.0223478;}
  else if (ht> 550 && ht<= 650 && met> 220 && met<= 230) {eff = 0.744409; errup = 0.0246551; errdown = 0.0246551;}
  else if (ht> 550 && ht<= 650 && met> 230 && met<= 240) {eff = 0.791304; errup = 0.0267957; errdown = 0.0267957;}
  else if (ht> 550 && ht<= 650 && met> 240 && met<= 250) {eff = 0.828125; errup = 0.0272272; errdown = 0.0272272;}
  else if (ht> 550 && ht<= 650 && met> 250 && met<= 275) {eff = 0.881459; errup = 0.0178212; errdown = 0.0178212;}
  else if (ht> 550 && ht<= 650 && met> 275 && met<= 300) {eff = 0.88587; errup = 0.023441; errdown = 0.023441;}
  else if (ht> 550 && ht<= 650 && met> 300 && met<= 400) {eff = 0.887029; errup = 0.014479; errdown = 0.014479;}
  else if (ht> 550 && ht<= 650 && met> 400 && met<= 9999) {eff = 0.977528; errup = 0.00473689; errdown = 0.00473689;}
  else if (ht> 650 && ht<= 800 && met> 150 && met<= 155) {eff = 0.312199; errup = 0.00737961; errdown = 0.00737961;}
  else if (ht> 650 && ht<= 800 && met> 155 && met<= 160) {eff = 0.355075; errup = 0.00835433; errdown = 0.00835433;}
  else if (ht> 650 && ht<= 800 && met> 160 && met<= 165) {eff = 0.397912; errup = 0.00962519; errdown = 0.00962519;}
  else if (ht> 650 && ht<= 800 && met> 165 && met<= 170) {eff = 0.446444; errup = 0.0104157; errdown = 0.0104157;}
  else if (ht> 650 && ht<= 800 && met> 170 && met<= 175) {eff = 0.505056; errup = 0.0115341; errdown = 0.0115341;}
  else if (ht> 650 && ht<= 800 && met> 175 && met<= 180) {eff = 0.540892; errup = 0.012404; errdown = 0.012404;}
  else if (ht> 650 && ht<= 800 && met> 180 && met<= 185) {eff = 0.593127; errup = 0.0128787; errdown = 0.0128787;}
  else if (ht> 650 && ht<= 800 && met> 185 && met<= 190) {eff = 0.609113; errup = 0.0137958; errdown = 0.0137958;}
  else if (ht> 650 && ht<= 800 && met> 190 && met<= 195) {eff = 0.628141; errup = 0.0153217; errdown = 0.0153217;}
  else if (ht> 650 && ht<= 800 && met> 195 && met<= 200) {eff = 0.692389; errup = 0.0150048; errdown = 0.0150048;}
  else if (ht> 650 && ht<= 800 && met> 200 && met<= 210) {eff = 0.703607; errup = 0.011694; errdown = 0.011694;}
  else if (ht> 650 && ht<= 800 && met> 210 && met<= 220) {eff = 0.790055; errup = 0.0114418; errdown = 0.0114418;}
  else if (ht> 650 && ht<= 800 && met> 220 && met<= 230) {eff = 0.805634; errup = 0.0121256; errdown = 0.0121256;}
  else if (ht> 650 && ht<= 800 && met> 230 && met<= 240) {eff = 0.862187; errup = 0.0116332; errdown = 0.0116332;}
  else if (ht> 650 && ht<= 800 && met> 240 && met<= 250) {eff = 0.868639; errup = 0.0116205; errdown = 0.0116205;}
  else if (ht> 650 && ht<= 800 && met> 250 && met<= 275) {eff = 0.878536; errup = 0.00769318; errdown = 0.00769318;}
  else if (ht> 650 && ht<= 800 && met> 275 && met<= 300) {eff = 0.922261; errup = 0.00711815; errdown = 0.00711815;}
  else if (ht> 650 && ht<= 800 && met> 300 && met<= 350) {eff = 0.934977; errup = 0.00508305; errdown = 0.00508305;}
  else if (ht> 650 && ht<= 800 && met> 350 && met<= 400) {eff = 0.95888; errup = 0.00484746; errdown = 0.00484746;}
  else if (ht> 650 && ht<= 800 && met> 400 && met<= 450) {eff = 0.969673; errup = 0.00484455; errdown = 0.00484455;}
  else if (ht> 650 && ht<= 800 && met> 450 && met<= 500) {eff = 0.979543; errup = 0.00491061; errdown = 0.00491061;}
  else if (ht> 650 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.986716; errup = 0.00311023; errdown = 0.00311023;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.539823; errup = 0.0034145; errdown = 0.0034145;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.571429; errup = 0.00358377; errdown = 0.00358377;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.61682; errup = 0.00383997; errdown = 0.00383997;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.65604; errup = 0.00399833; errdown = 0.00399833;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.675896; errup = 0.00424876; errdown = 0.00424876;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.719656; errup = 0.00434348; errdown = 0.00434348;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.735982; errup = 0.00455557; errdown = 0.00455557;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.765251; errup = 0.00464919; errdown = 0.00464919;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.773482; errup = 0.00482465; errdown = 0.00482465;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.806076; errup = 0.00486078; errdown = 0.00486078;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.833808; errup = 0.00351213; errdown = 0.00351213;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.858087; errup = 0.00366153; errdown = 0.00366153;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.874377; errup = 0.00379521; errdown = 0.00379521;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.891151; errup = 0.00394335; errdown = 0.00394335;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.909908; errup = 0.0039598; errdown = 0.0039598;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.920343; errup = 0.00270517; errdown = 0.00270517;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.940531; errup = 0.00282431; errdown = 0.00282431;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.943467; errup = 0.0025356; errdown = 0.0025356;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.952579; errup = 0.00332293; errdown = 0.00332293;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.968186; errup = 0.00382436; errdown = 0.00382436;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.969146; errup = 0.00486385; errdown = 0.00486385;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.987481; errup = 0.00248803; errdown = 0.00248803;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.432701; errup = 0.00130445; errdown = 0.00130445;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.468416; errup = 0.00141951; errdown = 0.00141951;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.500201; errup = 0.00153028; errdown = 0.00153028;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.534339; errup = 0.00164201; errdown = 0.00164201;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.565809; errup = 0.00176191; errdown = 0.00176191;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.597802; errup = 0.00187686; errdown = 0.00187686;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.624353; errup = 0.00198579; errdown = 0.00198579;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.657567; errup = 0.00209913; errdown = 0.00209913;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.683468; errup = 0.00220487; errdown = 0.00220487;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.707921; errup = 0.00231018; errdown = 0.00231018;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.736016; errup = 0.00174943; errdown = 0.00174943;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.77612; errup = 0.0018808; errdown = 0.0018808;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.81054; errup = 0.00199175; errdown = 0.00199175;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.831745; errup = 0.00215975; errdown = 0.00215975;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.855244; errup = 0.00226767; errdown = 0.00226767;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.882622; errup = 0.00156784; errdown = 0.00156784;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.90554; errup = 0.00178883; errdown = 0.00178883;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.927802; errup = 0.00151072; errdown = 0.00151072;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.941449; errup = 0.00195205; errdown = 0.00195205;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.944989; errup = 0.00260629; errdown = 0.00260629;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.954409; errup = 0.00317334; errdown = 0.00317334;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.960595; errup = 0.00225241; errdown = 0.00225241;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1el_trigeff2018("get_1el_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.0097629; errup = 0.00521891; errdown = 0.00359356;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.238559; errup = 0.0140446; errdown = 0.0135171;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 0 && met<= 50) {eff = 0.523018; errup = 0.00945585; errdown = 0.00947208;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.017301; errup = 0.00729105; errdown = 0.00535648;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.297573; errup = 0.0172358; errdown = 0.0166989;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.583857; errup = 0.0110342; errdown = 0.011117;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.260109; errup = 0.0153141; errdown = 0.0147704;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.444776; errup = 0.0162215; errdown = 0.0161081;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.650049; errup = 0.00869262; errdown = 0.00879159;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.387642; errup = 0.0066101; errdown = 0.00656955;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.596085; errup = 0.00792443; errdown = 0.00797407;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.751802; errup = 0.00372287; errdown = 0.00376028;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.386105; errup = 0.0171245; errdown = 0.0168562;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.643836; errup = 0.0259838; errdown = 0.0268136;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.754565; errup = 0.0118009; errdown = 0.0121843;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.761468; errup = 0.00860319; errdown = 0.00881662;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.800239; errup = 0.0140459; errdown = 0.0147947;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.87316; errup = 0.00563279; errdown = 0.0058494;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.00446429; errup = 0.0101903; errdown = 0.00369336;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.267442; errup = 0.0230545; errdown = 0.0219156;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.551826; errup = 0.0123281; errdown = 0.0123906;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.0408163; errup = 0.0195138; errdown = 0.0140021;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.307937; errup = 0.0282018; errdown = 0.0269105;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.622134; errup = 0.0139227; errdown = 0.0141213;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.33518; errup = 0.0266756; errdown = 0.0257133;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.524887; errup = 0.024796; errdown = 0.0249138;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.689604; errup = 0.010443; errdown = 0.0106352;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.490962; errup = 0.0112954; errdown = 0.0112864;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.663012; errup = 0.0115434; errdown = 0.0117363;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.785698; errup = 0.0044274; errdown = 0.00449406;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.470199; errup = 0.0304274; errdown = 0.0302191;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.636364; errup = 0.0394636; errdown = 0.0412433;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.793506; errup = 0.0148321; errdown = 0.0156277;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.790724; errup = 0.0138995; errdown = 0.014584;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.872; errup = 0.0175109; errdown = 0.0196448;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.897516; errup = 0.00642761; errdown = 0.00679249;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.0361446; errup = 0.0339112; errdown = 0.0196074;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.278607; errup = 0.0352149; errdown = 0.0328462;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.65252; errup = 0.0144558; errdown = 0.0147338;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.045977; errup = 0.0348613; errdown = 0.0218756;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.409639; errup = 0.0416348; errdown = 0.0404647;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.704675; errup = 0.0157244; errdown = 0.0162081;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.453988; errup = 0.0422296; errdown = 0.041623;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.61597; errup = 0.0313617; errdown = 0.0322973;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.763399; errup = 0.0110119; errdown = 0.0113664;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.626178; errup = 0.0160347; errdown = 0.0163076;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.76572; errup = 0.0137106; errdown = 0.0142695;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.825713; errup = 0.00481489; errdown = 0.00492038;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.610687; errup = 0.0453152; errdown = 0.0471475;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.775; errup = 0.0487807; errdown = 0.0565118;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.859083; errup = 0.0145352; errdown = 0.0158254;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.850498; errup = 0.0209401; errdown = 0.0234521;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.859459; errup = 0.0260778; errdown = 0.0303603;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.914036; errup = 0.00731012; errdown = 0.00789175;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.3125; errup = 0.0543931; errdown = 0.0500862;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.710778; errup = 0.0171191; errdown = 0.0177172;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0; errup = 0.0683597; errdown = 0;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.363636; errup = 0.0583279; errdown = 0.0549108;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.755848; errup = 0.0167559; errdown = 0.0175366;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.588235; errup = 0.0577422; errdown = 0.0600294;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.702479; errup = 0.0436783; errdown = 0.0473382;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.807843; errup = 0.012509; errdown = 0.0131362;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.705495; errup = 0.0219711; errdown = 0.0229195;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.808732; errup = 0.018281; errdown = 0.0196374;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.868353; errup = 0.0052035; errdown = 0.00537996;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.865385; errup = 0.0482807; errdown = 0.0649611;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.942308; errup = 0.0312343; errdown = 0.0529489;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.876147; errup = 0.0159935; errdown = 0.0178418;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.898438; errup = 0.0270367; errdown = 0.0341365;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.9; errup = 0.0338535; errdown = 0.0456142;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.933071; errup = 0.00788272; errdown = 0.00878826;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.1875; errup = 0.149399; errdown = 0.100212;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.583333; errup = 0.0918854; errdown = 0.0971962;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.799163; errup = 0.0186952; errdown = 0.0200167;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.238095; errup = 0.128988; errdown = 0.0987131;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.673077; errup = 0.0701338; errdown = 0.0776765;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.842784; errup = 0.0188114; errdown = 0.0206985;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.522727; errup = 0.0849595; errdown = 0.0861304;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.786885; errup = 0.0548156; errdown = 0.0655726;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.887324; errup = 0.0119924; errdown = 0.0131472;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.865546; errup = 0.0225156; errdown = 0.0258705;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.885621; errup = 0.0184512; errdown = 0.0211925;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.908734; errup = 0.00523456; errdown = 0.00551032;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.857143; errup = 0.0670785; errdown = 0.0986253;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.910448; errup = 0.0176206; errdown = 0.0210035;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.915254; errup = 0.0361154; errdown = 0.0532642;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.940299; errup = 0.0283546; errdown = 0.0446914;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.953344; errup = 0.00833639; errdown = 0.00987704;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.5; errup = 0.195182; errdown = 0.195182;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.871875; errup = 0.0189792; errdown = 0.021493;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.17521; errdown = 0.221361;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.772727; errup = 0.0944237; errdown = 0.12452;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.9375; errup = 0.014321; errdown = 0.0177218;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.764706; errup = 0.108959; errdown = 0.147312;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.913043; errup = 0.0411491; errdown = 0.0634308;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.924335; errup = 0.0120474; errdown = 0.013929;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.894309; errup = 0.028102; errdown = 0.0354144;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.94; errup = 0.0168114; errdown = 0.0218501;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.952295; errup = 0.00453493; errdown = 0.00496221;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.142229;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.969543; errup = 0.0120157; errdown = 0.0177485;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.969697; errup = 0.0250817; errdown = 0.0662602;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.972973; errup = 0.0223689; errdown = 0.0594217;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.970297; errup = 0.00752031; errdown = 0.00962462;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.666667; errup = 0.277375; errdown = 0.414535;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.8; errup = 0.0835235; errdown = 0.112668;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.952381; errup = 0.0116187; errdown = 0.0146509;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.8; errup = 0.0835235; errdown = 0.112668;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.983444; errup = 0.00713312; errdown = 0.0110449;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.9; errup = 0.0643201; errdown = 0.116971;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.918919; errup = 0.0438002; errdown = 0.0726265;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.965649; errup = 0.00793934; errdown = 0.00992767;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.957983; errup = 0.0180304; errdown = 0.027424;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.960199; errup = 0.0136568; errdown = 0.0190433;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.984898; errup = 0.00245839; errdown = 0.00288698;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.168149;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.95; errup = 0.0413995; errdown = 0.105764;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.975845; errup = 0.0103945; errdown = 0.0160098;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.994718; errup = 0.00287314; errdown = 0.0051109;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.15411;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0879414;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.991758; errup = 0.00448218; errdown = 0.00795189;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.264229;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.947368; errup = 0.0435805; errdown = 0.110836;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.988201; errup = 0.00563868; errdown = 0.00923117;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.184992;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.979592; errup = 0.0168888; errdown = 0.0453679;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.995208; errup = 0.00260705; errdown = 0.00463962;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0200277;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0104058;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.996068; errup = 0.00111791; errdown = 0.00148999;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.308024;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 0.993266; errup = 0.00434837; errdown = 0.00881245;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.102638;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0271039;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.99867; errup = 0.00110009; errdown = 0.00305118;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_1mu_trigeff2018("get_1mu_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig()), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.0509554; errup = 0.0103616; errdown = 0.00880306;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.406528; errup = 0.0160299; errdown = 0.0158392;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 50) {eff = 0.920266; errup = 0.00495963; errdown = 0.00524822;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.296203; errup = 0.024754; errdown = 0.0236669;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.57037; errup = 0.019671; errdown = 0.0198874;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.918153; errup = 0.00631369; errdown = 0.00677068;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.45394; errup = 0.0171845; errdown = 0.017079;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.734914; errup = 0.0120119; errdown = 0.0123595;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.94219; errup = 0.00339353; errdown = 0.0035852;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.904807; errup = 0.00908061; errdown = 0.00988317;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.939845; errup = 0.00607334; errdown = 0.00667259;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.953995; errup = 0.00248988; errdown = 0.00262123;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.0708333; errup = 0.0206706; errdown = 0.0166547;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.462222; errup = 0.0246733; errdown = 0.0244977;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.905621; errup = 0.00668024; errdown = 0.00711475;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.358382; errup = 0.0401117; errdown = 0.0383512;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.66548; errup = 0.0292488; errdown = 0.0304982;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.924117; errup = 0.00730077; errdown = 0.00797078;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.5629; errup = 0.0238039; errdown = 0.0240843;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.770013; errup = 0.015271; errdown = 0.0159857;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.950936; errup = 0.0036475; errdown = 0.00391321;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.923333; errup = 0.010941; errdown = 0.0124589;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.946875; errup = 0.00726638; errdown = 0.00826383;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.95772; errup = 0.00278904; errdown = 0.00297018;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.105263; errup = 0.0415156; errdown = 0.0318852;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.44898; errup = 0.0339661; errdown = 0.033524;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.924188; errup = 0.0071534; errdown = 0.0077968;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.410256; errup = 0.0629437; errdown = 0.0603993;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.71978; errup = 0.0346421; errdown = 0.0372475;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.93621; errup = 0.00752273; errdown = 0.00839199;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.681159; errup = 0.0291135; errdown = 0.0305079;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.833676; errup = 0.0171643; errdown = 0.0186162;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.953016; errup = 0.00387315; errdown = 0.00418804;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.918539; errup = 0.0146216; errdown = 0.0171965;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.956954; errup = 0.00826908; errdown = 0.00993237;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.965588; errup = 0.00275628; errdown = 0.00297726;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.111111; errup = 0.0791717; errdown = 0.0524058;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.519685; errup = 0.0478966; errdown = 0.0482324;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.938907; errup = 0.00787908; errdown = 0.00888368;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.68; errup = 0.102272; errdown = 0.119276;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.788462; errup = 0.0415895; errdown = 0.0477639;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.947566; errup = 0.00790478; errdown = 0.00911009;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.703125; errup = 0.042389; errdown = 0.0458523;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.882759; errup = 0.0191698; errdown = 0.0220444;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.965564; errup = 0.00386237; errdown = 0.00430287;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.976303; errup = 0.0101982; errdown = 0.0157124;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.971354; errup = 0.00844407; errdown = 0.01127;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.974014; errup = 0.00275157; errdown = 0.00305009;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.1; errup = 0.116971; errdown = 0.0643201;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.627119; errup = 0.0683561; errdown = 0.0731904;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.973875; errup = 0.00605301; errdown = 0.00759168;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.75; errup = 0.115499; errdown = 0.153966;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.869565; errup = 0.0504938; errdown = 0.0697763;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.97757; errup = 0.00634716; errdown = 0.00839086;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.923077; errup = 0.032824; errdown = 0.0486878;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.923858; errup = 0.01902; errdown = 0.0238853;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.980064; errup = 0.00353709; errdown = 0.0042104;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.963964; errup = 0.0171678; errdown = 0.0275761;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.995305; errup = 0.00388411; errdown = 0.0107125;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.984236; errup = 0.00249933; errdown = 0.00292235;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.875; errup = 0.103637; errdown = 0.23225;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.723404; errup = 0.0696186; errdown = 0.0805167;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.987923; errup = 0.00520699; errdown = 0.00808753;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.941176; errup = 0.048713; errdown = 0.12258;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.962963; errup = 0.0306592; errdown = 0.080075;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.980978; errup = 0.00698888; errdown = 0.0100953;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.957447; errup = 0.0202555; errdown = 0.0323679;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.981584; errup = 0.00406048; errdown = 0.00504708;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.967742; errup = 0.0208084; errdown = 0.0409676;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.986667; errup = 0.00860748; errdown = 0.0173148;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.992337; errup = 0.00201909; errdown = 0.00263025;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.9; errup = 0.082873; errdown = 0.194135;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.9; errup = 0.0539222; errdown = 0.0877974;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.989035; errup = 0.00472821; errdown = 0.00734954;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.857143; errup = 0.11848; errdown = 0.257124;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0659133;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.994819; errup = 0.00334597; errdown = 0.00679283;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.971429; errup = 0.0236478; errdown = 0.0626552;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.98913; errup = 0.00899357; errdown = 0.0245495;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.989756; errup = 0.00279678; errdown = 0.00367595;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.967213; errup = 0.0211491; errdown = 0.0416131;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0133482;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.99708; errup = 0.00115753; errdown = 0.0017398;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0769247;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00350724;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0802771;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.995204; errup = 0.00309728; errdown = 0.00629067;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0636358;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.016449;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.998025; errup = 0.00107467; errdown = 0.00191738;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0439098;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.01428;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998493; errup = 0.000720896; errdown = 0.00118964;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2el_trigeff2018("get_2el_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.897959; errup = 0.0249664; errdown = 0.0249664;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.950413; errup = 0.0139551; errdown = 0.0139551;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.968085; errup = 0.0104672; errdown = 0.0104672;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.95203; errup = 0.0129816; errdown = 0.0129816;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.92926; errup = 0.0145385; errdown = 0.0145385;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.953846; errup = 0.0130124; errdown = 0.0130124;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.962963; errup = 0.0114932; errdown = 0.0114932;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.960937; errup = 0.012109; errdown = 0.012109;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.952381; errup = 0.0140117; errdown = 0.0140117;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.961353; errup = 0.0133973; errdown = 0.0133973;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.975248; errup = 0.0109318; errdown = 0.0109318;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.978022; errup = 0.0108676; errdown = 0.0108676;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.979899; errup = 0.00994873; errdown = 0.00994873;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.987578; errup = 0.00872921; errdown = 0.00872921;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.996761; errup = 0.000784379; errdown = 0.000784379;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_2mu_trigeff2018("get_2mu_trigeff2018", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.984496; errup = 0.00769161; errdown = 0.00769161;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.960106; errup = 0.0100929; errdown = 0.0100929;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.988916; errup = 0.000322656; errdown = 0.000322656;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

const NamedFunc old_lepton_trig("old_lepton_trig",  [](const Baby &b) -> NamedFunc::ScalarType{
  return b.HLT_Ele27_WPTight_Gsf() || b.HLT_Ele35_WPTight_Gsf() || b.HLT_Ele115_CaloIdVT_GsfTrkIdT() || b.HLT_IsoMu24() || b.HLT_IsoMu27() || b.HLT_Mu50();
});

namespace{
  bool single_thread = false;
  string year_string = "2016";
  bool do_sr = true;
  bool do_qcd = true;
  bool do_ttbar = true;
  bool do_zll = true;
  bool unblind = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  //------------------------------------------------------------------------------------
  //                                 plot opts
  //------------------------------------------------------------------------------------

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_norm_nooverflow = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Overflow(OverflowType::none);
  PlotOpt lin_norm_nooverflow_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Overflow(OverflowType::none).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_lin_nooverflow = {lin_norm_nooverflow};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_lin_nooverflow = {lin_norm_nooverflow_data};
  if (unblind) plt_log = {log_norm_data};
  vector<PlotOpt> plt_lin_mc = {lin_norm};
  vector<PlotOpt> plt_log_mc = {log_norm};

  //------------------------------------------------------------------------------------
  //                                 samples
  //------------------------------------------------------------------------------------

  // Set options
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_klamath/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  //string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  //string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  //string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";
  string ttbar_mc_skim_folder = "mc/skim_higlep1T/";
  string zll_mc_skim_folder = "mc/skim_higlep2T/";
  string qcd_mc_skim_folder = "mc/skim_met150/";
  string met150_mc_skim_folder = "mc/skim_met150/";

  string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_klamath/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  //string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  //string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  //string qcd_data_skim_folder = "data/merged_higdata_higqcd/";
  string ttbar_data_skim_folder = "data/skim_higlep1T/";
  string zll_data_skim_folder = "data/skim_higlep2T/";
  string qcd_data_skim_folder = "data/skim_met150/";
  //string met150_data_skim_folder = "data/skim_met150/";

  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_klamath/";
  string search_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
  //string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";
  //string zll_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep2T/";
  //string qcd_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higqcd/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D/skim_higlep1T/";
  string zll_sig_skim_folder = "SMS-TChiHH_2D/skim_higlep2T/";
  string qcd_sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";
  string met150_sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  //years = {2016, 2017, 2018};
  //years = {2016};
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  // set qcd procs
  vector<shared_ptr<Process> > qcd_procs;
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["zjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["wjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["single_t"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["qcd"]),"stitch")); 
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    qcd_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
		    (met_trigger)));
  }

  // set ttbar procs
  vector<shared_ptr<Process> > ttbar_procs;
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
		    (el_trigger || mu_trigger || met_trigger)));
  }
  
  // set zll procs
  vector<shared_ptr<Process> > zll_procs;
  zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["zjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["wjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["single_t"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["qcd"]),"stitch")); 
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["other"]),"stitch"));
  if (unblind) {
    zll_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    (el_trigger || mu_trigger)));
  }

  // set sr procs
  vector<shared_ptr<Process> > sr_procs;
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder,mctags["zjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder,mctags["wjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["single_t"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["qcd"]),"stitch")); 
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["other"]),"stitch"));

  //set qcd procs btag
  vector<shared_ptr<Process> > qcd_procs_btag;
  qcd_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  qcd_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  qcd_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["all"]),"stitch&&(nbm==2)"));
  vector<shared_ptr<Process> > qcd_data_procs_btag;
  if (unblind) {
    qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("0b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
                    "(nbm==0)" && (met_trigger)
                    ));
    qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("1b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
                    "(nbm==1)" && (met_trigger)
                    ));
  } else {
    qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                    attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
    qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                    attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  }

  //set ttbar procs btag
  vector<shared_ptr<Process> > ttbar_procs_btag;
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  vector<shared_ptr<Process> > ttbar_data_procs_btag;
  if (unblind) {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("2b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt==2&&nbm==2)" && (el_trigger || mu_trigger || met_trigger)
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("3b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm==3&&nbl==3)" && (el_trigger || mu_trigger || met_trigger)
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("4b Data", Process::Type::data, kGreen,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm>=3&&nbl>=4)" && (el_trigger || mu_trigger || met_trigger)
                    ));
  } else {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  }

  //set zll procs btag
  vector<shared_ptr<Process> > zll_procs_btag;
  zll_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  zll_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  vector<shared_ptr<Process> > zll_data_procs_btag;
  if (unblind) {
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("0b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    "(nbm==0)" && (el_trigger || mu_trigger)
                    ));
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("1b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    "(nbm==1)" && (el_trigger || mu_trigger)
                    ));
  } else {
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                    attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  }

  //set sr procs btag
  vector<shared_ptr<Process> > sr_procs_btag;
  sr_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  sr_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, met150_mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  const NamedFunc trig_weights("trig_weights", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()/1000==0) return 1.0;
    float eff = 1.;
    if(b.nvlep()==0){ // search MC sample and qcd MC control sample
      if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
        if (b.SampleType()==2016) eff = get_0l_fakemet_trigeff2016.GetVector(b)[0];
        else if (b.SampleType()==2017) eff = get_0l_fakemet_trigeff2017.GetVector(b)[0];
        else if (b.SampleType()==2018) eff = get_0l_fakemet_trigeff2018.GetVector(b)[0];
      } else { // TRUE MET
        if (b.SampleType()==2016) eff = get_0l_trigeff2016.GetVector(b)[0];
        else if (b.SampleType()==2017) eff = get_0l_trigeff2017.GetVector(b)[0];
        else if (b.SampleType()==2018) eff = get_0l_trigeff2018.GetVector(b)[0];
      }
    } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
      if (b.SampleType()==2016) eff = get_1el_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_1el_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_1el_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
      if (b.SampleType()==2016) eff = get_1mu_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_1mu_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_1mu_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
      if (b.SampleType()==2016) eff = get_2el_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_2el_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_2el_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
      if (b.SampleType()==2016) eff = get_2mu_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_2mu_trigeff2017.GetVector(b)[0];
      else if (b.SampleType()==2018) eff = get_2mu_trigeff2018.GetVector(b)[0];
    }
    return eff;
  });

  const NamedFunc pass_filters_nolowneutraljet("pass_filters_nolowneutraljet", [](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    //if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //only for fastsim
    if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if ((abs(b.SampleType())==2018 && b.type()/1000 == 0) && !Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  const NamedFunc pass_filters_noecalnoisejet("pass_filters_noecalnoisejet", [](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    //if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    //if ((b.type()/1000 == 106)  && !b.pass_jets()) return false; //only for fastsim
    if (!b.pass_jets()) return false; //only for fastsim
    //if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if ((abs(b.SampleType())==2018 && b.type()/1000 == 0) && !Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  NamedFunc weight = "weight"*w_years*hem_weight*trig_weights;

  NamedFunc sr_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && !jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc zll_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);
  //dilepton mass cut in skim
  
  NamedFunc qcd_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc qcd_baseline_nolowneutraljet = pass_filters_nolowneutraljet && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc qcd_baseline_noecalnoisejet = pass_filters_noecalnoisejet && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);


  //------------------------------------------------------------------------------------
  //                                     plots
  //------------------------------------------------------------------------------------

  PlotMaker pm;

  // SR plots
  if (do_sr) {
    pm.Push<Hist1D>(Axis(20, 150, 800, "met", "p_{T}^{miss} [GeV]", {200,300,400}),
      sr_baseline,
      sr_procs, plt_log_mc).Weight(weight).Tag("FixName:validate_trig_sr_met_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, "ht", "H_{T} [GeV]", {}),
      sr_baseline,
      sr_procs, plt_lin_mc).Weight(weight).Tag("FixName:validate_trig_sr_ht_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(3, 1.5, 4.5, jetid_nb, "N_{b}", {2.5,3.5}),
      sr_baseline,
      sr_procs, plt_log_mc).Weight(weight).Tag("FixName:validate_trig_sr_nb_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, jetid_hig_cand_am, "<m_{bb}> [GeV]", {100,140}),
      sr_baseline,
      sr_procs, plt_lin_mc).Weight(weight).Tag("FixName:validate_trig_sr_am_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 100, jetid_hig_cand_dm, "#Delta m [GeV]", {40}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && !jetid_low_dphi_met &&
      (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      sr_procs, plt_lin_mc).Weight(weight).Tag("FixName:validate_trig_sr_dm_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && !jetid_low_dphi_met &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200),
      sr_procs, plt_lin_mc).Weight(weight).Tag("FixName:validate_trig_sr_drmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      sr_baseline,
      sr_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_sr_mc_shapes_"+year_string).LuminosityTag(total_luminosity_string);
  }

  // ttbar plots
  if (do_ttbar) {
    pm.Push<Hist1D>(Axis(20, 0, 800, "met", "p_{T}^{miss} [GeV]", {200,300,400}),
      ttbar_baseline,
      ttbar_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_ttbar_met_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 20, 300, lead_signal_lepton_pt, "leading p_{Tl} [GeV]", {30}),
      pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_ttbar_ptl_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, "ht", "H_{T} [GeV]", {}),
      ttbar_baseline,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_ttbar_ht_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(3, 1.5, 4.5, jetid_nb, "N_{b}", {2.5,3.5}),
      ttbar_baseline,
      ttbar_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_ttbar_nb_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, jetid_hig_cand_am, "<m_{bb}> [GeV]", {100,140}),
      ttbar_baseline,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_ttbar_am_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 100, jetid_hig_cand_dm, "#Delta m [GeV]", {40}),
      pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
      (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_ttbar_dm_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
      pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200),
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_ttbar_drmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_baseline,
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_ttbar_mc_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      ttbar_baseline,
      ttbar_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_ttbar_data_shapes_"+year_string).LuminosityTag(total_luminosity_string);
  }

  // QCD plots
  if (do_qcd) {
    pm.Push<Hist1D>(Axis(20, 150, 800, "met", "p_{T}^{miss} [GeV]", {200,300,400}),
      qcd_baseline,
      qcd_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_qcd_met_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1400, "ht", "H_{T} [GeV]", {}),
      qcd_baseline,
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_ht_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1400, "ht", "H_{T} [GeV]", {}),
      qcd_baseline_nolowneutraljet,
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_ht_nolowneutraljet_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1400, "ht", "H_{T} [GeV]", {}),
      qcd_baseline_noecalnoisejet,
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_ht_noecalnoisejet_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, -0.5, 4.5, jetid_nb, "N_{b}", {1.5,2.5,3.5}),
      qcd_baseline,
      qcd_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_qcd_nb_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, jetid_hig_cand_am, "<m_{bb}> [GeV]", {100,140}),
      qcd_baseline,
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_am_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 100, jetid_hig_cand_dm, "#Delta m [GeV]", {40}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
      (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
      (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_dm_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
      (jetid_njet>=4) && (jetid_njet<=5) && jetid_low_dphi_met &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200),
      qcd_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_qcd_drmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      qcd_baseline,
      qcd_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_qcd_mc_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      qcd_baseline,
      qcd_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_qcd_data_shapes_"+year_string).LuminosityTag(total_luminosity_string);
  }

  // Zll plots
  if (do_zll) {
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "nlep==2",
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow1_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2",
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow2_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    //  (jetid_njet>=4) && (jetid_njet<=5),
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow3_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    //  (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40),
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow4_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    //  (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    //  (jetid_hig_cand_dm<40),
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow5_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    //  (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    //  (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200),
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow6_"+year_string).LuminosityTag(total_luminosity_string);
    //pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
    //  pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2" &&
    //  (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    //  (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
    //  zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_cutflow7_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 400, "ll_pt[0]", "p_{Tll} [GeV]", {150,200,300,400}),
      zll_baseline,
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_met_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 20, 300, lead_signal_lepton_pt, "leading p_{Tl} [GeV]", {40}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met<50&&nlep==2" &&
      (jetid_njet>=4) && (jetid_njet<=5) &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_ptl_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 1000, "ht", "H_{T} [GeV]", {}),
      zll_baseline,
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_ht_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, -0.5, 4.5, jetid_nb, "N_{b}", {1.5,2.5,3.5}),
      zll_baseline,
      zll_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_zll_nb_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nbm", "N_{bm}", {1.5,2.5,3.5}),
      zll_baseline,
      zll_procs, plt_log).Weight(weight).Tag("FixName:validate_trig_zll_nbm_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, jetid_hig_cand_am, "<m_{bb}> [GeV]", {100,140}),
      zll_baseline,
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_am_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 100, jetid_hig_cand_dm, "#Delta m [GeV]", {40}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met<50&&nlep==2" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
      (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2),
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_dm_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {1.1,2.2}),
      pass_filters && "met/met_calo<2&&met/mht<2&&met<50&&nlep==2" &&
      (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
      (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200),
      zll_procs, plt_lin).Weight(weight).Tag("FixName:validate_trig_zll_drmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      zll_baseline,
      zll_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_zll_mc_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      zll_baseline,
      zll_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:validate_trig_zll_data_shapes_"+year_string).LuminosityTag(total_luminosity_string);
  }
  
  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"year", required_argument, 0, 0},
      {"unblind", no_argument, 0, 0},
      {"nosr", no_argument, 0, 0},
      {"noqcd", no_argument, 0, 0},
      {"nottbar", no_argument, 0, 0},
      {"nozll", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "year"){
        year_string = optarg;
      } else if (optname == "unblind") {
        unblind = true;
      } else if (optname == "nosr") {
        do_sr = false;
      } else if (optname == "noqcd") {
        do_qcd = false;
      } else if (optname == "nottbar") {
        do_ttbar = false;
      } else if (optname == "nozll") {
        do_zll = false;
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
