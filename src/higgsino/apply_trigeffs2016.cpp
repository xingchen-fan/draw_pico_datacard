#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/apply_trigeffs2016.hpp"

namespace Higfuncs{

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

}
