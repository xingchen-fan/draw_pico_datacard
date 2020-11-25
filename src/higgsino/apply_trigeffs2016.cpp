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


const NamedFunc get_0l_trigeff2016_v0("get_0l_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.421454; errup = 0.0122918; errdown = 0.0121974;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.511723; errup = 0.00553945; errdown = 0.0055423;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.524272; errup = 0.0229051; errdown = 0.0230034;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.502959; errup = 0.0412691; errdown = 0.0413067;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.575342; errup = 0.0630717; errdown = 0.0653595;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.478716; errup = 0.0135556; errdown = 0.0135252;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.58184; errup = 0.00573815; errdown = 0.00576013;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.523909; errup = 0.0237337; errdown = 0.0238375;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.564706; errup = 0.0404303; errdown = 0.0412473;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.5; errup = 0.0653624; errdown = 0.0653624;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.537445; errup = 0.0151993; errdown = 0.0152672;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.641288; errup = 0.00593655; errdown = 0.00597961;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.629712; errup = 0.023519; errdown = 0.0241214;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.580645; errup = 0.0474312; errdown = 0.0488429;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.682927; errup = 0.0786376; errdown = 0.0888868;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.556346; errup = 0.0176245; errdown = 0.0177627;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.699494; errup = 0.00600351; errdown = 0.00607174;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.699752; errup = 0.0235238; errdown = 0.0245674;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.59854; errup = 0.0445755; errdown = 0.0461321;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.680851; errup = 0.0733593; errdown = 0.0821322;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.638268; errup = 0.0184422; errdown = 0.0188431;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.749666; errup = 0.00603513; errdown = 0.00613206;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.733918; errup = 0.0246003; errdown = 0.0260493;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.72807; errup = 0.0436702; errdown = 0.0480675;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.574468; errup = 0.0799773; errdown = 0.083562;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.654701; errup = 0.02022; errdown = 0.0207713;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.784363; errup = 0.00597389; errdown = 0.00609421;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.75641; errup = 0.0249947; errdown = 0.0267432;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.734513; errup = 0.0434929; errdown = 0.0480594;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.634615; errup = 0.0726695; errdown = 0.0785167;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.741722; errup = 0.0210861; errdown = 0.022209;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.822244; errup = 0.00596603; errdown = 0.00612396;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.819672; errup = 0.0225045; errdown = 0.0247522;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.798319; errup = 0.0380713; errdown = 0.0436274;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.702703; errup = 0.0810427; errdown = 0.0937327;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.645533; errup = 0.0266391; errdown = 0.0275234;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.853995; errup = 0.00578957; errdown = 0.00598182;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.832753; errup = 0.0224984; errdown = 0.0249952;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.826531; errup = 0.0394452; errdown = 0.0469475;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.771429; errup = 0.0749545; errdown = 0.0932188;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.714286; errup = 0.0269016; errdown = 0.0284131;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.874115; errup = 0.00585518; errdown = 0.00609152;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.849817; errup = 0.022045; errdown = 0.0248186;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.818182; errup = 0.0425143; errdown = 0.0506834;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.880952; errup = 0.0504413; errdown = 0.0725148;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.75; errup = 0.0307629; errdown = 0.0333008;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.903775; errup = 0.00541855; errdown = 0.00569667;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.876068; errup = 0.0218991; errdown = 0.0254228;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.835443; errup = 0.042987; errdown = 0.0526851;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.870968; errup = 0.0607061; errdown = 0.0903257;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.794118; errup = 0.0236918; errdown = 0.0257459;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.926773; errup = 0.00363813; errdown = 0.00380788;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.915074; errup = 0.0129577; errdown = 0.0148636;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.9; errup = 0.0232988; errdown = 0.0285863;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.830769; errup = 0.0480266; errdown = 0.0597643;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.806122; errup = 0.0290347; errdown = 0.0324296;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.955567; errup = 0.00326322; errdown = 0.00349935;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.945026; errup = 0.0116963; errdown = 0.0142838;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.938356; errup = 0.0198781; errdown = 0.0268833;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.941176; errup = 0.0318434; errdown = 0.0539242;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.860656; errup = 0.0320534; errdown = 0.0387031;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.96675; errup = 0.00316544; errdown = 0.00346976;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.967262; errup = 0.00964044; errdown = 0.0128421;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.949153; errup = 0.0199859; errdown = 0.0291361;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.934426; errup = 0.0311194; errdown = 0.0488167;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.844444; errup = 0.0392546; errdown = 0.0479576;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.976041; errup = 0.00303115; errdown = 0.00342929;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.964968; errup = 0.0103099; errdown = 0.0137191;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.953488; errup = 0.0221284; errdown = 0.0352492;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.0449824;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.857143; errup = 0.0511428; errdown = 0.0684449;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.982009; errup = 0.00296575; errdown = 0.00348796;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.978102; errup = 0.00865233; errdown = 0.0128502;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0180628;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.964286; errup = 0.0295635; errdown = 0.077387;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.873239; errup = 0.0402311; errdown = 0.0525007;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.987697; errup = 0.00186062; errdown = 0.00216111;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.986667; errup = 0.00490435; errdown = 0.0071078;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.995238; errup = 0.00393961; errdown = 0.0108643;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0239329;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.870968; errup = 0.0607061; errdown = 0.0903257;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.992795; errup = 0.00183562; errdown = 0.00237029;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.99505; errup = 0.00319693; errdown = 0.00649192;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.984848; errup = 0.0097805; errdown = 0.0196341;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.983871; errup = 0.0133466; errdown = 0.0361115;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 9999) {eff = 0.846154; errup = 0.0721243; errdown = 0.105033;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.995964; errup = 0.00139558; errdown = 0.00198449;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.998106; errup = 0.00156681; errdown = 0.00434157;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.983784; errup = 0.00881271; errdown = 0.0155222;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.987179; errup = 0.0106082; errdown = 0.0288622;}
  return eff;
});

const NamedFunc get_0l_fakemet_trigeff2016_v0("get_0l_fakemet_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.34375; errup = 0.0551838; errdown = 0.0516032;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.329126; errup = 0.00748338; errdown = 0.00739907;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.340784; errup = 0.00565815; errdown = 0.00561368;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.318462; errup = 0.00311091; errdown = 0.00309492;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.267794; errup = 0.00247084; errdown = 0.00245656;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.30303; errup = 0.066999; errdown = 0.0602987;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.361364; errup = 0.00827805; errdown = 0.00819792;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.377637; errup = 0.00635196; errdown = 0.00631077;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.359295; errup = 0.00349402; errdown = 0.00347927;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.30763; errup = 0.00282814; errdown = 0.00281386;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.28; errup = 0.0777781; errdown = 0.067729;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.4; errup = 0.0093527; errdown = 0.00928171;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.418022; errup = 0.00701682; errdown = 0.00698425;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.407438; errup = 0.00392148; errdown = 0.00390982;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.346154; errup = 0.00319023; errdown = 0.00317655;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.186047; errup = 0.0788698; errdown = 0.0616867;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.458898; errup = 0.0104849; errdown = 0.0104493;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.456691; errup = 0.00778; errdown = 0.00775924;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.436683; errup = 0.00430337; errdown = 0.00429393;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.379798; errup = 0.00353981; errdown = 0.00352716;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.277778; errup = 0.0944452; errdown = 0.0800801;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.47762; errup = 0.0115283; errdown = 0.011505;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.496503; errup = 0.00850186; errdown = 0.00849987;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.470105; errup = 0.00468488; errdown = 0.00467967;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.421681; errup = 0.00393157; errdown = 0.00392175;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.324324; errup = 0.0946276; errdown = 0.0836685;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.508585; errup = 0.0124515; errdown = 0.0124619;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.532835; errup = 0.00914594; errdown = 0.00916767;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.511022; errup = 0.00506406; errdown = 0.0050663;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.45314; errup = 0.00425117; errdown = 0.0042444;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.423077; errup = 0.116981; errdown = 0.110076;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.535191; errup = 0.0138396; errdown = 0.0138925;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.562858; errup = 0.00986103; errdown = 0.00991;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.54331; errup = 0.0054092; errdown = 0.00541933;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.490859; errup = 0.00467494; errdown = 0.00467336;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.578947; errup = 0.1303; errdown = 0.140168;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.557904; errup = 0.0148052; errdown = 0.0149059;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.595641; errup = 0.0105275; errdown = 0.0106144;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.57839; errup = 0.00581006; errdown = 0.00583159;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.526404; errup = 0.00499295; errdown = 0.00499819;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.416667; errup = 0.183083; errdown = 0.166178;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.543564; errup = 0.0161155; errdown = 0.0162044;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.616804; errup = 0.0108963; errdown = 0.0110123;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.605413; errup = 0.0062108; errdown = 0.00624459;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.566151; errup = 0.00533269; errdown = 0.00534789;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.6; errup = 0.145412; errdown = 0.161475;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.579744; errup = 0.0173172; errdown = 0.0175091;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.632; errup = 0.0117334; errdown = 0.0118879;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.628879; errup = 0.00647922; errdown = 0.00652528;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.593062; errup = 0.0056954; errdown = 0.00572023;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.666667; errup = 0.135514; errdown = 0.162494;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.607942; errup = 0.0133945; errdown = 0.0135546;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.646354; errup = 0.00896477; errdown = 0.0090669;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.666143; errup = 0.00486205; errdown = 0.00489725;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.627819; errup = 0.00436025; errdown = 0.00438096;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.684211; errup = 0.117339; errdown = 0.140625;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.659864; errup = 0.0150872; errdown = 0.015408;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.66796; errup = 0.00991916; errdown = 0.0100672;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.711698; errup = 0.00522735; errdown = 0.00528364;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.673255; errup = 0.00471847; errdown = 0.00475342;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.375; errup = 0.156049; errdown = 0.137247;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.666667; errup = 0.0172876; errdown = 0.0177309;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.687773; errup = 0.0109918; errdown = 0.0112019;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.734524; errup = 0.00580289; errdown = 0.00588383;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.704943; errup = 0.00519146; errdown = 0.00524446;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.454545; errup = 0.189662; errdown = 0.179582;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.709898; errup = 0.0192106; errdown = 0.0199588;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.692726; errup = 0.0122306; errdown = 0.0124999;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.760484; errup = 0.00621321; errdown = 0.00632373;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.754891; errup = 0.00546306; errdown = 0.00554532;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.9; errup = 0.082873; errdown = 0.194135;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.696121; errup = 0.0219634; errdown = 0.0228502;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.71073; errup = 0.0135217; errdown = 0.0138949;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.79273; errup = 0.00655972; errdown = 0.00671362;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.774911; errup = 0.00590379; errdown = 0.00601387;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.71466; errup = 0.01669; errdown = 0.0172739;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.709596; errup = 0.0103435; errdown = 0.0105604;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.807188; errup = 0.0047099; errdown = 0.00479779;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.815172; errup = 0.00401727; errdown = 0.00408503;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0; errup = 0.308024; errdown = 0;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.694118; errup = 0.0230215; errdown = 0.0239817;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.722388; errup = 0.0124316; errdown = 0.0127735;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.845351; errup = 0.00547786; errdown = 0.00563772;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.851237; errup = 0.00456296; errdown = 0.00467933;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.777778; errup = 0.142118; errdown = 0.221429;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.735802; errup = 0.0225042; errdown = 0.0237324;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.716895; errup = 0.0126343; errdown = 0.0129743;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.862974; errup = 0.00495149; errdown = 0.00510348;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.874064; errup = 0.00403698; errdown = 0.00414878;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1; errup = 0; errdown = 0.458642;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.709677; errup = 0.0380961; errdown = 0.0410288;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.748673; errup = 0.0186589; errdown = 0.0195808;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.889896; errup = 0.00670688; errdown = 0.00707213;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.908459; errup = 0.00503657; errdown = 0.00529077;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.78481; errup = 0.0482156; errdown = 0.0563445;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.792531; errup = 0.0268436; errdown = 0.0294573;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.909605; errup = 0.00886479; errdown = 0.00967703;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.933579; errup = 0.00620482; errdown = 0.00676451;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.5; errup = 0.417248; errdown = 0.417248;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.772727; errup = 0.0944237; errdown = 0.12452;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.87037; errup = 0.0329774; errdown = 0.0407707;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.94084; errup = 0.0103533; errdown = 0.0121897;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.941714; errup = 0.00795617; errdown = 0.009038;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 9.53547e-322; errup = 8.74496e-322; errdown = 8.74496e-322;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.714286; errup = 0.182129; errdown = 0.259938;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.971014; errup = 0.0187; errdown = 0.0369543;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.961165; errup = 0.0109425; errdown = 0.0143599;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.948919; errup = 0.00978575; errdown = 0.0117251;}
  return eff;
});

const NamedFunc get_1el_trigeff2016_v0("get_1el_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 0 && met<= 110) {eff = 0.0582445; errup = 0.00513861; errdown = 0.00476062;}
  else if (el_pt> 25 && el_pt<= 30 && met> 0 && met<= 110) {eff = 0.157254; errup = 0.00903781; errdown = 0.00864611;}
  else if (el_pt> 30 && el_pt<= 110 && met> 0 && met<= 110) {eff = 0.530553; errup = 0.00429568; errdown = 0.00430018;}
  else if (el_pt> 110 && el_pt<= 120 && met> 0 && met<= 110) {eff = 0.650194; errup = 0.0151726; errdown = 0.015473;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 110) {eff = 0.89841; errup = 0.0059834; errdown = 0.00630236;}
  else if (el_pt> 20 && el_pt<= 25 && met> 110 && met<= 120) {eff = 0.149123; errup = 0.0411162; errdown = 0.0341911;}
  else if (el_pt> 25 && el_pt<= 30 && met> 110 && met<= 120) {eff = 0.338028; errup = 0.0652582; errdown = 0.0601709;}
  else if (el_pt> 30 && el_pt<= 110 && met> 110 && met<= 120) {eff = 0.672389; errup = 0.0181985; errdown = 0.0187112;}
  else if (el_pt> 110 && el_pt<= 120 && met> 110 && met<= 120) {eff = 0.686275; errup = 0.06986; errdown = 0.078157;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 110 && met<= 120) {eff = 0.969231; errup = 0.0146686; errdown = 0.02366;}
  else if (el_pt> 20 && el_pt<= 25 && met> 120 && met<= 130) {eff = 0.2; errup = 0.049646; errdown = 0.0425693;}
  else if (el_pt> 25 && el_pt<= 30 && met> 120 && met<= 130) {eff = 0.314286; errup = 0.0651353; errdown = 0.0592025;}
  else if (el_pt> 30 && el_pt<= 110 && met> 120 && met<= 130) {eff = 0.684039; errup = 0.0192486; errdown = 0.0198735;}
  else if (el_pt> 110 && el_pt<= 120 && met> 120 && met<= 130) {eff = 0.886364; errup = 0.0481931; errdown = 0.0695663;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 120 && met<= 130) {eff = 0.95082; errup = 0.0193366; errdown = 0.02822;}
  else if (el_pt> 20 && el_pt<= 25 && met> 130 && met<= 140) {eff = 0.425287; errup = 0.0593282; errdown = 0.0574387;}
  else if (el_pt> 25 && el_pt<= 30 && met> 130 && met<= 140) {eff = 0.513514; errup = 0.0641717; errdown = 0.0645758;}
  else if (el_pt> 30 && el_pt<= 110 && met> 130 && met<= 140) {eff = 0.723894; errup = 0.0192634; errdown = 0.0200923;}
  else if (el_pt> 110 && el_pt<= 120 && met> 130 && met<= 140) {eff = 0.909091; errup = 0.0429965; errdown = 0.0660622;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 130 && met<= 140) {eff = 0.892157; errup = 0.0311385; errdown = 0.0399966;}
  else if (el_pt> 20 && el_pt<= 25 && met> 140 && met<= 150) {eff = 0.453333; errup = 0.0643662; errdown = 0.0629896;}
  else if (el_pt> 25 && el_pt<= 30 && met> 140 && met<= 150) {eff = 0.464286; errup = 0.0755319; errdown = 0.0741026;}
  else if (el_pt> 30 && el_pt<= 110 && met> 140 && met<= 150) {eff = 0.765531; errup = 0.0193923; errdown = 0.0205114;}
  else if (el_pt> 110 && el_pt<= 120 && met> 140 && met<= 150) {eff = 0.844444; errup = 0.0555301; errdown = 0.0737002;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 140 && met<= 150) {eff = 0.947917; errup = 0.0223143; errdown = 0.0336975;}
  else if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 160) {eff = 0.606557; errup = 0.0681557; errdown = 0.0720644;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 160) {eff = 0.731707; errup = 0.0738471; errdown = 0.0869218;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 160) {eff = 0.816934; errup = 0.0188601; errdown = 0.0203981;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 160) {eff = 0.765957; errup = 0.0651312; errdown = 0.0782067;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 160) {eff = 0.987013; errup = 0.010746; errdown = 0.029229;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 170) {eff = 0.555556; errup = 0.0827174; errdown = 0.0855148;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 170) {eff = 0.75; errup = 0.0815095; errdown = 0.0999181;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 170) {eff = 0.836364; errup = 0.0192059; errdown = 0.0210707;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 170) {eff = 0.965517; errup = 0.0285434; errdown = 0.0748731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 170) {eff = 1; errup = 0; errdown = 0.0287997;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 180) {eff = 0.625; errup = 0.0842352; errdown = 0.0913836;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 180) {eff = 0.729167; errup = 0.0683073; errdown = 0.0792509;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 180) {eff = 0.869333; errup = 0.0176707; errdown = 0.019788;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 180) {eff = 0.92; errup = 0.051501; errdown = 0.0959187;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 180) {eff = 0.947368; errup = 0.02502; errdown = 0.0396606;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 190) {eff = 0.771429; errup = 0.0749545; errdown = 0.0932188;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 190) {eff = 0.897143; errup = 0.0164322; errdown = 0.018897;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 190) {eff = 1; errup = 0; errdown = 0.0738409;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 190) {eff = 0.984848; errup = 0.0125375; errdown = 0.0339781;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 200) {eff = 0.759259; errup = 0.0613529; errdown = 0.0723448;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 200) {eff = 0.75; errup = 0.0767393; errdown = 0.0929856;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 200) {eff = 0.907534; errup = 0.0171339; errdown = 0.0202008;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 200) {eff = 0.851852; errup = 0.0695101; errdown = 0.101731;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 200) {eff = 0.969231; errup = 0.0198493; errdown = 0.0391458;}
  else if (el_pt> 20 && el_pt<= 25 && met> 200 && met<= 210) {eff = 0.857143; errup = 0.0551652; errdown = 0.0755709;}
  else if (el_pt> 25 && el_pt<= 30 && met> 200 && met<= 210) {eff = 0.916667; errup = 0.0450072; errdown = 0.0744674;}
  else if (el_pt> 30 && el_pt<= 110 && met> 200 && met<= 210) {eff = 0.94386; errup = 0.0136633; errdown = 0.0171713;}
  else if (el_pt> 110 && el_pt<= 120 && met> 200 && met<= 210) {eff = 0.9375; errup = 0.051761; errdown = 0.129429;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 210) {eff = 1; errup = 0; errdown = 0.0292574;}
  else if (el_pt> 20 && el_pt<= 25 && met> 210 && met<= 9999) {eff = 0.918919; errup = 0.026025; errdown = 0.034831;}
  else if (el_pt> 25 && el_pt<= 30 && met> 210 && met<= 9999) {eff = 0.865169; errup = 0.0369968; errdown = 0.0464051;}
  else if (el_pt> 30 && el_pt<= 110 && met> 210 && met<= 9999) {eff = 0.965271; errup = 0.00585514; errdown = 0.0068881;}
  else if (el_pt> 110 && el_pt<= 120 && met> 210 && met<= 9999) {eff = 0.989796; errup = 0.00844284; errdown = 0.0230719;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 210 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00822172;}
  return eff;
});

const NamedFunc get_1mu_trigeff2016_v0("get_1mu_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 110) {eff = 0.115848; errup = 0.00702761; errdown = 0.00668609;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 110) {eff = 0.50953; errup = 0.0129078; errdown = 0.0129202;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 110) {eff = 0.646394; errup = 0.00726734; errdown = 0.00733457;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 110) {eff = 0.893396; errup = 0.00388525; errdown = 0.00401146;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 110 && met<= 120) {eff = 0.343284; errup = 0.0458476; errdown = 0.043308;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 110 && met<= 120) {eff = 0.709091; errup = 0.0455501; errdown = 0.049726;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 110 && met<= 120) {eff = 0.771331; errup = 0.0252171; errdown = 0.0271938;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 110 && met<= 120) {eff = 0.938636; errup = 0.0114961; errdown = 0.0136851;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 120 && met<= 130) {eff = 0.408333; errup = 0.0496489; errdown = 0.0479882;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 120 && met<= 130) {eff = 0.757895; errup = 0.0459589; errdown = 0.0519967;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 120 && met<= 130) {eff = 0.825758; errup = 0.0238784; errdown = 0.0265402;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 120 && met<= 130) {eff = 0.930946; errup = 0.0129022; errdown = 0.0153227;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 130 && met<= 140) {eff = 0.406977; errup = 0.0596105; errdown = 0.0572279;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 130 && met<= 140) {eff = 0.80303; errup = 0.0509296; errdown = 0.0614258;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 130 && met<= 140) {eff = 0.88; errup = 0.0220141; errdown = 0.0257281;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 130 && met<= 140) {eff = 0.950685; errup = 0.011346; errdown = 0.0141095;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 140 && met<= 150) {eff = 0.559524; errup = 0.0589468; errdown = 0.0605076;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 140 && met<= 150) {eff = 0.814286; errup = 0.0481926; errdown = 0.0584532;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 140 && met<= 150) {eff = 0.873786; errup = 0.0235445; errdown = 0.027542;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 140 && met<= 150) {eff = 0.951429; errup = 0.0114981; errdown = 0.0143934;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 160) {eff = 0.597015; errup = 0.0652245; errdown = 0.0684491;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 160) {eff = 0.90566; errup = 0.0401398; errdown = 0.0587823;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 160) {eff = 0.867299; errup = 0.0237886; errdown = 0.027612;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 160) {eff = 0.968254; errup = 0.00979478; errdown = 0.0132384;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 170) {eff = 0.716418; errup = 0.0584326; errdown = 0.065696;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 170) {eff = 0.93617; errup = 0.0345373; errdown = 0.0582116;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 170) {eff = 0.955128; errup = 0.0164018; errdown = 0.0233331;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 170) {eff = 0.9699; errup = 0.00977778; errdown = 0.0134446;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 180) {eff = 0.784615; errup = 0.0532994; errdown = 0.0632744;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 180) {eff = 0.872727; errup = 0.0457212; errdown = 0.061809;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 180) {eff = 0.951515; errup = 0.016605; errdown = 0.0230425;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 180) {eff = 0.977358; errup = 0.00894499; errdown = 0.0132786;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 190) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 190) {eff = 0.983051; errup = 0.0140254; errdown = 0.037896;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 190) {eff = 0.94702; errup = 0.0181267; errdown = 0.0250908;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 190) {eff = 0.987069; errup = 0.00702944; errdown = 0.0124182;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 200) {eff = 0.863636; errup = 0.0527265; errdown = 0.0725604;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0449824;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 200) {eff = 0.934211; errup = 0.0201296; errdown = 0.0267482;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 200) {eff = 0.988142; errup = 0.00644659; errdown = 0.0113996;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 210) {eff = 0.853659; errup = 0.0564709; errdown = 0.0771695;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 210) {eff = 0.916667; errup = 0.0450072; errdown = 0.0744674;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 210) {eff = 0.974359; errup = 0.0139228; errdown = 0.024313;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 210) {eff = 0.979381; errup = 0.00984218; errdown = 0.0160023;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 210 && met<= 9999) {eff = 0.92562; errup = 0.0239126; errdown = 0.0321199;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 210 && met<= 9999) {eff = 0.955752; errup = 0.018981; errdown = 0.0288242;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 210 && met<= 9999) {eff = 0.976959; errup = 0.00712399; errdown = 0.00966982;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.993455; errup = 0.00282399; errdown = 0.00440304;}
  return eff;
});

const NamedFunc get_2el_trigeff2016_v0("get_2el_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.819209; errup = 0.0289267; errdown = 0.0289267;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.887324; errup = 0.0216654; errdown = 0.0216654;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.865741; errup = 0.0231974; errdown = 0.0231974;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.880597; errup = 0.0198075; errdown = 0.0198075;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.906383; errup = 0.019002; errdown = 0.019002;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.851351; errup = 0.0238758; errdown = 0.0238758;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.885; errup = 0.0225583; errdown = 0.0225583;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.920455; errup = 0.0203964; errdown = 0.0203964;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.917197; errup = 0.021994; errdown = 0.021994;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.902439; errup = 0.0267544; errdown = 0.0267544;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.936508; errup = 0.0217235; errdown = 0.0217235;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.936508; errup = 0.0217235; errdown = 0.0217235;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.895238; errup = 0.0298866; errdown = 0.0298866;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.988395; errup = 0.00210652; errdown = 0.00210652;}
  return eff;
});

const NamedFunc get_2mu_trigeff2016_v0("get_2mu_trigeff2016_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.933594; errup = 0.0155619; errdown = 0.0155619;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.968288; errup = 0.00805725; errdown = 0.00805725;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.972572; errup = 0.000505438; errdown = 0.000505438;}
  return eff;
});

}


