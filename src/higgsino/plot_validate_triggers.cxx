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
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.616091; errup = 0.040626; errdown = 0.0408214;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.630731; errup = 0.0417236; errdown = 0.0417503;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.654605; errup = 0.0436736; errdown = 0.0437819;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.658147; errup = 0.0459552; errdown = 0.0471318;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.62381; errup = 0.0528005; errdown = 0.0532081;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.640449; errup = 0.0678307; errdown = 0.0858621;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.55665; errup = 0.0409086; errdown = 0.0400387;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.672944; errup = 0.0442055; errdown = 0.0442663;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.693822; errup = 0.0455462; errdown = 0.0455777;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.691021; errup = 0.0458807; errdown = 0.0459657;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.621951; errup = 0.0445279; errdown = 0.045123;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.711628; errup = 0.0552757; errdown = 0.0574471;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.579545; errup = 0.0677244; errdown = 0.0841137;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.62419; errup = 0.032886; errdown = 0.0300531;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.717157; errup = 0.029367; errdown = 0.0293312;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.733951; errup = 0.0298227; errdown = 0.0299625;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.740271; errup = 0.0306371; errdown = 0.0310111;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.726916; errup = 0.0337854; errdown = 0.0349358;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.702247; errup = 0.0442877; errdown = 0.0442432;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.703704; errup = 0.0714117; errdown = 0.0801084;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.674352; errup = 0.0361548; errdown = 0.0315178;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.768562; errup = 0.0310229; errdown = 0.0310815;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.771899; errup = 0.0313409; errdown = 0.0312548;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.771282; errup = 0.0317878; errdown = 0.0320403;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.774947; errup = 0.0349064; errdown = 0.0363371;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.710983; errup = 0.0446227; errdown = 0.0444359;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.655738; errup = 0.0699558; errdown = 0.0795429;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.784946; errup = 0.0386148; errdown = 0.0348829;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.80334; errup = 0.0323132; errdown = 0.0322761;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.818095; errup = 0.0327434; errdown = 0.0328414;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.827869; errup = 0.0334729; errdown = 0.0339392;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.808786; errup = 0.0363654; errdown = 0.0373441;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.798658; errup = 0.0450919; errdown = 0.0464521;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.740741; errup = 0.0689347; errdown = 0.0805693;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.765625; errup = 0.0394371; errdown = 0.0342832;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.829294; errup = 0.0332461; errdown = 0.0331725;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.859719; errup = 0.0338653; errdown = 0.0342819;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.868649; errup = 0.0344868; errdown = 0.0353221;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.813299; errup = 0.0363376; errdown = 0.0374796;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.858156; errup = 0.0438374; errdown = 0.0478987;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.703704; errup = 0.0714117; errdown = 0.0801084;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.784884; errup = 0.0374544; errdown = 0.0268282;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.867628; errup = 0.0242716; errdown = 0.0241506;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.880851; errup = 0.0237118; errdown = 0.024501;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.868759; errup = 0.0245513; errdown = 0.0252806;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.860724; errup = 0.0278056; errdown = 0.0302826;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.837209; errup = 0.0389778; errdown = 0.0409844;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.76; errup = 0.0663161; errdown = 0.0778814;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.759124; errup = 0.0421401; errdown = 0.0263964;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.89931; errup = 0.0243848; errdown = 0.0248113;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.897311; errup = 0.0240833; errdown = 0.0248437;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.904271; errup = 0.024543; errdown = 0.0259902;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.880117; errup = 0.0276601; errdown = 0.030603;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.857143; errup = 0.0387544; errdown = 0.0412218;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.826923; errup = 0.0578011; errdown = 0.0782744;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.873874; errup = 0.0383798; errdown = 0.0283735;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.890302; errup = 0.0248328; errdown = 0.0246228;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.914483; errup = 0.0243483; errdown = 0.0252029;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.932874; errup = 0.0247376; errdown = 0.026568;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.886435; errup = 0.0279373; errdown = 0.0307082;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.85567; errup = 0.0419407; errdown = 0.0412041;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.84; errup = 0.0570871; errdown = 0.0783549;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.777778; errup = 0.0456141; errdown = 0.0267084;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.934489; errup = 0.0249785; errdown = 0.0255515;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.949025; errup = 0.0243608; errdown = 0.0259304;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.926621; errup = 0.0247748; errdown = 0.0264412;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.903226; errup = 0.0281424; errdown = 0.0309896;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.844444; errup = 0.0441926; errdown = 0.0410701;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.880952; errup = 0.0547064; errdown = 0.0786143;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.844961; errup = 0.0492065; errdown = 0.0414398;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.942308; errup = 0.0417775; errdown = 0.0427968;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.94065; errup = 0.0415114; errdown = 0.0427708;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.952336; errup = 0.0419764; errdown = 0.0438442;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.952813; errup = 0.0424617; errdown = 0.0470105;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.9; errup = 0.0449583; errdown = 0.0530117;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.836066; errup = 0.0609471; errdown = 0.0840053;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.932584; errup = 0.0484419; errdown = 0.0448612;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.958926; errup = 0.0410742; errdown = 0.043491;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.967811; errup = 0.0321888; errdown = 0.0439045;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.975419; errup = 0.024581; errdown = 0.0447959;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.960784; errup = 0.0392157; errdown = 0.0473171;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.945652; errup = 0.0444227; errdown = 0.0544975;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.940299; errup = 0.0498009; errdown = 0.0860691;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.955556; errup = 0.0444444; errdown = 0.0457685;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.970588; errup = 0.0294118; errdown = 0.0439787;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.984869; errup = 0.0151307; errdown = 0.0446178;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.983247; errup = 0.0167526; errdown = 0.0451191;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.968668; errup = 0.0313316; errdown = 0.0476209;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.966216; errup = 0.0337838; errdown = 0.0551771;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.970149; errup = 0.0298507; errdown = 0.0866949;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.928571; errup = 0.0471241; errdown = 0.0216438;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.973134; errup = 0.0138291; errdown = 0.0162242;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.984686; errup = 0.0118413; errdown = 0.0164212;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.98609; errup = 0.0117738; errdown = 0.0179173;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.982143; errup = 0.0129215; errdown = 0.0246179;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.964286; errup = 0.0200601; errdown = 0.0372502;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 1; errup = 0; errdown = 0.0765059;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0220268;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.985915; errup = 0.0132899; errdown = 0.0163176;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.992537; errup = 0.00746269; errdown = 0.0164785;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.983986; errup = 0.0120326; errdown = 0.0179032;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.983221; errup = 0.0130238; errdown = 0.0246232;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 1; errup = 0; errdown = 0.0373644;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.969697; errup = 0.0272627; errdown = 0.0764585;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0209918;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.994065; errup = 0.00593472; errdown = 0.0149743;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.990919; errup = 0.00908059; errdown = 0.0150809;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 1; errup = 0; errdown = 0.0167291;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.994774; errup = 0.00522648; errdown = 0.0237701;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.995708; errup = 0.00429185; errdown = 0.0367548;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.98913; errup = 0.0108696; errdown = 0.0762035;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0209918;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.992188; errup = 0.0078125; errdown = 0.0149647;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.998155; errup = 0.00184502; errdown = 0.0151176;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.997222; errup = 0.00277778; errdown = 0.0167163;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0237869;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.993464; errup = 0.00653595; errdown = 0.0367502;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.986486; errup = 0.0135135; errdown = 0.0762009;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0192981;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0125259;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.997773; errup = 0.00222717; errdown = 0.0126706;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0145472;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0223064;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.987124; errup = 0.00757714; errdown = 0.0358206;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.990291; errup = 0.00854428; errdown = 0.0757643;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2016("get_0l_fakemet_trigeff2016", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 600 && met> 150 && met<= 155) {eff = 0.431034; errup = 0.010393; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.40057; errup = 0.00491586; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.390395; errup = 0.0024779; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.317039; errup = 0.00267592; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 155 && met<= 160) {eff = 0.458289; errup = 0.0118091; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.43192; errup = 0.00545693; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.430988; errup = 0.0027588; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.354034; errup = 0.00302543; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 160 && met<= 165) {eff = 0.512071; errup = 0.0129047; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.474905; errup = 0.00599445; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.477876; errup = 0.00303967; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.391472; errup = 0.00337327; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 165 && met<= 170) {eff = 0.535201; errup = 0.0140684; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.501234; errup = 0.00649546; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.512433; errup = 0.00332853; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.422406; errup = 0.00372064; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 170 && met<= 175) {eff = 0.574391; errup = 0.0152219; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.526943; errup = 0.00706778; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.550285; errup = 0.00358547; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.457542; errup = 0.00408883; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 175 && met<= 180) {eff = 0.559499; errup = 0.0164885; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.563263; errup = 0.00755281; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.588906; errup = 0.00385581; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.489956; errup = 0.0043993; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 180 && met<= 185) {eff = 0.627717; errup = 0.0183057; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.568481; errup = 0.00806442; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.617598; errup = 0.00408758; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.521893; errup = 0.00480928; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 185 && met<= 190) {eff = 0.627422; errup = 0.0191984; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.6; errup = 0.00847349; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.648006; errup = 0.00434604; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.557526; errup = 0.00511848; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 190 && met<= 195) {eff = 0.598182; errup = 0.0216079; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.614921; errup = 0.00902727; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.67869; errup = 0.00457369; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.592459; errup = 0.00544189; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 195 && met<= 200) {eff = 0.654397; errup = 0.0221754; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.617969; errup = 0.00958198; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.701713; errup = 0.00476542; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.611142; errup = 0.00583631; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 200 && met<= 210) {eff = 0.663265; errup = 0.0172894; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.647235; errup = 0.00733313; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.730722; errup = 0.00358954; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.641662; errup = 0.00444985; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 210 && met<= 220) {eff = 0.71777; errup = 0.0192453; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.661012; errup = 0.00819274; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.773799; errup = 0.00379672; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.686066; errup = 0.00481554; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 220 && met<= 230) {eff = 0.684327; errup = 0.0224902; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.672107; errup = 0.00915992; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.799541; errup = 0.00411108; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.711314; errup = 0.0053228; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 230 && met<= 240) {eff = 0.71547; errup = 0.0244348; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.696579; errup = 0.00982075; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.822887; errup = 0.00437569; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.757922; errup = 0.00563096; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 240 && met<= 250) {eff = 0.653285; errup = 0.0299272; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.708986; errup = 0.0110992; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.841388; errup = 0.00467685; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.783099; errup = 0.00605753; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 250 && met<= 275) {eff = 0.648649; errup = 0.02246; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.719601; errup = 0.00827798; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.854577; errup = 0.003362; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.815239; errup = 0.00417164; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 275 && met<= 300) {eff = 0.618375; errup = 0.0301455; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.735676; errup = 0.0103896; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.887328; errup = 0.00386041; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.850877; errup = 0.00474297; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 300 && met<= 350) {eff = 0.604895; errup = 0.0302126; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.752329; errup = 0.0102341; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.90065; errup = 0.00349259; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.874438; errup = 0.0042177; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 350 && met<= 400) {eff = 0.572727; errup = 0.0507375; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.802198; errup = 0.0150091; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.924596; errup = 0.00458628; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.910239; errup = 0.00523751; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 400 && met<= 450) {eff = 0.72; errup = 0.067729; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.85; errup = 0.0203299; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.943602; errup = 0.00570312; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.938462; errup = 0.00638369; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 450 && met<= 500) {eff = 0.875; errup = 0.103637; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.89781; errup = 0.0262134; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.963095; errup = 0.00650946; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.948784; errup = 0.0079152; errdown = 0.00371192;}
  else if (ht> 0 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.333333; errup = 0.414535; errdown = 0.0103338;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.976471; errup = 0.0151834; errdown = 0.00544063;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.978678; errup = 0.00659505; errdown = 0.00303804;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.954348; errup = 0.00974214; errdown = 0.00371192;}
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
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.197796; errup = 0.00954788; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.522956; errup = 0.0172376; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.592268; errup = 0.00609976; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.343518; errup = 0.00991985; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.632466; errup = 0.0136146; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.676785; errup = 0.00477725; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.51545; errup = 0.00513898; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.734544; errup = 0.0056898; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.756135; errup = 0.00212987; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.76007; errup = 0.0128345; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.85582; errup = 0.0129388; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.860034; errup = 0.00541667; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.851656; errup = 0.00672062; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.884747; errup = 0.00665368; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.90492; errup = 0.00282916; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.207113; errup = 0.0202127; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.564516; errup = 0.0268312; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.564735; errup = 0.00789192; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.371345; errup = 0.0279443; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.664093; errup = 0.0305449; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.651934; errup = 0.00869838; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.523629; errup = 0.0225912; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.722334; errup = 0.0206066; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.728115; errup = 0.00644023; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.656417; errup = 0.0101746; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.792651; errup = 0.00856531; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.802124; errup = 0.00282306; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.829493; errup = 0.0261432; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.852459; errup = 0.0231571; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.879075; errup = 0.00753224; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.866873; errup = 0.013536; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.90665; errup = 0.0104923; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.91716; errup = 0.00409646; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.283465; errup = 0.0456264; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.75; errup = 0.0330874; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.664792; errup = 0.00995901; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.578512; errup = 0.0481009; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.741722; errup = 0.0370708; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.735936; errup = 0.0107074; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.703883; errup = 0.033095; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.814545; errup = 0.0239949; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.78274; errup = 0.00762671; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.784091; errup = 0.0140963; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.871389; errup = 0.0103228; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.849313; errup = 0.00324759; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.870968; errup = 0.0354651; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.966942; errup = 0.015755; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.91099; errup = 0.00864424; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.893333; errup = 0.0208625; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.926509; errup = 0.013463; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.944102; errup = 0.00452597; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.5; errup = 0.0780973; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.739726; errup = 0.0541361; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.746437; errup = 0.0120965; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.54; errup = 0.0786763; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.794872; errup = 0.0475795; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.819343; errup = 0.0117716; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.843478; errup = 0.0347781; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.908397; errup = 0.0254458; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.842222; errup = 0.00867254; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.864662; errup = 0.0173909; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.926448; errup = 0.0103967; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.89891; errup = 0.00351003; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.966667; errup = 0.0275914; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.942529; errup = 0.0246014; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.9377; errup = 0.00971033; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.966292; errup = 0.0182896; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.950495; errup = 0.0152082; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.965256; errup = 0.00460922; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.851852; errup = 0.0695101; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.936508; errup = 0.0301398; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.853598; errup = 0.0126086; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.909091; errup = 0.0490632; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 1; errup = 0; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.871069; errup = 0.0134541; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.883333; errup = 0.0420083; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.929825; errup = 0.0332829; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.876818; errup = 0.00970247; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.920354; errup = 0.0181519; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.948276; errup = 0.0103093; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.935926; errup = 0.00348682; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.952381; errup = 0.0394264; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.981481; errup = 0.0153245; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.971074; errup = 0.00757748; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.981481; errup = 0.0153245; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.951613; errup = 0.0190274; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.97007; errup = 0.00505453; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.857143; errup = 0.0917089; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.925703; errup = 0.0118357; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.932367; errup = 0.0124155; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.810811; errup = 0.0670146; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.980769; errup = 0.0159141; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.950556; errup = 0.00764631; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.972603; errup = 0.0130667; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.977564; errup = 0.00823774; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.960628; errup = 0.00329173; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.975385; errup = 0.00847391; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.990654; errup = 0.00773258; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988032; errup = 0.00390372; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 200) {eff = 0.945455; errup = 0.0116067; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 200) {eff = 0.954545; errup = 0.0376328; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 200) {eff = 0.971963; errup = 0.00911194; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 200) {eff = 0.975332; errup = 0.0067088; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 200) {eff = 0.990826; errup = 0.00759067; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 200) {eff = 0.985782; errup = 0.00772816; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 200) {eff = 0.977743; errup = 0.00291454; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 200) {eff = 0.99569; errup = 0.00356599; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 200) {eff = 0.990307; errup = 0.00383826; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 200 && met<= 250) {eff = 0.985; errup = 0.00593413; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 200 && met<= 250) {eff = 0.983562; errup = 0.00650149; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 200 && met<= 250) {eff = 0.99572; errup = 0.00232823; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 200 && met<= 250) {eff = 0.995935; errup = 0.00336304; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 200 && met<= 250) {eff = 0.991002; errup = 0.0016287; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 200 && met<= 250) {eff = 0.99449; errup = 0.00355793; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 250) {eff = 0.99639; errup = 0.00196412; errdown = 0.00541802;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00563482;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0133946;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 250 && met<= 9999) {eff = 0.995833; errup = 0.00344712; errdown = 0.00507683;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0271513;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0318921;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 250 && met<= 9999) {eff = 0.986667; errup = 0.00724787; errdown = 0.00879899;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.035219;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0264544;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 250 && met<= 9999) {eff = 0.992629; errup = 0.00400894; errdown = 0.00782077;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0193496;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0118313;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 250 && met<= 9999) {eff = 0.997955; errup = 0.000978452; errdown = 0.00361948;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 250 && met<= 9999) {eff = 1.01909e-312; errup = 0; errdown = 0.101134;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0412969;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 250 && met<= 9999) {eff = 0.995215; errup = 0.00395846; errdown = 0.00978541;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0472931;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0211614;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00541802;}
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
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.433962; errup = 0.0171325; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.867868; errup = 0.0132858; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.908151; errup = 0.00408863; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.631429; errup = 0.0117375; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.929134; errup = 0.00660721; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.935745; errup = 0.0022462; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.938377; errup = 0.00596572; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.962165; errup = 0.00418374; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.961258; errup = 0.00142423; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.379603; errup = 0.0275657; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.91954; errup = 0.0147031; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.898394; errup = 0.00492584; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.595238; errup = 0.0324355; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.91791; errup = 0.0169135; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.928719; errup = 0.00479442; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.734667; errup = 0.0164532; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.927536; errup = 0.00906659; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.94032; errup = 0.0026961; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.948035; errup = 0.00793019; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.965458; errup = 0.00537137; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.963843; errup = 0.00173444; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.553333; errup = 0.0433853; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.930233; errup = 0.0194967; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.915444; errup = 0.00580904; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.675214; errup = 0.0457257; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.951515; errup = 0.016605; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.94313; errup = 0.00521512; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.848397; errup = 0.019717; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.973843; errup = 0.00711092; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.950236; errup = 0.00288025; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.946067; errup = 0.0107415; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.980488; errup = 0.00479992; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.967188; errup = 0.00193449; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.671233; errup = 0.0587767; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.914634; errup = 0.0309441; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.933731; errup = 0.00643465; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.775862; errup = 0.0574412; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.948454; errup = 0.0220862; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.958838; errup = 0.00565631; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.918478; errup = 0.0203311; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.947368; errup = 0.0132382; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.966378; errup = 0.00292606; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.983264; errup = 0.00799298; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.985102; errup = 0.00513915; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.97674; errup = 0.00198722; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.828571; errup = 0.0658118; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.955556; errup = 0.0286549; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.970874; errup = 0.00551699; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.973684; errup = 0.02178; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.964286; errup = 0.0230346; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.976055; errup = 0.00514512; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.94382; errup = 0.0240536; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.987952; errup = 0.00777825; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.983293; errup = 0.00258303; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.987261; errup = 0.0082239; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.988604; errup = 0.00544618; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.982685; errup = 0.00206541; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.987118; errup = 0.0044459; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.978873; errup = 0.00598041; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.952381; errup = 0.0258048; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.987744; errup = 0.00258316; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988764; errup = 0.00929678; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.986425; errup = 0.0073789; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.994496; errup = 0.00135936; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 200) {eff = 0.962963; errup = 0.0306592; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 200) {eff = 0.981735; errup = 0.00629624; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 200) {eff = 0.997283; errup = 0.00224806; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 200) {eff = 0.988095; errup = 0.00985028; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 200) {eff = 0.996252; errup = 0.00161803; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 200) {eff = 0.978723; errup = 0.0176077; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 200) {eff = 1; errup = 0; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 200) {eff = 0.99587; errup = 0.00134959; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 250) {eff = 0.990566; errup = 0.004069; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 250) {eff = 0.993789; errup = 0.00337848; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 250) {eff = 0.998163; errup = 0.000999657; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 250) {eff = 1; errup = 0; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 250) {eff = 0.998949; errup = 0.000572023; errdown = 0.00215996;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0119246;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0104239;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00347436;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0332379;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0203713;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00510027;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0218975;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00928315;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 250 && met<= 9999) {eff = 0.998957; errup = 0.00086263; errdown = 0.00304225;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0130352;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00726487;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 250 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00215996;}
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
  else if (ht> 250 && ht<= 350 && met> 150 && met<= 155) {eff = 0.247156; errup = 0.0263155; errdown = 0.0266406;}
  else if (ht> 350 && ht<= 450 && met> 150 && met<= 155) {eff = 0.313043; errup = 0.0331531; errdown = 0.033967;}
  else if (ht> 450 && ht<= 600 && met> 150 && met<= 155) {eff = 0.359584; errup = 0.0378066; errdown = 0.0386964;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.350832; errup = 0.0388109; errdown = 0.0418421;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.412698; errup = 0.0524354; errdown = 0.0569822;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.571429; errup = 0.0793748; errdown = 0.0842014;}
  else if (ht> 0 && ht<= 250 && met> 155 && met<= 160) {eff = 0.197836; errup = 0.0257522; errdown = 0.0235749;}
  else if (ht> 250 && ht<= 350 && met> 155 && met<= 160) {eff = 0.291738; errup = 0.0307667; errdown = 0.0306785;}
  else if (ht> 350 && ht<= 450 && met> 155 && met<= 160) {eff = 0.372247; errup = 0.0389705; errdown = 0.0393015;}
  else if (ht> 450 && ht<= 600 && met> 155 && met<= 160) {eff = 0.412712; errup = 0.0431359; errdown = 0.043505;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.457831; errup = 0.0492426; errdown = 0.0508362;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.424; errup = 0.0534294; errdown = 0.057776;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.42268; errup = 0.0695934; errdown = 0.0752684;}
  else if (ht> 0 && ht<= 250 && met> 160 && met<= 165) {eff = 0.208633; errup = 0.0216225; errdown = 0.0172532;}
  else if (ht> 250 && ht<= 350 && met> 160 && met<= 165) {eff = 0.350923; errup = 0.0223027; errdown = 0.0213953;}
  else if (ht> 350 && ht<= 450 && met> 160 && met<= 165) {eff = 0.427288; errup = 0.0266864; errdown = 0.0266198;}
  else if (ht> 450 && ht<= 600 && met> 160 && met<= 165) {eff = 0.451982; errup = 0.0281213; errdown = 0.028461;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.508527; errup = 0.0335266; errdown = 0.0356607;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.423645; errup = 0.0435078; errdown = 0.0457888;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.5; errup = 0.0635302; errdown = 0.0680247;}
  else if (ht> 0 && ht<= 250 && met> 165 && met<= 170) {eff = 0.291572; errup = 0.0278471; errdown = 0.0202741;}
  else if (ht> 250 && ht<= 350 && met> 165 && met<= 170) {eff = 0.376812; errup = 0.0243114; errdown = 0.0225664;}
  else if (ht> 350 && ht<= 450 && met> 165 && met<= 170) {eff = 0.488679; errup = 0.0300519; errdown = 0.0293645;}
  else if (ht> 450 && ht<= 600 && met> 165 && met<= 170) {eff = 0.516588; errup = 0.0313049; errdown = 0.0313213;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.515094; errup = 0.0351539; errdown = 0.0359173;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.593909; errup = 0.0482373; errdown = 0.0506957;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.551282; errup = 0.0679688; errdown = 0.0690991;}
  else if (ht> 0 && ht<= 250 && met> 170 && met<= 175) {eff = 0.29316; errup = 0.0321667; errdown = 0.0203366;}
  else if (ht> 250 && ht<= 350 && met> 170 && met<= 175) {eff = 0.449911; errup = 0.028093; errdown = 0.0259689;}
  else if (ht> 350 && ht<= 450 && met> 170 && met<= 175) {eff = 0.562955; errup = 0.0337113; errdown = 0.0327983;}
  else if (ht> 450 && ht<= 600 && met> 170 && met<= 175) {eff = 0.583333; errup = 0.0345929; errdown = 0.034375;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.553942; errup = 0.037316; errdown = 0.0374639;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.670103; errup = 0.049731; errdown = 0.0532282;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.662162; errup = 0.0682863; errdown = 0.0717102;}
  else if (ht> 0 && ht<= 250 && met> 175 && met<= 180) {eff = 0.33617; errup = 0.0379333; errdown = 0.0220804;}
  else if (ht> 250 && ht<= 350 && met> 175 && met<= 180) {eff = 0.493603; errup = 0.0308202; errdown = 0.0280542;}
  else if (ht> 350 && ht<= 450 && met> 175 && met<= 180) {eff = 0.575588; errup = 0.0345602; errdown = 0.0333921;}
  else if (ht> 450 && ht<= 600 && met> 175 && met<= 180) {eff = 0.616114; errup = 0.0365028; errdown = 0.0359038;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.691837; errup = 0.0420474; errdown = 0.0432779;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.664706; errup = 0.051484; errdown = 0.0530429;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.671875; errup = 0.0720834; errdown = 0.0719567;}
  else if (ht> 0 && ht<= 250 && met> 180 && met<= 185) {eff = 0.424893; errup = 0.0489984; errdown = 0.036998;}
  else if (ht> 250 && ht<= 350 && met> 180 && met<= 185) {eff = 0.52524; errup = 0.0462387; errdown = 0.0440442;}
  else if (ht> 350 && ht<= 450 && met> 180 && met<= 185) {eff = 0.598726; errup = 0.051834; errdown = 0.0507232;}
  else if (ht> 450 && ht<= 600 && met> 180 && met<= 185) {eff = 0.672794; errup = 0.0571534; errdown = 0.0568839;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.659955; errup = 0.0583602; errdown = 0.0586192;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.654545; errup = 0.0658697; errdown = 0.0665619;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.764706; errup = 0.0822803; errdown = 0.088313;}
  else if (ht> 0 && ht<= 250 && met> 185 && met<= 190) {eff = 0.384181; errup = 0.0507023; errdown = 0.0339368;}
  else if (ht> 250 && ht<= 350 && met> 185 && met<= 190) {eff = 0.568249; errup = 0.0501677; errdown = 0.0474334;}
  else if (ht> 350 && ht<= 450 && met> 185 && met<= 190) {eff = 0.672603; errup = 0.0574368; errdown = 0.0564977;}
  else if (ht> 450 && ht<= 600 && met> 185 && met<= 190) {eff = 0.682065; errup = 0.0581009; errdown = 0.0576071;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.75; errup = 0.0648599; errdown = 0.0653694;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.683544; errup = 0.0677075; errdown = 0.0684565;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.71875; errup = 0.0834579; errdown = 0.0857308;}
  else if (ht> 0 && ht<= 250 && met> 190 && met<= 195) {eff = 0.397163; errup = 0.0556659; errdown = 0.0349081;}
  else if (ht> 250 && ht<= 350 && met> 190 && met<= 195) {eff = 0.601351; errup = 0.0530601; errdown = 0.0500515;}
  else if (ht> 350 && ht<= 450 && met> 190 && met<= 195) {eff = 0.702663; errup = 0.0598218; errdown = 0.0588602;}
  else if (ht> 450 && ht<= 600 && met> 190 && met<= 195) {eff = 0.741214; errup = 0.0627835; errdown = 0.0622378;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.744318; errup = 0.06499; errdown = 0.0649399;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.738255; errup = 0.0707083; errdown = 0.0721047;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.72; errup = 0.0894686; errdown = 0.0857999;}
  else if (ht> 0 && ht<= 250 && met> 195 && met<= 200) {eff = 0.552632; errup = 0.0673381; errdown = 0.0468189;}
  else if (ht> 250 && ht<= 350 && met> 195 && met<= 200) {eff = 0.629699; errup = 0.0555027; errdown = 0.0522992;}
  else if (ht> 350 && ht<= 450 && met> 195 && met<= 200) {eff = 0.7552; errup = 0.0637816; errdown = 0.0630034;}
  else if (ht> 450 && ht<= 600 && met> 195 && met<= 200) {eff = 0.790774; errup = 0.0663718; errdown = 0.0661367;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.784741; errup = 0.0673963; errdown = 0.0680048;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.761905; errup = 0.0733868; errdown = 0.0737086;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.857143; errup = 0.0844708; errdown = 0.0937416;}
  else if (ht> 0 && ht<= 250 && met> 200 && met<= 210) {eff = 0.535484; errup = 0.046646; errdown = 0.0226955;}
  else if (ht> 250 && ht<= 350 && met> 200 && met<= 210) {eff = 0.702381; errup = 0.0289533; errdown = 0.0264565;}
  else if (ht> 350 && ht<= 450 && met> 200 && met<= 210) {eff = 0.798263; errup = 0.0301237; errdown = 0.0309384;}
  else if (ht> 450 && ht<= 600 && met> 200 && met<= 210) {eff = 0.823843; errup = 0.0304719; errdown = 0.0323719;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.82148; errup = 0.0317993; errdown = 0.0368323;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.887446; errup = 0.0370015; errdown = 0.0502992;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.789474; errup = 0.0512002; errdown = 0.0683763;}
  else if (ht> 0 && ht<= 250 && met> 210 && met<= 220) {eff = 0.607143; errup = 0.0611004; errdown = 0.0247211;}
  else if (ht> 250 && ht<= 350 && met> 210 && met<= 220) {eff = 0.769231; errup = 0.0322748; errdown = 0.0285542;}
  else if (ht> 350 && ht<= 450 && met> 210 && met<= 220) {eff = 0.845422; errup = 0.0315886; errdown = 0.0323745;}
  else if (ht> 450 && ht<= 600 && met> 210 && met<= 220) {eff = 0.891209; errup = 0.0322527; errdown = 0.0344016;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.864407; errup = 0.0332178; errdown = 0.0379673;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.921951; errup = 0.0367854; errdown = 0.051022;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.9; errup = 0.0474975; errdown = 0.0699597;}
  else if (ht> 0 && ht<= 250 && met> 220 && met<= 230) {eff = 0.673913; errup = 0.0782297; errdown = 0.0266734;}
  else if (ht> 250 && ht<= 350 && met> 220 && met<= 230) {eff = 0.823666; errup = 0.0338487; errdown = 0.030283;}
  else if (ht> 350 && ht<= 450 && met> 220 && met<= 230) {eff = 0.883436; errup = 0.0328181; errdown = 0.033544;}
  else if (ht> 450 && ht<= 600 && met> 220 && met<= 230) {eff = 0.902622; errup = 0.0326708; errdown = 0.0347489;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.929314; errup = 0.0339327; errdown = 0.0397251;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.914634; errup = 0.0382833; errdown = 0.0508673;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.940476; errup = 0.0410651; errdown = 0.0705815;}
  else if (ht> 0 && ht<= 250 && met> 230 && met<= 240) {eff = 0.73913; errup = 0.103; errdown = 0.0338081;}
  else if (ht> 250 && ht<= 350 && met> 230 && met<= 240) {eff = 0.865517; errup = 0.0416831; errdown = 0.037993;}
  else if (ht> 350 && ht<= 450 && met> 230 && met<= 240) {eff = 0.933712; errup = 0.040711; errdown = 0.041815;}
  else if (ht> 450 && ht<= 600 && met> 230 && met<= 240) {eff = 0.92515; errup = 0.040197; errdown = 0.0419819;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.942356; errup = 0.0412866; errdown = 0.0461783;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.942197; errup = 0.0433697; errdown = 0.0563299;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.981481; errup = 0.0185185; errdown = 0.0751297;}
  else if (ht> 0 && ht<= 250 && met> 240 && met<= 250) {eff = 0.631579; errup = 0.127314; errdown = 0.0297117;}
  else if (ht> 250 && ht<= 350 && met> 240 && met<= 250) {eff = 0.89083; errup = 0.0428662; errdown = 0.039012;}
  else if (ht> 350 && ht<= 450 && met> 240 && met<= 250) {eff = 0.946602; errup = 0.0412909; errdown = 0.0423234;}
  else if (ht> 450 && ht<= 600 && met> 240 && met<= 250) {eff = 0.948998; errup = 0.0409667; errdown = 0.0429112;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.971722; errup = 0.0282776; errdown = 0.0472403;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.972222; errup = 0.0277778; errdown = 0.0572234;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.969697; errup = 0.030303; errdown = 0.0748591;}
  else if (ht> 0 && ht<= 250 && met> 250 && met<= 275) {eff = 0.8125; errup = 0.10127; errdown = 0.0197939;}
  else if (ht> 250 && ht<= 350 && met> 250 && met<= 275) {eff = 0.934783; errup = 0.0224715; errdown = 0.0200817;}
  else if (ht> 350 && ht<= 450 && met> 250 && met<= 275) {eff = 0.959212; errup = 0.018765; errdown = 0.0225093;}
  else if (ht> 450 && ht<= 600 && met> 250 && met<= 275) {eff = 0.970817; errup = 0.0182131; errdown = 0.0235792;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.976934; errup = 0.0183933; errdown = 0.0295477;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.987755; errup = 0.0122449; errdown = 0.0438301;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.967742; errup = 0.0232094; errdown = 0.065165;}
  else if (ht> 0 && ht<= 250 && met> 275 && met<= 300) {eff = 1; errup = 0; errdown = 0.0223939;}
  else if (ht> 250 && ht<= 350 && met> 275 && met<= 300) {eff = 0.95614; errup = 0.0254774; errdown = 0.0204037;}
  else if (ht> 350 && ht<= 450 && met> 275 && met<= 300) {eff = 0.985366; errup = 0.0146341; errdown = 0.022871;}
  else if (ht> 450 && ht<= 600 && met> 275 && met<= 300) {eff = 0.989101; errup = 0.0108992; errdown = 0.0238232;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.994286; errup = 0.00571429; errdown = 0.0297339;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.994764; errup = 0.0052356; errdown = 0.0438813;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.980198; errup = 0.019802; errdown = 0.065225;}
  else if (ht> 0 && ht<= 250 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0378196;}
  else if (ht> 250 && ht<= 350 && met> 300 && met<= 9999) {eff = 0.962963; errup = 0.037037; errdown = 0.0358028;}
  else if (ht> 350 && ht<= 450 && met> 300 && met<= 9999) {eff = 0.990132; errup = 0.00986842; errdown = 0.037904;}
  else if (ht> 450 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.993781; errup = 0.00621891; errdown = 0.0385727;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.995399; errup = 0.00460123; errdown = 0.0424868;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0534581;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.99115; errup = 0.00884956; errdown = 0.0719287;}
  std::vector<double> ret = {eff, errup, errdown};
  return ret;
});

const NamedFunc get_0l_fakemet_trigeff2017("get_0l_fakemet_trigeff2017", [](const Baby &b) -> NamedFunc::VectorType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 600 && met> 150 && met<= 155) {eff = 0.179933; errup = 0.00556135; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.30717; errup = 0.00629238; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.489496; errup = 0.00372032; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.390291; errup = 0.00251209; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 155 && met<= 160) {eff = 0.212876; errup = 0.0068006; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.348242; errup = 0.00725148; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.522675; errup = 0.00398564; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.421634; errup = 0.00275205; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 160 && met<= 165) {eff = 0.231979; errup = 0.00789299; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.39797; errup = 0.00815974; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.568861; errup = 0.00424739; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.455842; errup = 0.0030279; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 165 && met<= 170) {eff = 0.297178; errup = 0.00990641; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.42671; errup = 0.00911225; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.592914; errup = 0.00454107; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.489976; errup = 0.00329713; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 170 && met<= 175) {eff = 0.329863; errup = 0.0113723; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.450157; errup = 0.0100702; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.635472; errup = 0.00473045; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.531632; errup = 0.00358512; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 175 && met<= 180) {eff = 0.399045; errup = 0.0131971; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.499315; errup = 0.0109079; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.674262; errup = 0.0049525; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.557535; errup = 0.00385331; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 180 && met<= 185) {eff = 0.425022; errup = 0.0152308; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.548507; errup = 0.0117264; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.686402; errup = 0.00521825; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.592875; errup = 0.00415504; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 185 && met<= 190) {eff = 0.47245; errup = 0.0177025; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.590937; errup = 0.0124981; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.726592; errup = 0.00538806; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.628726; errup = 0.0044328; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 190 && met<= 195) {eff = 0.496231; errup = 0.0183426; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.629434; errup = 0.0135397; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.737255; errup = 0.00555267; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.656775; errup = 0.00468535; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 195 && met<= 200) {eff = 0.557166; errup = 0.020625; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.691865; errup = 0.0134786; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.775837; errup = 0.00563502; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.6825; errup = 0.00497995; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 200 && met<= 210) {eff = 0.613588; errup = 0.0162625; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.714643; errup = 0.0102329; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.802235; errup = 0.00412988; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.721521; errup = 0.00372205; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 210 && met<= 220) {eff = 0.676248; errup = 0.0186648; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.767995; errup = 0.0108422; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.83659; errup = 0.00425486; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.771127; errup = 0.00401206; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 220 && met<= 230) {eff = 0.753715; errup = 0.0203287; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.815106; errup = 0.0110231; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.849854; errup = 0.00444347; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.798438; errup = 0.0043543; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 230 && met<= 240) {eff = 0.777457; errup = 0.0229277; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.838428; errup = 0.011004; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.880527; errup = 0.00456845; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.840774; errup = 0.00448606; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 240 && met<= 250) {eff = 0.808118; errup = 0.0245104; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.882955; errup = 0.0109464; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.893832; errup = 0.00468695; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.866433; errup = 0.00464701; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 250 && met<= 275) {eff = 0.898851; errup = 0.0146179; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.896425; errup = 0.00703359; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.916677; errup = 0.00307309; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.894316; errup = 0.00314429; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 275 && met<= 300) {eff = 0.909091; errup = 0.0191145; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.929204; errup = 0.00672708; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.936391; errup = 0.00324456; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.91712; errup = 0.0035512; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 300 && met<= 350) {eff = 0.962025; errup = 0.0123134; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.939075; errup = 0.00513974; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.941998; errup = 0.00287257; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.941336; errup = 0.00288951; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 350 && met<= 400) {eff = 0.962617; errup = 0.0178065; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.958088; errup = 0.00514037; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.952878; errup = 0.00377841; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.954771; errup = 0.00366682; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 400 && met<= 450) {eff = 0.954128; errup = 0.0196725; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.976703; errup = 0.004505; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.963603; errup = 0.00465925; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.970104; errup = 0.00420908; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 450 && met<= 500) {eff = 0.971631; errup = 0.0135284; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.989157; errup = 0.00353776; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.954495; errup = 0.00696122; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.964162; errup = 0.00632368; errdown = 0.00329626;}
  else if (ht> 0 && ht<= 600 && met> 500 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00543262;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.964353; errup = 0.00801863; errdown = 0.00718289;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.981785; errup = 0.0056382; errdown = 0.00425746;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.976471; errup = 0.0066564; errdown = 0.00329626;}
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
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.0580311; errup = 0.00853854; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.273902; errup = 0.0169699; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 0 && met<= 50) {eff = 0.512061; errup = 0.0119679; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.239788; errup = 0.0123203; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.432836; errup = 0.0161873; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 0 && met<= 50) {eff = 0.610241; errup = 0.00967965; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.368855; errup = 0.0057901; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.58404; errup = 0.00788592; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 0 && met<= 50) {eff = 0.718103; errup = 0.00425673; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.358079; errup = 0.0165383; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.574176; errup = 0.0270377; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 0 && met<= 50) {eff = 0.721596; errup = 0.01343; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.727555; errup = 0.00887959; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.776404; errup = 0.014199; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 50) {eff = 0.830705; errup = 0.00707675; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.0128205; errup = 0.00857983; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.251185; errup = 0.0228961; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 50 && met<= 75) {eff = 0.515088; errup = 0.0135744; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.07; errup = 0.0180052; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.359756; errup = 0.0284325; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 50 && met<= 75) {eff = 0.58209; errup = 0.0154441; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.313953; errup = 0.0239721; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.564612; errup = 0.0229432; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 50 && met<= 75) {eff = 0.668443; errup = 0.0116455; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.478528; errup = 0.0103138; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.690422; errup = 0.0112097; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 50 && met<= 75) {eff = 0.76283; errup = 0.00509721; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.404682; errup = 0.0303426; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.722973; errup = 0.0383999; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 50 && met<= 75) {eff = 0.766862; errup = 0.0165056; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.73375; errup = 0.0159396; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.814136; errup = 0.0203184; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 50 && met<= 75) {eff = 0.843819; errup = 0.00860726; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.0526316; errup = 0.0271984; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.295082; errup = 0.0375578; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 75 && met<= 100) {eff = 0.577465; errup = 0.0167048; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.196429; errup = 0.0448753; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.467105; errup = 0.0438353; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 75 && met<= 100) {eff = 0.620134; errup = 0.0182735; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.495098; errup = 0.037384; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.59375; errup = 0.03452; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 75 && met<= 100) {eff = 0.738619; errup = 0.0125052; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.584366; errup = 0.0156025; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.737811; errup = 0.0147447; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 75 && met<= 100) {eff = 0.796101; errup = 0.00566351; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.607843; errup = 0.0518085; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.77551; errup = 0.0439117; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 75 && met<= 100) {eff = 0.801205; errup = 0.0182366; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.767516; errup = 0.0244905; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.791878; errup = 0.0297937; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 75 && met<= 100) {eff = 0.863636; errup = 0.00987432; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.0338983; errup = 0.0429671; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.357895; errup = 0.0557905; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 100 && met<= 125) {eff = 0.63606; errup = 0.0202397; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.230769; errup = 0.0733859; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.532468; errup = 0.0624335; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 100 && met<= 125) {eff = 0.701613; errup = 0.0211067; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.645161; errup = 0.0455298; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.700787; errup = 0.0426733; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 100 && met<= 125) {eff = 0.807339; errup = 0.0135565; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.716075; errup = 0.0211522; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.797227; errup = 0.0170527; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 100 && met<= 125) {eff = 0.84023; errup = 0.0061067; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.785714; errup = 0.0664487; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.891892; errup = 0.0510106; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 100 && met<= 125) {eff = 0.827922; errup = 0.0219648; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.841667; errup = 0.0342095; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.893443; errup = 0.0283253; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 100 && met<= 125) {eff = 0.905502; errup = 0.0102026; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.2; errup = 0.130119; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.604167; errup = 0.0776099; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 125 && met<= 150) {eff = 0.801508; errup = 0.0204268; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.222222; errup = 0.109068; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.575; errup = 0.0872414; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 125 && met<= 150) {eff = 0.788406; errup = 0.0225238; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.745098; errup = 0.0646476; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.813333; errup = 0.0466366; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 125 && met<= 150) {eff = 0.85; errup = 0.0150446; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.815686; errup = 0.0248726; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.876254; errup = 0.0193393; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 125 && met<= 150) {eff = 0.895428; errup = 0.00633502; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.106107; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.892857; errup = 0.0577336; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 125 && met<= 150) {eff = 0.885965; errup = 0.0213678; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.982456; errup = 0.0145177; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.970588; errup = 0.0189746; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 125 && met<= 150) {eff = 0.93956; errup = 0.0102474; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.461538; errup = 0.172004; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.612903; errup = 0.0974957; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 150 && met<= 175) {eff = 0.865546; errup = 0.0225156; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.454545; errup = 0.189662; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.666667; errup = 0.106107; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 150 && met<= 175) {eff = 0.883178; errup = 0.0223037; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.75; errup = 0.0943497; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.906977; errup = 0.0439838; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 150 && met<= 175) {eff = 0.914286; errup = 0.0132412; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.873418; errup = 0.026946; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929648; errup = 0.0182155; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 150 && met<= 175) {eff = 0.929748; errup = 0.00601652; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.809524; errup = 0.0888157; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 150 && met<= 175) {eff = 0.907801; errup = 0.0246102; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.942857; errup = 0.0368224; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.975; errup = 0.0206905; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 175) {eff = 0.951351; errup = 0.011195; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.8; errup = 0.106751; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0536391; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 175 && met<= 215) {eff = 0.955645; errup = 0.0130229; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.769231; errup = 0.122762; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.869565; errup = 0.0701228; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 175 && met<= 215) {eff = 0.920833; errup = 0.0175659; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.944444; errup = 0.046004; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.97619; errup = 0.0197048; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 175 && met<= 215) {eff = 0.961039; errup = 0.00899222; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.950413; errup = 0.0194949; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.978022; errup = 0.0104893; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 175 && met<= 215) {eff = 0.969503; errup = 0.00381731; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.916667; errup = 0.0690403; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 175 && met<= 215) {eff = 0.964706; errup = 0.013912; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.966667; errup = 0.0275914; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 215) {eff = 0.980815; errup = 0.00661202; errdown = 0.0139251;}
  else if (ht> 0 && ht<= 400 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0020946;}
  else if (ht> 400 && ht<= 600 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.913043; errup = 0.0559625; errdown = 0.0133199;}
  else if (ht> 600 && ht<= 9999 && el_pt> 20 && el_pt<= 25 && met> 215 && met<= 9999) {eff = 0.981651; errup = 0.00725455; errdown = 0.0101585;}
  else if (ht> 0 && ht<= 400 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0148206;}
  else if (ht> 400 && ht<= 600 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0275299;}
  else if (ht> 600 && ht<= 9999 && el_pt> 25 && el_pt<= 30 && met> 215 && met<= 9999) {eff = 0.984674; errup = 0.00732055; errdown = 0.0156019;}
  else if (ht> 0 && ht<= 400 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0373327;}
  else if (ht> 400 && ht<= 600 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0354114;}
  else if (ht> 600 && ht<= 9999 && el_pt> 30 && el_pt<= 40 && met> 215 && met<= 9999) {eff = 0.997938; errup = 0.00170573; errdown = 0.0128916;}
  else if (ht> 0 && ht<= 400 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.990099; errup = 0.00819202; errdown = 0.0220989;}
  else if (ht> 400 && ht<= 600 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0181352;}
  else if (ht> 600 && ht<= 9999 && el_pt> 40 && el_pt<= 110 && met> 215 && met<= 9999) {eff = 0.992439; errup = 0.00171528; errdown = 0.00629736;}
  else if (ht> 0 && ht<= 400 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.122517;}
  else if (ht> 400 && ht<= 600 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0933526;}
  else if (ht> 600 && ht<= 9999 && el_pt> 110 && el_pt<= 120 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0250959;}
  else if (ht> 0 && ht<= 400 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0704444;}
  else if (ht> 400 && ht<= 600 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0551516;}
  else if (ht> 600 && ht<= 9999 && el_pt> 120 && el_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.998397; errup = 0.00132575; errdown = 0.0139251;}
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
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.271298; errup = 0.0170541; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.550075; errup = 0.0198882; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 50) {eff = 0.843672; errup = 0.0091337; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.431312; errup = 0.012669; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.704082; errup = 0.0120989; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 50) {eff = 0.886586; errup = 0.00491944; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.887264; errup = 0.00786673; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.92204; errup = 0.00652681; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 50) {eff = 0.936299; errup = 0.00313807; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.0614251; errup = 0.0143316; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.414583; errup = 0.0236888; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 50 && met<= 75) {eff = 0.814067; errup = 0.00972747; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.325843; errup = 0.0311932; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.614198; errup = 0.0281699; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 50 && met<= 75) {eff = 0.857381; errup = 0.0102039; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.543956; errup = 0.0190702; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.780458; errup = 0.0146207; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 50 && met<= 75) {eff = 0.893736; errup = 0.0056122; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.921429; errup = 0.0093482; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.952077; errup = 0.00699194; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 50 && met<= 75) {eff = 0.940981; errup = 0.0035208; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.15; errup = 0.0337209; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.45098; errup = 0.0332565; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 75 && met<= 100) {eff = 0.842715; errup = 0.0105919; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.520325; errup = 0.0487155; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.647668; errup = 0.0360589; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 75 && met<= 100) {eff = 0.861354; errup = 0.0115486; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.666667; errup = 0.025072; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.806922; errup = 0.0171572; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 75 && met<= 100) {eff = 0.903462; errup = 0.00599439; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.941414; errup = 0.0106018; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.954397; errup = 0.00843586; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 75 && met<= 100) {eff = 0.95625; errup = 0.00337964; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.14; errup = 0.0672437; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.556452; errup = 0.0479412; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 100 && met<= 125) {eff = 0.868726; errup = 0.0122566; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.619048; errup = 0.0664109; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.816514; errup = 0.0382778; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 100 && met<= 125) {eff = 0.919672; errup = 0.0110898; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.73262; errup = 0.0336118; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.855596; errup = 0.0215151; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 100 && met<= 125) {eff = 0.931906; errup = 0.00593101; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.956364; errup = 0.0122799; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.972569; errup = 0.00808855; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 100 && met<= 125) {eff = 0.968242; errup = 0.00341417; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.423077; errup = 0.116981; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.666667; errup = 0.0654592; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 125 && met<= 150) {eff = 0.912046; errup = 0.0124958; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.724138; errup = 0.089299; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.850746; errup = 0.044684; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 125 && met<= 150) {eff = 0.943262; errup = 0.0112896; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.94186; errup = 0.0248848; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.928571; errup = 0.0191788; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 125 && met<= 150) {eff = 0.956884; errup = 0.00560047; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.980392; errup = 0.0106526; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.973913; errup = 0.0102998; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 125 && met<= 150) {eff = 0.981845; errup = 0.00295218; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.533333; errup = 0.152935; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.882353; errup = 0.0554381; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 175) {eff = 0.929412; errup = 0.01398; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.769231; errup = 0.122762; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.904762; errup = 0.0450174; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 175) {eff = 0.981308; errup = 0.0073897; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 1; errup = 0; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.94898; errup = 0.0218627; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 175) {eff = 0.973348; errup = 0.00524904; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.988636; errup = 0.00940245; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.993711; errup = 0.0052034; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 175) {eff = 0.985469; errup = 0.00306041; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.888889; errup = 0.0598486; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 215) {eff = 0.98374; errup = 0.00643122; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.947368; errup = 0.0435805; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 215) {eff = 0.988473; errup = 0.00550887; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 215) {eff = 0.990955; errup = 0.00295229; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 1; errup = 0; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.993711; errup = 0.0052034; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 215) {eff = 0.995131; errup = 0.00168324; errdown = 0.00378781;}
  else if (ht> 0 && ht<= 400 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00767985;}
  else if (ht> 400 && ht<= 600 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0147661;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 20 && mu_pt<= 25 && met> 215 && met<= 9999) {eff = 0.992788; errup = 0.00392226; errdown = 0.00832374;}
  else if (ht> 0 && ht<= 400 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0298068;}
  else if (ht> 400 && ht<= 600 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0289136;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 25 && mu_pt<= 30 && met> 215 && met<= 9999) {eff = 0.99505; errup = 0.00319693; errdown = 0.0108243;}
  else if (ht> 0 && ht<= 400 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0260004;}
  else if (ht> 400 && ht<= 600 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0183345;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 30 && mu_pt<= 50 && met> 215 && met<= 9999) {eff = 0.994958; errup = 0.00199818; errdown = 0.00633423;}
  else if (ht> 0 && ht<= 400 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.01608;}
  else if (ht> 400 && ht<= 600 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 1; errup = 0; errdown = 0.0108016;}
  else if (ht> 600 && ht<= 9999 && mu_pt> 50 && mu_pt<= 9999 && met> 215 && met<= 9999) {eff = 0.997174; errup = 0.00112047; errdown = 0.00378781;}
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


using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

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
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  //string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  //string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  //string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";
  string ttbar_mc_skim_folder = "mc/skim_higlep1T/";
  string zll_mc_skim_folder = "mc/skim_higlep2T/";
  string qcd_mc_skim_folder = "mc/skim_met150/";
  string met150_mc_skim_folder = "mc/skim_met150/";

  string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  //string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  //string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  //string qcd_data_skim_folder = "data/merged_higdata_higqcd/";
  string ttbar_data_skim_folder = "data/skim_higlep1T/";
  string zll_data_skim_folder = "data/skim_higlep2T/";
  string qcd_data_skim_folder = "data/skim_met150/";
  //string met150_data_skim_folder = "data/skim_met150/";

  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
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
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),"1"));
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
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),"1"));
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
                    "1"));
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
                    "(nbm==0)"
                    ));
    qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("1b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
                    "(nbm==1)"
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
                    "(nbt==2&&nbm==2)"
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("3b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm==3&&nbl==3)"
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("4b Data", Process::Type::data, kGreen,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm>=3&&nbl>=4)"
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
                    "(nbm==0)"
                    ));
    zll_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("1b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                    "(nbm==1)"
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
    if (b.type()/1000==0 || b.SampleType() == 2018) return 1.0;
    float eff = 1.;
    if(b.nvlep()==0){ // search MC sample and qcd MC control sample
      if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
        if (b.SampleType()==2016) eff = get_0l_fakemet_trigeff2016.GetVector(b)[0];
        else if (b.SampleType()==2017) eff = get_0l_fakemet_trigeff2017.GetVector(b)[0];
        //else if (b.SampleType()==2018) eff = get_0l_fakemet_trigeff2018.GetVector(b)[0];
      } else { // TRUE MET
        if (b.SampleType()==2016) eff = get_0l_trigeff2016.GetVector(b)[0];
        else if (b.SampleType()==2017) eff = get_0l_trigeff2017.GetVector(b)[0];
        //else if (b.SampleType()==2018) eff = get_0l_trigeff2018.GetVector(b)[0];
      }
    } else if (b.nel()==1 && b.nmu()==0) { // 1 lepton MC control sample
      if (b.SampleType()==2016) eff = get_1el_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_1el_trigeff2017.GetVector(b)[0];
      //else if (b.SampleType()==2018) eff = get_1el_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==0 && b.nmu()==1) { // 1 lepton MC control sample
      if (b.SampleType()==2016) eff = get_1mu_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_1mu_trigeff2017.GetVector(b)[0];
      //else if (b.SampleType()==2018) eff = get_1mu_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==2 && b.nmu()==0) { // 2 lepton MC control sample
      if (b.SampleType()==2016) eff = get_2el_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_2el_trigeff2017.GetVector(b)[0];
      //else if (b.SampleType()==2018) eff = get_2el_trigeff2018.GetVector(b)[0];
    } else if (b.nel()==0 && b.nmu()==2) { // 2 lepton MC control sample
      if (b.SampleType()==2016) eff = get_2mu_trigeff2016.GetVector(b)[0];
      else if (b.SampleType()==2017) eff = get_2mu_trigeff2017.GetVector(b)[0];
      //else if (b.SampleType()==2018) eff = get_2mu_trigeff2018.GetVector(b)[0];
    }
    return eff;
  });

  NamedFunc weight = "weight"*w_years*hem_weight*trig_weights;

  NamedFunc sr_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && !jetid_low_dphi_met &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (jetid_nb>=2) && (lead_signal_lepton_pt>30) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);

  NamedFunc zll_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met<50&&nlep==2" &&
    (jetid_njet>=4) && (jetid_njet<=5) && (lead_signal_lepton_pt>40) &&
    (jetid_hig_cand_dm<40) && (jetid_hig_cand_am<200) && (jetid_hig_cand_drmax<2.2);
  //dilepton mass cut in skim
  
  NamedFunc qcd_baseline = pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
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
      sr_procs, plt_lin_mc).Weight(weight).Tag("FixName:validate_trig_sr_met_"+year_string).LuminosityTag(total_luminosity_string);
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
