#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/apply_trigeffs2017.hpp"

namespace Higfuncs{

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

}
