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



const NamedFunc get_0l_trigeff2017_v0("get_0l_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.118202; errup = 0.00647114; errdown = 0.00618617;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.228561; errup = 0.00414684; errdown = 0.00409514;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.362438; errup = 0.0205543; errdown = 0.0200824;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.351759; errup = 0.0370682; errdown = 0.0354728;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.493333; errup = 0.0640194; errdown = 0.0638228;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.140786; errup = 0.00792831; errdown = 0.00758193;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.283364; errup = 0.0047502; errdown = 0.00470307;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.395731; errup = 0.0207976; errdown = 0.0204413;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.42328; errup = 0.038911; errdown = 0.0380424;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.493151; errup = 0.0649718; errdown = 0.0647641;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.166766; errup = 0.00959853; errdown = 0.00918874;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.324521; errup = 0.00523779; errdown = 0.00519473;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.458; errup = 0.0233451; errdown = 0.0231698;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.477157; errup = 0.038145; errdown = 0.0378973;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.38961; errup = 0.0632202; errdown = 0.0600452;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.205904; errup = 0.0115753; errdown = 0.0111268;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.390841; errup = 0.00582333; errdown = 0.00579277;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.5; errup = 0.0234834; errdown = 0.0234834;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.535032; errup = 0.042604; errdown = 0.043084;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.588235; errup = 0.0650517; errdown = 0.0679384;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.241224; errup = 0.0135216; errdown = 0.0130397;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.432834; errup = 0.00624755; errdown = 0.00622653;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.581395; errup = 0.0247251; errdown = 0.0251216;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.620915; errup = 0.0414935; errdown = 0.0431987;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.716981; errup = 0.0659706; errdown = 0.0752878;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.26533; errup = 0.0160304; errdown = 0.0154558;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.480098; errup = 0.00670629; errdown = 0.00669922;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.593301; errup = 0.0249575; errdown = 0.0254254;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.649007; errup = 0.0409427; errdown = 0.0430755;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.550725; errup = 0.0657743; errdown = 0.0674068;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.341679; errup = 0.0191662; errdown = 0.0186814;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.539748; errup = 0.00708134; errdown = 0.0070972;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.654731; errup = 0.0248756; errdown = 0.0257078;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.71875; errup = 0.0416165; errdown = 0.0453516;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.692308; errup = 0.0686825; errdown = 0.0770809;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.369732; errup = 0.0223259; errdown = 0.0218048;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.596124; errup = 0.00745613; errdown = 0.00750012;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.701333; errup = 0.0243702; errdown = 0.0255026;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.691176; errup = 0.0416075; errdown = 0.0446653;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.7; errup = 0.0780739; errdown = 0.0896007;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.383966; errup = 0.0236247; errdown = 0.0231127;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.633548; errup = 0.0078334; errdown = 0.00790347;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.667665; errup = 0.0267017; errdown = 0.027762;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.766129; errup = 0.0395363; errdown = 0.0442593;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.733333; errup = 0.0702056; errdown = 0.082142;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.430267; errup = 0.0286282; errdown = 0.0281924;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.681906; errup = 0.00783883; errdown = 0.00794142;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.775801; errup = 0.0255714; errdown = 0.0276701;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.704918; errup = 0.043372; errdown = 0.0470452;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.777778; errup = 0.0651643; errdown = 0.0795029;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.51153; errup = 0.0238841; errdown = 0.0239346;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.7266; errup = 0.00589238; errdown = 0.0059715;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.818681; errup = 0.0167789; errdown = 0.0180097;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.813084; errup = 0.027352; errdown = 0.0305252;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.818182; errup = 0.0454972; errdown = 0.0549007;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.545763; errup = 0.0304653; errdown = 0.030793;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.797219; errup = 0.00592227; errdown = 0.00605172;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.864198; errup = 0.0157643; errdown = 0.0173603;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.882716; errup = 0.0256997; errdown = 0.0309806;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.861111; errup = 0.041698; errdown = 0.0533274;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.585859; errup = 0.0369679; errdown = 0.0378947;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.841494; errup = 0.0061156; errdown = 0.0063088;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.893665; errup = 0.0148345; errdown = 0.0167507;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.948529; errup = 0.0187887; errdown = 0.0266233;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.892308; errup = 0.0388522; errdown = 0.0531815;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.740157; errup = 0.0406252; errdown = 0.0447687;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.885054; errup = 0.00587368; errdown = 0.00613912;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.929775; errup = 0.0136293; errdown = 0.0162897;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.920863; errup = 0.0230263; errdown = 0.0299998;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.942308; errup = 0.0312343; errdown = 0.0529489;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.7625; errup = 0.0498414; errdown = 0.0572049;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.921169; errup = 0.00569931; errdown = 0.00608683;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.959752; errup = 0.0109001; errdown = 0.0141436;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.972727; errup = 0.0148066; errdown = 0.0258178;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.953488; errup = 0.0299851; errdown = 0.0580755;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.82; errup = 0.0396683; errdown = 0.0468539;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.938905; errup = 0.00390786; errdown = 0.00414787;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.964561; errup = 0.00725305; errdown = 0.00883537;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.969828; errup = 0.0110614; errdown = 0.0158737;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.951456; errup = 0.0208096; errdown = 0.0315047;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.868421; errup = 0.0556306; errdown = 0.0792159;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.965547; errup = 0.00381497; errdown = 0.00424421;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.970954; errup = 0.00760867; errdown = 0.0098252;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.96988; errup = 0.0129497; errdown = 0.0198622;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.986486; errup = 0.0111817; errdown = 0.0303874;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 9999) {eff = 0.878788; errup = 0.0570897; errdown = 0.0855134;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 9999) {eff = 0.968519; errup = 0.00376154; errdown = 0.004222;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 9999) {eff = 0.968439; errup = 0.00710846; errdown = 0.0088446;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 9999) {eff = 0.981308; errup = 0.00892453; errdown = 0.0145322;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 9999) {eff = 0.928571; errup = 0.0259675; errdown = 0.0363535;}
  return eff;
});

const NamedFunc get_0l_fakemet_trigeff2017_v0("get_0l_fakemet_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.123596; errup = 0.0452794; errdown = 0.0355329;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.188815; errup = 0.00513021; errdown = 0.00502708;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.279651; errup = 0.00605781; errdown = 0.00597971;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.31749; errup = 0.00544657; errdown = 0.00539765;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.37227; errup = 0.00307408; errdown = 0.00306385;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.15; errup = 0.060943; errdown = 0.0473246;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.238992; errup = 0.00625925; errdown = 0.00615056;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.319574; errup = 0.00694655; errdown = 0.00686858;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.350866; errup = 0.00620533; errdown = 0.00615597;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.412794; errup = 0.00341814; errdown = 0.00340981;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.266667; errup = 0.082142; errdown = 0.0702056;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.266804; errup = 0.0072869; errdown = 0.00716463;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.361294; errup = 0.0079572; errdown = 0.00788304;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.378864; errup = 0.00686246; errdown = 0.00681498;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.443561; errup = 0.00374757; errdown = 0.00374121;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.104167; errup = 0.0643268; errdown = 0.044248;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.309732; errup = 0.00870177; errdown = 0.00857157;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.404604; errup = 0.00890297; errdown = 0.00884172;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.43023; errup = 0.00759906; errdown = 0.00756681;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.477556; errup = 0.00411404; errdown = 0.00411102;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.371429; errup = 0.0986578; errdown = 0.0901926;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.352278; errup = 0.0102298; errdown = 0.0100988;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.423592; errup = 0.00988964; errdown = 0.00982997;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.463069; errup = 0.00837741; errdown = 0.00835695;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.507831; errup = 0.00445136; errdown = 0.00445259;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.166667; errup = 0.112297; errdown = 0.0779889;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.38976; errup = 0.0114647; errdown = 0.0113469;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.486221; errup = 0.0110393; errdown = 0.0110261;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.49855; errup = 0.009136; errdown = 0.00913505;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.546107; errup = 0.00487945; errdown = 0.00488823;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.277778; errup = 0.144222; errdown = 0.114258;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.425186; errup = 0.0132496; errdown = 0.0131457;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.496301; errup = 0.0122113; errdown = 0.012207;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.540927; errup = 0.00996673; errdown = 0.00999895;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.574945; errup = 0.00522855; errdown = 0.00524519;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.235294; errup = 0.147312; errdown = 0.108959;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.493966; errup = 0.0145821; errdown = 0.0145721;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.561616; errup = 0.0131657; errdown = 0.0132508;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.575424; errup = 0.0103516; errdown = 0.0104168;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.623677; errup = 0.0055526; errdown = 0.0055849;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.272727; errup = 0.196072; errdown = 0.144396;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.507331; errup = 0.0161042; errdown = 0.016119;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.597656; errup = 0.0140129; errdown = 0.0141698;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.579246; errup = 0.0114304; errdown = 0.0115141;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.642391; errup = 0.0059573; errdown = 0.00600106;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.25; errup = 0.239567; errdown = 0.159659;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.539424; errup = 0.0181966; errdown = 0.0182987;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.653775; errup = 0.0147448; errdown = 0.0150369;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.623623; errup = 0.01221; errdown = 0.0123651;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.67077; errup = 0.00624379; errdown = 0.00630382;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.571429; errup = 0.222488; errdown = 0.247841;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.601674; errup = 0.0144865; errdown = 0.0146616;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.664822; errup = 0.0119066; errdown = 0.0121146;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.64353; errup = 0.00906688; errdown = 0.00916893;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.694319; errup = 0.00476164; errdown = 0.00480306;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.4; errup = 0.303366; errdown = 0.253348;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.661078; errup = 0.0167699; errdown = 0.0171696;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.753149; errup = 0.0126899; errdown = 0.013129;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.693149; errup = 0.0102007; errdown = 0.0103888;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.747028; errup = 0.0051753; errdown = 0.00524533;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.2; errup = 0.324251; errdown = 0.166039;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.721501; errup = 0.0174053; errdown = 0.018071;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.786091; errup = 0.0135196; errdown = 0.014146;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.708743; errup = 0.0110955; errdown = 0.0113436;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.778038; errup = 0.00561447; errdown = 0.00571617;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.5; errup = 0.417248; errdown = 0.417248;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.781053; errup = 0.019387; errdown = 0.0206344;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.792402; errup = 0.0151964; errdown = 0.0160253;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.739066; errup = 0.0121406; errdown = 0.0125059;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.806997; errup = 0.00609081; errdown = 0.00623776;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.75; errup = 0.207731; errdown = 0.368402;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.809249; errup = 0.0215917; errdown = 0.0234981;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.834813; errup = 0.0159028; errdown = 0.0171583;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.760676; errup = 0.0129268; errdown = 0.0134066;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.849689; errup = 0.00619514; errdown = 0.00640745;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.444444; errup = 0.213463; errdown = 0.198267;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.847122; errup = 0.0154938; errdown = 0.0168146;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.865031; errup = 0.011045; errdown = 0.0118259;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.772081; errup = 0.0095614; errdown = 0.00984496;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.854076; errup = 0.0045077; errdown = 0.00462405;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0; errup = 0.458642; errdown = 0;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.892857; errup = 0.0178531; errdown = 0.0206337;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.879159; errup = 0.0138072; errdown = 0.015219;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.799681; errup = 0.0114553; errdown = 0.0119503;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.892884; errup = 0.0049459; errdown = 0.00514997;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 1; errup = 0; errdown = 0.841345;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.929515; errup = 0.0170808; errdown = 0.0213439;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.873702; errup = 0.0139923; errdown = 0.015365;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.773749; errup = 0.0114217; errdown = 0.0118313;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.911984; errup = 0.0042738; errdown = 0.00446466;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 1.80842e-308; errup = 9.35748e-307; errdown = 1.69349e-306;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.923077; errup = 0.0279331; errdown = 0.038974;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.853061; errup = 0.0230683; errdown = 0.0262025;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.790164; errup = 0.0167972; errdown = 0.0177949;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.929591; errup = 0.00560601; errdown = 0.00603208;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 9.53547e-322; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.947368; errup = 0.0339218; errdown = 0.0652364;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.834862; errup = 0.0365832; errdown = 0.0434785;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.798535; errup = 0.0248984; errdown = 0.0272452;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.945255; errup = 0.00689888; errdown = 0.00776494;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 1; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.818182; errup = 0.116508; errdown = 0.191402;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.888889; errup = 0.0400561; errdown = 0.0547111;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.780488; errup = 0.0334179; errdown = 0.0371466;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.949275; errup = 0.00936678; errdown = 0.0111507;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 1; errup = -1; errdown = -1;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 1; errup = 0; errdown = 0.601684;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.807692; errup = 0.0804259; errdown = 0.109172;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.8625; errup = 0.0393773; errdown = 0.0498223;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.957377; errup = 0.0115359; errdown = 0.0149532;}
  return eff;
});

const NamedFunc get_1el_trigeff2017_v0("get_1el_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 0 && met<= 110) {eff = 0.0619935; errup = 0.00524257; errdown = 0.00487288;}
  else if (el_pt> 25 && el_pt<= 30 && met> 0 && met<= 110) {eff = 0.118907; errup = 0.00796992; errdown = 0.00754721;}
  else if (el_pt> 30 && el_pt<= 110 && met> 0 && met<= 110) {eff = 0.552149; errup = 0.00415803; errdown = 0.00416527;}
  else if (el_pt> 110 && el_pt<= 120 && met> 0 && met<= 110) {eff = 0.585569; errup = 0.0153569; errdown = 0.01552;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 0 && met<= 110) {eff = 0.82829; errup = 0.0072596; errdown = 0.00750497;}
  else if (el_pt> 20 && el_pt<= 25 && met> 110 && met<= 120) {eff = 0.213793; errup = 0.0395552; errdown = 0.0352269;}
  else if (el_pt> 25 && el_pt<= 30 && met> 110 && met<= 120) {eff = 0.260504; errup = 0.0464748; errdown = 0.0420534;}
  else if (el_pt> 30 && el_pt<= 110 && met> 110 && met<= 120) {eff = 0.676737; errup = 0.0151551; errdown = 0.0155229;}
  else if (el_pt> 110 && el_pt<= 120 && met> 110 && met<= 120) {eff = 0.701493; errup = 0.0595131; errdown = 0.0662634;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 110 && met<= 120) {eff = 0.917647; errup = 0.0212488; errdown = 0.0268265;}
  else if (el_pt> 20 && el_pt<= 25 && met> 120 && met<= 130) {eff = 0.364407; errup = 0.0495787; errdown = 0.0470745;}
  else if (el_pt> 25 && el_pt<= 30 && met> 120 && met<= 130) {eff = 0.350515; errup = 0.0550056; errdown = 0.0516183;}
  else if (el_pt> 30 && el_pt<= 110 && met> 120 && met<= 130) {eff = 0.738876; errup = 0.0153181; errdown = 0.0158989;}
  else if (el_pt> 110 && el_pt<= 120 && met> 120 && met<= 130) {eff = 0.746269; errup = 0.056044; errdown = 0.0643443;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 120 && met<= 130) {eff = 0.913907; errup = 0.0230206; errdown = 0.0292709;}
  else if (el_pt> 20 && el_pt<= 25 && met> 130 && met<= 140) {eff = 0.359551; errup = 0.0578891; errdown = 0.0544101;}
  else if (el_pt> 25 && el_pt<= 30 && met> 130 && met<= 140) {eff = 0.464286; errup = 0.0604404; errdown = 0.0595044;}
  else if (el_pt> 30 && el_pt<= 110 && met> 130 && met<= 140) {eff = 0.734637; errup = 0.0168473; errdown = 0.0175299;}
  else if (el_pt> 110 && el_pt<= 120 && met> 130 && met<= 140) {eff = 0.710145; errup = 0.0579991; errdown = 0.0648324;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 130 && met<= 140) {eff = 0.93985; errup = 0.0205473; errdown = 0.0283269;}
  else if (el_pt> 20 && el_pt<= 25 && met> 140 && met<= 150) {eff = 0.639535; errup = 0.0554437; errdown = 0.0590257;}
  else if (el_pt> 25 && el_pt<= 30 && met> 140 && met<= 150) {eff = 0.6; errup = 0.0661652; errdown = 0.0695965;}
  else if (el_pt> 30 && el_pt<= 110 && met> 140 && met<= 150) {eff = 0.78254; errup = 0.0167488; errdown = 0.0176882;}
  else if (el_pt> 110 && el_pt<= 120 && met> 140 && met<= 150) {eff = 0.822222; errup = 0.0590681; errdown = 0.0759182;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 140 && met<= 150) {eff = 0.919643; errup = 0.0257972; errdown = 0.0345396;}
  else if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 160) {eff = 0.573529; errup = 0.0655803; errdown = 0.0679846;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 160) {eff = 0.563636; errup = 0.0738937; errdown = 0.0764906;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 160) {eff = 0.82459; errup = 0.0156505; errdown = 0.0167703;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 160) {eff = 0.893617; errup = 0.0451725; errdown = 0.0655621;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 160) {eff = 0.963303; errup = 0.0174813; errdown = 0.028065;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 170) {eff = 0.607143; errup = 0.0713147; errdown = 0.0756132;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 170) {eff = 0.538462; errup = 0.0770645; errdown = 0.078728;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 170) {eff = 0.873646; errup = 0.0142978; errdown = 0.0157315;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 170) {eff = 0.827586; errup = 0.0723692; errdown = 0.0998191;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 170) {eff = 0.906667; errup = 0.0337754; errdown = 0.04665;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 180) {eff = 0.696429; errup = 0.0657537; errdown = 0.0736929;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 180) {eff = 0.734694; errup = 0.0670433; errdown = 0.0780196;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 180) {eff = 0.859922; errup = 0.0155312; errdown = 0.0170189;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 180) {eff = 0.911765; errup = 0.0476324; errdown = 0.0784416;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 180) {eff = 0.974026; errup = 0.0167592; errdown = 0.0332329;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 190) {eff = 0.816327; errup = 0.0574355; errdown = 0.0725146;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 190) {eff = 0.837209; errup = 0.0580173; errdown = 0.0766334;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 190) {eff = 0.919725; errup = 0.0131195; errdown = 0.0152117;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 190) {eff = 0.913043; errup = 0.0559625; errdown = 0.103371;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 190) {eff = 0.985294; errup = 0.0121686; errdown = 0.0330032;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 200) {eff = 0.767442; errup = 0.0679853; errdown = 0.0824369;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 200) {eff = 0.791667; errup = 0.0613414; errdown = 0.0754508;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 200) {eff = 0.936471; errup = 0.0118929; errdown = 0.0141481;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 200) {eff = 1; errup = 0; errdown = 0.0263287;}
  else if (el_pt> 20 && el_pt<= 25 && met> 200 && met<= 210) {eff = 0.878049; errup = 0.0516458; errdown = 0.0740833;}
  else if (el_pt> 25 && el_pt<= 30 && met> 200 && met<= 210) {eff = 0.810811; errup = 0.0670146; errdown = 0.0869581;}
  else if (el_pt> 30 && el_pt<= 110 && met> 200 && met<= 210) {eff = 0.960526; errup = 0.00996633; errdown = 0.0127051;}
  else if (el_pt> 110 && el_pt<= 120 && met> 200 && met<= 210) {eff = 0.931034; errup = 0.0444185; errdown = 0.0838105;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 200 && met<= 210) {eff = 0.958904; errup = 0.0222832; errdown = 0.0383633;}
  else if (el_pt> 20 && el_pt<= 25 && met> 210 && met<= 9999) {eff = 0.935714; errup = 0.0207171; errdown = 0.0279785;}
  else if (el_pt> 25 && el_pt<= 30 && met> 210 && met<= 9999) {eff = 0.938053; errup = 0.0225646; errdown = 0.0317723;}
  else if (el_pt> 30 && el_pt<= 110 && met> 210 && met<= 9999) {eff = 0.965116; errup = 0.00529352; errdown = 0.00612651;}
  else if (el_pt> 110 && el_pt<= 120 && met> 210 && met<= 9999) {eff = 0.976744; errup = 0.015007; errdown = 0.0298505;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 210 && met<= 9999) {eff = 1; errup = 0; errdown = 0.00719369;}
  return eff;
});

const NamedFunc get_1mu_trigeff2017_v0("get_1mu_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 0 && met<= 110) {eff = 0.140487; errup = 0.0069928; errdown = 0.00672041;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 0 && met<= 110) {eff = 0.509969; errup = 0.0117043; errdown = 0.011715;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 0 && met<= 110) {eff = 0.664176; errup = 0.00668724; errdown = 0.00675275;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 0 && met<= 110) {eff = 0.915724; errup = 0.00327313; errdown = 0.00339022;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 110 && met<= 120) {eff = 0.397661; errup = 0.0408696; errdown = 0.039584;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 110 && met<= 120) {eff = 0.709924; errup = 0.0415572; errdown = 0.0450534;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 110 && met<= 120) {eff = 0.833773; errup = 0.0194869; errdown = 0.0213668;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 110 && met<= 120) {eff = 0.966499; errup = 0.00735202; errdown = 0.009091;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 120 && met<= 130) {eff = 0.447853; errup = 0.042223; errdown = 0.0415354;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 120 && met<= 130) {eff = 0.723881; errup = 0.0403806; errdown = 0.0440273;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 120 && met<= 130) {eff = 0.86921; errup = 0.0178717; errdown = 0.0200359;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 120 && met<= 130) {eff = 0.958084; errup = 0.00895556; errdown = 0.0109858;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 130 && met<= 140) {eff = 0.575; errup = 0.0483942; errdown = 0.0497523;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 130 && met<= 140) {eff = 0.762887; errup = 0.0451032; errdown = 0.0511298;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 130 && met<= 140) {eff = 0.852778; errup = 0.0189996; errdown = 0.0211017;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 130 && met<= 140) {eff = 0.969072; errup = 0.00782771; errdown = 0.0100131;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 140 && met<= 150) {eff = 0.554545; errup = 0.051138; errdown = 0.0522185;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 140 && met<= 150) {eff = 0.768116; errup = 0.0532749; errdown = 0.0620693;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 140 && met<= 150) {eff = 0.88755; errup = 0.0203114; errdown = 0.023728;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 140 && met<= 150) {eff = 0.973214; errup = 0.0075712; errdown = 0.00998961;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 160) {eff = 0.666667; errup = 0.0443388; errdown = 0.0472195;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 160) {eff = 0.866667; errup = 0.0366017; errdown = 0.0459433;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 160) {eff = 0.945607; errup = 0.0146741; errdown = 0.0189235;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 160) {eff = 0.99; errup = 0.00477986; errdown = 0.00783614;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 170) {eff = 0.771429; errup = 0.0525767; errdown = 0.0613516;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 170) {eff = 0.851852; errup = 0.0404926; errdown = 0.0504559;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 170) {eff = 0.92766; errup = 0.0170011; errdown = 0.0210882;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 170) {eff = 0.972637; errup = 0.00806857; errdown = 0.0107752;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 180) {eff = 0.802469; errup = 0.0459391; errdown = 0.0543801;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 180) {eff = 0.888889; errup = 0.0354029; errdown = 0.0466071;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 180) {eff = 0.963303; errup = 0.0126003; errdown = 0.0176003;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 180) {eff = 0.993671; errup = 0.00408698; errdown = 0.0082865;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 190) {eff = 0.969697; errup = 0.0195489; errdown = 0.0385739;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 190) {eff = 0.940887; errup = 0.0165669; errdown = 0.0215411;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 190) {eff = 0.983165; errup = 0.00725289; errdown = 0.0112281;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 200) {eff = 0.857143; errup = 0.0478753; errdown = 0.0628782;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 200) {eff = 0.94; errup = 0.0324767; errdown = 0.054936;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 200) {eff = 0.975309; errup = 0.0117802; errdown = 0.0190923;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 200) {eff = 0.992647; errup = 0.00474792; errdown = 0.00961549;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 200 && met<= 210) {eff = 0.865385; errup = 0.0482807; errdown = 0.0649611;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 200 && met<= 210) {eff = 1; errup = 0; errdown = 0.0297298;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 200 && met<= 210) {eff = 0.992857; errup = 0.00590966; errdown = 0.0162324;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 200 && met<= 210) {eff = 0.985816; errup = 0.00677638; errdown = 0.0110731;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 210 && met<= 9999) {eff = 0.938547; errup = 0.0179637; errdown = 0.0236054;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 210 && met<= 9999) {eff = 0.97191; errup = 0.0120805; errdown = 0.0185555;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 210 && met<= 9999) {eff = 0.983264; errup = 0.00577123; errdown = 0.00815023;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 210 && met<= 9999) {eff = 0.994152; errup = 0.0025237; errdown = 0.00393672;}
  return eff;
});

const NamedFunc get_2el_trigeff2017_v0("get_2el_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig());
  errup+=errdown; //suppress unused warning
  if (el_pt> 40 && el_pt<= 45) {eff = 0.846154; errup = 0.0353795; errdown = 0.0353795;}
  else if (el_pt> 45 && el_pt<= 50) {eff = 0.875; errup = 0.0268248; errdown = 0.0268248;}
  else if (el_pt> 50 && el_pt<= 55) {eff = 0.888889; errup = 0.0240328; errdown = 0.0240328;}
  else if (el_pt> 55 && el_pt<= 60) {eff = 0.911765; errup = 0.0217539; errdown = 0.0217539;}
  else if (el_pt> 60 && el_pt<= 65) {eff = 0.872222; errup = 0.0248831; errdown = 0.0248831;}
  else if (el_pt> 65 && el_pt<= 70) {eff = 0.862857; errup = 0.0260038; errdown = 0.0260038;}
  else if (el_pt> 70 && el_pt<= 75) {eff = 0.897436; errup = 0.0242905; errdown = 0.0242905;}
  else if (el_pt> 75 && el_pt<= 80) {eff = 0.926667; errup = 0.0212847; errdown = 0.0212847;}
  else if (el_pt> 80 && el_pt<= 85) {eff = 0.881481; errup = 0.0278184; errdown = 0.0278184;}
  else if (el_pt> 85 && el_pt<= 90) {eff = 0.900763; errup = 0.026122; errdown = 0.026122;}
  else if (el_pt> 90 && el_pt<= 95) {eff = 0.925532; errup = 0.027078; errdown = 0.027078;}
  else if (el_pt> 95 && el_pt<= 100) {eff = 0.897727; errup = 0.0323006; errdown = 0.0323006;}
  else if (el_pt> 100 && el_pt<= 105) {eff = 0.93617; errup = 0.025213; errdown = 0.025213;}
  else if (el_pt> 105 && el_pt<= 110) {eff = 0.880435; errup = 0.0338265; errdown = 0.0338265;}
  else if (el_pt> 110 && el_pt<= 9999) {eff = 0.975464; errup = 0.00289636; errdown = 0.00289636;}
  return eff;
});

const NamedFunc get_2mu_trigeff2017_v0("get_2mu_trigeff2017_v0", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.932039; errup = 0.0175353; errdown = 0.0175353;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.952381; errup = 0.011271; errdown = 0.011271;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.977828; errup = 0.000485423; errdown = 0.000485423;}
  return eff;
});
}
