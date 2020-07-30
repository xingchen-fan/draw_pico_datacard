#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/apply_trigeffs2017.hpp"

namespace Higfuncs{

const NamedFunc get_0l_trigeff2017("get_0l_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_0l_fakemet_trigeff2017("get_0l_fakemet_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), ht = b.ht();
  errup+=errdown; //suppress unused warning
  if (ht> 0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.810127; errup = 0.0457539; errdown = 0.054656;}
  else if (ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.642169; errup = 0.00917933; errdown = 0.00928275;}
  else if (ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.554865; errup = 0.0164403; errdown = 0.0165575;}
  else if (ht> 800 && ht<= 1000 && met> 150 && met<= 155) {eff = 0.571858; errup = 0.00432187; errdown = 0.00433277;}
  else if (ht> 1000 && ht<= 9999 && met> 150 && met<= 155) {eff = 0.402514; errup = 0.00273356; errdown = 0.00272755;}
  else if (ht> 0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.684932; errup = 0.0579627; errdown = 0.0636162;}
  else if (ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.654285; errup = 0.0092453; errdown = 0.00936105;}
  else if (ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.55787; errup = 0.0173951; errdown = 0.0175336;}
  else if (ht> 800 && ht<= 1000 && met> 155 && met<= 160) {eff = 0.606505; errup = 0.00452896; errdown = 0.00454717;}
  else if (ht> 1000 && ht<= 9999 && met> 155 && met<= 160) {eff = 0.435122; errup = 0.00300178; errdown = 0.00299706;}
  else if (ht> 0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.671642; errup = 0.0614633; errdown = 0.0671964;}
  else if (ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.645223; errup = 0.00952407; errdown = 0.00963822;}
  else if (ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.583123; errup = 0.0180071; errdown = 0.0182239;}
  else if (ht> 800 && ht<= 1000 && met> 160 && met<= 165) {eff = 0.633412; errup = 0.00480826; errdown = 0.0048347;}
  else if (ht> 1000 && ht<= 9999 && met> 160 && met<= 165) {eff = 0.469185; errup = 0.00328289; errdown = 0.00328025;}
  else if (ht> 0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.868852; errup = 0.0440748; errdown = 0.0582986;}
  else if (ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.667719; errup = 0.00948519; errdown = 0.00962032;}
  else if (ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.6603; errup = 0.0179362; errdown = 0.0183903;}
  else if (ht> 800 && ht<= 1000 && met> 165 && met<= 170) {eff = 0.662715; errup = 0.00501358; errdown = 0.00505004;}
  else if (ht> 1000 && ht<= 9999 && met> 165 && met<= 170) {eff = 0.506302; errup = 0.00356749; errdown = 0.00356812;}
  else if (ht> 0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.782609; errup = 0.0518935; errdown = 0.06119;}
  else if (ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.675754; errup = 0.00964388; errdown = 0.00979218;}
  else if (ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.69375; errup = 0.0186721; errdown = 0.0193029;}
  else if (ht> 800 && ht<= 1000 && met> 170 && met<= 175) {eff = 0.68951; errup = 0.00518397; errdown = 0.00523141;}
  else if (ht> 1000 && ht<= 9999 && met> 170 && met<= 175) {eff = 0.54457; errup = 0.00386342; errdown = 0.00386875;}
  else if (ht> 0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.853933; errup = 0.0383743; errdown = 0.0474531;}
  else if (ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.67605; errup = 0.00972527; errdown = 0.00987642;}
  else if (ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.683849; errup = 0.0197869; errdown = 0.0204461;}
  else if (ht> 800 && ht<= 1000 && met> 175 && met<= 180) {eff = 0.728483; errup = 0.00536027; errdown = 0.00542658;}
  else if (ht> 1000 && ht<= 9999 && met> 175 && met<= 180) {eff = 0.568572; errup = 0.00418586; errdown = 0.0041956;}
  else if (ht> 0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.702128; errup = 0.0498352; errdown = 0.0545876;}
  else if (ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.653027; errup = 0.00986813; errdown = 0.00999863;}
  else if (ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.736243; errup = 0.0196582; errdown = 0.0205979;}
  else if (ht> 800 && ht<= 1000 && met> 180 && met<= 185) {eff = 0.744226; errup = 0.00547086; errdown = 0.00554765;}
  else if (ht> 1000 && ht<= 9999 && met> 180 && met<= 185) {eff = 0.600231; errup = 0.00447722; errdown = 0.00449387;}
  else if (ht> 0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.6875; errup = 0.0550977; errdown = 0.0603075;}
  else if (ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.677798; errup = 0.0100687; errdown = 0.0102328;}
  else if (ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.718504; errup = 0.0204684; errdown = 0.02137;}
  else if (ht> 800 && ht<= 1000 && met> 185 && met<= 190) {eff = 0.758144; errup = 0.00563256; errdown = 0.00572195;}
  else if (ht> 1000 && ht<= 9999 && met> 185 && met<= 190) {eff = 0.641484; errup = 0.00471161; errdown = 0.0047388;}
  else if (ht> 0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.709302; errup = 0.0517623; errdown = 0.0571663;}
  else if (ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.670208; errup = 0.0102524; errdown = 0.0104132;}
  else if (ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.753086; errup = 0.020024; errdown = 0.0211189;}
  else if (ht> 800 && ht<= 1000 && met> 190 && met<= 195) {eff = 0.775367; errup = 0.00587796; errdown = 0.00598742;}
  else if (ht> 1000 && ht<= 9999 && met> 190 && met<= 195) {eff = 0.660188; errup = 0.00501455; errdown = 0.00505032;}
  else if (ht> 0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.723404; errup = 0.0485602; errdown = 0.0538247;}
  else if (ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.666186; errup = 0.0104988; errdown = 0.0106623;}
  else if (ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.768473; errup = 0.0214445; errdown = 0.0228427;}
  else if (ht> 800 && ht<= 1000 && met> 195 && met<= 200) {eff = 0.80348; errup = 0.00582893; errdown = 0.00596008;}
  else if (ht> 1000 && ht<= 9999 && met> 195 && met<= 200) {eff = 0.693627; errup = 0.00527469; errdown = 0.00532525;}
  else if (ht> 0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.724719; errup = 0.0348293; errdown = 0.0375563;}
  else if (ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.681522; errup = 0.00760242; errdown = 0.00769866;}
  else if (ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.807829; errup = 0.0137771; errdown = 0.0145388;}
  else if (ht> 800 && ht<= 1000 && met> 200 && met<= 210) {eff = 0.820866; errup = 0.00431055; errdown = 0.00439199;}
  else if (ht> 1000 && ht<= 9999 && met> 200 && met<= 210) {eff = 0.718784; errup = 0.00398585; errdown = 0.0040202;}
  else if (ht> 0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.707447; errup = 0.0345607; errdown = 0.0369364;}
  else if (ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.65846; errup = 0.00811979; errdown = 0.00821212;}
  else if (ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.782369; errup = 0.0155887; errdown = 0.0164009;}
  else if (ht> 800 && ht<= 1000 && met> 210 && met<= 220) {eff = 0.845687; errup = 0.00447447; errdown = 0.00458126;}
  else if (ht> 1000 && ht<= 9999 && met> 210 && met<= 220) {eff = 0.766566; errup = 0.00433146; errdown = 0.00438741;}
  else if (ht> 0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.68617; errup = 0.0353444; errdown = 0.0374734;}
  else if (ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.661848; errup = 0.00849065; errdown = 0.00859426;}
  else if (ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.824757; errup = 0.014394; errdown = 0.0153409;}
  else if (ht> 800 && ht<= 1000 && met> 220 && met<= 230) {eff = 0.855788; errup = 0.00472357; errdown = 0.0048533;}
  else if (ht> 1000 && ht<= 9999 && met> 220 && met<= 230) {eff = 0.7932; errup = 0.00470347; errdown = 0.00478276;}
  else if (ht> 0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.675393; errup = 0.0354172; errdown = 0.0373891;}
  else if (ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.648879; errup = 0.00912516; errdown = 0.00923315;}
  else if (ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.85381; errup = 0.014125; errdown = 0.0152858;}
  else if (ht> 800 && ht<= 1000 && met> 230 && met<= 240) {eff = 0.876031; errup = 0.00481756; errdown = 0.00498015;}
  else if (ht> 1000 && ht<= 9999 && met> 230 && met<= 240) {eff = 0.823779; errup = 0.00499071; errdown = 0.00510241;}
  else if (ht> 0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.60804; errup = 0.0364418; errdown = 0.0376028;}
  else if (ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.648638; errup = 0.00969323; errdown = 0.0098148;}
  else if (ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.866667; errup = 0.0137177; errdown = 0.0149488;}
  else if (ht> 800 && ht<= 1000 && met> 240 && met<= 250) {eff = 0.893153; errup = 0.0049167; errdown = 0.00511894;}
  else if (ht> 1000 && ht<= 9999 && met> 240 && met<= 250) {eff = 0.85932; errup = 0.00505085; errdown = 0.00520394;}
  else if (ht> 0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.639405; errup = 0.0213408; errdown = 0.0218816;}
  else if (ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.620462; errup = 0.00683414; errdown = 0.00688158;}
  else if (ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.885174; errup = 0.00866595; errdown = 0.00924923;}
  else if (ht> 800 && ht<= 1000 && met> 250 && met<= 275) {eff = 0.90347; errup = 0.00338453; errdown = 0.00349182;}
  else if (ht> 1000 && ht<= 9999 && met> 250 && met<= 275) {eff = 0.87082; errup = 0.00364741; errdown = 0.00373583;}
  else if (ht> 0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.65896; errup = 0.0214318; errdown = 0.0220719;}
  else if (ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.593524; errup = 0.00814207; errdown = 0.00819295;}
  else if (ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.894431; errup = 0.00893204; errdown = 0.00961823;}
  else if (ht> 800 && ht<= 1000 && met> 275 && met<= 300) {eff = 0.921142; errup = 0.00370607; errdown = 0.00386822;}
  else if (ht> 1000 && ht<= 9999 && met> 275 && met<= 300) {eff = 0.902105; errup = 0.00400252; errdown = 0.00415049;}
  else if (ht> 0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.739726; errup = 0.0165732; errdown = 0.0172571;}
  else if (ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.570788; errup = 0.00715092; errdown = 0.00718017;}
  else if (ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.906707; errup = 0.00657081; errdown = 0.00699672;}
  else if (ht> 800 && ht<= 1000 && met> 300 && met<= 350) {eff = 0.920547; errup = 0.00342041; errdown = 0.00355717;}
  else if (ht> 1000 && ht<= 9999 && met> 300 && met<= 350) {eff = 0.915116; errup = 0.00360477; errdown = 0.00374586;}
  else if (ht> 0 && ht<= 200 && met> 350 && met<= 400) {eff = 0.824042; errup = 0.0161611; errdown = 0.0173508;}
  else if (ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.540258; errup = 0.0101823; errdown = 0.0102153;}
  else if (ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.922977; errup = 0.00682144; errdown = 0.00739492;}
  else if (ht> 800 && ht<= 1000 && met> 350 && met<= 400) {eff = 0.931567; errup = 0.00450118; errdown = 0.0047829;}
  else if (ht> 1000 && ht<= 9999 && met> 350 && met<= 400) {eff = 0.929544; errup = 0.00474075; errdown = 0.00504365;}
  else if (ht> 0 && ht<= 200 && met> 400 && met<= 450) {eff = 0.883382; errup = 0.0175727; errdown = 0.0199919;}
  else if (ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.54553; errup = 0.0147641; errdown = 0.0148423;}
  else if (ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.941818; errup = 0.00708915; errdown = 0.0079438;}
  else if (ht> 800 && ht<= 1000 && met> 400 && met<= 450) {eff = 0.946429; errup = 0.005608; errdown = 0.00618791;}
  else if (ht> 1000 && ht<= 9999 && met> 400 && met<= 450) {eff = 0.943038; errup = 0.00585416; errdown = 0.00644522;}
  else if (ht> 0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.861842; errup = 0.0285803; errdown = 0.0338786;}
  else if (ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.57594; errup = 0.0197822; errdown = 0.0200194;}
  else if (ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.959524; errup = 0.00681006; errdown = 0.00799924;}
  else if (ht> 800 && ht<= 1000 && met> 450 && met<= 500) {eff = 0.923326; errup = 0.00880117; errdown = 0.00977137;}
  else if (ht> 1000 && ht<= 9999 && met> 450 && met<= 500) {eff = 0.94075; errup = 0.00824823; errdown = 0.0093919;}
  else if (ht> 0 && ht<= 200 && met> 500 && met<= 9999) {eff = 0.906977; errup = 0.0315422; errdown = 0.0426732;}
  else if (ht> 200 && ht<= 600 && met> 500 && met<= 9999) {eff = 0.626263; errup = 0.0292556; errdown = 0.0301549;}
  else if (ht> 600 && ht<= 800 && met> 500 && met<= 9999) {eff = 0.946524; errup = 0.00953202; errdown = 0.0112716;}
  else if (ht> 800 && ht<= 1000 && met> 500 && met<= 9999) {eff = 0.95539; errup = 0.00891233; errdown = 0.0107801;}
  else if (ht> 1000 && ht<= 9999 && met> 500 && met<= 9999) {eff = 0.956427; errup = 0.00953231; errdown = 0.0117457;}
  return eff;
});

const NamedFunc get_1el_trigeff2017("get_1el_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), el_pt = b.el_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (el_pt> 20 && el_pt<= 25 && met> 150 && met<= 110) {eff = 0.0619935; errup = 0.00524257; errdown = 0.00487288;}
  else if (el_pt> 25 && el_pt<= 30 && met> 150 && met<= 110) {eff = 0.118907; errup = 0.00796992; errdown = 0.00754721;}
  else if (el_pt> 30 && el_pt<= 110 && met> 150 && met<= 110) {eff = 0.552149; errup = 0.00415803; errdown = 0.00416527;}
  else if (el_pt> 110 && el_pt<= 120 && met> 150 && met<= 110) {eff = 0.585569; errup = 0.0153569; errdown = 0.01552;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 150 && met<= 110) {eff = 0.82829; errup = 0.0072596; errdown = 0.00750497;}
  else if (el_pt> 20 && el_pt<= 25 && met> 155 && met<= 120) {eff = 0.213793; errup = 0.0395552; errdown = 0.0352269;}
  else if (el_pt> 25 && el_pt<= 30 && met> 155 && met<= 120) {eff = 0.260504; errup = 0.0464748; errdown = 0.0420534;}
  else if (el_pt> 30 && el_pt<= 110 && met> 155 && met<= 120) {eff = 0.676737; errup = 0.0151551; errdown = 0.0155229;}
  else if (el_pt> 110 && el_pt<= 120 && met> 155 && met<= 120) {eff = 0.701493; errup = 0.0595131; errdown = 0.0662634;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 155 && met<= 120) {eff = 0.917647; errup = 0.0212488; errdown = 0.0268265;}
  else if (el_pt> 20 && el_pt<= 25 && met> 160 && met<= 130) {eff = 0.364407; errup = 0.0495787; errdown = 0.0470745;}
  else if (el_pt> 25 && el_pt<= 30 && met> 160 && met<= 130) {eff = 0.350515; errup = 0.0550056; errdown = 0.0516183;}
  else if (el_pt> 30 && el_pt<= 110 && met> 160 && met<= 130) {eff = 0.738876; errup = 0.0153181; errdown = 0.0158989;}
  else if (el_pt> 110 && el_pt<= 120 && met> 160 && met<= 130) {eff = 0.746269; errup = 0.056044; errdown = 0.0643443;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 160 && met<= 130) {eff = 0.913907; errup = 0.0230206; errdown = 0.0292709;}
  else if (el_pt> 20 && el_pt<= 25 && met> 165 && met<= 140) {eff = 0.359551; errup = 0.0578891; errdown = 0.0544101;}
  else if (el_pt> 25 && el_pt<= 30 && met> 165 && met<= 140) {eff = 0.464286; errup = 0.0604404; errdown = 0.0595044;}
  else if (el_pt> 30 && el_pt<= 110 && met> 165 && met<= 140) {eff = 0.734637; errup = 0.0168473; errdown = 0.0175299;}
  else if (el_pt> 110 && el_pt<= 120 && met> 165 && met<= 140) {eff = 0.710145; errup = 0.0579991; errdown = 0.0648324;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 165 && met<= 140) {eff = 0.93985; errup = 0.0205473; errdown = 0.0283269;}
  else if (el_pt> 20 && el_pt<= 25 && met> 170 && met<= 150) {eff = 0.639535; errup = 0.0554437; errdown = 0.0590257;}
  else if (el_pt> 25 && el_pt<= 30 && met> 170 && met<= 150) {eff = 0.6; errup = 0.0661652; errdown = 0.0695965;}
  else if (el_pt> 30 && el_pt<= 110 && met> 170 && met<= 150) {eff = 0.78254; errup = 0.0167488; errdown = 0.0176882;}
  else if (el_pt> 110 && el_pt<= 120 && met> 170 && met<= 150) {eff = 0.822222; errup = 0.0590681; errdown = 0.0759182;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 170 && met<= 150) {eff = 0.919643; errup = 0.0257972; errdown = 0.0345396;}
  else if (el_pt> 20 && el_pt<= 25 && met> 175 && met<= 160) {eff = 0.573529; errup = 0.0655803; errdown = 0.0679846;}
  else if (el_pt> 25 && el_pt<= 30 && met> 175 && met<= 160) {eff = 0.563636; errup = 0.0738937; errdown = 0.0764906;}
  else if (el_pt> 30 && el_pt<= 110 && met> 175 && met<= 160) {eff = 0.82459; errup = 0.0156505; errdown = 0.0167703;}
  else if (el_pt> 110 && el_pt<= 120 && met> 175 && met<= 160) {eff = 0.893617; errup = 0.0451725; errdown = 0.0655621;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 175 && met<= 160) {eff = 0.963303; errup = 0.0174813; errdown = 0.028065;}
  else if (el_pt> 20 && el_pt<= 25 && met> 180 && met<= 170) {eff = 0.607143; errup = 0.0713147; errdown = 0.0756132;}
  else if (el_pt> 25 && el_pt<= 30 && met> 180 && met<= 170) {eff = 0.538462; errup = 0.0770645; errdown = 0.078728;}
  else if (el_pt> 30 && el_pt<= 110 && met> 180 && met<= 170) {eff = 0.873646; errup = 0.0142978; errdown = 0.0157315;}
  else if (el_pt> 110 && el_pt<= 120 && met> 180 && met<= 170) {eff = 0.827586; errup = 0.0723692; errdown = 0.0998191;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 180 && met<= 170) {eff = 0.906667; errup = 0.0337754; errdown = 0.04665;}
  else if (el_pt> 20 && el_pt<= 25 && met> 185 && met<= 180) {eff = 0.696429; errup = 0.0657537; errdown = 0.0736929;}
  else if (el_pt> 25 && el_pt<= 30 && met> 185 && met<= 180) {eff = 0.734694; errup = 0.0670433; errdown = 0.0780196;}
  else if (el_pt> 30 && el_pt<= 110 && met> 185 && met<= 180) {eff = 0.859922; errup = 0.0155312; errdown = 0.0170189;}
  else if (el_pt> 110 && el_pt<= 120 && met> 185 && met<= 180) {eff = 0.911765; errup = 0.0476324; errdown = 0.0784416;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 185 && met<= 180) {eff = 0.974026; errup = 0.0167592; errdown = 0.0332329;}
  else if (el_pt> 20 && el_pt<= 25 && met> 190 && met<= 190) {eff = 0.816327; errup = 0.0574355; errdown = 0.0725146;}
  else if (el_pt> 25 && el_pt<= 30 && met> 190 && met<= 190) {eff = 0.837209; errup = 0.0580173; errdown = 0.0766334;}
  else if (el_pt> 30 && el_pt<= 110 && met> 190 && met<= 190) {eff = 0.919725; errup = 0.0131195; errdown = 0.0152117;}
  else if (el_pt> 110 && el_pt<= 120 && met> 190 && met<= 190) {eff = 0.913043; errup = 0.0559625; errdown = 0.103371;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 190 && met<= 190) {eff = 0.985294; errup = 0.0121686; errdown = 0.0330032;}
  else if (el_pt> 20 && el_pt<= 25 && met> 195 && met<= 200) {eff = 0.767442; errup = 0.0679853; errdown = 0.0824369;}
  else if (el_pt> 25 && el_pt<= 30 && met> 195 && met<= 200) {eff = 0.791667; errup = 0.0613414; errdown = 0.0754508;}
  else if (el_pt> 30 && el_pt<= 110 && met> 195 && met<= 200) {eff = 0.936471; errup = 0.0118929; errdown = 0.0141481;}
  else if (el_pt> 110 && el_pt<= 120 && met> 195 && met<= 200) {eff = 1; errup = 0; errdown = 0.0576587;}
  else if (el_pt> 120 && el_pt<= 9999 && met> 195 && met<= 200) {eff = 1; errup = 0; errdown = 0.0263287;}
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

const NamedFunc get_1mu_trigeff2017("get_1mu_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., met = b.met(), mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 20 && mu_pt<= 25 && met> 150 && met<= 110) {eff = 0.140487; errup = 0.0069928; errdown = 0.00672041;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 150 && met<= 110) {eff = 0.509969; errup = 0.0117043; errdown = 0.011715;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 150 && met<= 110) {eff = 0.664176; errup = 0.00668724; errdown = 0.00675275;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 150 && met<= 110) {eff = 0.915724; errup = 0.00327313; errdown = 0.00339022;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 155 && met<= 120) {eff = 0.397661; errup = 0.0408696; errdown = 0.039584;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 155 && met<= 120) {eff = 0.709924; errup = 0.0415572; errdown = 0.0450534;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 155 && met<= 120) {eff = 0.833773; errup = 0.0194869; errdown = 0.0213668;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 155 && met<= 120) {eff = 0.966499; errup = 0.00735202; errdown = 0.009091;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 160 && met<= 130) {eff = 0.447853; errup = 0.042223; errdown = 0.0415354;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 160 && met<= 130) {eff = 0.723881; errup = 0.0403806; errdown = 0.0440273;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 160 && met<= 130) {eff = 0.86921; errup = 0.0178717; errdown = 0.0200359;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 160 && met<= 130) {eff = 0.958084; errup = 0.00895556; errdown = 0.0109858;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 165 && met<= 140) {eff = 0.575; errup = 0.0483942; errdown = 0.0497523;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 165 && met<= 140) {eff = 0.762887; errup = 0.0451032; errdown = 0.0511298;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 165 && met<= 140) {eff = 0.852778; errup = 0.0189996; errdown = 0.0211017;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 165 && met<= 140) {eff = 0.969072; errup = 0.00782771; errdown = 0.0100131;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 170 && met<= 150) {eff = 0.554545; errup = 0.051138; errdown = 0.0522185;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 170 && met<= 150) {eff = 0.768116; errup = 0.0532749; errdown = 0.0620693;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 170 && met<= 150) {eff = 0.88755; errup = 0.0203114; errdown = 0.023728;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 170 && met<= 150) {eff = 0.973214; errup = 0.0075712; errdown = 0.00998961;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 175 && met<= 160) {eff = 0.666667; errup = 0.0443388; errdown = 0.0472195;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 175 && met<= 160) {eff = 0.866667; errup = 0.0366017; errdown = 0.0459433;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 175 && met<= 160) {eff = 0.945607; errup = 0.0146741; errdown = 0.0189235;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 175 && met<= 160) {eff = 0.99; errup = 0.00477986; errdown = 0.00783614;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 180 && met<= 170) {eff = 0.771429; errup = 0.0525767; errdown = 0.0613516;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 180 && met<= 170) {eff = 0.851852; errup = 0.0404926; errdown = 0.0504559;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 180 && met<= 170) {eff = 0.92766; errup = 0.0170011; errdown = 0.0210882;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 180 && met<= 170) {eff = 0.972637; errup = 0.00806857; errdown = 0.0107752;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 185 && met<= 180) {eff = 0.802469; errup = 0.0459391; errdown = 0.0543801;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 185 && met<= 180) {eff = 0.888889; errup = 0.0354029; errdown = 0.0466071;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 185 && met<= 180) {eff = 0.963303; errup = 0.0126003; errdown = 0.0176003;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 185 && met<= 180) {eff = 0.993671; errup = 0.00408698; errdown = 0.0082865;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 190 && met<= 190) {eff = 0.714286; errup = 0.0643308; errdown = 0.0730111;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 190 && met<= 190) {eff = 0.969697; errup = 0.0195489; errdown = 0.0385739;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 190 && met<= 190) {eff = 0.940887; errup = 0.0165669; errdown = 0.0215411;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 190 && met<= 190) {eff = 0.983165; errup = 0.00725289; errdown = 0.0112281;}
  else if (mu_pt> 20 && mu_pt<= 25 && met> 195 && met<= 200) {eff = 0.857143; errup = 0.0478753; errdown = 0.0628782;}
  else if (mu_pt> 25 && mu_pt<= 30 && met> 195 && met<= 200) {eff = 0.94; errup = 0.0324767; errdown = 0.054936;}
  else if (mu_pt> 30 && mu_pt<= 50 && met> 195 && met<= 200) {eff = 0.975309; errup = 0.0117802; errdown = 0.0190923;}
  else if (mu_pt> 50 && mu_pt<= 9999 && met> 195 && met<= 200) {eff = 0.992647; errup = 0.00474792; errdown = 0.00961549;}
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

const NamedFunc get_2el_trigeff2017("get_2el_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., el_pt = b.el_pt()->at(0); //assumes first lepton is signal
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

const NamedFunc get_2mu_trigeff2017("get_2mu_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = b.mu_pt()->at(0); //assumes first lepton is signal
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.932039; errup = 0.0175353; errdown = 0.0175353;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.952381; errup = 0.011271; errdown = 0.011271;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.977828; errup = 0.000485423; errdown = 0.000485423;}
  return eff;
});

}
