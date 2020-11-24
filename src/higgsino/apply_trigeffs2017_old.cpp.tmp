#include <vector>
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "higgsino/hig_utilities.hpp"
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

const NamedFunc get_1el_trigeff2017("get_1el_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_1mu_trigeff2017("get_1mu_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_2el_trigeff2017("get_2el_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
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

const NamedFunc get_2mu_trigeff2017("get_2mu_trigeff2017", [](const Baby &b) -> NamedFunc::ScalarType{
  float errup=0., errdown=0.; // Not used, but for reference
  float eff = 1., mu_pt = HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig());
  errup+=errdown; //suppress unused warning
  if (mu_pt> 40 && mu_pt<= 45) {eff = 0.932039; errup = 0.0175353; errdown = 0.0175353;}
  else if (mu_pt> 45 && mu_pt<= 50) {eff = 0.952381; errup = 0.011271; errdown = 0.011271;}
  else if (mu_pt> 50 && mu_pt<= 9999) {eff = 0.977828; errup = 0.000485423; errdown = 0.000485423;}
  return eff;
});

}
