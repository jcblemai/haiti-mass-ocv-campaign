 double eff_v_2d(double t_since_vacc, int scenario) {
      return 0;
  double eff_v_2d = 0.0;
  switch(scenario){
    case 1:
      if      (t_since_vacc <=   1./12) eff_v_2d =  0.76              ;
      else if (t_since_vacc <=   2./12) eff_v_2d =  0.753527484533759;
      else if (t_since_vacc <=   3./12) eff_v_2d =  0.746961516042262;
      else if (t_since_vacc <=   4./12) eff_v_2d =  0.740300745209625;
      else if (t_since_vacc <=   5./12) eff_v_2d =  0.733543803237948;
      else if (t_since_vacc <=   6./12) eff_v_2d =  0.726689301566023;
      else if (t_since_vacc <=   7./12) eff_v_2d =  0.719735831583988;
      else if (t_since_vacc <=   8./12) eff_v_2d =  0.712681964343848;
      else if (t_since_vacc <=   9./12) eff_v_2d =  0.705526250265831;
      else if (t_since_vacc <=  10./12) eff_v_2d =  0.698267218840495;
      else if (t_since_vacc <=  11./12) eff_v_2d =  0.690903378326535;
      else if (t_since_vacc <=  12./12) eff_v_2d =  0.683433215444227;
      else if (t_since_vacc <=  13./12) eff_v_2d =  0.675855195064453;
      else if (t_since_vacc <=  14./12) eff_v_2d =  0.668167759893222;
      else if (t_since_vacc <=  15./12) eff_v_2d =  0.660369330151649;
      else if (t_since_vacc <=  16./12) eff_v_2d =  0.652458303251305;
      else if (t_since_vacc <=  17./12) eff_v_2d =  0.644433053464886;
      else if (t_since_vacc <=  18./12) eff_v_2d =  0.636291931592122;
      else if (t_since_vacc <=  19./12) eff_v_2d =  0.628033264620864;
      else if (t_since_vacc <=  20./12) eff_v_2d =  0.619655355383277;
      else if (t_since_vacc <=  21./12) eff_v_2d =  0.61115648220707 ;
      else if (t_since_vacc <=  22./12) eff_v_2d =  0.602534898561692;
      else if (t_since_vacc <=  23./12) eff_v_2d =  0.593788832699414;
      else if (t_since_vacc <=  24./12) eff_v_2d =  0.584916487291234;
      else if (t_since_vacc <=  25./12) eff_v_2d =  0.575916039057525;
      else if (t_since_vacc <=  26./12) eff_v_2d =  0.566785638393345;
      else if (t_since_vacc <=  27./12) eff_v_2d =  0.557523408988343;
      else if (t_since_vacc <=  28./12) eff_v_2d =  0.548127447441173;
      else if (t_since_vacc <=  29./12) eff_v_2d =  0.538595822868346;
      else if (t_since_vacc <=  30./12) eff_v_2d =  0.528926576507423;
      else if (t_since_vacc <=  31./12) eff_v_2d =  0.519117721314497;
      else if (t_since_vacc <=  32./12) eff_v_2d =  0.509167241555842;
      else if (t_since_vacc <=  33./12) eff_v_2d =  0.499073092393685;
      else if (t_since_vacc <=  34./12) eff_v_2d =  0.488833199465984;
      else if (t_since_vacc <=  35./12) eff_v_2d =  0.478445458460148;
      else if (t_since_vacc <=  36./12) eff_v_2d =  0.467907734680592;
      else if (t_since_vacc <=  37./12) eff_v_2d =  0.457217862610059;
      else if (t_since_vacc <=  38./12) eff_v_2d =  0.446373645464601;
      else if (t_since_vacc <=  39./12) eff_v_2d =  0.435372854742138;
      else if (t_since_vacc <=  40./12) eff_v_2d =  0.424213229764494;
      else if (t_since_vacc <=  41./12) eff_v_2d =  0.412892477212831;
      else if (t_since_vacc <=  42./12) eff_v_2d =  0.401408270656362;
      else if (t_since_vacc <=  43./12) eff_v_2d =  0.38975825007427  ;
      else if (t_since_vacc <=  44./12) eff_v_2d =  0.377940021370718;
      else if (t_since_vacc <=  45./12) eff_v_2d =  0.365951155882864;
      else if (t_since_vacc <=  46./12) eff_v_2d =  0.353789189881759;
      else if (t_since_vacc <=  47./12) eff_v_2d =  0.341451624066056;
      else if (t_since_vacc <=  48./12) eff_v_2d =  0.328935923048392;
      else if (t_since_vacc <=  49./12) eff_v_2d =  0.316239514834368;
      else if (t_since_vacc <=  50./12) eff_v_2d =  0.303359790293999;
      else if (t_since_vacc <=  51./12) eff_v_2d =  0.290294102625533;
      else if (t_since_vacc <=  52./12) eff_v_2d =  0.27703976681153  ;
      else if (t_since_vacc <=  53./12) eff_v_2d =  0.263594059067087;
      else if (t_since_vacc <=  54./12) eff_v_2d =  0.249954216280098;
      else if (t_since_vacc <=  55./12) eff_v_2d =  0.236117435443426;
      else if (t_since_vacc <=  56./12) eff_v_2d =  0.222080873078887;
      else if (t_since_vacc <=  57./12) eff_v_2d =  0.207841644652907;
      else if (t_since_vacc <=  58./12) eff_v_2d =  0.193396823983748;
      else if (t_since_vacc <=  59./12) eff_v_2d =  0.178743442640173;
      else if (t_since_vacc <=  60./12) eff_v_2d = 0.163878489331427;
      else if (t_since_vacc <=  61./12) eff_v_2d =  0.148798909288418;
    break;
    case 2:
      eff_v_2d = 0.76;
    break;
    case 3:
      if      (t_since_vacc <=   1./12) eff_v_2d = 0.62899141753524; 
      else if (t_since_vacc <=   2./12) eff_v_2d = 0.62363463243243;
      else if (t_since_vacc <=   3./12) eff_v_2d = 0.61820050371012;
      else if (t_since_vacc <=   4./12) eff_v_2d = 0.61268791464710;
      else if (t_since_vacc <=   5./12) eff_v_2d = 0.60709573239845;
      else if (t_since_vacc <=   6./12) eff_v_2d = 0.60142280776277;
      else if (t_since_vacc <=   7./12) eff_v_2d = 0.59566797494594;
      else if (t_since_vacc <=   8./12) eff_v_2d = 0.58983005132162;
      else if (t_since_vacc <=   9./12) eff_v_2d = 0.58390783718819;
      else if (t_since_vacc <=  10./12) eff_v_2d = 0.57790011552220;
      else if (t_since_vacc <=  11./12) eff_v_2d = 0.57180565172828;
      else if (t_since_vacc <=  12./12) eff_v_2d = 0.56562319338543;
      else if (t_since_vacc <=  13./12) eff_v_2d = 0.55935146998966;
      else if (t_since_vacc <=  14./12) eff_v_2d = 0.55298919269287;
      else if (t_since_vacc <=  15./12) eff_v_2d = 0.54653505403800;
      else if (t_since_vacc <=  16./12) eff_v_2d = 0.53998772769036;
      else if (t_since_vacc <=  17./12) eff_v_2d = 0.53334586816505;
      else if (t_since_vacc <=  18./12) eff_v_2d = 0.52660811055048;
      else if (t_since_vacc <=  19./12) eff_v_2d = 0.51977307022784;
      else if (t_since_vacc <=  20./12) eff_v_2d = 0.51283934258662;
      else if (t_since_vacc <=  21./12) eff_v_2d = 0.50580550273589;
      else if (t_since_vacc <=  22./12) eff_v_2d = 0.49867010521154;
      else if (t_since_vacc <=  23./12) eff_v_2d = 0.49143168367921;
      else if (t_since_vacc <=  24./12) eff_v_2d = 0.48408875063295;
      else if (t_since_vacc <=  25./12) eff_v_2d = 0.47663979708957;
      else if (t_since_vacc <=  26./12) eff_v_2d = 0.46908329227848;
      else if (t_since_vacc <=  27./12) eff_v_2d = 0.46141768332718;
      else if (t_since_vacc <=  28./12) eff_v_2d = 0.45364139494210;
      else if (t_since_vacc <=  29./12) eff_v_2d = 0.44575282908489;
      else if (t_since_vacc <=  30./12) eff_v_2d = 0.43775036464403;
      else if (t_since_vacc <=  31./12) eff_v_2d = 0.42963235710167;
      else if (t_since_vacc <=  32./12) eff_v_2d = 0.42139713819568;
      else if (t_since_vacc <=  33./12) eff_v_2d = 0.41304301557684;
      else if (t_since_vacc <=  34./12) eff_v_2d = 0.40456827246104;
      else if (t_since_vacc <=  35./12) eff_v_2d = 0.39597116727650;
      else if (t_since_vacc <=  36./12) eff_v_2d = 0.38724993330585;
      else if (t_since_vacc <=  37./12) eff_v_2d = 0.37840277832307;
      else if (t_since_vacc <=  38./12) eff_v_2d = 0.36942788422520;
      else if (t_since_vacc <=  39./12) eff_v_2d = 0.36032340665871;
      else if (t_since_vacc <=  40./12) eff_v_2d = 0.35108747464049;
      else if (t_since_vacc <=  41./12) eff_v_2d = 0.34171819017333;
      else if (t_since_vacc <=  42./12) eff_v_2d = 0.33221362785594;
      else if (t_since_vacc <=  43./12) eff_v_2d = 0.32257183448719;
      else if (t_since_vacc <=  44./12) eff_v_2d = 0.31279082866482;
      else if (t_since_vacc <=  45./12) eff_v_2d = 0.30286860037818;
      else if (t_since_vacc <=  46./12) eff_v_2d = 0.29280311059522;
      else if (t_since_vacc <=  47./12) eff_v_2d = 0.28259229084344;
      else if (t_since_vacc <=  48./12) eff_v_2d = 0.27223404278483;
      else if (t_since_vacc <=  49./12) eff_v_2d = 0.26172623778464;
      else if (t_since_vacc <=  50./12) eff_v_2d = 0.25106671647396;
      else if (t_since_vacc <=  51./12) eff_v_2d = 0.24025328830599;
      else if (t_since_vacc <=  52./12) eff_v_2d = 0.22928373110581;
      else if (t_since_vacc <=  53./12) eff_v_2d = 0.21815579061378;
      else if (t_since_vacc <=  54./12) eff_v_2d = 0.20686718002227;
      else if (t_since_vacc <=  55./12) eff_v_2d = 0.19541557950571;
      else if (t_since_vacc <=  56./12) eff_v_2d = 0.18379863574388;
      else if (t_since_vacc <=  57./12) eff_v_2d = 0.17201396143827;
      else if (t_since_vacc <=  58./12) eff_v_2d = 0.16005913482151;
      else if (t_since_vacc <=  59./12) eff_v_2d = 0.14793169915969;
      else if (t_since_vacc <=  60./12) eff_v_2d = 0.13562916224751;      
      else if (t_since_vacc <=  61./12) eff_v_2d = 0.12314899589607;
    break;
  }
  return eff_v_2d;
}