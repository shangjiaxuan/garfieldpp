#include <cmath>
#include <iostream>

#include <TMath.h>

#include "Garfield/Random.hh"

namespace {

double denlan(const double v) {
  const double p1[5] = {0.4259894875, -0.1249762550, 0.03984243700,
                        -0.006298287635, 0.001511162253};
  const double q1[5] = {1.0, -0.3388260629, 0.09594393323, -0.01608042283,
                        0.003778942063};

  const double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518,
                        -0.001394989411, 0.0001283617211};
  const double q2[5] = {1.0, 0.7428795082, 0.3153932961, 0.06694219548,
                        0.008790609714};

  const double p3[5] = {0.1788544503, 0.09359161662, 0.006325387654,
                        0.00006611667319, -0.000002031049101};
  const double q3[5] = {1.0, 0.6097809921, 0.2560616665, 0.04746722384,
                        0.006957301675};

  const double p4[5] = {0.9874054407, 118.6723273, 849.2794360, -743.7792444,
                        427.0262186};
  const double q4[5] = {1.0, 106.8615961, 337.6496214, 2016.712389,
                        1597.063511};

  const double p5[5] = {1.003675074, 167.5702434, 4789.711289, 21217.86767,
                        -22324.94910};
  const double q5[5] = {1.0, 156.9424537, 3745.310488, 9834.698876,
                        66924.28357};

  const double p6[5] = {1.000827619, 664.9143136, 62972.92665, 475554.6998,
                        -5743609.109};
  const double q6[5] = {1.0, 651.4101098, 56974.73333, 165917.4725,
                        -2815759.939};

  const double a1[3] = {0.04166666667, -0.01996527778, 0.02709538966};

  const double a2[2] = {-1.845568670, -4.284640743};

  if (v < -5.5) {
    const double u = std::exp(v + 1.0);
    if (u < 1e-10) return 0.0;
    const double ue = std::exp(-1 / u);
    const double us = std::sqrt(u);
    return 0.3989422803 * (ue / us) *
           (1 + (a1[0] + (a1[1] + a1[2] * u) * u) * u);
  } else if (v < -1) {
    const double u = std::exp(-v - 1);
    return std::exp(-u) * std::sqrt(u) *
           (p1[0] + (p1[1] + (p1[2] + (p1[3] + p1[4] * v) * v) * v) * v) /
           (q1[0] + (q1[1] + (q1[2] + (q1[3] + q1[4] * v) * v) * v) * v);
  } else if (v < 1) {
    return (p2[0] + (p2[1] + (p2[2] + (p2[3] + p2[4] * v) * v) * v) * v) /
           (q2[0] + (q2[1] + (q2[2] + (q2[3] + q2[4] * v) * v) * v) * v);
  } else if (v < 5) {
    return (p3[0] + (p3[1] + (p3[2] + (p3[3] + p3[4] * v) * v) * v) * v) /
           (q3[0] + (q3[1] + (q3[2] + (q3[3] + q3[4] * v) * v) * v) * v);
  } else if (v < 12) {
    const double u = 1. / v;
    return u * u *
           (p4[0] + (p4[1] + (p4[2] + (p4[3] + p4[4] * u) * u) * u) * u) /
           (q4[0] + (q4[1] + (q4[2] + (q4[3] + q4[4] * u) * u) * u) * u);
  } else if (v < 50) {
    const double u = 1. / v;
    return u * u *
           (p5[0] + (p5[1] + (p5[2] + (p5[3] + p5[4] * u) * u) * u) * u) /
           (q5[0] + (q5[1] + (q5[2] + (q5[3] + q5[4] * u) * u) * u) * u);
  } else if (v < 300) {
    const double u = 1. / v;
    return u * u *
           (p6[0] + (p6[1] + (p6[2] + (p6[3] + p6[4] * u) * u) * u) * u) /
           (q6[0] + (q6[1] + (q6[2] + (q6[3] + q6[4] * u) * u) * u) * u);
  } else {
    const double u = 1. / (v - v * std::log(v) / (v + 1));
    return u * u * (1 + (a2[0] + a2[1] * u) * u);
  }
}
}  // namespace
namespace Garfield {

double RndmLandau() {
  constexpr double f[] = {
      0,         0,         0,         0,         0,         -2.244733,
      -2.204365, -2.168163, -2.135219, -2.104898, -2.076740, -2.050397,
      -2.025605, -2.002150, -1.979866, -1.958612, -1.938275, -1.918760,
      -1.899984, -1.881879, -1.864385, -1.847451, -1.831030, -1.815083,
      -1.799574, -1.784473, -1.769751, -1.755383, -1.741346, -1.727620,
      -1.714187, -1.701029, -1.688130, -1.675477, -1.663057, -1.650858,
      -1.638868, -1.627078, -1.615477, -1.604058, -1.592811, -1.581729,
      -1.570806, -1.560034, -1.549407, -1.538919, -1.528565, -1.518339,
      -1.508237, -1.498254, -1.488386, -1.478628, -1.468976, -1.459428,
      -1.449979, -1.440626, -1.431365, -1.422195, -1.413111, -1.404112,
      -1.395194, -1.386356, -1.377594, -1.368906, -1.360291, -1.351746,
      -1.343269, -1.334859, -1.326512, -1.318229, -1.310006, -1.301843,
      -1.293737, -1.285688, -1.277693, -1.269752, -1.261863, -1.254024,
      -1.246235, -1.238494, -1.230800, -1.223153, -1.215550, -1.207990,
      -1.200474, -1.192999, -1.185566, -1.178172, -1.170817, -1.163500,
      -1.156220, -1.148977, -1.141770, -1.134598, -1.127459, -1.120354,
      -1.113282, -1.106242, -1.099233, -1.092255, -1.085306, -1.078388,
      -1.071498, -1.064636, -1.057802, -1.050996, -1.044215, -1.037461,
      -1.030733, -1.024029, -1.017350, -1.010695, -1.004064, -.997456,
      -.990871,  -.984308,  -.977767,  -.971247,  -.964749,  -.958271,
      -.951813,  -.945375,  -.938957,  -.932558,  -.926178,  -.919816,
      -.913472,  -.907146,  -.900838,  -.894547,  -.888272,  -.882014,
      -.875773,  -.869547,  -.863337,  -.857142,  -.850963,  -.844798,
      -.838648,  -.832512,  -.826390,  -.820282,  -.814187,  -.808106,
      -.802038,  -.795982,  -.789940,  -.783909,  -.777891,  -.771884,
      -.765889,  -.759906,  -.753934,  -.747973,  -.742023,  -.736084,
      -.730155,  -.724237,  -.718328,  -.712429,  -.706541,  -.700661,
      -.694791,  -.688931,  -.683079,  -.677236,  -.671402,  -.665576,
      -.659759,  -.653950,  -.648149,  -.642356,  -.636570,  -.630793,
      -.625022,  -.619259,  -.613503,  -.607754,  -.602012,  -.596276,
      -.590548,  -.584825,  -.579109,  -.573399,  -.567695,  -.561997,
      -.556305,  -.550618,  -.544937,  -.539262,  -.533592,  -.527926,
      -.522266,  -.516611,  -.510961,  -.505315,  -.499674,  -.494037,
      -.488405,  -.482777,  -.477153,  -.471533,  -.465917,  -.460305,
      -.454697,  -.449092,  -.443491,  -.437893,  -.432299,  -.426707,
      -.421119,  -.415534,  -.409951,  -.404372,  -.398795,  -.393221,
      -.387649,  -.382080,  -.376513,  -.370949,  -.365387,  -.359826,
      -.354268,  -.348712,  -.343157,  -.337604,  -.332053,  -.326503,
      -.320955,  -.315408,  -.309863,  -.304318,  -.298775,  -.293233,
      -.287692,  -.282152,  -.276613,  -.271074,  -.265536,  -.259999,
      -.254462,  -.248926,  -.243389,  -.237854,  -.232318,  -.226783,
      -.221247,  -.215712,  -.210176,  -.204641,  -.199105,  -.193568,
      -.188032,  -.182495,  -.176957,  -.171419,  -.165880,  -.160341,
      -.154800,  -.149259,  -.143717,  -.138173,  -.132629,  -.127083,
      -.121537,  -.115989,  -.110439,  -.104889,  -.099336,  -.093782,
      -.088227,  -.082670,  -.077111,  -.071550,  -.065987,  -.060423,
      -.054856,  -.049288,  -.043717,  -.038144,  -.032569,  -.026991,
      -.021411,  -.015828,  -.010243,  -.004656,  .000934,   .006527,
      .012123,   .017722,   .023323,   .028928,   .034535,   .040146,
      .045759,   .051376,   .056997,   .062620,   .068247,   .073877,
      .079511,   .085149,   .090790,   .096435,   .102083,   .107736,
      .113392,   .119052,   .124716,   .130385,   .136057,   .141734,
      .147414,   .153100,   .158789,   .164483,   .170181,   .175884,
      .181592,   .187304,   .193021,   .198743,   .204469,   .210201,
      .215937,   .221678,   .227425,   .233177,   .238933,   .244696,
      .250463,   .256236,   .262014,   .267798,   .273587,   .279382,
      .285183,   .290989,   .296801,   .302619,   .308443,   .314273,
      .320109,   .325951,   .331799,   .337654,   .343515,   .349382,
      .355255,   .361135,   .367022,   .372915,   .378815,   .384721,
      .390634,   .396554,   .402481,   .408415,   .414356,   .420304,
      .426260,   .432222,   .438192,   .444169,   .450153,   .456145,
      .462144,   .468151,   .474166,   .480188,   .486218,   .492256,
      .498302,   .504356,   .510418,   .516488,   .522566,   .528653,
      .534747,   .540850,   .546962,   .553082,   .559210,   .565347,
      .571493,   .577648,   .583811,   .589983,   .596164,   .602355,
      .608554,   .614762,   .620980,   .627207,   .633444,   .639689,
      .645945,   .652210,   .658484,   .664768,   .671062,   .677366,
      .683680,   .690004,   .696338,   .702682,   .709036,   .715400,
      .721775,   .728160,   .734556,   .740963,   .747379,   .753807,
      .760246,   .766695,   .773155,   .779627,   .786109,   .792603,
      .799107,   .805624,   .812151,   .818690,   .825241,   .831803,
      .838377,   .844962,   .851560,   .858170,   .864791,   .871425,
      .878071,   .884729,   .891399,   .898082,   .904778,   .911486,
      .918206,   .924940,   .931686,   .938446,   .945218,   .952003,
      .958802,   .965614,   .972439,   .979278,   .986130,   .992996,
      .999875,   1.006769,  1.013676,  1.020597,  1.027533,  1.034482,
      1.041446,  1.048424,  1.055417,  1.062424,  1.069446,  1.076482,
      1.083534,  1.090600,  1.097681,  1.104778,  1.111889,  1.119016,
      1.126159,  1.133316,  1.140490,  1.147679,  1.154884,  1.162105,
      1.169342,  1.176595,  1.183864,  1.191149,  1.198451,  1.205770,
      1.213105,  1.220457,  1.227826,  1.235211,  1.242614,  1.250034,
      1.257471,  1.264926,  1.272398,  1.279888,  1.287395,  1.294921,
      1.302464,  1.310026,  1.317605,  1.325203,  1.332819,  1.340454,
      1.348108,  1.355780,  1.363472,  1.371182,  1.378912,  1.386660,
      1.394429,  1.402216,  1.410024,  1.417851,  1.425698,  1.433565,
      1.441453,  1.449360,  1.457288,  1.465237,  1.473206,  1.481196,
      1.489208,  1.497240,  1.505293,  1.513368,  1.521465,  1.529583,
      1.537723,  1.545885,  1.554068,  1.562275,  1.570503,  1.578754,
      1.587028,  1.595325,  1.603644,  1.611987,  1.620353,  1.628743,
      1.637156,  1.645593,  1.654053,  1.662538,  1.671047,  1.679581,
      1.688139,  1.696721,  1.705329,  1.713961,  1.722619,  1.731303,
      1.740011,  1.748746,  1.757506,  1.766293,  1.775106,  1.783945,
      1.792810,  1.801703,  1.810623,  1.819569,  1.828543,  1.837545,
      1.846574,  1.855631,  1.864717,  1.873830,  1.882972,  1.892143,
      1.901343,  1.910572,  1.919830,  1.929117,  1.938434,  1.947781,
      1.957158,  1.966566,  1.976004,  1.985473,  1.994972,  2.004503,
      2.014065,  2.023659,  2.033285,  2.042943,  2.052633,  2.062355,
      2.072110,  2.081899,  2.091720,  2.101575,  2.111464,  2.121386,
      2.131343,  2.141334,  2.151360,  2.161421,  2.171517,  2.181648,
      2.191815,  2.202018,  2.212257,  2.222533,  2.232845,  2.243195,
      2.253582,  2.264006,  2.274468,  2.284968,  2.295507,  2.306084,
      2.316701,  2.327356,  2.338051,  2.348786,  2.359562,  2.370377,
      2.381234,  2.392131,  2.403070,  2.414051,  2.425073,  2.436138,
      2.447246,  2.458397,  2.469591,  2.480828,  2.492110,  2.503436,
      2.514807,  2.526222,  2.537684,  2.549190,  2.560743,  2.572343,
      2.583989,  2.595682,  2.607423,  2.619212,  2.631050,  2.642936,
      2.654871,  2.666855,  2.678890,  2.690975,  2.703110,  2.715297,
      2.727535,  2.739825,  2.752168,  2.764563,  2.777012,  2.789514,
      2.802070,  2.814681,  2.827347,  2.840069,  2.852846,  2.865680,
      2.878570,  2.891518,  2.904524,  2.917588,  2.930712,  2.943894,
      2.957136,  2.970439,  2.983802,  2.997227,  3.010714,  3.024263,
      3.037875,  3.051551,  3.065290,  3.079095,  3.092965,  3.106900,
      3.120902,  3.134971,  3.149107,  3.163312,  3.177585,  3.191928,
      3.206340,  3.220824,  3.235378,  3.250005,  3.264704,  3.279477,
      3.294323,  3.309244,  3.324240,  3.339312,  3.354461,  3.369687,
      3.384992,  3.400375,  3.415838,  3.431381,  3.447005,  3.462711,
      3.478500,  3.494372,  3.510328,  3.526370,  3.542497,  3.558711,
      3.575012,  3.591402,  3.607881,  3.624450,  3.641111,  3.657863,
      3.674708,  3.691646,  3.708680,  3.725809,  3.743034,  3.760357,
      3.777779,  3.795300,  3.812921,  3.830645,  3.848470,  3.866400,
      3.884434,  3.902574,  3.920821,  3.939176,  3.957640,  3.976215,
      3.994901,  4.013699,  4.032612,  4.051639,  4.070783,  4.090045,
      4.109425,  4.128925,  4.148547,  4.168292,  4.188160,  4.208154,
      4.228275,  4.248524,  4.268903,  4.289413,  4.310056,  4.330832,
      4.351745,  4.372794,  4.393982,  4.415310,  4.436781,  4.458395,
      4.480154,  4.502060,  4.524114,  4.546319,  4.568676,  4.591187,
      4.613854,  4.636678,  4.659662,  4.682807,  4.706116,  4.729590,
      4.753231,  4.777041,  4.801024,  4.825179,  4.849511,  4.874020,
      4.898710,  4.923582,  4.948639,  4.973883,  4.999316,  5.024942,
      5.050761,  5.076778,  5.102993,  5.129411,  5.156034,  5.182864,
      5.209903,  5.237156,  5.264625,  5.292312,  5.320220,  5.348354,
      5.376714,  5.405306,  5.434131,  5.463193,  5.492496,  5.522042,
      5.551836,  5.581880,  5.612178,  5.642734,  5.673552,  5.704634,
      5.735986,  5.767610,  5.799512,  5.831694,  5.864161,  5.896918,
      5.929968,  5.963316,  5.996967,  6.030925,  6.065194,  6.099780,
      6.134687,  6.169921,  6.205486,  6.241387,  6.277630,  6.314220,
      6.351163,  6.388465,  6.426130,  6.464166,  6.502578,  6.541371,
      6.580553,  6.620130,  6.660109,  6.700495,  6.741297,  6.782520,
      6.824173,  6.866262,  6.908795,  6.951780,  6.995225,  7.039137,
      7.083525,  7.128398,  7.173764,  7.219632,  7.266011,  7.312910,
      7.360339,  7.408308,  7.456827,  7.505905,  7.555554,  7.605785,
      7.656608,  7.708035,  7.760077,  7.812747,  7.866057,  7.920019,
      7.974647,  8.029953,  8.085952,  8.142657,  8.200083,  8.258245,
      8.317158,  8.376837,  8.437300,  8.498562,  8.560641,  8.623554,
      8.687319,  8.751955,  8.817481,  8.883916,  8.951282,  9.019600,
      9.088889,  9.159174,  9.230477,  9.302822,  9.376233,  9.450735,
      9.526355,  9.603118,  9.681054,  9.760191,  9.840558,  9.922186,
      10.005107, 10.089353, 10.174959, 10.261958, 10.350389, 10.440287,
      10.531693, 10.624646, 10.719188, 10.815362, 10.913214, 11.012789,
      11.114137, 11.217307, 11.322352, 11.429325, 11.538283, 11.649285,
      11.762390, 11.877664, 11.995170, 12.114979, 12.237161, 12.361791,
      12.488946, 12.618708, 12.751161, 12.886394, 13.024498, 13.165570,
      13.309711, 13.457026, 13.607625, 13.761625, 13.919145, 14.080314,
      14.245263, 14.414134, 14.587072, 14.764233, 14.945778, 15.131877,
      15.322712, 15.518470, 15.719353, 15.925570, 16.137345, 16.354912,
      16.578520, 16.808433, 17.044929, 17.288305, 17.538873, 17.796967,
      18.062943, 18.337176, 18.620068, 18.912049, 19.213574, 19.525133,
      19.847249, 20.180480, 20.525429, 20.882738, 21.253102, 21.637266,
      22.036036, 22.450278, 22.880933, 23.329017, 23.795634, 24.281981,
      24.789364, 25.319207, 25.873062, 26.452634, 27.059789, 27.696581,
      28.365274, 29.068370, 29.808638, 30.589157, 31.413354, 32.285060,
      33.208568, 34.188705, 35.230920, 36.341388, 37.527131, 38.796172,
      40.157721, 41.622399, 43.202525, 44.912465, 46.769077, 48.792279,
      51.005773, 53.437996, 56.123356, 59.103894};

  const double x = RndmUniformPos();
  double u = 1000 * x;
  const int i = static_cast<int>(u);
  u = u - i;
  if (i >= 70 && i <= 800) {
    return f[i - 1] + u * (f[i] - f[i - 1]);
  } else if (i >= 7 && i <= 980) {
    return f[i - 1] +
           u * (f[i] - f[i - 1] -
                0.25 * (1 - u) * (f[i + 1] - f[i] - f[i - 1] + f[i - 2]));
  } else if (i < 7) {
    const double v = log(x);
    u = 1. / v;
    return ((0.99858950 + (3.45213058e1 + 1.70854528e1 * u) * u) /
            (1 + (3.41760202e1 + 4.01244582 * u) * u)) *
           (-log(-0.91893853 - v) - 1);
  } else {
    u = 1. - x;
    const double v = u * u;
    if (x <= 0.999) {
      return (1.00060006 + 2.63991156e2 * u + 4.37320068e3 * v) /
             ((1 + 2.57368075e2 * u + 3.41448018e3 * v) * u);
    } else {
      return (1.00001538 + 6.07514119e3 * u + 7.34266409e5 * v) /
             ((1 + 6.06511919e3 * u + 6.94021044e5 * v) * u);
    }
  }
}

double RndmVavilov(const double rkappa, const double beta2) {
  double ran = RndmUniform();
  constexpr double bkmnx1 = 0.02;
  constexpr double bkmny1 = 0.05;
  constexpr double bkmnx2 = 0.12;
  constexpr double bkmny2 = 0.05;
  constexpr double bkmnx3 = 0.22;
  constexpr double bkmny3 = 0.05;
  constexpr double bkmxx1 = 0.1;
  constexpr double bkmxy1 = 1;
  constexpr double bkmxx2 = 0.2;
  constexpr double bkmxy2 = 1;
  constexpr double bkmxx3 = 0.3;
  constexpr double bkmxy3 = 1;
  constexpr double fbkx1 = 2 / (bkmxx1 - bkmnx1);
  constexpr double fbkx2 = 2 / (bkmxx2 - bkmnx2);
  constexpr double fbkx3 = 2 / (bkmxx3 - bkmnx3);
  constexpr double fbky1 = 2 / (bkmxy1 - bkmny1);
  constexpr double fbky2 = 2 / (bkmxy2 - bkmny2);
  constexpr double fbky3 = 2 / (bkmxy3 - bkmny3);

  double ac[14] = {0};
  double hc[9] = {0};
  double h[9] = {0};
  double drk[5] = {0};
  double dsigm[5] = {0};
  double alfa[5] = {0};

  double fninv[] = {1, 0.5, 0.33333333, 0.25, 0.2};

  double edgec[] = {0.,            0.16666667e+0, 0.41666667e-1, 0.83333333e-2,
                    0.13888889e-1, 0.69444444e-2, 0.77160493e-3};

  double u1[] = {0.25850868e+0, 0.32477982e-1, -0.59020496e-2, 0.,  
                 0.24880692e-1, 0.47404356e-2, -0.74445130e-3, 0.73225731e-2,
                 0.,            0.11668284e-2, 0.,             -0.15727318e-2,
                 -0.11210142e-2};

  double u2[] = {0.43142611e+0, 0.40797543e-1, -0.91490215e-2, 0.,  
                 0.42127077e-1, 0.73167928e-2, -0.14026047e-2, 0.16195241e-1,
                 0.24714789e-2, 0.20751278e-2, 0.,             -0.25141668e-2,
                 -0.14064022e-2};

  double u3[] = {0.25225955e+0, 0.64820468e-1, -0.23615759e-1, 0.,
                 0.23834176e-1, 0.21624675e-2, -0.26865597e-2, -0.54891384e-2,
                 0.39800522e-2, 0.48447456e-2, -0.89439554e-2, -0.62756944e-2,
                 -0.24655436e-2};

  double u4[] = {0.12593231e+1,  -0.20374501e+0, 0.95055662e-1, -0.20771531e-1,
                 -0.46865180e-1, -0.77222986e-2, 0.32241039e-2, 0.89882920e-2,
                 -0.67167236e-2, -0.13049241e-1, 0.18786468e-1, 0.14484097e-1};

  double u5[] = {-0.24864376e-1, -0.10368495e-2, 0.14330117e-2, 0.20052730e-3,
                 0.18751903e-2,  0.12668869e-2,  0.48736023e-3, 0.34850854e-2,
                 0.,             -0.36597173e-3, 0.19372124e-2, 0.70761825e-3,
                 0.46898375e-3};

  double u6[] = {0.35855696e-1,  -0.27542114e-1, 0.12631023e-1, -0.30188807e-2,
                 -0.84479939e-3, 0.,             0.45675843e-3, -0.69836141e-2,
                 0.39876546e-2,  -0.36055679e-2, 0.,            0.15298434e-2,
                 0.19247256e-2};

  double u7[] = {0.10234691e+2,  -0.35619655e+1, 0.69387764e+0, -0.14047599e+0,
                 -0.19952390e+1, -0.45679694e+0, 0.,            0.50505298e+0};

  double u8[] = {0.21487518e+2,  -0.11825253e+2, 0.43133087e+1,  -0.14500543e+1,
                 -0.34343169e+1, -0.11063164e+1, -0.21000819e+0, 0.17891643e+1,
                 -0.89601916e+0, 0.39120793e+0,  0.73410606e+0,  0.,  
                 -0.32454506e+0};

  double v1[] = {0.27827257e+0,  -0.14227603e-2, 0.24848327e-2,  0.,  
                 0.45091424e-1,  0.80559636e-2,  -0.38974523e-2, 0.,  
                 -0.30634124e-2, 0.75633702e-3,  0.54730726e-2,  0.19792507e-2};

  double v2[] = {0.41421789e+0,  -0.30061649e-1, 0.52249697e-2,  0.,  
                 0.12693873e+0,  0.22999801e-1,  -0.86792801e-2, 0.31875584e-1,
                 -0.61757928e-2, 0.,             0.19716857e-1,  0.32596742e-2};

  double v3[] = {0.20191056e+0,  -0.46831422e-1, 0.96777473e-2,  -0.17995317e-2,
                 0.53921588e-1,  0.35068740e-2,  -0.12621494e-1, -0.54996531e-2,
                 -0.90029985e-2, 0.34958743e-2,  0.18513506e-1,  0.68332334e-2,
                 -0.12940502e-2};

  double v4[] = {0.13206081e+1,  0.10036618e+0,  -0.22015201e-1,
                 0.61667091e-2,  -0.14986093e+0, -0.12720568e-1,
                 0.24972042e-1,  -0.97751962e-2, 0.26087455e-1,
                 -0.11399062e-1, -0.48282515e-1, -0.98552378e-2};

  double v5[] = {0.16435243e-1,  0.36051400e-1, 0.23036520e-2,  -0.61666343e-3,
                 -0.10775802e-1, 0.51476061e-2, 0.56856517e-2,  -0.13438433e-1,
                 0.,             0.,            -0.25421507e-2, 0.20169108e-2,
                 -0.15144931e-2};

  double v6[] = {0.33432405e-1,  0.60583916e-2,  -0.23381379e-2, 0.83846081e-3,
                 -0.13346861e-1, -0.17402116e-2, 0.21052496e-2,  0.15528195e-2,
                 0.21900670e-2,  -0.13202847e-2, -0.45124157e-2, -0.15629454e-2,
                 0.22499176e-3};

  double v7[] = {0.54529572e+1,  -0.90906096e+0, 0.86122438e-1,  0.,  
                 -0.12218009e+1, -0.32324120e+0, -0.27373591e-1, 0.12173464e+0,
                 0.,             0.,             0.40917471e-1};

  double v8[] = {0.93841352e+1,  -0.16276904e+1, 0.16571423e+0,  0.,  
                 -0.18160479e+1, -0.50919193e+0, -0.51384654e-1, 0.21413992e+0,
                 0.,             0.,             0.66596366e-1};

  double w1[] = {0.29712951e+0, 0.97572934e-2, 0.,             -0.15291686e-2,
                 0.35707399e-1, 0.96221631e-2, -0.18402821e-2, -0.49821585e-2,
                 0.18831112e-2, 0.43541673e-2, 0.20301312e-2,  -0.18723311e-2,
                 -0.73403108e-3};

  double w2[] = {0.40882635e+0, 0.14474912e-1, 0.25023704e-2, -0.37707379e-2,
                 0.18719727e+0, 0.56954987e-1, 0.,            0.23020158e-1,
                 0.50574313e-2, 0.94550140e-2, 0.19300232e-1};

  double w3[] = {0.16861629e+0, 0.,            0.36317285e-2,  -0.43657818e-2,
                 0.30144338e-1, 0.13891826e-1, -0.58030495e-2, -0.38717547e-2,
                 0.85359607e-2, 0.14507659e-1, 0.82387775e-2,  -0.10116105e-1,
                 -0.55135670e-2};

  double w4[] = {0.13493891e+1,  -0.26863185e-2, -0.35216040e-2, 0.24434909e-1,
                 -0.83447911e-1, -0.48061360e-1, 0.76473951e-2,  0.24494430e-1,
                 -0.16209200e-1, -0.37768479e-1, -0.47890063e-1, 0.17778596e-1,
                 0.13179324e-1};

  double w5[] = {0.10264945e+0,  0.32738857e-1,  0.,             0.43608779e-2,
                 -0.43097757e-1, -0.22647176e-2, 0.94531290e-2,  -0.12442571e-1,
                 -0.32283517e-2, -0.75640352e-2, -0.88293329e-2, 0.52537299e-2,
                 0.13340546e-2};

  double w6[] = {0.29568177e-1,  -0.16300060e-2, -0.21119745e-3, 0.23599053e-2,
                 -0.48515387e-2, -0.40797531e-2, 0.40403265e-3,  0.18200105e-2,
                 -0.14346306e-2, -0.39165276e-2, -0.37432073e-2, 0.19950380e-2,
                 0.12222675e-2};

  double w8[] = {0.66184645e+1,  -0.73866379e+0, 0.44693973e-1,  0.,  
                 -0.14540925e+1, -0.39529833e+0, -0.44293243e-1, 0.88741049e-1};

  if (rkappa < 0.01 || rkappa > 12) {
    return 0.;
  }

  int itype = 0;
  int npt = 1;
  if (rkappa >= 0.29) {
    itype = 1;
    npt = 100;
    const double wk = 1. / sqrt(rkappa);
    ac[0] = (-0.032227 * beta2 - 0.074275) * rkappa +
            (0.24533 * beta2 + 0.070152) * wk + (-0.55610 * beta2 - 3.1579);
    ac[8] = (-0.013483 * beta2 - 0.048801) * rkappa +
            (-1.6921 * beta2 + 8.3656) * wk + (-0.73275 * beta2 - 3.5226);
    drk[0] = wk * wk;
    dsigm[0] = sqrt(rkappa / (1 - 0.5 * beta2));
    for (int j = 1; j <= 4; j++) {
      drk[j] = drk[0] * drk[j - 1];
      dsigm[j] = dsigm[0] * dsigm[j - 1];
      alfa[j] = (fninv[j - 1] - beta2 * fninv[j]) * drk[j - 1];
    }
    hc[0] = log(rkappa) + beta2 + 1. - Gamma;
    hc[1] = dsigm[0];
    hc[2] = alfa[2] * dsigm[2];
    hc[3] = (3 * alfa[1] * alfa[1] + alfa[3]) * dsigm[3] - 3;
    hc[4] = (10 * alfa[1] * alfa[2] + alfa[4]) * dsigm[4] - 10 * hc[2];
    hc[5] = hc[2] * hc[2];
    hc[6] = hc[2] * hc[3];
    hc[7] = hc[2] * hc[5];
    for (int j = 2; j <= 7; j++) {
      hc[j] = edgec[j - 1] * hc[j];
    }
    hc[8] = 0.39894228 * hc[1];
  } else if (rkappa >= 0.22) {
    itype = 2;
    npt = 150;
    const double x = 1 + (rkappa - bkmxx3) * fbkx3;
    const double y = 1 + (sqrt(beta2) - bkmxy3) * fbky3;
    const double xx = 2 * x;
    const double yy = 2 * y;
    const double x2 = xx * x - 1;
    const double x3 = xx * x2 - x;
    const double y2 = yy * y - 1;
    const double y3 = yy * y2 - y;
    const double xy = x * y;
    const double p2 = x2 * y;
    const double p3 = x3 * y;
    const double q2 = y2 * x;
    const double q3 = y3 * x;
    const double pq = x2 * y2;
    ac[1] = w1[0] + w1[1] * x + w1[3] * x3 + w1[4] * y + w1[5] * y2 +
            w1[6] * y3 + w1[7] * xy + w1[8] * p2 + w1[9] * p3 + w1[10] * q2 +
            w1[11] * q3 + w1[12] * pq;
    ac[2] = w2[0] + w2[1] * x + w2[2] * x2 + w2[3] * x3 + w2[4] * y +
            w2[5] * y2 + w2[7] * xy + w2[8] * p2 + w2[9] * p3 + w2[10] * q2;
    ac[3] = w3[0] + w3[2] * x2 + w3[3] * x3 + w3[4] * y + w3[5] * y2 +
            w3[6] * y3 + w3[7] * xy + w3[8] * p2 + w3[9] * p3 + w3[10] * q2 +
            w3[11] * q3 + w3[12] * pq;
    ac[4] = w4[0] + w4[1] * x + w4[2] * x2 + w4[3] * x3 + w4[4] * y +
            w4[5] * y2 + w4[6] * y3 + w4[7] * xy + w4[8] * p2 + w4[9] * p3 +
            w4[10] * q2 + w4[11] * q3 + w4[12] * pq;
    ac[5] = w5[0] + w5[1] * x + w5[3] * x3 + w5[4] * y + w5[5] * y2 +
            w5[6] * y3 + w5[7] * xy + w5[8] * p2 + w5[9] * p3 + w5[10] * q2 +
            w5[11] * q3 + w5[12] * pq;
    ac[6] = w6[0] + w6[1] * x + w6[2] * x2 + w6[3] * x3 + w6[4] * y +
            w6[5] * y2 + w6[6] * y3 + w6[7] * xy + w6[8] * p2 + w6[9] * p3 +
            w6[10] * q2 + w6[11] * q3 + w6[12] * pq;
    ac[8] = w8[0] + w8[1] * x + w8[2] * x2 + w8[4] * y + w8[5] * y2 +
            w8[6] * y3 + w8[7] * xy;
    ac[0] = -3.05;
  } else if (rkappa >= 0.12) {
    itype = 3;
    npt = 200;
    const double x = 1 + (rkappa - bkmxx2) * fbkx2;
    const double y = 1 + (sqrt(beta2) - bkmxy2) * fbky2;
    const double xx = 2 * x;
    const double yy = 2 * y;
    const double x2 = xx * x - 1;
    const double x3 = xx * x2 - x;
    const double y2 = yy * y - 1;
    const double y3 = yy * y2 - y;
    const double xy = x * y;
    const double p2 = x2 * y;
    const double p3 = x3 * y;
    const double q2 = y2 * x;
    const double q3 = y3 * x;
    const double pq = x2 * y2;
    ac[1] = v1[0] + v1[1] * x + v1[2] * x2 + v1[4] * y + v1[5] * y2 +
            v1[6] * y3 + v1[8] * p2 + v1[9] * p3 + v1[10] * q2 + v1[11] * q3;
    ac[2] = v2[0] + v2[1] * x + v2[2] * x2 + v2[4] * y + v2[5] * y2 +
            v2[6] * y3 + v2[7] * xy + v2[8] * p2 + v2[10] * q2 + v2[11] * q3;
    ac[3] = v3[0] + v3[1] * x + v3[2] * x2 + v3[3] * x3 + v3[4] * y +
            v3[5] * y2 + v3[6] * y3 + v3[7] * xy + v3[8] * p2 + v3[9] * p3 +
            v3[10] * q2 + v3[11] * q3 + v3[12] * pq;
    ac[4] = v4[0] + v4[1] * x + v4[2] * x2 + v4[3] * x3 + v4[4] * y +
            v4[5] * y2 + v4[6] * y3 + v4[7] * xy + v4[8] * p2 + v4[9] * p3 +
            v4[10] * q2 + v4[11] * q3;
    ac[5] = v5[0] + v5[1] * x + v5[2] * x2 + v5[3] * x3 + v5[4] * y +
            v5[5] * y2 + v5[6] * y3 + v5[7] * xy + v5[10] * q2 + v5[11] * q3 +
            v5[12] * pq;
    ac[6] = v6[0] + v6[1] * x + v6[2] * x2 + v6[3] * x3 + v6[4] * y +
            v6[5] * y2 + v6[6] * y3 + v6[7] * xy + v6[8] * p2 + v6[9] * p3 +
            v6[10] * q2 + v6[11] * q3 + v6[12] * pq;
    ac[7] = v7[0] + v7[1] * x + v7[2] * x2 + v7[4] * y + v7[5] * y2 +
            v7[6] * y3 + v7[7] * xy + v7[10] * q2;
    ac[8] = v8[0] + v8[1] * x + v8[2] * x2 + v8[4] * y + v8[5] * y2 +
            v8[6] * y3 + v8[7] * xy + v8[10] * q2;
    ac[0] = -3.04;
  } else {
    itype = 4;
    if (rkappa >= 0.02) {
      itype = 3;
    }
    npt = 200;
    const double x = 1 + (rkappa - bkmxx1) * fbkx1;
    const double y = 1 + (sqrt(beta2) - bkmxy1) * fbky1;
    const double xx = 2 * x;
    const double yy = 2 * y;
    const double x2 = xx * x - 1;
    const double x3 = xx * x2 - x;
    const double y2 = yy * y - 1;
    const double y3 = yy * y2 - y;
    const double xy = x * y;
    const double p2 = x2 * y;
    const double p3 = x3 * y;
    const double q2 = y2 * x;
    const double q3 = y3 * x;
    const double pq = x2 * y2;
    if (itype == 3) {
      ac[1] = u1[0] + u1[1] * x + u1[2] * x2 + u1[4] * y + u1[5] * y2 +
              u1[6] * y3 + u1[7] * xy + u1[9] * p3 + u1[11] * q3 + u1[12] * pq;
      ac[2] = u2[0] + u2[1] * x + u2[2] * x2 + u2[4] * y + u2[5] * y2 +
              u2[6] * y3 + u2[7] * xy + u2[8] * p2 + u2[9] * p3 + u2[11] * q3 +
              u2[12] * pq;
      ac[3] = u3[0] + u3[1] * x + u3[2] * x2 + u3[4] * y + u3[5] * y2 +
              u3[6] * y3 + u3[7] * xy + u3[8] * p2 + u3[9] * p3 + u3[10] * q2 +
              u3[11] * q3 + u3[12] * pq;
      ac[4] = u4[0] + u4[1] * x + u4[2] * x2 + u4[3] * x3 + u4[4] * y +
              u4[5] * y2 + u4[6] * y3 + u4[7] * xy + u4[8] * p2 + u4[9] * p3 +
              u4[10] * q2 + u4[11] * q3;
      ac[5] = u5[0] + u5[1] * x + u5[2] * x2 + u5[3] * x3 + u5[4] * y +
              u5[5] * y2 + u5[6] * y3 + u5[7] * xy + u5[9] * p3 + u5[10] * q2 +
              u5[11] * q3 + u5[12] * pq;
      ac[6] = u6[0] + u6[1] * x + u6[2] * x2 + u6[3] * x3 + u6[4] * y +
              u6[6] * y3 + u6[7] * xy + u6[8] * p2 + u6[9] * p3 + u6[11] * q3 +
              u6[12] * pq;
      ac[7] = u7[0] + u7[1] * x + u7[2] * x2 + u7[3] * x3 + u7[4] * y +
              u7[5] * y2 + u7[7] * xy;
    }
    ac[8] = u8[0] + u8[1] * x + u8[2] * x2 + u8[3] * x3 + u8[4] * y +
            u8[5] * y2 + u8[6] * y3 + u8[7] * xy + u8[8] * p2 + u8[9] * p3 +
            u8[10] * q2 + u8[12] * pq;
    ac[0] = -3.03;
  }
  ac[9] = (ac[8] - ac[0]) / npt;
  if (itype == 3) {
    const double x = (ac[7] - ac[8]) / (ac[7] * ac[8]);
    const double y = 1 / log(ac[8] / ac[7]);
    const double p2 = ac[7] * ac[7];
    ac[11] = p2 *
             (ac[1] * exp(-ac[2] * (ac[7] + ac[5] * p2) -
                          ac[3] * exp(-ac[4] * (ac[7] + ac[6] * p2))) -
              0.045 * y / ac[7]) /
             (1 + x * y * ac[7]);
    ac[12] = (0.045 + x * ac[11]) * y;
  }
  if (itype == 4) {
    ac[10] = 0.995 / TMath::LandauI(ac[8]);
  }

  const double t = 2 * ran / ac[9];
  double rlam = ac[0];
  double fl = 0.;
  double fu = 0.;
  double s = 0;
  for (int n = 1; n <= npt; n++) {
    rlam += ac[9];
    if (itype == 1) {
      double fn = 1;
      const double x = (rlam + hc[0]) * hc[1];
      h[0] = x;
      h[1] = x * x - 1;
      for (int k = 2; k <= 8; k++) {
        fn += 1;
        h[k] = x * h[k - 1] - fn * h[k - 2];
      }
      double y = 1 + hc[7] * h[8];
      for (int k = 2; k <= 6; k++) {
        y = y + hc[k] * h[k];
      }
      if (y < 0) {
        fu = 0;
      } else {
        fu = hc[8] * exp(-0.5 * x * x) * y;
      }
    } else if (itype == 2) {
      const double x = rlam * rlam;
      fu = ac[1] * exp(-ac[2] * (rlam + ac[5] * x) -
                       ac[3] * exp(-ac[4] * (rlam + ac[6] * x)));
    } else if (itype == 3) {
      if (rlam < ac[7]) {
        const double x = rlam * rlam;
        fu = ac[1] * exp(-ac[2] * (rlam + ac[5] * x) -
                         ac[3] * exp(-ac[4] * (rlam + ac[6] * x)));
      } else {
        const double x = 1 / rlam;
        fu = (ac[11] * x + ac[12]) * x;
      }
    } else {
      fu = ac[10] * denlan(rlam);
    }
    s += fl + fu;
    if (s > t) {
      break;
    }
    fl = fu;
  }

  const double s0 = s - fl - fu;
  double v = rlam - ac[9];
  if (s > s0) {
    v += ac[9] * (t - s0) / (s - s0);
  }
  return v;
}

double RndmHeedWF(const double w, const double f) {
  // RNDHWF - Generates random energies needed to create a single e- in
  //          a gas with asymptotic work function W and Fano factor F,
  //          according to Igor Smirnov's phenomenological model.
  constexpr double wref = 30.0;
  constexpr double fref = 0.174;
  // Check parameters
  if (w <= 0 || f < 0) {
    std::cerr << "RndmHeedWF: Work and/or Fano parameter out of range. "
              << "Returning 0.\n";
    return 0.;
  } else if (f == 0) {
    // Special case of F = 0
    return w;
  }
  // First generate a standardised (W = 30, F = 0.174) random energy.
  const double x = RndmUniform() * wref * 0.82174;
  double e;
  if (x < 0) {
    // E = 0 to w/2:     p = 0,       integral = 0
    std::cerr << "RndmHeedWF: Random number is below the applicable range. "
              << "Program error. Returning w/2.\n";
    e = 0.5 * wref;
  } else if (x < 0.5 * wref) {
    // E = w/2 to w:     p = 1,       integral = E-w/2
    e = 0.5 * wref + x;
  } else if (x < wref * 0.82174) {
    // E = w to 3.064 w: p = (w/E)^4, integral = w^4/3 (1/E^3 - 1/w^3)
    constexpr double w4 = wref * wref * wref * wref;
    e = pow(2 * w4 / (5 * wref - 6 * x), 1.0 / 3.0);
  } else {
    // E > 3.064 w:      p = 0,       integral = 0
    std::cerr << "RndmHeedWF: Random number is above applicable range. "
              << "Program error. Returning 3.064 w.\n";
    e = 3.064 * wref;
  }
  // Scale.
  const double sqf = sqrt(f / fref);
  return (w / wref) * sqf * e + w * (1. - sqf);
}
}
