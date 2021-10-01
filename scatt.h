/************************************************************************
Copyright [2021] [Christian B. Huebschle chuebsch@moliso.de]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
************************************************************************/
#ifndef SCATT_H
#define SCATT_H 1
#include <math.h>
#include <QString>
#include "molecule.h"

// lab,rad,mass,a1,b1,a2,b2,a3,b3,a4,b4,c,f'cu,f''cu,mucu,f'mo,f''mo,mumo,f'ag,f''ag,muag

struct Atomformfactor{
  char symb[3];
  double rad,mass;
  double ab4c[9];
  double fpcu,fppcu,mucu,fpmo,fppmo,mumo,fpag,fppag,muag;
};

Atomformfactor ff[94]= {//copy paste from shelxl.f 97
  {"H",   .32,   1.008,{ 0.49300, 10.51091,  0.32291, 26.12573,  0.14019,   3.14236, 0.04081,  57.79977,  0.00304},  0.0000,  0.0000,    .655,   0.000, 0.000,  .624,  0.000, 0.000,  .614 },
  {"HE", 1.50,    4.00,{ 0.87340,  9.10371,  0.63090,  3.35680,  0.31120,  22.92763, 0.17800,   0.98210,  0.00640},  0.0000,  0.0000,    1.94,   0.000, 0.000,  1.34,  0.000, 0.000,  1.28 },
  {"LI", 1.52,    6.94,{ 1.12820,  3.95460,  0.75080,  1.05240,  0.61750,  85.39058, 0.46530, 168.26120,  0.03770},  0.0008,  0.0003,    5.76, -0.0003, 0.0001,  2.28, -0.0004, 0.0000,  2.06 },
  {"BE", 1.11,    9.01,{ 1.59190, 43.64275,  1.12780,  1.86230,  0.53910, 103.48310, 0.70290,   0.54200,  0.03850},  0.0038,  0.0014,    16.6,  0.0005, 0.0002,  3.83,  0.0001, 0.0001,  3.13 },
  {"B",  0.82,   10.81,{ 2.05450, 23.21852,  1.33260,  1.02100,  1.09790,  60.34987, 0.70680,   0.14030, -0.19320},  0.0090,  0.0039,    41.5,  0.0013, 0.0007,  6.61,  0.0004, 0.0004,  4.79 },
  {"C",  0.77,   12.01,{ 2.31000, 20.84392,  1.02000, 10.20751,  1.58860,   0.56870, 0.86500,  51.65125,  0.21560},  0.0181,  0.0091,    89.9,  0.0033, 0.0016,  11.5,  0.0015, 0.0009,  7.45 },
  {"N",  0.70,   14.01,{12.21261,  0.00570,  3.13220,  9.89331,  2.01250,  28.99754, 1.16630,   0.58260,-11.52901},  0.0311,  0.0180,    173.,  0.0061, 0.0033,  19.6,  0.0030, 0.0019,  11.7 },
  {"O",  0.66,   16.00,{ 3.04850, 13.27711,  2.28680,  5.70111,  1.54630,   0.32390, 0.86700,  32.90894,  0.25080},  0.0492,  0.0322,    304.,  0.0106, 0.0060,  32.5,  0.0056, 0.0036,  18.2 },
  {"F",  0.64,   19.00,{ 3.53920, 10.28251,  2.64120,  4.29440,  1.51700,   0.26150, 1.02430,  26.14763,  0.27760},  0.0727,  0.0534,    498.,  0.0171, 0.0103,  51.5,  0.0096, 0.0061,  27.7 },
  {"NE", 1.50,   20.18,{ 3.95530,  8.40421,  3.11250,  3.42620,  1.45460,   0.23060, 1.12510,  21.71841,  0.35150},  0.1019,  0.0833,    768.,  0.0259, 0.0164,  78.6,  0.0152, 0.0098,  41.2 },
  {"NA", 1.86,   22.99,{ 4.76260,  3.28500,  3.17360,  8.84221,  1.26740,   0.31360, 1.11280, 129.42410,  0.67600},  0.1353,  0.1239,   1140.,  0.0362, 0.0249,  116.,  0.0218, 0.0150,  59.6 },
  {"MG", 1.60,   24.31,{ 5.42041,  2.82750,  2.17350, 79.26118,  1.22690,   0.38080, 2.30730,   7.19371,  0.85840},  0.1719,  0.1771,   1610.,  0.0486, 0.0363,  165.,  0.0298, 0.0220,  84.2 },
  {"AL", 1.25,   26.98,{ 6.42021,  3.03870,  1.90020,  0.74260,  1.59360,  31.54724, 1.96460,  85.08868,  1.11510},  0.2130,  0.2455,   2220.,  0.0645, 0.0514,  229.,  0.0406, 0.0313,  116. },
  {"SI", 1.17,   28.09,{ 6.29151,  2.43860,  3.03530, 32.33374,  1.98910,   0.67850, 1.54100,  81.69379,  1.14070},  0.2541,  0.3302,   2970.,  0.0817, 0.0704,  310.,  0.0522, 0.0431,  156. },
  {"P",  1.10,   30.97,{ 6.43451,  1.90670,  4.17910, 27.15704,  1.78000,   0.52600, 1.49080,  68.16457,  1.11490},  0.2955,  0.4335,   3880.,  0.1023, 0.0942,  410.,  0.0667, 0.0580,  206. },
  {"S",  1.03,   32.06,{ 6.90531,  1.46790,  5.20341, 22.21512,  1.43790,   0.25360, 1.58630,  56.17207,  0.86690},  0.3331,  0.5567,   4970.,  0.1246, 0.1234,  532.,  0.0826, 0.0763,  267. },
  {"CL", 0.99,   35.45,{11.46041,  0.01040,  7.19641,  1.16620,  6.25561,  18.51942, 1.64550,  47.77846, -9.55741},  0.3639,  0.7018,   6240.,  0.1484, 0.1585,  678.,  0.0998, 0.0984,  341. },
  {"AR", 1.50,   39.95,{ 7.48451,  0.90720,  6.77231, 14.84071,  0.65390,  43.89835, 1.64420,  33.39293,  1.44450},  0.3843,  0.8717,   7720.,  0.1743, 0.2003,  851.,  0.1191, 0.1249,  429. },
  {"K",  2.27,   39.10,{ 8.21861, 12.79491,  7.43981,  0.77480,  1.05190, 213.18720, 0.86590,  41.68416,  1.42280},  0.3868,  1.0657,   9400.,  0.2009, 0.2494, 1050.,  0.1399, 0.1562,  532. },
  {"CA", 1.97,   40.08,{ 8.62661, 10.44211,  7.38731,  0.65990,  1.58990,  85.74849, 1.02110, 178.43720,  1.37510},  0.3641,  1.2855,  11300.,  0.2262, 0.3064, 1290.,  0.1611, 0.1926,  652. },
  {"SC", 1.61,   44.96,{ 9.18901,  9.02131,  7.36791,  0.57290,  1.64090, 136.10810, 1.46800,  51.35315,  1.33290},  0.3119,  1.5331,  13500.,  0.2519, 0.3716, 1560.,  0.1829, 0.2348,  789. },
  {"TI", 1.45,   47.90,{ 9.75951,  7.85081,  7.35581,  0.50000,  1.69910,  35.63383, 1.90210, 116.10510,  1.28070},  0.2191,  1.8069,  15900.,  0.2776, 0.4457, 1860.,  0.2060, 0.2830,  947. },
  {"V",  1.31,   50.94,{10.29711,  6.86571,  7.35111,  0.43850,  2.07030,  26.89383, 2.05710, 102.47810,  1.21990},  0.0687,  2.1097,  18500.,  0.3005, 0.5294, 2200.,  0.2276, 0.3376, 1120. },
  {"CR", 1.24,   52.00,{10.64061,  6.10381,  7.35371,  0.39200,  3.32400,  20.26262, 1.49220,  98.73999,  1.18320}, -0.1635,  2.4439,  21300.,  0.3209, 0.6236, 2580.,  0.2496, 0.3992, 1330. },
  {"MN", 1.37,   54.94,{11.28191,  5.34091,  7.35731,  0.34320,  3.01930,  17.86742, 2.24410,  83.75438,  1.08960}, -0.5299,  2.8052,  24600.,  0.3368, 0.7283, 3020.,  0.2704, 0.4681, 1550. },
  {"FE", 1.24,   55.85,{11.76951,  4.76111,  7.35731,  0.30720,  3.52220,  15.35351, 2.30450,  76.88058,  1.03690}, -1.1336,  3.1974,  28000.,  0.3463, 0.8444, 3490.,  0.2886, 0.5448, 1800. },
  {"CO", 1.25,   58.93,{12.28411,  4.27910,  7.34091,  0.27840,  4.00340,  13.53591, 2.34880,  71.16927,  1.01180}, -2.3653,  3.6143,  31400.,  0.3494, 0.9721, 4010.,  0.3050, 0.6296, 2070. },
  {"NI", 1.25,   58.71,{12.83761,  3.87850,  7.29201,  0.25650,  4.44380,  12.17631, 2.38000,  66.34216,  1.03410}, -3.0029,  0.5091,   4760.,  0.3393, 1.1124, 4570.,  0.3147, 0.7232, 2380. },
  {"CU", 1.28,   63.54,{13.33801,  3.58280,  7.16761,  0.24700,  5.61581,  11.39661, 1.67350,  64.81267,  1.19100}, -1.9646,  0.5888,   5470.,  0.3201, 1.2651, 5180.,  0.3240, 0.8257, 2710. },
  {"ZN", 1.33,   65.37,{14.07431,  3.26550,  7.03181,  0.23330,  5.16521,  10.31631, 2.41000,  58.70976,  1.30410}, -1.5491,  0.6778,   6290.,  0.2839, 1.4301, 5860.,  0.3242, 0.9375, 3070. },
  {"GA", 1.26,   69.72,{15.23541,  3.06690,  6.70061,  0.24120,  4.35910,  10.78051, 2.96230,  61.41357,  1.71890}, -1.2846,  0.7763,   7190.,  0.2307, 1.6083, 6600.,  0.3179, 1.0589, 3460. },
  {"GE", 1.22,   72.59,{16.08162,  2.85090,  6.37471,  0.25160,  3.70680,  11.44681, 3.68300,  54.76256,  2.13130}, -1.0885,  0.8855,   8190.,  0.1547, 1.8001, 7380.,  0.3016, 1.1903, 3870. },
  {"AS", 1.21,   74.92,{16.67232,  2.63450,  6.07011,  0.26470,  3.43130,  12.94791, 4.27790,  47.79726,  2.53100}, -0.9300,  1.0051,   9290.,  0.0499, 2.0058, 8220.,  0.2758, 1.3314, 4330. },
  {"SE", 1.17,   78.96,{17.00063,  2.40980,  5.81961,  0.27260,  3.97310,  15.23721, 4.35430,  43.81635,  2.84090}, -0.7943,  1.1372,  10500., -0.0929, 2.2259, 9110.,  0.2367, 1.4831, 4820. },
  {"BR", 1.14,   79.91,{17.17892,  2.17230,  5.23581, 16.57962,  5.63771,   0.26090, 3.98510,  41.43285,  2.95570}, -0.6763,  1.2805,  11800., -0.2901, 2.4595,10000.,  0.1811, 1.6452, 5350. },
  {"KR", 1.50,   83.80,{17.35551,  1.93840,  6.72861, 16.56232,  5.54931,   0.22610, 3.53750,  39.39723,  2.82500}, -0.5657,  1.4385,  13200., -0.5574, 2.7079,11000.,  0.1067, 1.8192, 5920. },
  {"RB", 2.48,   85.47,{17.17842,  1.78880,  9.64351, 17.31512,  5.13990,   0.27480, 1.52920, 164.93420,  3.48730}, -0.4688,  1.6079,  14800., -0.9393, 2.9676,12100.,  0.0068, 2.0025, 6520. },
  {"SR", 2.15,   87.62,{17.56631,  1.55640,  9.81841, 14.09881,  5.42200,   0.16640, 2.66940, 132.37610,  2.50640}, -0.3528,  1.8200,  16500., -1.5307, 3.2498,13200., -0.1172, 2.2025, 7150. },
  {"Y",  1.78,   88.91,{17.77602,  1.40290, 10.29461, 12.80061,  5.72630,   0.12560, 3.26588, 104.35410,  1.91213}, -0.2670,  2.0244,  18300., -2.7962, 3.5667,14300., -0.2879, 2.4099, 7800. },
  {"ZR", 1.59,   91.22,{17.87653,  1.27618, 10.94801, 11.91601,  5.41733,   0.11762, 3.65721,  87.66278,  2.06929}, -0.1862,  2.2449,  20300., -2.9673, 0.5597, 2470., -0.5364, 2.6141, 8470. },
  {"NB", 1.43,   92.91,{17.61423,  1.18865, 12.01441, 11.76601,  4.04183,   0.20479, 3.53346,  69.79576,  3.75591}, -0.1121,  2.4826,  22300., -2.0727, 0.6215, 2730., -0.8282, 2.8404, 9220. },
  {"MO", 1.36,   95.94,{ 3.70250,  0.27720, 17.23563,  1.09580, 12.88761,  11.00401, 3.74290,  61.65846,  4.38750}, -0.0483,  2.7339,  24600., -1.6832, 0.6857, 3000., -1.2703, 3.0978, 11500. },
  {"TC", 1.35,   98.00,{19.13013,  0.86413, 11.09481,  8.14488,  4.64902,  21.57072, 2.71263,  86.84727,  5.40429},  0.0057,  3.0049,  27000., -1.4390, 0.7593, 3320., -2.0087, 3.3490, 10700. },
  {"RU", 1.33,  101.07,{19.26743,  0.80852, 12.91821,  8.43468,  4.86337,  24.79974, 1.56756,  94.29289,  5.37875},  0.0552,  3.2960,  29500., -1.2594, 0.8363, 3640., -5.3630, 3.6506, 1920. },
  {"RH", 1.35,  102.91,{19.29572,  0.75154, 14.35011,  8.21759,  4.73425,  25.87494, 1.28918,  98.60629,  5.32800},  0.0927,  3.6045,  32300., -1.1178, 0.9187, 3990., -2.5280, 0.5964, 2100. },
  {"PD", 1.38,  106.40,{19.33192,  0.69866, 15.50172,  7.98930,  5.29537,  25.20523, 0.60584,  76.89868,  5.26593},  0.1215,  3.9337,  35200., -0.9988, 1.0072, 4360., -1.9556, 0.6546, 2300. },
  {"AG", 1.44,  107.87,{19.28082,  0.64460, 16.68852,  7.47261,  4.80451,  24.66054, 1.04630,  99.81570,  5.17900},  0.1306,  4.2820,  38200., -0.8971, 1.1015, 4760., -1.6473, 0.7167, 2510. },
  {"CD", 1.49,  112.40,{19.22142,  0.59460, 17.64442,  6.90891,  4.46100,  24.70084, 1.60290,  87.48257,  5.06941},  0.1185,  4.6533,  41500., -0.8075, 1.2024, 5180., -1.4396, 0.7832, 2730. },
  {"IN", 1.44,  114.82,{19.16241,  0.54760, 18.55962,  6.37761,  4.29480,  25.84993, 2.03960,  92.80299,  4.93911},  0.0822,  5.0449,  45000., -0.7276, 1.3100, 5630., -1.2843, 0.8542, 2970. },
  {"SN", 1.40,  118.69,{19.18892,  5.83031, 19.10052,  0.50310,  4.45850,  26.89093, 2.46630,  83.95718,  4.78211},  0.0259,  5.4591,  48600., -0.6537, 1.4246, 6110., -1.1587, 0.9299, 3230. },
  {"SB", 1.41,  121.75,{19.64182,  5.30340, 19.04552,  0.46070,  5.03711,  27.90744, 2.68270,  75.28258,  4.59091}, -0.0562,  5.8946,  52500., -0.5866, 1.5461, 6620., -1.0547, 1.0104, 3500. },
  {"TE", 1.37,  127.60,{19.96442,  4.81742, 19.01382,  0.42089,  6.14488,  28.52844, 2.52390,  70.84036,  4.35200}, -0.1759,  6.3531,  56500., -0.5308, 1.6751, 7160., -0.9710, 1.0960, 3780. },
  {"I",  1.33,  126.90,{20.14722,  4.34700, 18.99492,  0.38140,  7.51381,  27.76604, 2.27350,  66.87767,  4.07120}, -0.3257,  6.8362,  60700., -0.4742, 1.8119, 7730., -0.8919, 1.1868, 4090. },
  {"XE", 1.50,  131.30,{20.29332,  3.92820, 19.02982,  0.34400,  8.97671,  26.46594, 1.99000,  64.26587,  3.71180}, -0.5179,  7.3500,  65200., -0.4205, 1.9578, 8340., -0.8200, 1.2838, 4410. },
  {"CS", 2.65,  132.91,{20.38922,  3.56900, 19.10622,  0.31070, 10.66201,  24.38794, 1.49530, 213.90420,  3.33520}, -0.7457,  7.9052,  70000., -0.3680, 2.1192, 8980., -0.7527, 1.3916, 4750. },
  {"BA", 2.17,  137.34,{20.33612,  3.21600, 19.29703,  0.27560, 10.88801,  20.20732, 2.69590, 167.20220,  2.77310}, -1.0456,  8.4617,  75000., -0.3244, 2.2819, 9650., -0.6940, 1.5004, 5110. },
  {"LA", 1.87,  138.91,{20.57802,  2.94817, 19.59901,  0.24448, 11.37271,  18.77261, 3.28719, 133.12410,  2.14678}, -1.4094,  9.0376,  80300., -0.2871, 2.4523,10400., -0.6411, 1.6148, 5490. },
  {"CE", 1.83,  140.12,{21.16711,  2.81219, 19.76952,  0.22684, 11.85131,  17.60832, 3.33049, 127.11310,  1.86264}, -1.8482,  9.6596,  85700., -0.2486, 2.6331,11100., -0.5890, 1.7358, 5880. },
  {"PR", 1.82,  140.91,{22.04402,  2.77393, 19.66972,  0.22209, 12.38561,  16.76692, 2.82428, 143.64410,  2.05830}, -2.4164, 10.2820,  91200., -0.2180, 2.8214,11900., -0.5424, 1.8624, 6300. },
  {"ND", 1.81,  144.24,{22.68452,  2.66248, 19.68472,  0.21063, 12.77401,  15.88502, 2.85137, 137.90310,  1.98486}, -3.1807, 10.9079,  96800., -0.1943, 3.0179,12700., -0.5012, 1.9950, 6740. },
  {"PM", 1.81,  147.00,{23.34052,  2.56270, 19.60953,  0.20209, 13.12351,  15.10091, 2.87516, 132.72110,  2.02876}, -4.0598, 11.5523, 102000., -0.1753, 3.2249,13500., -0.4626, 2.1347, 7200. },
  {"SM", 1.80,  150.35,{24.00424,  2.47274, 19.42583,  0.19645, 13.43961,  14.39961, 2.89604, 128.00710,  2.20963}, -5.3236, 12.2178, 108000., -0.1638, 3.4418,14400., -0.4287, 2.2815, 7680. },
  {"EU", 2.00,  151.96,{24.62744,  2.38790, 19.08862,  0.19420, 13.76031,  13.75461, 2.92270, 123.17410,  2.57450}, -8.9294, 11.1857, 110000., -0.1578, 3.6682,15400., -0.3977, 2.4351, 8190. },
  {"GD", 1.79,  157.25,{25.07094,  2.25341, 19.07982,  0.18195, 13.85181,  12.93311, 3.54545, 101.39810,  2.41960}, -8.8380, 11.9157, 105000., -0.1653, 3.9035,16300., -0.3741, 2.5954, 8720. },
  {"TB", 1.76,  158.92,{25.89763,  2.24256, 18.21852,  0.19614, 14.31671,  12.66481, 2.95354, 115.36210,  3.58324}, -9.1472,  9.1891,  84700., -0.1723, 4.1537,17400., -0.3496, 2.7654, 9270. },
  {"DY", 1.75,  162.50,{26.50703,  2.18020, 17.63832,  0.20217, 14.55962,  12.18991, 2.96577, 111.87410,  4.29728}, -9.8046,  9.8477,  97700., -0.1892, 4.4098,18400., -0.3302, 2.9404, 9850. },
  {"HO", 1.74,  164.93,{26.90494,  2.07051, 17.29402,  0.19794, 14.55831,  11.44071, 3.63837,  92.65669,  4.56797},-14.9734,  3.7046,  34700., -0.2175, 4.6783,19500., -0.3168, 3.1241, 10400. },
  {"ER", 1.73,  167.26,{27.65634,  2.07356, 16.42853,  0.22355, 14.97791,  11.36041, 2.98233, 105.70310,  5.92047}, -9.4367,  3.9380,  36700., -0.2586, 4.9576,20700., -0.3091, 3.3158, 11100. },
  {"TM", 1.72,  168.93,{28.18193,  2.02859, 15.88512,  0.23885, 15.15421,  10.99751, 2.98706, 102.96110,  6.75622}, -8.0393,  4.1821,  39300., -0.3139, 5.2483,21900., -0.3084, 3.5155, 11700. },
  {"YB", 1.94,  173.04,{28.66414,  1.98890, 15.43451,  0.25712, 15.30871,  10.66471, 2.98963, 100.41710,  7.56673}, -7.2108,  4.4329,  41000., -0.3850, 5.5486,23100., -0.3157, 3.7229, 12400. },
  {"LU", 1.72,  174.97,{28.94763,  1.90182, 15.22081,  9.98520, 15.10001,   0.26103, 3.71601,  84.32988,  7.97629}, -6.6179,  4.6937,  45000., -0.4720, 5.8584,24400., -0.3299, 3.9377, 13100. },
  {"HF", 1.56,  178.49,{29.14404,  1.83262, 15.17261,  9.59991, 14.75861,   0.27512, 4.30013,  72.02908,  8.58155}, -6.1794,  4.9776,  46000., -0.5830, 6.1852,25800., -0.3548, 4.1643, 13900. },
  {"TA", 1.43,  180.95,{29.20244,  1.77333, 15.22931,  9.37047, 14.51351,   0.29598, 4.76492,  63.36447,  9.24355}, -5.7959,  5.2718,  48500., -0.7052, 6.5227,27200., -0.3831, 4.3992, 14600. },
  {"W",  1.37,  183.85,{29.08183,  1.72029, 15.43001,  9.22591, 14.43271,   0.32170, 5.11983,  57.05606,  9.88751}, -5.4734,  5.5774,  51300., -0.8490, 6.8722,28600., -0.4201, 4.6430, 15400. },
  {"RE", 1.37,  186.20,{28.76213,  1.67191, 15.71892,  9.09228, 14.55641,   0.35050, 5.44174,  52.08615, 10.47201}, -5.2083,  5.8923,  57200., -1.0185, 7.2310,30100., -0.4693, 4.8944, 16200. },
  {"OS", 1.34,  190.20,{28.18944,  1.62903, 16.15501,  8.97949, 14.93051,   0.38266, 5.67590,  48.16475, 11.00051}, -4.9801,  6.2216,  58000., -1.2165, 7.6030,31600., -0.5280, 5.1558, 17100. },
  {"IR", 1.36,  192.20,{27.30493,  1.59279, 16.72961,  8.86554, 15.61152,   0.41792, 5.83378,  45.00114, 11.47221}, -4.7710,  6.5667,  62400., -1.4442, 7.9887,33100., -0.5977, 5.4269, 18000. },
  {"PT", 1.37,  195.09,{27.00594,  1.51293, 17.76392,  8.81175, 15.71312,   0.42459, 5.78371,  38.61034, 11.68831}, -4.5932,  6.9264,  63400., -1.7033, 8.3905,34800., -0.6812, 5.7081, 18900. },
  {"AU", 1.44,196.9665,{16.88193,  0.46110, 18.59132,  8.62161, 25.55824,   1.48260, 5.86001,  36.39563, 12.06581}, -4.4197,  7.2980,  66900., -2.0133, 8.8022,36500., -0.7638, 5.9978, 19900. },
  {"HG", 1.50,  200.59,{20.68092,  0.54500, 19.04172,  8.44841, 21.65752,   1.57290, 5.96761,  38.32463, 12.60891}, -4.2923,  7.6849,  66800., -2.3894, 9.2266,38200., -0.8801, 6.2989, 20900. },
  {"TL", 1.64,  204.37,{27.54463,  0.65515, 19.15842,  8.70752, 15.53802,   1.96347, 5.52594,  45.81496, 13.17461}, -4.1627,  8.0900,  75400., -2.8358, 9.6688,40100., -1.0117, 6.6090, 21900. },
  {"PB", 1.60,  207.19,{31.06174,  0.69020, 13.06371,  2.35760, 18.44202,   8.61801, 5.96961,  47.25795, 13.41181}, -4.0753,  8.5060,  79800., -3.3944,10.1111,41900., -1.1676, 6.9287, 22900. },
  {"BI", 1.60,  208.98,{33.36894,  0.70400, 12.95101,  2.92380, 16.58772,   8.79371, 6.46921,  48.00935, 13.57821}, -4.0111,  8.9310,  84300., -4.1077,10.2566,43800., -1.3494, 7.2566, 24000. },
  {"PO", 1.60,  210.00,{34.67264,  0.70100, 15.47331,  3.55078, 13.11381,   9.55643, 7.02589,  47.00455, 13.67701}, -3.9670,  9.3834,  88100., -5.1210,11.0496,45800., -1.5613, 7.5986, 25100. },
  {"AT", 1.60,  210.00,{35.31633,  0.68587, 19.02112,  3.97458,  9.49888,  11.38241, 7.42519,  45.47156, 13.71081}, -3.9588,  9.8433,  86500., -7.9122, 9.9777,40700., -1.8039, 7.9509, 26200. },
  {"RN", 1.80,  222.00,{35.56314,  0.66310, 21.28162,  4.06910,  8.00371,  14.04221, 7.44331,  44.24734, 13.69051}, -3.9487, 10.3181,  97200., -8.0659,10.4580,39800., -2.0847, 8.3112, 27300. },
  {"FR", 2.80,  223.00,{35.92993,  0.64645, 23.05472,  4.17619, 12.14391,  23.10522, 2.11253, 150.64510, 13.72471}, -3.9689, 10.8038, 102000., -7.2224, 7.7847,32200., -2.4129, 8.6839, 28500. },
  {"RA", 2.20,  226.00,{35.76303,  0.61634, 22.90642,  3.87135, 12.47391,  19.98872, 3.21097, 142.32510, 13.62111}, -4.0088, 11.2969, 102000., -6.7704, 8.1435,33000., -2.8081, 9.0614, 29800. },
  {"AC", 1.90,  227.00,{35.65973,  0.58909, 23.10323,  3.65155, 12.59771,  18.59901, 4.08655, 117.02010, 13.52661}, -4.0794, 11.7994, 143000., -6.8494, 8.5178,54000., -3.2784, 9.4502, 31100. },
  {"TH", 1.85,  232.04,{35.56453,  0.56336, 23.42192,  3.46204, 12.74731,  17.83092, 4.80704,  99.17230, 13.43141}, -4.1491, 12.3296, 118000., -7.2400, 8.8979,37000., -3.8533, 9.8403, 32300. },
  {"PA", 1.80,  231.00,{35.88474,  0.54775, 23.29482,  3.41519, 14.18911,  16.92352, 4.17287, 105.25110, 13.42871}, -4.2473, 12.8681, 106000., -8.0334, 9.2807,38700., -4.6067,10.2413, 34200. },
  {"U",  1.80,  238.03,{36.02284,  0.52930, 23.41283,  3.32530, 14.94911,  16.09273, 4.18800, 100.61310, 13.39661}, -4.3638, 13.4090, 112000., -9.6767, 9.6646,40300., -5.7225,10.6428, 35000. },
  {"NP", 1.80,  237.00,{36.18744,  0.51193, 23.59642,  3.25396, 15.64022,  15.36222, 4.18550,  97.49089, 13.35731}, -4.5053, 13.9666, 123000.,-11.4937, 4.1493,25700., -6.9995, 9.5876, 29900. },
  {"PU", 1.80,  242.00,{36.52544,  0.49938, 23.80832,  3.26371, 16.77072,  14.94551, 3.47947, 105.98010, 13.38121}, -4.6563, 14.3729, 113000., -9.4100, 4.3056,16200.,-13.5905, 6.9468, 22700.}
};
#define SCATTCU 3
#define SCATTAG 2
#define SCATTMO 1
class Scatt{
  typedef struct {
    int  ih,//!< h
         ik,//!< k
         il;//!< l 
    double fo,//!< F observed
          so,//!< \f$\sigma(observed)\f$
          si,
          fc,//!< F calculated
          phi, //!< \f$\varphi\f$ 
          dis; 
  } Rec;
  public: 
  Scatt(Molecule *m) { 
    mol=m; 
    f000=0.0;
    Fmax=0.0;
    resd=-1.0;
    //printf("Scatt(Molecule *m) %f\n",resd);
    //   qDebug()<<"12345"<<atomformfactor(12,0.0);
  }

  void setSCAT(QStringList tok){
    if (tok.size()<16)return;
    int an = mol->getOZ(tok.at(1));
    if (an<0)return;
    for (int i=0; i<9; i++) {
      ff[an].ab4c[i]=tok.at(2+i).toDouble();
      //    qDebug()<<i+2<<ff[an].ab4c[i];
    }
    ff[an].fpcu=ff[an].fpmo=ff[an].fpag    = tok.at(11).toDouble();
    ff[an].fppcu=ff[an].fppmo=ff[an].fppag = tok.at(12).toDouble();
    ff[an].mucu=ff[an].mumo=ff[an].muag    = tok.at(13).toDouble();
    ff[an].rad=tok.at(14).toDouble();
    ff[an].mass=tok.at(15).toDouble();
  }

  void setDISP(QStringList tok){
    if (tok.size()<4)return;
    int an = mol->getOZ(tok.at(1));
    ff[an].fpcu=ff[an].fpmo=ff[an].fpag    = tok.at(2).toDouble();
    if (tok.size()<4)return;
    ff[an].fppcu=ff[an].fppmo=ff[an].fppag = tok.at(3).toDouble();
    if (tok.size()<5)return;
    ff[an].mucu=ff[an].mumo=ff[an].muag    = tok.at(4).toDouble();

  }

  double fprime(int an,int radiation){
    switch (radiation){
      case SCATTCU: return ff[an].fpcu;
      case SCATTMO: return ff[an].fpmo;
      case SCATTAG: return ff[an].fpag;
      default: return ff[an].fpcu; 
    } 
  }
  double fdprime(int an,int radiation){
    switch (radiation){
      case SCATTCU: return ff[an].fppcu;
      case SCATTMO: return ff[an].fppmo;
      case SCATTAG: return ff[an].fppag;
      default: return ff[an].fppcu; 
    } 
  }
  double atomformfactor(int an, double s){
    double erg=0.0;
    for (int i=0; i<4;i++){
      erg+=ff[an].ab4c[2*i]*exp(-ff[an].ab4c[2*i+1]*s);
    }
    erg+=ff[an].ab4c[8];
    return erg; 
  }
  double FF(V3 hkl,int radiation,double &phang,double &I){
    const double twopi=M_PI*2.0;
    const double twopi2=M_PI*M_PI*2.0;
    phang=0;
    double A=0.,Aa=0.;
    double B=0.,Ba=0.;
    double s2=sintl2(hkl)*0.25;
    //printf("s2 %f %f\n",s2,resd);
    if ((s2< mol->hklShellLow)||(s2> mol->hklShellHig)){
      phang=0.0;
      //printf("XXX\n");
      I=-666.0;
      return 0.0;
    }
    if (s2>0.001)resd=fmax(s2,resd);
    QList<double> f0;
    f000=0.0;
    for (int i=0; i<isfac.size(); i++) f0.append(atomformfactor(isfac.at(i),s2));
    for (int i=0; i<xyz.size(); i++){
      double sof = sofs.at(i);
      int an= isfac.at(type.at(i));
      f000+=sof*(an+1);
      double p=static_cast<double>(twopi*(xyz.at(i)*hkl));
      double fp=fprime(an,radiation);
      double fpp=fdprime(an,radiation);
      //T = 8*(pi**2)*Uiso*sin(theta/lambda)**2    
      if (isos.at(i)) {
          sof*=exp(bet.at(i).m11*s2);
          //printf("beta =%f\n", bet.at(i).m11);
      }
      else {
        sof*=exp(-twopi2* ((hkl*bet.at(i))*hkl));
      }//  */
      double cs=cos(p),si=sin(p);
      //double cs=1,si=0.0;
      double a=(f0.at(type.at(i))+fp);
      A+= sof*(a*cs-fpp*si);
      B+= sof*(a*si+fpp*cs);
      Aa+=sof*(a*cs);
      Ba+=sof*(a*si);
    }
    I=A*A+B*B;
    double Ia=Aa*Aa+Ba*Ba;
    phang=atan2(Ba,Aa)*180.0/M_PI;
    phang=fmod(720.01+phang,360)-0.01;
    double fc=sqrt(Ia);
    if (s2>0.001) Fmax=fmax(Fmax,fc);
    return fc;
  }

  double F(int h, int k, int l,int radiation,double &phang,double &disp){
    const double twopi=M_PI*2.0;
    const double twopi2=M_PI*M_PI*2.0;
    phang=0;
    double A=0,Aa=0;
    f000=0.0;
    double B=0,Ba=0;
    Matrix N=Matrix(mol->cell.as,0,0,0,mol->cell.bs,0,0,0,mol->cell.cs),beta_,beta;//,ucc;
    double s2=sintl2(h,k,l)*0.25;
    //printf("Fs2 %f %f\n",s2,resd);
    if (s2>0.001)resd=fmin(0.5/sqrt(s2),resd);
    int ns = mol->cell.symmops.size();
    int na = mol->asymm.size();
    for (int s=0; s<ns; s++){
      for (int i=0; i<na; i++){
        int an = mol->asymm.at(i).an;
        if (an<0) continue;
        bool iso=(mol->asymm.at(i).isIso);
        if (!iso){
          beta=(N*mol->asymm.at(i).uf*N);
          mol->Usym(beta,mol->cell.symmops.at(s),beta_);
        }
        double sof = mol->asymm.at(i).sof;
        f000+=sof*(an+1);
        V3 uvw = mol->cell.symmops.at(s)*mol->asymm[i].frac+mol->cell.trans.at(s);

        V3 hkl(h,k,l);
        double p=twopi*uvw*hkl;
        double f0=atomformfactor(an,s2);
        double fp=fprime(an,radiation);
        double fpp=fdprime(an,radiation);
        //T = 8*(pi**2)*Uiso*sin(theta/lambda)**2
        if (iso) sof*=exp(-4*twopi2*mol->asymm.at(i).uf.m11*s2);
        else {
          sof*=exp(-twopi2* ((hkl*beta_)*hkl));

        }
        double cs=cos(p),si=sin(p);
        A+=sof*((f0+fp)*cs-fpp*si);
        B+=sof*((f0+fp)*si+fpp*cs);
        Aa+=sof*((f0+fp)*cs);
        Ba+=sof*((f0+fp)*si);
      }
    }
    double I=A*A+B*B;
    double Ia=Aa*Aa+Ba*Ba;
    phang=atan2(Ba,Aa)*180.0/M_PI;
    phang=fmod(720.01+phang,360)-0.01;
    disp=I;
    double fc=sqrt(Ia);
    if (s2>0.001) Fmax=fmax(Fmax,fc);
    return fc;
  }
  /*double F(int h, int k, int l,int radiation,double &phang,double &disp){
    const double twopi=M_PI*2.0;
    const double twopi2=M_PI*M_PI*2.0;
    phang=0;
    double A=0,Aa=0;
    f000=0.0;
    Matrix N=Matrix(mol->cell.as,0,0,0,mol->cell.bs,0,0,0,mol->cell.cs),beta_,beta,ucc;
    double s2=sintl2(h,k,l)*0.25;
    if (s2>0.001)resd=fmin(0.5/sqrt(s2),resd);
    int ns = mol->cell.symmops.size();
    int na = mol->asymm.size();
    double sofsum=0;
    double sofsum2=0;
    for (int s=0; s<ns; s++){
    for (int i=0; i<na; i++){
    int an = mol->asymm.at(i).an;
    if (an<0) continue;
    bool iso=(mol->asymm.at(i).isIso);
    if (!iso){
    beta=(N*mol->asymm.at(i).uf*N);
    mol->Usym(beta,mol->cell.symmops.at(s),beta_);
    Matrix uff;
    mol->Usym(mol->asymm.at(i).uf,mol->cell.symmops.at(s),uff);
    mol->Uf2Uo(uff,ucc);
    }
    double sof = mol->asymm.at(i).sof;
    f000+=sof*(an+1);
    sofsum+=sof;
    V3 uvw = mol->cell.symmops.at(s)*mol->asymm[i].frac+mol->cell.trans.at(s);
    V3 hkl(h,k,l);
    double p=twopi*uvw*hkl;
    double f0=atomformfactor(an,s2);
    double fp=fdprime(an,radiation);
    double fpp=fprime(an,radiation);
    if (iso) sof*=exp(-4*twopi2*mol->asymm.at(i).uf.m11*s2);
    else {
    sof*=exp(-twopi2* ((hkl*beta_)*hkl));
    }
    sofsum2+=sof;
    double cs=cos(p),si=sin(p);
    A+=sof*((f0+fp)*cs+fpp*si);
    B+=sof*((f0+fp)*si+fpp*cs);
    Aa+=sof*((f0+fp)*cs);
    Ba+=sof*((f0+fp)*si);
    }
    }
    double I=A*A+B*B;
    double Ia=Aa*Aa+Ba*Ba;
    phang=atan2(Ba,Aa)*180.0/M_PI;
    phang=fmod(720+phang,360);
    disp= (I-Ia)/Ia;
    if ((h==0)&&(k==0)&&(l==0)){ 
    if (I>0.001){
    printf("%4d%4d%4d %12.6f %12.6f %10.6f s2= %f s= %f d= %f ph%f %f %f%%\n",h,k,l,I,A,B,s2,sqrt(s2),0.5/sqrt(s2) ,phang,I,100*I/Ia);
    }
    else printf("%4d%4d%4d absent d=%f\n",h,k,l,0.5/sqrt(s2));

    }
    if ((h==0)&&(k==0)&&(l==0)) printf("%f %f\n",sofsum, sofsum2);
    disp=I;
    double fc=sqrt(Ia);
    if (s2>0.001) Fmax=fmax(Fmax,fc);
    return fc;
    }*/
  QString myfcf(QString fn){
    qDebug()<<fn;
    //  qDebug()<<mol->hklf;
    //  qDebug()<<mol->hklScale<<mol->hklSigmaScale;
    if ((mol->hklf<3)||(mol->hklf>4)) return "";
    //  const double twopi=M_PI*2.0;
    const double twopi2=M_PI*M_PI*2.0;
    //QTime zeit; zeit.start();

    Rec *lr=NULL,ob;
    QString base=fn, mfcf=fn, hklf=fn;
    hklf=hklf.replace(QRegExp("(.res$)|(.ins$)"),".hkl");
    mfcf=mfcf.replace(QRegExp("(.res$)|(.ins$)"),".fcf6");
    base=base.replace(QRegExp("(.res$)|(.ins$)"),"");
    base=base.section('/', -1);
    FILE *hkl=fopen(hklf.toStdString().c_str(),"rt");
    if (hkl==NULL) {fprintf(stderr,"can't open hkl\n"); return "";}
    int nr=0;
    char line[128];
    int wl=0;
    if (fabs(mol->cell.wave-0.71073)<0.01) wl=SCATTMO;
    else if (fabs(mol->cell.wave-1.5418)<0.01) wl=SCATTCU;
    else if (fabs(mol->cell.wave-0.56083)<0.01) wl=SCATTAG;
    else wl=-1;
    size_t nnr=0;
    mol->hklShellLow = pow(0.5/mol->hklShellLow,2);
    mol->hklShellHig = pow(0.5/mol->hklShellHig,2);
    //  qDebug()<<"@@"<<mol->hklShellLow<<mol->hklShellHig;
    while (!feof(hkl)){
      if (fgets(line,120,hkl)) {}
      nnr++;
    }
    rewind(hkl);
    lr=static_cast<Rec*>( malloc(sizeof(Rec)*nnr));
    FILE *fc6=fopen(mfcf.toStdString().c_str(),"wt");
    if (fc6==NULL) {
      fprintf(stderr,"can't open %s\n",mfcf.toStdString().c_str()); 
      qDebug()<<"can't open " << mfcf;
      return "";
    }

    //FILE *tst=NULL;//:qfopen("scale_test_123.456","wt");
    if (1) {
      Matrix N=Matrix(mol->cell.as,0,0,0,mol->cell.bs,0,0,0,mol->cell.cs),beta_,beta;//,ucc;
      type.clear();
      bet.clear();
      sofs.clear();
      xyz.clear();
      isfac.clear();
      isos.clear();
      int ns = mol->cell.symmops.size();
      int na = mol->asymm.size();
      for (int s=0; s<ns; s++){
        for (int i=0; i<na; i++){
          int an = mol->asymm.at(i).an;
          if (an<0) continue;
          if (!isfac.contains(an)) isfac.append(an);
          type.append(isfac.indexOf(an));
          bool iso=(mol->asymm.at(i).isIso);
          if (!iso){
            beta=(N*mol->asymm.at(i).uf*N);
            mol->Usym(beta,mol->cell.symmops.at(s),beta_);
          }else{beta_.m11=-4*twopi2*mol->asymm.at(i).uf.m11;}
          isos.append(iso);
          bet.append(beta_);
          double sof = mol->asymm.at(i).sof;
          sofs.append(sof);
          xyz.append(mol->cell.symmops.at(s)*mol->asymm[i].frac+mol->cell.trans.at(s));
          //printf("%s%d %d %f [%f %f %f] %f %f %f\n",mol->asymm.at(i).Label.toStdString().c_str(),s,isos.last(),sof,beta_.m11,beta_.m22,beta_.m33,beta_.m23,beta_.m13,beta_.m12);
        }
      }
    }




    //  printf("wave:%f %f %f %f %d %d wl=%d\n",mol->cell.wave,fabs(mol->cell.wave-0.71073),fabs(mol->cell.wave-1.5418),fabs(mol->cell.wave-0.56083), mol->cell.ns0, mol->cell.symmops.size(),wl);

    /*printf("hklf N=%d s:%g \n%g %g %g\n%g %g %g\n%g %g %g\nsm: %g \n",mol->hklf,mol->hklScale ,
        mol->hklMat.m11,mol->hklMat.m12,mol->hklMat.m13,
        mol->hklMat.m21,mol->hklMat.m22,mol->hklMat.m23,
        mol->hklMat.m31,mol->hklMat.m32,mol->hklMat.m33,
        mol->hklSigmaScale);// */
    while (!feof(hkl)){    
      if (fgets(line,120,hkl)) {}
      if (feof(hkl))continue;
      if (strlen(line)<28) continue;
      char chkl[13],cisig[17];
      for (int ci=0; ci<12; ci++) chkl[ci]=line[ci];
      chkl[12]='\0';
      for (int ci=0; ci<16; ci++) cisig[ci]=line[ci+12];
      cisig[16]='\0';
      sscanf(chkl,"%4d%4d%4d", &ob.ih, &ob.ik, &ob.il);
      sscanf(cisig,"%8lf%8lf", &ob.fo, &ob.si);
      if ((ob.ih==0)&&(ob.il==0)&&(ob.ik==0)) {
        //      printf("#[%d,%d,%d,%d]%s",nr,ob.ih,ob.ik,ob.il,line);

        break; 
      }
      V3 hkl=V3(ob.ih,ob.ik,ob.il)*mol->hklMat; 
      ob.ih=static_cast<int>(hkl.x);
      ob.ik=static_cast<int>(hkl.y);
      ob.il=static_cast<int>(hkl.z);
      ob.fo*=mol->hklScale;
      ob.si=(mol->hklScale*mol->hklSigmaScale*ob.si);
      if (mol->hklf==3){
        double xxx=fmax(0.01,ob.si);
        ob.si=2.0*xxx * fabs(fmax(fabs(ob.fo),xxx));
        ob.fo*=ob.fo;
      }
      //ob.fc=FF(hkl,wl,ob.phi,ob.dis);
      //    if ((ob.ih==1)&&(ob.ik==0)&&(ob.il==0)){ printf("%4d%4d%4d %18.9f %18.9f\n",ob.ih,ob.ik,ob.il,ob.fo,ob.si); }
      //   ob.fc=F(ob.ih,ob.ik,ob.il,wl,ob.phi,ob.dis);

      lr[nr]=ob;
      nr++;
    }
    //  qDebug()<<resd <<0.5/sqrt(resd);
    //resd=0.5f/sqrtf(resd);
    //printf("resd %f\n",resd);
    //  for (int iw=0; iw<xyz.size(); iw++){ qDebug()<<type.at(iw)<<sofs.at(iw)<<xyz.at(iw).x<<xyz.at(iw).y<<xyz.at(iw).z; }
    //if (tst!=NULL) fclose(tst);

    double myosf=0.0;
    for (int i=0; i<nr; i++){
      V3 uvw=V3(lr[i].ih,lr[i].ik,lr[i].il);
      V3 m=uvw;
      double p,q=lr[i].phi;
      lr[i].phi=fmod(720.0+q,360.0);
      for (int s=0; s<mol->cell.ns0;s++){
        double t=1.0;
        V3 nhkl=uvw*mol->cell.symmops.at(s);
        //      if((nl<0)||((nl==0)&&(nk<0))||((nl==0)&&(nk==0)&&(nh<0)))
        //if((nl<0)||((nl==0)&&(nk<0))||((nl==0)&&(nk==0)&&(nh<0)))
        if((nhkl.z<0)||((nhkl.z==0)&&(nhkl.y<0))||((nhkl.z==0)&&(nhkl.y==0)&&(nhkl.x<0))) {
          nhkl.x*=-1;
          nhkl.y*=-1;
          nhkl.z*=-1;
          t=-1.0;
        }
        if ((nhkl.z<m.z)||((nhkl.z==m.z)&&(nhkl.y<m.y))||((nhkl.z==m.z)&&(nhkl.y==m.y)&&(nhkl.x<=m.x))) continue;
        m=nhkl;
        //p=u*sy[9][k]+v*sy[10][k]+w*sy[11][k];
        p=uvw*mol->cell.trans.at(s);
        lr[i].phi=fmod(719.99+t*fmod(q-360.0*p,360.0),360.0)+0.01;
        /*      if (m==V3(1,1,0)) {printf("?%4d%4d%4d sym=%d %g \n %4g%4g%4g fo%12.6f si%12.6f fc%12.6f phi%12.3f %d %f %f (%g %g %g)\n",lr[i].ih,lr[i].ik,lr[i].il,s,t,
                nhkl.x,nhkl.y,nhkl.z,
                lr[i].fo,lr[i].si,lr[i].fc,lr[i].phi,i,q,p, mol->cell.trans.at(s).x,mol->cell.trans.at(s).y,mol->cell.trans.at(s).z);}*/
      }
      lr[i].ih=static_cast<int>(m.x);
      lr[i].ik=static_cast<int>(m.y);
      lr[i].il=static_cast<int>(m.z);
    }
    sorthkl(nr,lr);
    int n=-1;
    //  double sw=0.,swxy=0.0,swx=0.0,swy=0.0,swx2=0.0,sfo=0,sfoc=0;
    double swxy=0.0,swx2=0.0,sfo=0,sfoc=0;
    {int i=0;
      while(i<nr){
        double t=0.;
        double u=0.;
        double v=0.;
        double v2=0.;
        double ws=0;
        int m;
        int k=i;
        while ((i<nr)&&(lr[i].ih==lr[k].ih)&&(lr[i].ik==lr[k].ik)&&(lr[i].il==lr[k].il)) {
          t=t+1.;
          double w=fmax(4.0,lr[i].fo/lr[i].si)/lr[i].si;//1.0;//
          ws+=w; 
          u+=lr[i].fo*w;//*lr[i].fc*lr[i].fc/lr[i].dis;
          v+=1./(lr[i].si*lr[i].si);
          v2+=lr[i].si;
          // p=lr[i].phi;
          // ddi+=lr[i].dis;
          //         if (t>1) {printf("!%4d%4d%4d fo%12.6f sig%12.6f fc%12.6f phi%12.3f i%d  k%dw=%f dis%f\n",lr[i].ih,lr[i].ik,lr[i].il,lr[i].fo,lr[i].si,lr[i].fc,lr[i].phi,i,k,w,lr[i].dis);}
          i++;
        }
        m=n+1;
        double yy=(u/ws);
        lr[m].fo = sqrt(fmax(0.0,yy));
        lr[m].so=sqrt(lr[m].fo*lr[m].fo+sqrt(1./v))-lr[m].fo;
        lr[m].si=1.0/sqrt(v);
        n=m;
        lr[n].ih=lr[k].ih;
        lr[n].ik=lr[k].ik;
        lr[n].il=lr[k].il;
      }
      //printf("%d %d %d \n",nr,n,i);
    }
    n++;
    nr=n;
    for (int i=0; i<nr; i++){
      V3 chkl=V3(lr[i].ih,lr[i].ik,lr[i].il);
      lr[i].fc=0.0;
      lr[i].fc=FF(chkl,wl,lr[i].phi,lr[i].dis);
      lr[i].fo =lr[i].fo * (lr[i].fc*lr[i].fc/lr[i].dis);//correct for dispersion
      if (lr[i].dis<=-665.0) {
        lr[i].fc=0.0;
        continue;
      }
      if ((mol->exti!=-666.0)&&(mol->swat==-666.0)){
        lr[i].fc*=pow(1 + 0.001 * mol->exti * lr[i].fc * lr[i].fc  *mol->cell.wave * mol->cell.wave /sqrt(sintl2(chkl)),-0.25);
      }
      double xx=lr[i].fc*lr[i].fc;
      double yy=lr[i].fo*lr[i].fo;
      double v=1.0/(lr[i].si*lr[i].si);
      if (xx>0.001){
        swx2+=v*xx*xx;
        swxy+=v*yy*xx;
      }
    }
    double osf2=(mol->osf<=0.0)?swxy/swx2:mol->osf*mol->osf;//myosf;
    myosf=sqrt(osf2);
    // qDebug() << mol->osf << myosf << osf2 << swxy/swx2<<swxy <<swx2;
    //printf("%g %g %g\n",myosf,osf2,swxy/swx2);
    type.clear();
    bet.clear();
    sofs.clear();
    xyz.clear();
    isos.clear();  
    //  myosf=mol->osf;
    resd=0.5/sqrt(resd);
    //printf("???? %f \n",resd);
    fprintf(fc6,
        "#\n# h,k,l, Fo-squared, sigma(Fo-squared), Fc and phi(calc)\n#\ndata_%s\n_shelx_title '%s'\n_shelx_refln_list_code          6\n"
        "_shelx_F_calc_maximum      %6.2f\n_exptl_crystal_F_000       %6.2f\n_reflns_d_resolution_high  %6.4f\n\nloop_\n_space_group_symop_operation_xyz\n"
        ,base.toStdString().c_str()
        ,mol->titl.toStdString().c_str()
        ,Fmax
        ,f000
        ,resd
        );
    for (int s=0; s<mol->cell.symmops.size();s++) 
      fprintf(fc6,"'%s'\n",mol->symmcode2human(s).trimmed().toStdString().c_str());
    fprintf(fc6,
        "\n_cell_length_a   %8.4f\n_cell_length_b   %8.4f\n_cell_length_c   %8.4f\n_cell_angle_alpha %7.3f\n_cell_angle_beta  %7.3f\n_cell_angle_gamma %7.3f\n\n"
        "loop_\n_refln_index_h\n_refln_index_k\n_refln_index_l\n_refln_F_squared_meas\n_refln_F_squared_sigma\n_refln_F_calc\n_refln_phase_calc\n",
        mol->cell.a,mol->cell.b,mol->cell.c,mol->cell.al,mol->cell.be,mol->cell.ga
        );
    int igut=0;
    for (int i=0; i<nr; i++){
      /*
         fprintf(fc6,"%d %d %d %.2f %.2f %2f %.2f %.1f %f\n"
         ,lr[i].ih
         ,lr[i].ik
         ,lr[i].il
      //        ,(pow(lr[i].fo/mol->osf,2.0))/(1.0+lr[i].dis)
      ,pow(lr[i].fo/mol->osf,2.0)
      ,lr[i].si/(mol->osf*mol->osf)
      ,sqrt(lr[i].dis)
      ,lr[i].fc
      ,lr[i].phi
      ,lr[i].fc*lr[i].fc
      // ,lr[i].dis
      );// */

      double iosig = pow(lr[i].fo/myosf,2.0)/(lr[i].si/(osf2));
      if ((lr[i].fc>0.001)&& (iosig>-1.0)){
        igut++;
        sfo+=lr[i].fo/myosf;
        sfoc+=fabs(lr[i].fo/myosf-lr[i].fc);
      }//else{printf("notgut? %d %d %d %g %g\n",lr[i].ih,lr[i].ik,lr[i].il,iosig,lr[i].fc);}
      if (lr[i].fc>0.001)
        fprintf(fc6,"%d %d %d %.2f %.2f %.2f %.1f\n"
            ,lr[i].ih
            ,lr[i].ik
            ,lr[i].il
            //        ,(pow(lr[i].fo/mol->osf,2.0))/(1.0+lr[i].dis)
            ,pow(lr[i].fo/myosf,2.0)
            ,lr[i].si/(osf2)
            //,sqrt(lr[i].dis)
            ,lr[i].fc
            ,lr[i].phi
            // ,lr[i].fc*lr[i].fc
            // ,lr[i].dis
            );//  */
    }
    fclose(hkl);
    fprintf(fc6,"\n");
    fclose(fc6);
    printf("R1=%f %d %d %g %g %f\n",sfoc/sfo,igut,nr,sfoc,sfo,f000);

    //printf("=>=>%d %d\n",zeit.elapsed(),xyz.size());
    return mfcf;
  }
  private:
  Molecule *mol;
  double Fmax,f000,resd;
  QList<Matrix> bet;
  QList<V3> xyz;
  QList<double> sofs;
  QList<bool> isos;
  QList<int> type,isfac;

  double sintl2(int h,int k, int l){
    V3 H=V3(h,k,l);
    double s=H*mol->cell.Gi*H;
    return s;
  }
  double sintl2(V3 H){
    double s=H*mol->cell.Gi*H;
    return s;
  }
  void sorthkl(int nr, Rec r[]){
    /*! sorts the reflection list
    */
    Rec *hilf=static_cast<Rec*> (malloc(sizeof(Rec)*nr));
    if (hilf==NULL)return ; 
    int i,j,k,nj,ni,spalte;int index[4096];
    for (spalte=0; spalte<3; spalte++){
      j=-999999;
      k=999999;
      switch (spalte) {
        case 0: for (i=0; i<nr; i++){ j=(j<r[i].ih)?r[i].ih:j; k=(k>r[i].ih)?r[i].ih:k; } break;
        case 1: for (i=0; i<nr; i++){ j=(j<r[i].ik)?r[i].ik:j; k=(k>r[i].ik)?r[i].ik:k; } break;
        case 2: for (i=0; i<nr; i++){ j=(j<r[i].il)?r[i].il:j; k=(k>r[i].il)?r[i].il:k; } break;
      }
      nj=-k;
      ni=(nj+j+1);
      for (i=0; i<=ni; i++) index[i]=0;
      for (i=0; i<nr; i++){
        switch (spalte){
          case 0: j=r[i].ih+nj; break;
          case 1: j=r[i].ik+nj; break;
          case 2: j=r[i].il+nj; break;
        }
        index[j]++;/*brauch ich das? -->JA!*/
        hilf[i].ih=r[i].ih;
        hilf[i].ik=r[i].ik;
        hilf[i].il=r[i].il;
        hilf[i].fo=r[i].fo;
        hilf[i].so=r[i].so;
        hilf[i].si=r[i].si;
        hilf[i].fc=r[i].fc;
        hilf[i].dis=r[i].dis;
        hilf[i].phi=r[i].phi;
      }/*/4*/
      j=0;
      for (i=0; i<ni; i++){
        k=j;
        j+=index[i];
        index[i]=k;     
      }/*/5*/
      for (i=0; i<nr;i++){
        switch (spalte) {
          case 0: j=hilf[i].ih +nj;break;
          case 1: j=hilf[i].ik +nj;break;
          case 2: j=hilf[i].il +nj;break;
        }     
        index[j]++;   
        j=index[j]-1;
        r[j].ih=hilf[i].ih;
        r[j].ik=hilf[i].ik;
        r[j].il=hilf[i].il;
        r[j].fo=hilf[i].fo;
        r[j].so=hilf[i].so;
        r[j].si=hilf[i].si;
        r[j].fc=hilf[i].fc;
        r[j].dis=hilf[i].dis;
        r[j].phi=hilf[i].phi;
      }/*/6*/
    }/*/spalten*/
    free(hilf);
  }
};
#endif
