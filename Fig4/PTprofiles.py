#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# Kevin's pressure array, the same for tp_crp.dat and tp_orp.dat
# WASP12b, Stevenson et al 2014
pres_WASP12b = np.array([  1.00000000e+02,   8.49753000e+01,   7.22081000e+01,
         6.13591000e+01,   5.21401000e+01,   4.43062000e+01,
         3.76494000e+01,   3.19927000e+01,   2.71859000e+01,
         2.31013000e+01,   1.96304000e+01,   1.66810000e+01,
         1.41747000e+01,   1.20450000e+01,   1.02353000e+01,
         8.69749000e+00,   7.39072000e+00,   6.28029000e+00,
         5.33670000e+00,   4.53488000e+00,   3.85353000e+00,
         3.27455000e+00,   2.78256000e+00,   2.36449000e+00,
         2.00923000e+00,   1.70735000e+00,   1.45083000e+00,
         1.23285000e+00,   1.04762000e+00,   8.90215000e-01,
         7.56463000e-01,   6.42807000e-01,   5.46228000e-01,
         4.64159000e-01,   3.94421000e-01,   3.35160000e-01,
         2.84804000e-01,   2.42013000e-01,   2.05651000e-01,
         1.74753000e-01,   1.48497000e-01,   1.26186000e-01,
         1.07227000e-01,   9.11163000e-02,   7.74264000e-02,
         6.57933000e-02,   5.59081000e-02,   4.75081000e-02,
         4.03702000e-02,   3.43047000e-02,   2.91505000e-02,
         2.47708000e-02,   2.10490000e-02,   1.78865000e-02,
         1.51991000e-02,   1.29155000e-02,   1.09750000e-02,
         9.32603000e-03,   7.92483000e-03,   6.73415000e-03,
         5.72237000e-03,   4.86260000e-03,   4.13201000e-03,
         3.51119000e-03,   2.98365000e-03,   2.53536000e-03,
         2.15443000e-03,   1.83074000e-03,   1.55568000e-03,
         1.32194000e-03,   1.12332000e-03,   9.54548000e-04,
         8.11131000e-04,   6.89261000e-04,   5.85702000e-04,
         4.97702000e-04,   4.22924000e-04,   3.59381000e-04,
         3.05386000e-04,   2.59502000e-04,   2.20513000e-04,
         1.87382000e-04,   1.59228000e-04,   1.35305000e-04,
         1.14976000e-04,   9.77010000e-05,   8.30218000e-05,
         7.05480000e-05,   5.99484000e-05,   5.09414000e-05,
         4.32876000e-05,   3.67838000e-05,   3.12572000e-05,
         2.65609000e-05,   2.25702000e-05,   1.91791000e-05,
         1.62975000e-05,   1.38489000e-05,   1.17681000e-05,
         1.00000000e-05])


# Kevin tp_orp.dat C/O = 0.5 solar
# WASP12b, Stevenson et al 2014
temp_WASP12b_solar = np.array([ 2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.26,  2947.68,
        2946.6 ,  2945.04,  2943.03,  2940.59,  2937.75,  2934.53,
        2930.96,  2927.08,  2922.89,  2918.48,  2914.35,  2910.51,
        2906.96,  2903.69,  2900.71,  2898.01,  2895.59,  2891.16,
        2883.82,  2873.64,  2860.69,  2845.05,  2826.77,  2805.92,
        2782.58,  2756.81,  2728.67,  2698.24,  2667.89,  2638.57,
        2610.28,  2583.01,  2556.78,  2531.57,  2507.39,  2484.25,
        2462.13,  2441.04,  2420.97,  2401.94,  2383.94,  2366.96,
        2351.02,  2336.1 ,  2322.21,  2309.35,  2297.52,  2286.71,
        2276.94,  2268.2 ,  2260.48,  2253.79,  2248.13,  2238.36,
        2234.76,  2232.19,  2230.64,  2230.13])


# Kevin tp_crp.dat C/O=1.2
# WASP12b, Stevenson et al 2014
temp_WASP12b_swapped = np.array([ 2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.51,  2963.85,  2962.85,  2961.53,  2959.96,
        2957.23,  2949.72,  2937.53,  2920.75,  2899.48,  2873.8 ,
        2844.07,  2810.51,  2773.19,  2732.13,  2687.4 ,  2639.98,
        2593.58,  2548.2 ,  2503.84,  2460.5 ,  2418.18,  2376.87,
        2336.59,  2297.33,  2259.09,  2221.86,  2185.66,  2150.48,
        2116.31,  2083.17,  2051.05,  2019.94,  1989.86,  1960.79,
        1932.75,  1905.72,  1879.72,  1854.73,  1830.77,  1807.82,
        1785.9 ,  1764.99,  1745.1 ,  1726.24,  1708.39,  1691.56,
        1675.76,  1660.97,  1647.2 ,  1634.46,  1622.73,  1612.02,
        1602.33,  1593.66,  1586.01,  1579.39,  1573.78,  1564.09,
        1560.52,  1557.97,  1556.44,  1555.93])

# Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt
# WASP43b, Stevenson et al 2014
pres_WASP43b = np.array([  3.16227770e+01,   2.63026800e+01,   2.18776160e+01,
         1.81970090e+01,   1.51356120e+01,   1.25892540e+01,
         1.04712850e+01,   8.70963590e+00,   7.24435960e+00,
         6.02559590e+00,   5.01187230e+00,   4.16869380e+00,
         3.46736850e+00,   2.88403150e+00,   2.39883290e+00,
         1.99526230e+00,   1.65958690e+00,   1.38038430e+00,
         1.14815360e+00,   9.54992590e-01,   7.94328230e-01,
         6.60693450e-01,   5.49540870e-01,   4.57088190e-01,
         3.80189400e-01,   3.16227770e-01,   2.63026800e-01,
         2.18776160e-01,   1.81970090e-01,   1.51356120e-01,
         1.25892540e-01,   1.04712850e-01,   8.70963590e-02,
         7.24435960e-02,   6.02559590e-02,   5.01187230e-02,
         4.16869380e-02,   3.46736850e-02,   2.88403150e-02,
         2.39883290e-02,   1.99526230e-02,   1.65958690e-02,
         1.38038430e-02,   1.14815360e-02,   9.54992590e-03,
         7.94328230e-03,   6.60693450e-03,   5.49540870e-03,
         4.57088190e-03,   3.80189400e-03,   3.16227770e-03,
         2.63026800e-03,   2.18776160e-03,   1.81970090e-03,
         1.51356120e-03,   1.25892540e-03,   1.04712850e-03,
         8.70963590e-04,   7.24435960e-04,   6.02559590e-04,
         5.01187230e-04,   4.16869380e-04,   3.46736850e-04,
         2.88403150e-04,   2.39883290e-04,   1.99526230e-04,
         1.65958690e-04,   1.38038430e-04,   1.14815360e-04,
         9.54992590e-05,   7.94328230e-05,   6.60693450e-05,
         5.49540870e-05,   4.57088190e-05,   3.80189400e-05,
         3.16227770e-05,   2.63026800e-05,   2.18776160e-05,
         1.81970090e-05,   1.51356120e-05,   1.25892540e-05,
         1.04712850e-05, 1e-5])

# Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt
# WASP43b, Stevenson et al 2014
temp_WASP43b = np.array([ 1811.8938 ,  1810.9444 ,  1810.1535 ,  1809.4948 ,  1808.9463 ,
        1808.4898 ,  1808.1098 ,  1807.7936 ,  1807.5304 ,  1807.3114 ,
        1807.1291 ,  1806.9766 ,  1806.8464 ,  1806.7212 ,  1806.53   ,
        1806.1269 ,  1805.2849 ,  1803.7403 ,  1800.5841 ,  1794.9518 ,
        1786.6255 ,  1775.3705 ,  1761.0973 ,  1742.3631 ,  1719.6396 ,
        1694.0976 ,  1666.1517 ,  1636.2055 ,  1603.0265 ,  1567.6227 ,
        1531.3326 ,  1494.5529 ,  1457.6432 ,  1419.5923 ,  1381.4921 ,
        1344.4864 ,  1308.8483 ,  1274.7949 ,  1241.6997 ,  1210.3302 ,
        1181.2851 ,  1154.5844 ,  1130.2066 ,  1107.7735 ,  1087.5396 ,
        1069.5885 ,  1053.7529 ,  1039.8606 ,  1027.6493 ,  1017.057  ,
        1007.9573 ,  1000.1681 ,   993.52384,   987.85759,   983.05871,
         979.01311,   975.60886,   972.74937,   970.3489 ,   968.33983,
         966.66146,   965.26054,   964.09216,   963.11826,   962.30738,
         961.63264,   961.0714 ,   960.60473,   960.2169 ,   959.89468,
         959.627  ,   959.40467,   959.22003,   959.06677,   958.93955,
         958.83394,   958.74627,   958.6735 ,   958.61313,   958.56303,
         958.52146,   958.48696, 958.48])


# WASP-12b C/O = 0.5
plt.figure(1, figsize=(2, 3))
plt.semilogy(temp_WASP12b_solar, pres_WASP12b, color='k', linestyle='-', linewidth=2)
plt.xlim(2000, 3100)   
plt.xticks(np.arange(2000, 3100, 500))  
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Pressure [bar]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(max(pres_WASP12b), min(pres_WASP12b))
plt.savefig('WASP12b-Kevin-solar-PT.png')
plt.savefig('WASP12b-Kevin-solar-PT.ps')

# WASP-12b C/O swapped
plt.figure(2, figsize=(2, 3))
plt.semilogy(temp_WASP12b_swapped, pres_WASP12b, color='orange', linestyle='--', linewidth=3)
plt.xlim(1400, 3100)   
plt.xticks(np.arange(1500, 3100, 500))  
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Pressure [bar]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(max(pres_WASP12b), min(pres_WASP12b))
plt.savefig('WASP12b-Kevin-CO1.2-PT.png')
plt.savefig('WASP12b-Kevin-CO1.2-PT.ps')

# WASP-43b all metallicities
plt.figure(3, figsize=(2, 3))
plt.semilogy(temp_WASP43b, pres_WASP43b, color='r', linestyle='-', linewidth=2)
plt.xlim(700, 2000)   
plt.xticks(np.arange(1000, 2001, 500))   
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Pressure [bar]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(max(pres_WASP43b), min(pres_WASP43b))
plt.savefig('WASP43b-Kevin.png')
plt.savefig('WASP43b-Kevin-PT.ps')


