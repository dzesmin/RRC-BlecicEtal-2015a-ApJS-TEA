import numpy as np

##################################################################################################################
  ######################################  CEA-TEA comparison free energies CEA 
##################################################################################################################

species = np.array(['H', 'C', 'N', 'O', 'H2', 'CO', 'CH4', 'H2O', 'N2', 'NH3']) 

# read TEA theory document, section 7.1
coeff_order = np.array(['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'b1', 'b2'])

# CEA T range 200 - 1000 K, taken from thermo.inp
coeff = \
np.array([[0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 2.547370801e+04,  -4.466828530e-01], \
[6.495031470e+02, -9.649010860e-01, 2.504675479e+00, -1.281448025e-05, 1.980133654e-08, -1.606144025e-11, 5.314483411e-15, 8.545763110e+04, 4.747924288e+00],\
[0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 5.610463780e+04, 4.193905036e+00], \
[-7.953611300e+03, 1.607177787e+02, 1.966226438e+00, 1.013670310e-03, -1.110415423e-06, 6.517507500e-10, -1.584779251e-13, 2.840362437e+04, 8.404241820e+00], \
[4.078323210e+04, -8.009186040e+02, 8.214702010e+00, -1.269714457e-02, 1.753605076e-05, -1.202860270e-08, 3.368093490e-12, 2.682484665e+03, -3.043788844e+01],\
[1.489045326e+04, -2.922285939e+02, 5.724527170e+00, -8.176235030e-03, 1.456903469e-05, -1.087746302e-08, 3.027941827e-12, -1.303131878e+04, -7.859241350e+00], \
[-1.766850998e+05, 2.786181020e+03, -1.202577850e+01, 3.917619290e-02, -3.619054430e-05, 2.026853043e-08, -4.976705490e-12, -2.331314360e+04, 8.904322750e+01], \
[-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03, -7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04, 1.724205775e+01], \
[2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03, 1.384646189e-05, -9.625793620e-09, 2.519705809e-12, 7.108460860e+02, -1.076003744e+01], \
[-7.681226150e+04, 1.270951578e+03, -3.893229130e+00, 2.145988418e-02, -2.183766703e-05, 1.317385706e-08, -3.332322060e-12, -1.264886413e+04, 4.366014588e+01]]) 

# divide temperature range
T = np.linspace(100.0, 900.0, 9)

# allocate array
g_RT_CEA = np.empty([len(T), len(species)])

# read TEA theory document and refer to CEA theory document
# S_R = -a1*T**(-2)/2.0 - a2*T**(-1)+ a3*np.log(T) + a4*T + a5*(T**2)/2.0 + a6*(T**3)/3.0 + a7*(T**4)/4.0 + b2
# H_RT = -a1*T**(-2) + a2*T**(-1)*np.log(T)+ a3 + a4*T/2.0 + a5*(T**2)/3.0 + a6*(T**3)/4.0 + a7*(T**4)/5.0 + b1/T
# g_RT_CEA = -S_R + H_RT

for i in np.arange(len(T)):
    for j in np.arange(len(species)):
            S_R = -coeff[j,0]*T[i]**(-2)/2.0 - coeff[j,1]*T[i]**(-1)+ coeff[j,2]*np.log(T[i]) + coeff[j,3]*T[i] + coeff[j,4]*(T[i]**2)/2.0 + coeff[j,5]*(T[i]**3)/3.0 + coeff[j,6]*(T[i]**4)/4.0 + coeff[j,8]
            H_RT = -coeff[j,0]*T[i]**(-2) + coeff[j,1]*T[i]**(-1)*np.log(T[i])+ coeff[j,2] + coeff[j,3]*T[i]/2.0 + coeff[j,4]*(T[i]**2)/3.0 + coeff[j,5]*(T[i]**3)/4.0 + coeff[j,6]*(T[i]**4)/5.0 + coeff[j,7]/T[i]
            g_RT_CEA[i, j] = -S_R + H_RT
print species
print g_RT_CEA


# CEA T range 1000 - 6000 K, taken from thermo.inp
coeff2 = \
np.array([[ 6.078774250e+01, -1.819354417e-01, 2.500211817e+00, -1.226512864e-07, 3.732876330e-11, -5.687744560e-15, 3.410210197e-19, 2.547486398e+04, -4.481917770e-01], \
[-1.289136472e+05, 1.719528572e+02, 2.646044387e+00, -3.353068950e-04, 1.742092740e-07, -2.902817829e-11, 1.642182385e-15, 8.410597850e+04, 4.130047418e+00], \
[ 8.876501380e+04, -1.071231500e+02, 2.362188287e+00, 2.916720081e-04, -1.729515100e-07, 4.012657880e-11, -2.677227571e-15, 5.697351330e+04, 4.865231506e+00], \
[2.619020262e+05, -7.298722030e+02, 3.317177270e+00, -4.281334360e-04, 1.036104594e-07, -9.438304330e-12, 2.725038297e-16, 3.392428060e+04, -6.679585350e-01], \
[5.608128010e+05, -8.371504740e+02, 2.975364532e+00, 1.252249124e-03, -3.740716190e-07, 5.936625200e-11, -3.606994100e-15, 5.339824410e+03, -2.202774769e+00], \
[4.619197250e+05, -1.944704863e+03, 5.916714180e+00, -5.664282830e-04, 1.398814540e-07, -1.787680361e-11, 9.620935570e-16, -2.466261084e+03, -1.387413108e+01], \
[3.730042760e+06, -1.383501485e+04, 2.049107091e+01, -1.961974759e-03, 4.727313040e-07, -3.728814690e-11, 1.623737207e-15, 7.532066910e+04, -1.219124889e+02], \
[1.034972096e+06, -2.412698562e+03, 4.646110780e+00, 2.291998307e-03, -6.836830480e-07, 9.426468930e-11, -4.822380530e-15, -1.384286509e+04, -7.978148510e+00], \
[5.877124060e+05, -2.239249073e+03, 6.066949220e+00, -6.139685500e-04, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 1.283210415e+04, -1.586640027e+01], \
[2.452389535e+06, -8.040894240e+03, 1.271346201e+01, -3.980186580e-04, 3.552502750e-08, 2.530923570e-12, -3.322700530e-16, 4.386191960e+04, -6.462330602e+01]])

# divide temperature range
T = np.linspace(1000.0, 3000.0, 21)

# allocate array
g_RT_CEA2 = np.empty([len(T), len(species)])

for i in np.arange(len(T)):
    for j in np.arange(len(species)):
            S_R = -coeff2[j,0]*T[i]**(-2)/2.0 - coeff2[j,1]*T[i]**(-1)+ coeff2[j,2]*np.log(T[i]) + coeff2[j,3]*T[i] + coeff2[j,4]*(T[i]**2)/2.0 + coeff2[j,5]*(T[i]**3)/3.0 + coeff2[j,6]*(T[i]**4)/4.0 + coeff2[j,8]
            H_RT = -coeff2[j,0]*T[i]**(-2) + coeff2[j,1]*T[i]**(-1)*np.log(T[i])+ coeff2[j,2] + coeff2[j,3]*T[i]/2.0 + coeff2[j,4]*(T[i]**2)/3.0 + coeff2[j,5]*(T[i]**3)/4.0 + coeff2[j,6]*(T[i]**4)/5.0 + coeff2[j,7]/T[i]
            g_RT_CEA2[i, j] = -S_R + H_RT
print species
print g_RT_CEA2

