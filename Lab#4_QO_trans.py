import numpy as np
import scipy.optimize as opt
from scipy import interpolate
import matplotlib.pyplot as plt
from peakdetect import peakdetect

h_bar = 1.055e-34       # Planck constant (J*s)
e_vacc = 8.854e-12      # Permittivity in vaccum (F/m)
e_SiO2 = 3.9 * e_vacc   # Dielectric const for SiO2
q_e = 1.60218e-19       # Electron charge (C)


def f(x, a, b):
    return a * x + b


# Prepare parameters in the loop


file_list = ('T_33A_LHe.csv', 'T_30A_LHe.csv', 'T_27A_LHe.csv', 'T_25A_LHe.csv', 'T_LN2.csv', 'T_RoomTemp.csv')
label_list = ('LHe, 33A', 'LHe, 30A', 'LHe, 27A', 'LHe, 25A', 'LN2', 'Room Temp')
color_list = ((0, 0, 1), (0, 0.3, 1), (0, 0.5, 1), (0, 0.4, 0), (1, 0.4, 0), (1, 0, 0))


I_list = [33, 30, 27, 25]
LIdxStart_list = [0, 1, 1, 1]

V_T = []
errV_T = []
V_T_def2 = []
errV_T_def2 = []
t_ox = []
V_G = [[], [], [], []]

for i in [0,1,2,3,5]:

# Import transconductance data
    # Select data in one cycle

    Data = np.loadtxt(file_list[i], delimiter=',', skiprows=16, usecols=(1, 3))

    minVGidx = Data[:, 0].argmin()
    Data = Data[minVGidx::]
    minVGidx = Data[::-1, 0].argmin()
    Data = Data[:(len(Data) - minVGidx):]

# Averaging data with same V_gate value
    idx = np.argsort(Data[:, 0])
    sorted = Data[idx, 0]
    V1, idx_start, count = np.unique(sorted, return_counts=True, return_index=True)
    res = np.split(idx, idx_start[1:])
    V2 = np.ones(len(res))
    for (i_uniq_value, j) in zip(res, range(0, len(res))):
        unique_value = np.mean(Data[i_uniq_value, 1])
        V2[j] = unique_value
    Data = np.transpose([V1, V2])

    N = len(Data)
    Idx = np.linspace(0, N, N)

    x = Idx
    y = Data[:, 1]

    plt.figure('Original')
    plt.plot(Data[:, 0], y, '--', markersize=3, color=color_list[i], label=label_list[i])
    plt.xlabel(r'$V_G$(V)')
    plt.ylabel(r'dG/dV')
    plt.legend()
plt.show()
'''
# Find local minimum between 2 peaks using peakdetect.py
    # Delete other minimum outside 2 peaks
    if i == 1:
        continue
    (locMax, locMin) = peakdetect(y, lookahead=3)
    locMax = np.array(locMax)
    locMin = np.array(locMin)

    sortedMaxIdx = np.argsort(locMax[:, 1])
    idxStart = locMax[sortedMaxIdx[-1], 0]

# Find V_T(def2)
    ThV_def2 = Data[int(idxStart), 0]
    errThV_def2 = 0.2

    V_T_def2 = np.append(V_T_def2, ThV_def2)
    errV_T_def2 = np.append(errV_T_def2, errThV_def2)

    condition = (locMin[:, 0] > idxStart)
    locMin = locMin[condition]
    #plt.figure()
    #plt.plot(locMin[:, 0], locMin[:, 1], 'r+')
    #plt.plot(locMax[:, 0], locMax[:, 1], 'b+')

# Smooth signal using interpolate spline
    # Set smoothing factor to 0.3

    x2 = x[::1]
    y2 = y[::1]
    spl = interpolate.UnivariateSpline(x2, y2)
    spl.set_smoothing_factor(0.05)
    ySmooth = spl(x2)

    plt.figure('Transconductance Curve')
    plt.plot(Data[:, 0], ySmooth, 'k--', linewidth=1)
    plt.plot(Data[:, 0], y, '.', markersize=5, color=color_list[i], label=label_list[i])
    plt.xlabel(r'$V_G$(V)')
    plt.ylabel(r'dG/dV')
    plt.legend()


    plt.plot([Data[int(x), 0] for x in locMin[:, 0]], locMin[:, 1], 'r+', markersize=10)

# Find corresponding gate voltage using index
    # Deal with asymmetric peakdetect result

    idxGV = [int(x) for x in locMin[:, 0]]

    gateV = np.sort(Data[idxGV, 0])
    errGV = 0.2
    print(gateV)

    V_G[i] = np.append(V_G[i], gateV)
    print('\nGate Voltage is {}'.format(gateV))


    LandauIdx = np.arange(LIdxStart_list[i], len(gateV) + LIdxStart_list[i])

# Fit gate voltage points to a straight line
    # and plot four lines

    LIdx_range = np.arange(-1, 10)
    popt, pcov = opt.curve_fit(f, LandauIdx, gateV)
    perr = np.sqrt(np.diag(pcov))
    print('Fitting param {}, param error {}'.format(popt, perr))

    plt.figure('QO')
    plt.plot(LIdx_range, f(LIdx_range, *popt), color=color_list[i],
             label=label_list[i])
    plt.errorbar(LandauIdx, gateV, yerr=errGV, fmt='k.', ecolor=color_list[i], capsize=5)
    plt.axis((-1, 8, 0, 25))
    plt.grid('on')
    plt.xlabel(r'Landau Indices')
    plt.ylabel(r'$V_G$(V)')
    plt.legend(loc='lower right')

# Calculate Threshold voltage V_T
    # V_T = f(-1/2)

    ThV = f(-1/2, popt[0], popt[1])
    errThV = perr[0]/2 + 0.2
    V_T = np.append(V_T, ThV)
    errV_T = np.append(errV_T, errThV)

# Calculate thickness t_ox
    # t_ox = slope * pi * h_bar * epsilon / e^2 / B

    B = 1737 * I_list[i] * 1e-4
    print(B)
    tox = popt[0] * np.pi * h_bar * e_SiO2 / (q_e**2 * B)
    t_ox = np.append(t_ox, tox)


print('\nThreshold voltage', V_T)
print('Error in threshold voltage', errV_T)
print('\nThreshold voltage(def2)', V_T_def2)
print('Error in threshold voltage(def2)', errV_T_def2)

sd_VT = np.std(V_T)
V_T = np.mean(V_T)
errV_T = np.mean(errV_T)
V_T_def2 = np.mean(V_T_def2)
errV_T_def2 = np.mean(errV_T_def2)

print('\nThickness is', t_ox)
sd_t_ox = np.std(t_ox)
t_ox = np.mean(t_ox)

print('\nThreshold voltage is {:5.3f} +/- {:5.3f} V'.format(V_T, sd_VT+errV_T))
print('\nThreshold voltage(def2) is {:5.3f} +/- {:5.3f} V'.format(V_T_def2, errV_T_def2))
print('Thickness is {:5.3f} +/- {:5.3f} nm'.format(t_ox * 1e9, sd_t_ox * 1e9))  # Convert m to nm


# Save gate voltage to a file
transfile = open('transfile_v2.csv', 'w')

for i in range(0, 4):
    for j in V_G[i]:
        transfile.write(str(j)+',')
    transfile.write('\n')
transfile.write('\n')
# Calculate number of electrons in all energy levels per unit N
    # and plot N(V_gate)

for i in (0,2,3):
    N = [e_SiO2 * (x - V_T) / (q_e * t_ox) * 1e-4 for x in V_G[i]]    # Convert m^-2 to cm^-2
    for j in N:
        transfile.write(str(j)+',')

    transfile.write('\n')
    #plt.figure('Electron Density N')
    #plt.plot(V_G[i], N, 'o', markersize=5, color=color_list[i], label=label_list[i])
    #plt.xlabel(r'$V_G$(V)')
    #plt.ylabel(r'2DEG density($cm^{-2}$)')
    #plt.legend(loc='lower right')

transfile.close()

plt.show()
'''