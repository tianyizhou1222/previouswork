import numpy as np
import xlrd
from scipy import interpolate
import scipy.optimize as opt
import matplotlib.pyplot as plt


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


e_vacc = 8.854e-12      # Permittivity in vaccum (F/m)
e_SiO2 = 3.9 * e_vacc   # Dielectric const for SiO2
q_e = 1.60218e-19       # Electron charge (C)
m_e = 9.10938e-31       # Electron mass (kg)
m_eff = 0.19 * m_e      # Effective mass (kg)

L = 1e-5                # Length of substrate (m)
W = 5.72e-4             # Width of substrate (m)

t_ox = 281.589e-9        # Thickness of dioxide (m)
V_T = 4.185             # Threshold voltage (V)

# Prepare parameters in the loop

file_list = ('C_33A_LHe.csv', 'C_30A_LHe.csv', 'C_27A_LHe.csv', 'C_25A_LHe.csv', 'C_LN2.csv', 'C_RoomTemp.csv')
label_list = ('LHe, 33A', 'LHe, 30A', 'LHe, 27A', 'LHe, 25A', 'LN2', 'Room Temp')
color_list = ((0, 0, 1), (0, 0.3, 1), (0, 0.5, 1), (0, 0.4, 0), (1, 0.4, 0), (1, 0, 0))

smooth1_list = [0.2, 0.5, 0.5, 0.5, 0.005, 0.001]

book1 = xlrd.open_workbook('transfile.xlsx')
sheet1 = book1.sheet_by_index(0)

fig, ax1 = plt.subplots()

for i in range(0, 6):

# Import Calibration and Conductance data

    Cali_data = np.loadtxt('Calibration.csv', delimiter=',', skiprows=10, usecols=(1, i+2))
    Data = np.loadtxt(file_list[i], delimiter=',', skiprows=16, usecols=(1, 3))
    VarG = Cali_data[::-1, 0]
    MeasV = Cali_data[::-1, 1]

# Plot Calibration curve: condunctance v.s. V_gate

    SatVm = 11

    spl = interpolate.UnivariateSpline(MeasV, VarG)
    Vm_space = np.linspace(0.2, 11.2, 1000)
    spl.set_smoothing_factor(0.00000001)
    G_spline = spl(Vm_space)

    for Vm in Vm_space:
        if Vm > SatVm:
            idx_sat = list(Vm_space).index(Vm)
            break
    G_sat = G_spline[idx_sat - 1]
    G_spline[idx_sat:] = [G_sat] * len(G_spline[idx_sat:])

    plt.figure('Calibration')
    plt.plot(Vm_space, G_spline, color=color_list[i], label=label_list[i])
    plt.plot(MeasV, VarG, 'k.', markersize=2)
    plt.xlabel(r'Measured V(V)')
    plt.ylabel(r'Conductance($\Omega^{-1}$)')


# Plot original data points

    plt.figure('Scope data')
    plt.plot(Data[:, 1], Data[:, 0], '.', color=color_list[i],  markersize=1)
    plt.xlabel(r'Measured V(V)')
    plt.ylabel(r'$V_G$(V)')

# Plot averaged curve
    # Select data with mono-increasing V_measure
    #  by deleting starting zero values
    #  and cut off by saturation V_measure
    # Then using spline interpolation to get smooth curve

    idx = np.argsort(Data[:, 0])
    sorted = Data[idx, 0]
    V1, idx_start, count = np.unique(sorted, return_counts=True, return_index=True)
    res = np.split(idx, idx_start[1:])
    V2 = np.ones(len(res))
    for (i_uniq_value, j) in zip(res, range(0, len(res))):
        unique_value = np.mean(Data[i_uniq_value, 1])
        V2[j] = unique_value

    print('V1 is', V1, '\nV2 is', V2)

    V1V2 = np.transpose([V1, V2])
    condition = ((V1V2[:, 1] - V1V2[0, 1]) > 0.05)
    V1V2 = V1V2[condition]

    idx = np.argsort(V1V2[:, 1])
    sorted = V1V2[idx, 1]
    V2, idx_start, count = np.unique(sorted, return_counts=True, return_index=True)
    res = np.split(idx, idx_start[1:])
    V1 = np.ones(len(res))
    for (i_uniq_value, j) in zip(res, range(0, len(res))):
        unique_value = np.mean(V1V2[i_uniq_value, 0])
        V1[j] = unique_value
    V1V2 = np.transpose([V1, V2])

    V1V2lt10 = V1V2[np.where(V1V2[:, 1] < 10)]
    V1V2gt10 = V1V2[np.where(V1V2[:, 1] > 10)]
    print(V1V2gt10)
    V1V2 = np.append(V1V2lt10, V1V2gt10[::6], axis=0)
    print(V1V2)

    spl = interpolate.UnivariateSpline(V1V2[:, 1], V1V2[:, 0])
    spl.set_smoothing_factor(smooth1_list[i])
    VG_spline = spl(Vm_space)

# Plot V_gate v.s. V_measure

    plt.figure('Averaged scope curve')
    plt.plot(Vm_space, VG_spline, color=color_list[i], label=label_list[i])
    plt.plot(V1V2[:, 1], V1V2[:, 0], 'k+', markersize=5)
    plt.xlabel(r'Measured V(V)')
    plt.ylabel(r'$V_G$(V)')
    plt.legend()



# Plot Conductance Curve: conductance v.s. V_gate
    plt.figure('VGspline-Gspline')
    plt.plot(VG_spline[:idx_sat:], G_spline[:idx_sat:], color=color_list[i], label=label_list[i])
    plt.xlabel(r'$V_G$(V)')
    plt.ylabel(r'Conductance($\Omega^{-1}$)')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0, ymax=0.005)
    plt.legend()

    if i == 0:
        continue


    VG_spline = np.delete(VG_spline, list(range(idx_sat-1, idx_sat+9)))
    G_spline = np.delete(G_spline, list(range(idx_sat-1, idx_sat+9)))

    VG = np.append(VG_spline[:980:12], VG_spline[984::4])
    G = np.append(G_spline[:980:12], G_spline[984::4])

    VG_append = np.linspace(VG[-1]+2, 25, int((25-VG_spline[idx_sat-1])/2))
    G_append = [G_sat] * len(VG_append)

    VG = np.append(VG, VG_append)
    G = np.append(G, G_append)

    print(VG)
    VG_space = np.linspace(min(VG), max(VG), 1000)
    spl = interpolate.UnivariateSpline(VG, G)
    spl.set_smoothing_factor(0)
    Cond = spl(VG_space)

    plt.figure('Conductance Curve')
    plt.plot(VG_space, Cond, color=color_list[i], label=label_list[i])
    plt.plot(VG, G, 'k.', markersize=2)
    plt.xlabel(r'$V_G$(V)')
    plt.ylabel(r'Conductance($\Omega^{-1}$)')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0, ymax=0.006)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()

# Plot Conductivity v.s. VG

    sigma = L / W * Cond
    plt.figure('Conductivity')
    plt.plot(VG_space, sigma, color=color_list[i], label=label_list[i])
    plt.xlabel(r'$V_G$(V)')
    plt.ylabel(r'Single Layer Conductivity(S)')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()

    if i == 0:
        continue
    if i == 4:
        V_T = 2.76
    if i == 5:
        V_T = 1.92

# Calculate N
    N = [e_SiO2 * (x - V_T) / (q_e * t_ox) * 1e-4 for x in VG_space]    # Convert m^-2 to cm^-2


# Plot mobility v.s. VG


    #VG_discrete = sheet1.row_values(i)
    #N_discrete = sheet1.row_values(i+4)
    #VG_discrete = list(filter(None, VG_discrete))
    #N_discrete = list(filter(None, N_discrete))
    #print(VG_discrete)
    #nearest_idx = [find_nearest(VG_space, x) for x in VG_discrete]

    #miu_discrete = [sigma[x1] / q_e / x2 for (x1, x2) in zip(nearest_idx, N_discrete)]  # in unit cm^2/(Vs)
    miu = [x1 / q_e / x2 for (x1, x2) in zip(sigma, N)]
    #fig, ax1 = plt.subplots()
    plt.grid('on')
    #plt.figure('Mobility')
    #plt.plot(VG_discrete, miu_discrete, 'D', color=color_list[i])
    ax1.plot(VG_space, miu, color=color_list[i], label=label_list[i])
    ax1.set_xlabel(r'$V_G$(V)')
    ax1.set_ylabel(r'2DEG mobility($cm^2$/(Vs)')
    ax1.set_xlim([0, 26])
    ax1.set_ylim([0, 1.1e3])
    #ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1.grid('on')
    ax1.legend()

# Plot mean free time tau v.s. VG
    #tau_discrete = [m_eff * x * 1e-4 / q_e for x in miu_discrete]  # Convert miu back to m^2/(Vs)
    tau = [m_eff * x * 1e-4 / q_e for x in miu]                    # Convert s to ns
    tau_lim = m_eff * 1.1e6 * 1e-4 * 1e9 / q_e
    #plt.figure('Mean free time')
    #plt.plot(VG_discrete, tau_discrete, 'D', color=color_list[i])
    ax2 = ax1.twinx()
    ax2.plot(VG_space, [x*1e9 for x in tau], color=color_list[i], label=label_list[i])
    ax2.set_ylabel(r'2DEG mean free time(ns)')
    ax2.set_ylim([0, tau_lim])
    #ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.legend()

# Plot mobility v.s. density
    plt.figure('Mobility v.s. Density')
    plt.grid('on')
    #plt.plot(N_discrete, miu_discrete, 'D', color=color_list[i])
    plt.plot(N, miu, color=color_list[i], label=label_list[i])
    plt.xlabel(r'2DEG density($cm^{-2}$)')
    plt.ylabel(r'2DEG mobility($cm^2$/(Vs)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()

    ax3 = ax1.twiny()
    N_lim_min = e_SiO2 * (0 - V_T) / (q_e * t_ox) * 1e-4 * 1e-12      # gate limit_min set to be 0V
    N_lim_max = e_SiO2 * (26 - V_T) / (q_e * t_ox) * 1e-4 * 1e-12      # gate limit_max set to be 26V
    ax3.plot([x*1e-12 for x in N], miu, color=color_list[i], label=label_list[i])
    ax3.set_xlabel(r'2DEG density($\times 10^{12} cm^{-2}$)')
    ax3.set_xlim([N_lim_min, N_lim_max])
    #ax3.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

plt.show()
