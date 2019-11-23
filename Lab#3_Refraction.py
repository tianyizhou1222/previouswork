from math import *
import xlrd
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit

# Import data from xlsx file

book1 = xlrd.open_workbook('Refraction2.xlsx')
sheet1 = book1.sheet_by_index(0)

book2 = xlrd.open_workbook('Data2.xlsx')
sheet2 = book2.sheet_by_index(0)

lambda_legend_list = sheet1.col_values(0, 3, 8)
lambda_list = [650, 600, 550, 500, 450]
color_list = ((1, 0, 0), (1, 0.75, 0), (0.64, 1, 0), (0, 1, 0.57), (0, 0.27, 1))

# Calculate index of refraction
delta_m = sheet1.col_values(4, 3, 8)
sd_angle = 1.45e-4          # Error of angle measured in rad
n_air = 1.00028             # Index of refraction of air in STP
n_prism = [sin(x/2+pi/8)*n_air/sin(pi/8) for x in delta_m]
sd_n_prism = [n_air*cos(x/2+pi/8)/2/sin(pi/8)*sd_angle for x in delta_m]
print(r'Refractive index is', n_prism)
print(r'Refractive index error is', sd_n_prism)


# Calculate critical angle and needed angle of divided circle


def r2dm(x):    # Rad to Degree Minute transform
    return [int(x/pi*180), (x/pi*180-int(x/pi*180))*60]


alpha_c = [asin(n_air/x) for x in n_prism]
alpha_c_dm = [r2dm(x) for x in alpha_c]
sd_alpha_c = [n_air/x/sqrt(x**2-n_air**2)*0.0002*180/pi for x in n_prism]


def cal_theta(x1, x2):
    return asin(x1/n_air*sin(pi/4-x2))


theta_theory = [cal_theta(x1, x2) for (x1, x2) in zip(n_prism, alpha_c)]
theta_theory_dm = [r2dm(x) for x in theta_theory]
print(theta_theory_dm)


def cal_exp_alpha(t, n):
    if t > 0:
        return pi / 4 - asin(n_air * sin(t) / n)
    else:
        return pi / 4 + asin(n_air * sin(t) / n)


# Plot dispersion relation


c = 3e8                   # Speed of light (m/s)
h_bar = 6.58212e-16       # Planck constant (eV*s)
omega_p_theory = 8.9      # eV/h_bar

epsilon_d = 1.00059     # Dielectric const for air
epsilon_0 = 8.85419e-12 # Vacuum permittivity (F/m)
q_e = 1.60218e-19       # Electron charge (C)
m_eff = 9.10938e-31     # Effective mass of electron in metal (assume = electron mass)


omega = [2*pi*c/x/10**-9 for x in lambda_list]

exp_theta_r = sheet1.col_values(11, 3, 8)
exp_alpha = [cal_exp_alpha(t, n) for (t, n) in zip(exp_theta_r, n_prism)]
exp_alpha_dm = [r2dm(x) for x in exp_alpha]
sd_exp_alpha = [sqrt((sin(t)**2*0.0002**2/n**2/(n**2-sin(t)**2))+
                     (cos(t)**2*sd_angle**2/(n**2-sin(t)**2)))
                for (t, n) in zip(exp_theta_r, n_prism)]
print(r'Critical angle is', alpha_c_dm)
print(r'Critical angle error (in minute) is', sd_alpha_c)
print(r'Resonance angle is', exp_alpha_dm)
print(r'Resonance angle error is', sd_exp_alpha)


k_space = np.linspace(0, 4e7)
omega_space = np.linspace(0, 9e15)
omega_air = k_space*c/n_air
omega_prism = [k_space*c/x for x in n_prism]
k_air = [x/c*n_air for x in omega]
k_sp = [x1/c*x2*sin(x3) for (x1, x2, x3) in zip(omega, n_prism, exp_alpha)]
k_theory = [i/c*sqrt((i**2-(omega_p_theory/h_bar)**2)/
                     (2*i**2-(omega_p_theory/h_bar)**2)) for i in omega_space]
omega_asp = omega_p_theory/h_bar/sqrt(2)

sd_k_sp = [x1/c*sqrt((sin(x3)*0.0002)**2+(cos(x3)*x2*1e-4)**2)
           for (x1, x2, x3) in zip(omega, n_prism, exp_alpha)]
print(r'k_sp error is', sd_k_sp)

# Curve-fit of dispersion relation


# def disp(x, a):
#    return np.sqrt((2*(c*x)**2 - np.sqrt(4*(c*x)**4 + a**4) + a**2) / 2 )


# def disp(x, a):
#    return a * np.sqrt((1 - (x*c)**2) / (1 - 2 * (x*c)**2))

def disp(x, a):
    return a * np.sqrt(1 + (c**2 / ((1 / x**2) - 2*c**2)))


popt, pcov = curve_fit(disp, np.array(k_sp), np.array(omega))
perr = np.sqrt(np.diag(pcov))
print(popt, perr)


plt.figure(0)
plt.plot(disp(k_space, popt))

plt.figure(1)
plt.plot(k_space, omega_air, 'k')
for i in range(0, 5):
    plt.plot(k_space, omega_prism[i], color=color_list[i])
plt.plot(k_theory, omega_space, 'r')
plt.plot(k_sp, omega, 'k+')
# plt.errorbar(x=k_sp, y=omega, xerr=np.array(sd_k_sp))
plt.axhline(y=omega_asp, color='r', linestyle='--')
plt.xlabel(r'k(m^-1)')
plt.ylabel(r'$\omega$(s^-1)')


# Calculate electron density of silver


def cal_n(npr, alpha, omega):
    return (1 / (npr**2 * sin(alpha)**2 - 1) + 2) * omega**2 \
           * (epsilon_0 * m_eff / q_e**2) / 1e6


density = [cal_n(x1, x2, x3) for (x1, x2, x3) in zip(n_prism, exp_alpha, omega)]
print(r'Electron density of silver is', density)


# Plot absorption spectrum near resonance
for i in range(0,3):
    # Import raw data from book2.sheet2
    angle_raw = sheet2.col_values(0+5*i, 4, 17)
    TM_raw = sheet2.col_values(3+5*i, 4, 17)
    TE_raw = sheet2.col_values(4+5*i, 4, 17)

    # Prepare data for plot
    # 12669 is calibration angle in minutes
    # Calculate angle alpha with function cal_exp_alpha
    # Remove background signal (32.9 counts) and calibration (5 counts)
    #   from raw data
    theta = [(12669 - x - 45 * 60) / 60 / 180 * pi for x in angle_raw]      # in rad
    alpha = [cal_exp_alpha(t, n_prism[0]) / pi * 180 for t in theta]    # in degree
    TM = [x - 32.9 for x in TM_raw]
    TE = [x - 32.9 for x in TE_raw]
    Intensity = [x1 / x2 for (x1, x2) in zip(TM, TE)]

    # Fitting to data points using spline interpolation
    # Set smoothing factor to 0
    spl = interpolate.UnivariateSpline(alpha, Intensity)
    alpha_space = np.linspace(min(alpha), max(alpha), 1000)
    spl.set_smoothing_factor(0)
    Intensity_spline = spl(alpha_space)

    plt.figure(2)
    #plt.plot(alpha, Intensity, color=color_list[i], linewidth=0.9)
    plt.plot(alpha_space, Intensity_spline, color=color_list[i], linewidth=0.9)
    plt.axvline(x=alpha_c[i]/pi*180, color=color_list[i], linestyle='--', linewidth=0.7)
    plt.plot(alpha, Intensity, 'k.', markersize=2)

plt.xlabel(r'Angle of incidence $\alpha$(degrees)')
plt.ylabel(r'Relative intensity I(TM)/I(TE)')
# plt.legend(lambda_legend_list, loc='best')
plt.show()

