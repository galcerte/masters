#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.optimize import fsolve
from scipy.constants import pi
from scipy import integrate

# This script has the contents of:
# - sec23.py (lepton flavor violation graphs in section 2.3)
# - sec51.py (5D active neutrino mass matrix convergence, sum, effective couplings and mass splitting with mode in section 5.1)
# - ScalarMassSpectrum.py (normalized mass splitting among eta components in section 4.2, difference among wavefunctions)
# in order to clean up the code used in this thesis and have it all be in one place. Graphs from
# Chapter 2 have not been included yet. This script produces all graphs found in Chapter 4 and
# Chapter 5.

# We will use natural units as usual, ħ = c = 1

# The curvature and the value of krc are taken assuming gravitationally interacting
# dark matter (DM), something which requires a rather high DM mass (> 1 TeV, around 4-11 TeV)
# and given LHC bounds on the KK-graviton, a smoothing of the hierarchy problem is implied
# rather than a complete solution, so Lambda ~ 10 TeV.

# AdS curvature (GeV)
k = 4.01577*(10**17)
# k times r_c (dimensionless)
krc = 32.7961/pi
# Compactification radius (GeV⁻¹)
rc = krc/k
# Planck scale (GeV)
MPl = 1.22*(10**19)
# Reduced Planck scale (GeV)
MPlr = MPl/np.sqrt(8*pi)
# Fundamental scale of gravity (GeV)
M5 = np.cbrt( (k*(MPlr**2))/(1 - np.exp(-2*krc*pi)) )
# Here we initialize where Yukawa interaction parameters will be stored
h = np.ones((3, 3))
# Scalar potential parameters (adimensional)
lambda3 = 1
lambda4 = 1
lambda5 = 0.5
# Bulk mass m2 of eta scalars (GeV)
m2 = 3*(10**17)
# SM Higgs VEV (GeV)
vev = 246.21965079413732
# Order of Bessel functions
v = np.sqrt((m2/k)**2 + 4)
# Masses of fermionic singlets, in GeV
M = np.array([10000, 15000, 20000])

# We first define the corrections with respect to the (brane mass)-less terms,
warp = np.exp(-krc*pi)
correta = ((lambda3*(vev**2))/(4*k*M5))*np.exp(2*krc*pi)
corrr = (((lambda3 + lambda4 + lambda5)*(vev**2))/(4*k*M5))*np.exp(2*krc*pi)
corri = (((lambda3 + lambda4 - lambda5)*(vev**2))/(4*k*M5))*np.exp(2*krc*pi)

# We wish to find the roots of the expressions below,
# which'll tell us the mass spectrums of the bulk scalars we consider.
# Define the expression whose roots we want to find:
def masseta(x):
    return ((2 + correta)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + correta)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    # Here we leave the simplified expression as well.
    #return (2 + correta)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def massr(x):
    return ((2 + corrr)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + corrr)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    #return (2 + corrr)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def massi(x):
    return ((2 + corri)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + corri)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    #return (2 + corri)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def mass2(x):
    return (2*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - (2*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
def massless(x):
    return (2*sp.jv(2, x) + x*sp.jvp(2, x, 1))*(2*sp.yv(2, warp*x) + warp*x*sp.yvp(2, warp*x, 1)) - (2*sp.yv(2, x) + x*sp.yvp(2, x, 1))*(2*sp.jv(2, warp*x) + warp*x*sp.jvp(2, warp*x, 1))

# Charged eta boson
# Here we have our guesses for the roots we want to find,
x_init_guesseta = [4.66874849, 7.72492277, 10.7, 13.88, 16.988, 20.103, 23.22, 26.35, 29.5, 32.61, 35.7,
                   38.9, 42.01, 45.1, 48.28, 51.42, 54.56, 57.7, 60.84, 64, 67.11, 70.25, 73, 76]
# We'll use fsolve(), a numerical solver, to find the roots
x_rootseta = fsolve(masseta, x_init_guesseta)

# Real, neutral eta boson
x_init_guessr = [4.97, 8.07, 11.13, 14.19, 17.27, 20.35, 23.45, 26.55, 29.66, 32.78, 35.9, 39.02,
                 42.1, 45.28, 48.4, 51.42, 54.67, 57.8, 60.84, 64, 67.11, 70.25, 73, 76]
x_rootsr = fsolve(massr, x_init_guessr)

# Imaginary, neutral eta boson
x_init_guessi = [4.97, 8.07, 11.13, 14.19, 17.27, 20.35, 23.45, 26.55, 29.66, 32.78, 35.9, 39.02,
                 42.1, 45.28, 48.4, 51.42, 54.67, 57.8, 60.84, 64, 67.11, 70.25, 73, 76]
x_rootsi = fsolve(massi, x_init_guessi)

x_init_guess = [3.96089237, 7.16961637, 10.33728283, 13.49278243, 16.64302665, 19.79052262, 22.93640086, 26.08124643, 29.22539253, 32.36904292, 35.51232917, 38.65534017, 41.79813792, 44.94076717, 48.08326097, 51.22564422, 54.36793611, 57.51015159, 60, 63, 66, 69, 72, 76]
x_roots2 = fsolve(mass2, x_init_guess)
x_rootsmassless = fsolve(massless, x_init_guess)
num_roots = np.shape(x_init_guess)[0]

# Turns out some results can be bogus because the root solver above
# is not finding the right root, it's solving the ones farther away,
# like 79 instead of 72-73 and such.

# Masses of the modes (GeV)
maetan       = k*x_rootseta*warp
marn         = k*x_rootsr*warp
main         = k*x_rootsi*warp
m2n          = k*x_roots2*warp
masslessspec = k*x_rootsmassless*warp

# Why is the spectrum of m_{\eta,n} and m_{I,n} identical???
# Oh silly me. if all lambdas are 1, the correction is identical between both.

# Now that we know the mass spectrum of the eta doublet, it is time to evaluate
# their wavefunctions at the IR brane.
ind = np.linspace(0, num_roots - 2, num_roots - 1, dtype=int)
# Mass difference among modes
mnsp = np.zeros(num_roots - 1)
# Coefficients. If x_n << the warp factor, then we can take bnv = 0.
betanv = - (2*sp.jv(v, maetan/k) + (maetan/k)*sp.jvp(v, maetan/k, 1))/(2*sp.yv(v, maetan/k) + (maetan/k)*sp.yvp(v, maetan/k, 1))
brnv   = - (2*sp.jv(v, marn/k)   + (marn/k)*sp.jvp(v, marn/k, 1))/(2*sp.yv(v, marn/k)       + (marn/k)*sp.yvp(v, marn/k, 1))
binv   = - (2*sp.jv(v, main/k)   + (main/k)*sp.jvp(v, main/k, 1))/(2*sp.yv(v, main/k)       + (main/k)*sp.yvp(v, main/k, 1))

# Normalization constants
netan = np.zeros(num_roots)
nrn   = np.zeros(num_roots)
nin   = np.zeros(num_roots)

index = np.linspace(0, num_roots - 1, num_roots, dtype=int)

# Initializing wavefunctions and where mass matrix'll go
fetan   = np.zeros(num_roots)
frn     = np.zeros(num_roots)
fin     = np.zeros(num_roots)
massmat = np.zeros((num_roots, num_roots))

for i in index:
    # We first determine the normalization constant for each eta field,
    def wavefuncetan(phi):
        return np.exp(2*krc*phi)*((sp.jv(v, (maetan[i]/k)*np.exp(krc*phi)) + betanv[i]*sp.yv(v, (maetan[i]/k)*np.exp(krc*phi)))**2)
        #return np.exp(2*krc*np.abs(phi))*(sp.jv(v, (mn[i]/k)*np.exp(krc*np.abs(phi)))**2)
    def wavefuncrn(phi):
        return np.exp(2*krc*phi)*((sp.jv(v, (marn[i]/k)*np.exp(krc*np.abs(phi))) + brnv[i]*sp.yv(v, (marn[i]/k)*np.exp(krc*np.abs(phi))))**2)
    def wavefuncin(phi):
        return np.exp(2*krc*np.abs(phi))*((sp.jv(v, (main[i]/k)*np.exp(krc*np.abs(phi))) + binv[i]*sp.yv(v, (main[i]/k)*np.exp(krc*np.abs(phi))))**2)
    netan[i] = np.sqrt(2*integrate.quad(wavefuncetan, 0, pi)[0])
    nrn[i] = np.sqrt(2*integrate.quad(wavefuncrn, 0, pi)[0])
    nin[i] = np.sqrt(2*integrate.quad(wavefuncin, 0, pi)[0])
    
    fetan[i] = (1/netan[i])*(sp.jv(v, (maetan[i]/k)*np.exp(krc*pi)) + betanv[i]*sp.yv(v, (maetan[i]/k)*np.exp(krc*pi)))
    frn[i] = (1/nrn[i])*(sp.jv(v, (marn[i]/k)*np.exp(krc*pi)) + brnv[i]*sp.yv(v, (marn[i]/k)*np.exp(krc*pi)))
    fin[i] = (1/nin[i])*(sp.jv(v, (main[i]/k)*np.exp(krc*pi)) + binv[i]*sp.yv(v, (main[i]/k)*np.exp(krc*pi)))
    #fetan[i] = (1/netan[i])*(sp.jv(v, (metan[i]/k)*np.exp(krc*pi)))

modes = np.linspace(0, num_roots - 1, num_roots, dtype=int)

plt.subplots_adjust(left = 0.15, right = 0.96, top = 0.95, bottom = 0.11)
plt.plot(modes, np.abs(maetan - marn)/(maetan + marn), '.', label=r'$|m_{\eta,n} - m_{R,n}|/(m_{\eta,n} + m_{R,n})$')
plt.plot(modes, np.abs(maetan - main)/(maetan + main), '.', label=r'$|m_{\eta,n} - m_{I,n}|/(m_{\eta,n} + m_{R,n})$')
plt.plot(modes, np.abs(marn - main)/(marn + main), '.', label=r'$|m_{R,n} - m_{I,n}|/(m_{R,n} + m_{I,n})$')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$\frac{|m_i - m_j|}{m_i + m_j}$', fontsize=20)
plt.legend()
plt.grid(True)
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=13)
plt.savefig('Figure41a.png', dpi = 300, pad_inches = 1)
plt.clf()

plt.subplots_adjust(left = 0.13, right = 0.96, top = 0.95, bottom = 0.11)
plt.plot(modes, np.abs(m2n - marn)/(m2n + marn), '.', label=r'$|m_{2,n} - m_{R,n}|/(m_{2,n} + m_{R,n})$')
plt.plot(modes, np.abs(m2n - main)/(m2n + main), '.', label=r'$|m_{2,n} - m_{I,n}|/(m_{2,n} + m_{I,n})$')
plt.plot(modes, np.abs(m2n - maetan)/(m2n + maetan), '.', label=r'$|m_{2,n} - m_{\eta,n}|/(m_{2,n} + m_{\eta,n})$')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$\frac{|m_i - m_j|}{m_i + m_j}$', fontsize=20)
plt.legend()
plt.grid(True)
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=13)
plt.savefig('Figure41b.png', dpi = 300, pad_inches = 1)
plt.clf()

# Needed for correct spacing of plots to image border
plt.subplots_adjust(left = 0.16, right = 0.98, top = 0.95, bottom = 0.11)

plt.plot(modes, (fin - frn)/(10**(-17)), '.')
# Axis labels
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$f^{(n)}_{\eta^0_I} - f^{(n)}_{\eta^0_R}$ at $\varphi = \pi$ ($\cdot 10^{-17}$)', fontsize=13)
# Render grid or not
plt.grid(True)
# First array sets the value at which ticks should be put and the second one sets their labels
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=12)
# Save the resulting figure as an image
plt.savefig('Figure52.png', dpi = 300, pad_inches = 1)
# Clear the figure, in order to allow other plots to be made
plt.clf()

modes = np.linspace(1, num_roots, num_roots, dtype=int)
modes0 = np.linspace(0, num_roots - 1, num_roots, dtype=int)
#plt.plot(modes, - sumand, '.', label="Mass matrix sumand")
#plt.plot(modes, sumand2, '.')
#print(sumand)

index1 = np.linspace(0, 2, 3, dtype=int)
Nwf = np.sqrt(krc/(2*(np.exp(krc*pi) - 1)))
warp32 = np.exp(1.5*krc*pi)
M5rc = M5*rc
pi162 = 16*(pi**2)
s = 0
heffr = []
heffi = []
for i in modes0:
    heffr.append(np.zeros((3, 3)))
    heffi.append(np.zeros((3, 3)))
    for j in index1:
        for k in index1:
            heffr[i][j, k] = (h[j, k]*warp32)/(M5rc)*frn[i]*Nwf
            heffi[i][j, k] = (h[j, k]*warp32)/(M5rc)*fin[i]*Nwf
Cijs = np.zeros(num_roots)
for n in modes0:
    Cijs[n] = heffr[n][0, 0]*heffr[n][0, 0]*((marn[n]**2)/(marn[n]**2 - M[0]**2))*np.log(marn[n]**2/M[0]**2) - heffi[n][0, 0]*heffi[n][0, 0]*((main[n]**2)/(main[n]**2 - M[0]**2))*np.log(main[n]**2/M[0]**2)
print(((marn**2)/(marn**2 - M[0]**2))*np.log(marn**2/M[0]**2))
print(((main**2)/(main**2 - M[0]**2))*np.log(main**2/M[0]**2))
print((((warp32)/(M5rc)*frn*Nwf)**2)*((marn**2)/(marn**2 - M[0]**2))*np.log(marn**2/M[0]**2) - (((warp32)/(M5rc)*fin*Nwf)**2)*((main**2)/(main**2 - M[0]**2))*np.log(main**2/M[0]**2))
plt.plot(modes, np.abs(Cijs[0])/(modes**2), '.', label=r"$C_{ij,s}(1)/n^2$")
plt.plot(modes, np.abs(Cijs[0])/modes, '.', label=r"$C_{ij,s}(1)/n$")
plt.plot(modes, np.abs(Cijs), '.', label=r"$C_{ij,s}(n)$")
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'Sequence value', fontsize=14)
plt.legend()
plt.grid(True)
plt.xticks(modes, np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=13)
plt.savefig('Figure53.png', dpi = 300, pad_inches = 1)
plt.clf()

indexlambda = np.linspace(0, 4, 5, dtype=int)
lam5 = []
l5 = np.array([1, 0.9, 0.8, 0.7, 0.6])
for i in indexlambda:
    corrr = (((lambda3 + lambda4 + l5[i])*(vev**2))/(4*k*M5))*np.exp(2*krc*pi)
    corri = (((lambda3 + lambda4 - l5[i])*(vev**2))/(4*k*M5))*np.exp(2*krc*pi)
    ass = k*fsolve(massr, x_init_guessr)*np.exp(-krc*pi)
    dicks = k*fsolve(massi, x_init_guessi)*np.exp(-krc*pi)
    lam5.append(np.abs(ass**2 - dicks**2))

plt.plot(modes, lam5[0]/(vev**2), '.', label=r'$\lambda_5 = 1.0$')
plt.plot(modes, lam5[1]/(vev**2), '.', label=r'$\lambda_5 = 0.9$')
plt.plot(modes, lam5[2]/(vev**2), '.', label=r'$\lambda_5 = 0.8$')
plt.plot(modes, lam5[3]/(vev**2), '.', label=r'$\lambda_5 = 0.7$')
plt.plot(modes, lam5[4]/(vev**2), '.', label=r'$\lambda_5 = 0.6$')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$|m^2_R - m^2_I|/v^2$', fontsize=13)
plt.grid(True)
plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0)
plt.xlim([-0.5, num_roots - 1 + 0.5])
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=12)
plt.savefig('Figure54a.png', dpi = 300, pad_inches = 1)
plt.clf()

plt.subplots_adjust(left = 0.15, right = 0.96, top = 0.95, bottom = 0.11)
plt.plot(modes, (np.abs(marn**2 - main**2))/(246**2), '.')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$|m^2_R - m^2_I|/v^2$', fontsize=13)
plt.grid(True)
plt.xlim([-0.5, num_roots - 1 + 0.5])
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=12)
plt.savefig('Figure54b.png', dpi = 300, pad_inches = 1)

plt.subplots_adjust(left = 0.14, right = 0.96, top = 0.95, bottom = 0.11)
M5rc = M5*rc
Nwf = np.sqrt(krc/(2*(np.exp(krc*pi) - 1)))
plt.plot(modes, np.abs(((Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.', label=r'$h_{ij} = 1.0$')
plt.plot(modes, np.abs(((0.9*Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.', label=r'$h_{ij} = 0.9$')
plt.plot(modes, np.abs(((0.8*Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.', label=r'$h_{ij} = 0.8$')
plt.plot(modes, np.abs(((0.7*Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.', label=r'$h_{ij} = 0.7$')
plt.plot(modes, np.abs(((0.6*Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.', label=r'$h_{ij} = 0.6$')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$|h^{\mathrm{eff}, R}_{ij,n0}|$', fontsize=13)
# Plotting limits on each axis
plt.xlim([-0.5, num_roots - 0.5])
# Call this function in order to show the legend. These arguments allow for the legend to be
# placed in the upper right corner of the graph.
plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0)
plt.grid(True)
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=12)
plt.savefig('Figure55a.png', dpi = 300, pad_inches = 1)
plt.clf()

plt.subplots_adjust(left = 0.16, right = 0.96, top = 0.95, bottom = 0.11)
plt.plot(modes, np.abs(((Nwf*frn*np.exp((3/2)*krc*pi))/(M5rc))), '.')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$|h^{\mathrm{eff}, R}_{ij,n0}|$', fontsize=13)
plt.xlim([-0.5, num_roots - 0.5])
plt.grid(True)
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=12)
plt.savefig('Figure55b.png', dpi = 300, pad_inches = 1)
plt.clf()
# Needed for correct spacing of plots to image border
plt.subplots_adjust(left = 0.14, right = 0.96, top = 0.95, bottom = 0.11)

# Let's first deal with the case of the bulk Majorana particle.
def calcmassmat(h):
    massmat = np.zeros((num_roots, num_roots))
    heffr = []
    heffi = []
    for i in modes0:
        heffr.append(np.zeros((3, 3)))
        heffi.append(np.zeros((3, 3)))
        for j in index1:
            for k in index1:
                heffr[i][j, k] = (h[j, k]*warp32)/(M5rc)*frn[i]*Nwf
                heffi[i][j, k] = (h[j, k]*warp32)/(M5rc)*fin[i]*Nwf
    massmat = np.zeros((3, 3))
    for i in index1:
        for j in index1:
            sumsings = 0
            for s in index1:
                summodes = 0
                for n in modes0:
                    summodes = summodes + heffr[n][i, s]*heffr[n][j, s]*((marn[n]**2)/(marn[n]**2 - M[s]**2))*np.log(marn[n]**2/M[s]**2) - heffi[n][i, s]*heffi[n][j, s]*((main[n]**2)/(main[n]**2 - M[s]**2))*np.log(main[n]**2/M[s]**2)
                sumsings = sumsings + ((M[s])/(pi162))*summodes
            massmat[i, j] = sumsings
    return massmat

# We now plot the sum of various Cijs(n) at different n compared to the 4D case
index5 = np.linspace(0, 4, 5, dtype=int)
index10 = np.linspace(0, 9, 10, dtype=int)
Cijs5 = 0
Cijs10 = 0
for i in index5:
    Cijs5 = Cijs5 + Cijs[i]
for i in index10:
    Cijs10 = Cijs10 + Cijs[i]

plt.subplots_adjust(left = 0.15, right = 0.96, top = 0.95, bottom = 0.11)
m24D = marn[0]**2 - (lambda3 + lambda4 + lambda5)*((vev**2)/2)
mar2 = m24D + (lambda3 + lambda4 + lambda5)*((vev**2)/2)
mai2 = m24D + (lambda3 + lambda4 - lambda5)*((vev**2)/2)
Cc = ((mar2)/(mar2 - M[0]**2))*np.log(mar2/M[0]**2) - ((mai2)/(mai2 - M[0]**2))*np.log(mai2/M[0]**2)

plt.bar([0, 1, 2, 3], [Cijs5/(10**(-4)), Cijs10/(10**(-4)), np.sum(Cijs)/(10**(-4)), Cc/(10**(-4))])
plt.ylabel(r'Sum of coefficients ($\cdot 10^{-4}$)', fontsize=14)
plt.grid(True)
plt.xticks([0, 1, 2, 3], ['n = 5', 'n = 10', 'n = 24', '4D case'], fontsize=15)
plt.yticks(fontsize=13)
plt.savefig('Figure56.png', dpi = 300, pad_inches = 1)
plt.clf()
