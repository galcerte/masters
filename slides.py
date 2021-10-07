#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.optimize import fsolve
from scipy.optimize import brentq
from scipy import constants
from scipy import integrate
from scipy import linalg

pi = constants.pi
# AdS curvature, GeV
k = 4.01577*(10**17)
# k times r_c
krc = 32.7961/pi
# Compactification radius (GeV⁻¹)
rc = krc/k
# Planck scale (GeV)
MPl = 1.22*(10**19)
# Reduced Planck scale (GeV)
MPlr = MPl/np.sqrt(8*pi)
# Fundamental scale of gravity (GeV)
M5 = np.cbrt( (k*(MPlr**2))/(1 - np.exp(-2*krc*pi)) )
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

# We first define the corrections with respect to the (brane mass)-less terms,
warp = np.exp(-krc*pi)
correta = ((lambda3*(vev**2))/(4*k*M5))*np.exp(2*krc*constants.pi)
corrr = (((lambda3 + lambda4 + lambda5)*(vev**2))/(4*k*M5))*np.exp(2*krc*constants.pi)
corri = (((lambda3 + lambda4 - lambda5)*(vev**2))/(4*k*M5))*np.exp(2*krc*constants.pi)

# We wish to find the roots of the expressions below,
# which'll tell us the mass spectrums of the bulk scalars we consider.
# Define the expression whose roots we want to find:
def masseta(x, v):
    return ((2 + correta)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + correta)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    #return (2 + correta)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def massr(x, v):
    return ((2 + corrr)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + corrr)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    #return (2 + corrr)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def massi(x, v):
    return ((2 + corri)*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - ((2 + corri)*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
    #return (2 + corri)*sp.jv(v, x) + x*sp.jvp(v, x, 1)
def mass2(x, v):
    return (2*sp.jv(v, x) + x*sp.jvp(v, x, 1))*(2*sp.yv(v, warp*x) + warp*x*sp.yvp(v, warp*x, 1)) - (2*sp.yv(v, x) + x*sp.yvp(v, x, 1))*(2*sp.jv(v, warp*x) + warp*x*sp.jvp(v, warp*x, 1))
def massless(x):
    return (2*sp.jv(2, x) + x*sp.jvp(2, x, 1))*(2*sp.yv(2, warp*x) + warp*x*sp.yvp(2, warp*x, 1)) - (2*sp.yv(2, x) + x*sp.yvp(2, x, 1))*(2*sp.jv(2, warp*x) + warp*x*sp.jvp(2, warp*x, 1))

# Use the numerical solver to find the roots
x_init_guesseta = [4.66874849, 7.72492277, 10.7, 13.88, 16.988, 20.103, 23.22, 26.35, 29.5, 32.61, 35.7, 38.9, 42.01, 45.1, 48.28, 51.42, 54.56, 57.7, 60.84, 64, 67.11, 70.25, 73, 76]
x_rootseta = fsolve(masseta, x_init_guesseta, args=v)

x_init_guessr = [4.97, 8.07, 11.13, 14.19, 17.27, 20.35, 23.45, 26.55, 29.66, 32.78, 35.9, 39.02, 42.1, 45.28, 48.4, 51.42, 54.67, 57.8, 60.84, 64, 67.11, 70.25, 73, 76]
x_rootsr = fsolve(massr, x_init_guessr, args=v)

x_init_guessi = [4.97, 8.07, 11.13, 14.19, 17.27, 20.35, 23.45, 26.55, 29.66, 32.78, 35.9, 39.02, 42.1, 45.28, 48.4, 51.42, 54.67, 57.8, 60.84, 64, 67.11, 70.25, 73, 76]
x_rootsi = fsolve(massi, x_init_guessi, args=v)

x_init_guess = [3.96089237, 7.16961637, 10.33728283, 13.49278243, 16.64302665, 19.79052262, 22.93640086, 26.08124643, 29.22539253, 32.36904292, 35.51232917, 38.65534017, 41.79813792, 44.94076717, 48.08326097, 51.22564422, 54.36793611, 57.51015159, 60, 63, 66, 69, 72, 76]
x_roots2 = fsolve(mass2, x_init_guess, args=v)
x_rootsmassless = fsolve(massless, x_init_guess)
num_roots = np.shape(x_init_guess)[0]

# Turns out some results can be bogus because the root solver above
# is not finding the right root, it's solving the ones farther away
# 79 instead of 72-73 and such.

# Masses of the modes are in GeV
maetan = k*x_rootseta*warp
marn = k*x_rootsr*warp
main = k*x_rootsi*warp
m2n = k*x_roots2*warp
masslessspec = k*x_rootsmassless*warp

# Now that we know the mass spectrum of the eta doublet, it is time to evaluate
# their wavefunctions at the IR brane.
ind = np.linspace(0, num_roots - 2, num_roots - 1, dtype=int)
# Coefficients. If x_n << the warp factor, then we can take bnv = 0.
betanv = - (2*sp.jv(v, maetan/k) + (maetan/k)*sp.jvp(v, maetan/k, 1))/(2*sp.yv(v, maetan/k) + (maetan/k)*sp.yvp(v, maetan/k, 1))
brnv = - (2*sp.jv(v, marn/k) + (marn/k)*sp.jvp(v, marn/k, 1))/(2*sp.yv(v, marn/k) + (marn/k)*sp.yvp(v, marn/k, 1))
binv = - (2*sp.jv(v, main/k) + (main/k)*sp.jvp(v, main/k, 1))/(2*sp.yv(v, main/k) + (main/k)*sp.yvp(v, main/k, 1))

# Normalization constants
netan = np.zeros(num_roots)
nrn = np.zeros(num_roots)
nin = np.zeros(num_roots)

index = np.linspace(0, num_roots - 1, num_roots, dtype=int)

# Initializing wavefunctions and where mass matrix'll go
fetan = np.zeros(num_roots)
frn = np.zeros(num_roots)
fin = np.zeros(num_roots)
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

plt.subplots_adjust(left = 0.15, right = 0.96, top = 0.95, bottom = 0.11)
plt.plot(index, (marn - main)/(vev**2), '.')
plt.xlabel('KK-mode number', fontsize=13)
plt.ylabel(r'$|m^2_R - m^2_I|/v^2$', fontsize=20)
#plt.legend()
plt.grid(True)
plt.xticks(np.linspace(0, num_roots - 1, num_roots, dtype=int), np.linspace(0, num_roots - 1, num_roots, dtype=int), fontsize=11)
plt.yticks(fontsize=13)
plt.savefig('splittingm2.png', dpi = 300, pad_inches = 1)
plt.clf()

# Index keeping track of how many values of m_2 we're going to be taking a look at.
m2_times = np.linspace(0, 4, 5, dtype=int)
m2 = np.array([0.2, 0.5, 1, 2, 5])*4.01577*(10**17)
v = np.sqrt((m2/k)**2 + 4)

phi = np.linspace(0.001, pi, 200)
index_phi = np.linspace(0, 199, 200, dtype=int)
plt.subplots_adjust(left = 0.15, right = 0.96, top = 0.95, bottom = 0.11)

def frnf(nrn, brnv, metarn, phir, vr):
    return (1/nrn)*(sp.jv(vr, (metarn/k)*np.exp(krc*phir)) + brnv*sp.yv(vr, (metarn/k)*np.exp(krc*phir)))
def finf(nin, binv, metain, phir, vr):
    return (1/nin)*(sp.jv(vr, (metain/k)*np.exp(krc*phir)) + binv*sp.yv(vr, (metain/k)*np.exp(krc*phir)))

mrnfm = np.zeros(5)
minfm = np.zeros(5)
brv = np.zeros(5)
biv = np.zeros(5)
nr = np.zeros(5)
ni = np.zeros(5)

def wavefuncrn2(phir, marn, brnv, vr):
    return np.exp(2*krc*phir)*((sp.jv(vr, (marn/k)*np.exp(krc*phir)) + brnv*sp.yv(vr, (marn/k)*np.exp(krc*phir)))**2)
def wavefuncin2(phir, main, binv, vr):
    return np.exp(2*krc*phir)*((sp.jv(vr, (main/k)*np.exp(krc*phir)) + binv*sp.yv(vr, (main/k)*np.exp(krc*phir)))**2)

for i in m2_times:
    x_rootsr2 = fsolve(massr, 3, args=v[i])
    x_rootsi2 = fsolve(massi, 3, args=v[i])
    mrnfm[i] = k*x_rootsr2*warp
    minfm[i] = k*x_rootsi2*warp

    brv[i] = - (2*sp.jv(v[i], mrnfm[i]/k) + (mrnfm[i]/k)*sp.jvp(v[i], mrnfm[i]/k, 1))/(2*sp.yv(v[i], mrnfm[i]/k) + (mrnfm[i]/k)*sp.yvp(v[i], mrnfm[i]/k, 1))
    biv[i] = - (2*sp.jv(v[i], minfm[i]/k) + (minfm[i]/k)*sp.jvp(v[i], minfm[i]/k, 1))/(2*sp.yv(v[i], minfm[i]/k) + (minfm[i]/k)*sp.yvp(v[i], minfm[i]/k, 1))

    nr[i] = np.sqrt(2*integrate.quad(wavefuncrn2, 0, pi, args=(mrnfm[i], brv[i], v[i]))[0])
    ni[i] = np.sqrt(2*integrate.quad(wavefuncin2, 0, pi, args=(minfm[i], biv[i], v[i]))[0])

frnf1 = np.zeros(200)
frnf2 = np.zeros(200)
frnf3 = np.zeros(200)
frnf4 = np.zeros(200)
frnf5 = np.zeros(200)

finf1 = np.zeros(200)
finf2 = np.zeros(200)
finf3 = np.zeros(200)
finf4 = np.zeros(200)
finf5 = np.zeros(200)
for i in index_phi:
    frnf1[i] = frnf(nr[0], brv[0], mrnfm[0], phi[i], v[0])
for i in index_phi:
    frnf2[i] = frnf(nr[1], brv[1], mrnfm[1], phi[i], v[1])
for i in index_phi:
    frnf3[i] = frnf(nr[2], brv[2], mrnfm[2], phi[i], v[2])
for i in index_phi:
    frnf4[i] = frnf(nr[3], brv[3], mrnfm[3], phi[i], v[3])
for i in index_phi:
    frnf5[i] = frnf(nr[4], brv[4], mrnfm[4], phi[i], v[4])

for i in index_phi:
    finf1[i] = finf(ni[0], biv[0], minfm[0], phi[i], v[0])
for i in index_phi:
    finf2[i] = finf(ni[1], biv[1], minfm[1], phi[i], v[1])
for i in index_phi:
    finf3[i] = finf(ni[2], biv[2], minfm[2], phi[i], v[2])
for i in index_phi:
    finf4[i] = finf(ni[3], biv[3], minfm[3], phi[i], v[3])
for i in index_phi:
    finf5[i] = finf(ni[4], biv[4], minfm[4], phi[i], v[4])

plt.plot(phi, frnf1, '-', label=r'm/k = 0.2')
plt.plot(phi, frnf2, '-', label=r'm/k = 0.5')
plt.plot(phi, frnf3, '-', label=r'm/k = 1.0')
plt.plot(phi, frnf4, '-', label=r'm/k = 2.0')
plt.plot(phi, frnf5, '-', label=r'm/k = 5.0')
plt.xlabel(r'$\varphi$', fontsize=13)
plt.ylabel(r'Wavefunction value', fontsize=20)
plt.legend()
plt.grid(True)
plt.xticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi], ['0: UV-brane', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$: IR-brane'])
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig('profilesm2R.png', dpi = 300, pad_inches = 1)
plt.clf()

plt.plot(phi, finf1, '-', label=r'm/k = 0.2')
plt.plot(phi, finf2, '-', label=r'm/k = 0.5')
plt.plot(phi, finf3, '-', label=r'm/k = 1.0')
plt.plot(phi, finf4, '-', label=r'm/k = 2.0')
plt.plot(phi, finf5, '-', label=r'm/k = 5.0')
plt.xlabel(r'$\varphi$', fontsize=13)
plt.ylabel(r'Wavefunction value', fontsize=20)
plt.legend()
plt.grid(True)
plt.xticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi], ['0: UV-brane', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$: IR-brane'])
plt.yticks(fontsize=13)
plt.yscale('log')
plt.savefig('profilesm2I.png', dpi = 300, pad_inches = 1)
plt.clf()