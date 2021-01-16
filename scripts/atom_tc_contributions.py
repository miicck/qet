import sys
import numpy as np
import matplotlib.pyplot as plt
from qet import parser
from qet.calculations import tc_allen_dynes
from qet.constants import RY_TO_CMM1

modes      = parser.phonon_interp_modes("ph_interp.modes")
a2f        = parser.a2f_dos_out("a2F.dos10")
a2f_freqs  = a2f["frequencies"]
mode_freqs = modes["frequencies"]
mode_evecs = modes["eigenvectors"]
mode_qpts  = modes["q-points"] 
mode_count = int(len(mode_freqs)/len(mode_qpts))
atom_count = int(mode_count/3)
CMM_SCALE  = 10

def tc_from_a2f(a2f): return tc_allen_dynes(a2f_freqs, a2f, mu_stars=[0.125])[0.125]

# Work out the Tc contribution from each frequency
tc_contribs = []
tc_with_all = tc_from_a2f(a2f["a2f"])
for i in range(len(a2f_freqs)):
    a2fi = list(a2f["a2f"])
    a2fi[i] = 0
    tc_i = tc_with_all - tc_from_a2f(a2fi)
    if tc_i < 0: tc_i = 0
    tc_contribs.append(tc_i)

def norm(e): return (e[0]*e[0] + e[1]*e[1] + e[2]*e[2])**0.5

# Work out which atoms move at each frequency
atom_amps = []

for i in range(len(a2f_freqs)):
    wi = a2f_freqs[i] * RY_TO_CMM1
    evec_i = None

    for j in range(len(mode_freqs)):
        wj = mode_freqs[j]
        exp = np.exp(-(wi-wj)**2.0/(CMM_SCALE**2.0))

        evec = [norm(e)*exp for e in mode_evecs[j]]
        if evec_i is None: evec_i = evec
        else: evec_i = [a+b for a,b in zip(evec_i, evec)]

    evec_i = [e/len(mode_freqs) for e in evec_i]
    atom_amps.append(evec_i)

# Work out Tc contribution from atom moves
atom_tcs = []
for i in range(atom_count):
    
    tc_av = 0
    for j in range(len(a2f_freqs)):
        tc_av += tc_contribs[j] * atom_amps[j][i]

    print(tc_av)
    atom_tcs.append(tc_av)

plt.subplot(411)
plt.plot(a2f_freqs, a2f["a2f"])
plt.xlabel("$\omega$")
plt.ylabel("$\\alpha^2F(\\omega)$")

plt.subplot(412)
plt.plot(a2f_freqs, tc_contribs)
plt.xlabel("$\omega$")
plt.ylabel("$T_c$ contribution")

plt.subplot(413)
for i in range(0, atom_count):
    plt.plot(a2f_freqs, [aa[i] for aa in atom_amps], label="Atom {0} movement".format(i))
plt.xlabel("$\omega$")
plt.legend()

plt.subplot(414)
plt.plot(atom_tcs)
plt.xlabel("Atom number")
plt.ylabel("$T_c$ contributions")
plt.show()
