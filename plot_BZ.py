import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

infile_name = sys.argv[1]

infile = h5py.File(infile_name, 'r')

print(list(infile.keys()))
key0 = 'Gk'

print("plotting key :", key0)

data = infile[key0][()]

infile.close()

print('data.shape = {}'.format(data.shape))

#momenta for which to plot data as function of frequencies
kx = 0
ky = 0
Nw = data.shape[0]
nu_data = data[:, ky, kx]
nu_arr = np.linspace(- Nw + 1, Nw - 1, Nw, endpoint=True)

plt.plot(nu_arr, nu_data.real, marker='x', linestyle='None')
plt.plot(nu_arr, nu_data.imag, marker='x', linestyle='None')
plt.show()




