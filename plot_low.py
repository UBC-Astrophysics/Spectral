import numpy as np
from sys import argv
import matplotlib.pyplot as plt

data=np.loadtxt(argv[1],unpack=True)
x=data[0]
y=np.transpose(data[3:-1])
plt.ylim(300,800)
plt.plot(x*4.9e9,-91.16/y,'k,')
plt.xlabel('Magnetic Field [Gauss]')
plt.ylabel('Wavelength [nm]')
plt.show()
