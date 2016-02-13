import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.animation as animation

#
# p6.m - variable coefficient wave equation
#   ported from p6.py from Trefethen Spectral Methods in MATLAB
# Grid, variable coefficient, and initial data:
#
N=128
h=2*math.pi/N
x = h*np.arange(1,N+1)
t = 0
dt = h/4

c = .2 + np.sin(x-1)**2;
v = np.exp(-100*(x-1)**2);
vold = np.exp(-100*(x-.2*dt-1)**2);

# Time-stepping by leap frog formula:
tmax = 8
tplot = .15

plotgap = round(tplot/dt)
nplots = round(tmax/tplot)
dt = tplot/plotgap

data = [v]
tdata=[t]

sampfreq=1j*np.array(range(0,N/2)+[0]+range(-N/2+1,0))

for i in range(1,int(nplots+1)):
    for n in range(1,int(plotgap+1)):
        t = t+dt
        v_hat = np.fft.fft(v)
        w_hat = sampfreq*v_hat
        w = np.real(np.fft.ifft(w_hat))                
        vnew = vold - 2*dt*c*w
        vold = v
        v = vnew
    data.append(v)
    tdata.append(t)

fig, ax = plt.subplots()
line, = ax.plot(x, data[0])

def animate(i):
#    line.set_ydata(np.convolve(wdata[i],w,mode='same'))  # update the data
    line.set_ydata(data[i])
    return line,


# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, len(data)), init_func=init,
                              interval=25, blit=True)
plt.show()





