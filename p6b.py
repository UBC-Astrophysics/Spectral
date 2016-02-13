import numpy as np
import matplotlib.pyplot as plt
import math
#
# p6b.py - non-linear wave equation
#
omegap2=0.01
#omegap2=0.0
A=0.006
#A=0.00
b0=0;
amp=0.5
k=1
resist=0.00

# Grid and initial data:
N = 1024
h = 2*math.pi/N
x = h*np.arange(1,N+1)
t = 0
dt = h/4

v = amp*np.sin(k*x)
vold = amp*np.sin(k*(x+dt))

# v = np.exp(-100*(x-1)**2)
# vold = np.exp(-100*(x-dt-1)**2)


# Time-stepping by leap frog formula:
tmax = 256
tplot = .15

plotgap = round(tplot/dt)
nplots = round(tmax/tplot)
dt = tplot/plotgap
data = [v]
wdata = []
tdata=[t]
bedata=[np.mean((v*k)**2)]
eedata=[np.mean(((v-vold)/dt)**2)]
cedata=[np.mean(omegap2*v**2)]
qedata=[-4*A*np.mean( ((v*k)**2+((v-vold)/dt)**2)*((v*k)**2-((v-vold)/dt)**2))]

sampfreq=1j*np.array(range(0,N/2)+[0]+range(-N/2+1,0))

for i in range(1,int(nplots+1)):
    for n in range(1,int(plotgap+1)):
        t = t+dt
        dvdt = (v-vold)/dt                              # E = -diff(Ay,t)
        v_hat = np.fft.fft(v)
#       w_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* v_hat;
        w_hat = sampfreq*v_hat
        w = np.real(np.fft.ifft(w_hat))                     # B = diff(Ay,x)
        b = w + b0
#       w2_hat = 1i*[0:N/2-1 0 -N/2+1:-1] .* w_hat;     
#        w2_hat = (range(0,N/2)+[0]+range(-N/2+1,0)) * w_hat
        w2_hat = sampfreq*w_hat
#        w2 = np.real(np.fft.ifft(w2_hat))                  # diff(Ay,x$2)
        w2 = np.real(np.fft.ifft(w2_hat))                       # diff(Ay,x$2)
        n1factor=-1+A*(24*b*b-8*dvdt*dvdt)
        n2factor=1+A*(24*dvdt*dvdt-8*b*b)
        vnew = 2*v - vold + dt**2*(-n1factor*w2-omegap2*v-resist*dvdt)/n2factor
        dvdt = (vnew-vold)/(2*dt);                      # E = -diff(Ay,t)
        n1factor=-1+A*(24*b*b-8*dvdt*dvdt)
        n2factor=1+A*(24*dvdt*dvdt-8*b*b)
        vnew = 2*v - vold + dt**2*(-n1factor*w2-omegap2*v-resist*dvdt)/n2factor
        vold = v
        v = vnew
#    plt.plot(x,w)
#    plt.savefig('Plot_%03d.png' % i)
#    plt.close()
    wdata.append(w)
    data.append(v)
    tdata.append(t)
    bedata.append(np.mean(w*w))
    eedata.append(np.mean(dvdt*dvdt))
    cedata.append(np.mean(omegap2*v*v))
    qedata.append(-4*A*np.mean( (w*w+dvdt*dvdt)*(b*b-dvdt*dvdt)))

#  waterfall(x,tdata,data), view(10,70), colormap(1e-6*[1 1 1]);
#  axis([0 2*pi 0 tmax -amp amp]), ylabel t, zlabel u, grid off
etotal = bedata+eedata+cedata+qedata
# rat=etotal/etotal(1)-1
np.save("wdata.npy",wdata)
np.save("vdata.npy",data)
# plot(x,ifft(1i*[0:N/2-1 0 -N/2+1:-1].*fft(data(1,:))),x,ifft(1i*[0:N/2-1 0 -N/2+1:-1].*fft(data(length(tdata),:))))







