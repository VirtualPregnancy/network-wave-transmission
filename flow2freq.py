# Testing For Fourier Transforms
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

###
# TEST CASE FOR HUMAN
###
HR = 72. #bpm
T = 60./HR
flow = np.array(
    [210.653207, 219.1012794, 229.4833496, 241.2798739, 253.8980065, 266.7824325, 279.494959, 291.7446355, 303.3671806,
     314.2686824, 324.3586639, 333.4986557, 341.4846111, 348.0680046, 353.0062406, 356.122971, 357.356294, 356.7780408,
     354.5781526, 351.0204085, 346.3851163, 340.917732, 334.7989511, 328.1432577, 321.0225826, 313.5033642, 305.6817134,
     297.7034503, 289.7623095, 282.0779147, 274.8620812, 268.2852683, 262.4537651, 257.4033666, 253.1090406,
     249.5048545, 246.5060522, 244.026171, 241.9856176, 240.312427, 238.9391125, 237.8003591, 236.8347073, 235.9903076,
     235.2318802, 234.5446444, 233.9318808, 233.4055654, 232.9729094, 232.6240469, 232.3263263, 232.0284089,
     231.6734783, 231.2169114, 230.6414723, 229.9635654, 229.2272942, 228.4878526, 227.7902138, 227.1513291,
     226.5530555, 225.9490447, 225.2833288, 224.5135091, 223.6292916, 222.6586299, 221.658467, 220.6932494, 219.8095451,
     219.0170911, 218.284448, 217.5518172, 216.7566025, 215.8617342, 214.8749984, 213.8505564, 212.8706223, 212.01336,
     211.3192789, 210.7700827, 210.2900419, 209.7716461, 209.1175185, 208.2831246])
# comes roughly from Mo et al. Ultrasound Med Biol 14(5)1988 (fourier series converted to time series)
flow = flow * 1000. / 60.  # scale to mm^3/s
###
# RAT
###
HR = 408
T = 60./HR
flow = np.array([105.3258471,145.8323486,180.4213225,188.4143689,206.8027905,216.2868497,231.4115793,240.4401918,
                 238.6985516,235.8919197,233.2664555,230.7410468,228.215638,209.8801124,192.6008121,183.4611276,
                 174.6004657,168.7419038,159.5639903,151.5794089,139.1027223,133.281779,128.7315909,128.7315909,
                 126.9838225,125.2311349,123.3248408,121.321373,120.3723967,120.3723967,120.3723967,120.3723967,
                 119.6903991,117.7270627,115.7637264,114.1094796,112.9883184,111.8671572,110.5030612,107.3363274,
                 106.997686])

#donvert to ml/min
flow = flow * 60./1000.

# CALCULATE THE COEFFICIENTS FOR EACH HARMONIC
Y = fft(flow)

#print(Y)

a0 = 0
for i in range(0,len(flow)):
    a0 = a0 + (1./len(flow))*flow[i] #basically the mean value of flow

print(a0)

no_possible_freqs = int(np.ceil(len(flow)/2.))

#print(no_possible_freqs)
a = np.zeros(no_possible_freqs-1)
b = np.zeros(no_possible_freqs-1)
amp = np.zeros(no_possible_freqs-1)
phase = np.zeros(no_possible_freqs-1)
for n in range (1,no_possible_freqs):
    a[n-1] = (2./(len(flow)))*Y[n].real
    b[n-1] = -(2./(len(flow)))*Y[n].imag

    amp[n-1] = (a[n-1]**2.+b[n-1]**2.)**(1./2.)
    phase[n-1] = -1.*np.arctan2(b[n-1],a[n-1])

#print(a)
#print(b)

print(amp)
print(phase)#*180./np.pi)

flow_time = np.zeros(40) + a0
time = np.zeros(40)
for n in range(0,20):
    omega = 2.*(n+1)*np.pi/T
    for i in range(0,40):
        time[i] = i*T/40.
        flow_time[i]=flow_time[i] + amp[n]*np.cos(omega*time[i] + phase[n])

#print(len(flow))
#plt.plot(flow_time)
#plt.plot(flow)
#plt.show()

#print(np.max(flow),np.min(flow))
#print(np.max(flow_time)-np.max(flow), np.min(flow_time)-np.min(flow))




