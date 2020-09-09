import numpy as np
#Distance from 'inlet' to insonation site
Dist2Inson=90 #mm
#Blood viscosity (Pa.s)
mu=3.4e-3
#Blood density (g/mm3)
rho=1.05e-03
#Vascular elastance parameters
E=1.50e+06 #Pa
h=0.1 #no units

#Definition of geometry
#Generation |Number of vessels at this level | Vessel Radius (mm) | Vessel length (mm) | 
#To eliminate anastomoses simply give them zero length (note you also need to do this in baseline case)
vessels = np.array([(1, 1, 2.0, 100.0,'Uterine'),(2, 1, 2.0, 18.0,'Arcuate'),(3, 50, 0.2, 6.0,'Radial'),(4, 50, 0.2, 9.0,'Anastomose')],
                  dtype=[('generation', 'i4'),('number', 'i4'),('radius', 'f8'),('length', 'f8'),('vessel_type', 'U10')])
#spirals and IVS are defined by  resistance (Pa.s/mm) and compliance (/Pa) [R|C|0=off 1=on]
#To remove these from the model (eg post partum) set third parameter to be zero, otherwise set as 1 
SA_IVS = np.array([1.6,1e-8,1]);

## DEFINE INCIDENT WAVEFORM
#incident harmonics
#Heart rate in beats per minute
HeartRate=72
#steady flow component (baseline) in ml/min - will be scaled with resistance
SteadyFlow=249.0
#Number of flow harmonics
NHar=10

#First ten harmonics of incident waveform [[omega],[A_n],[Phi_n]]. See Mo et al. A transmission line modelling approach to the interpretation of uterine doppler waveforms. Ultrasound in Medicine and Biology, 1988. 14(5): p. 365-376
IWavHar=np.array([[1.0*HeartRate/60.0, 2.0*HeartRate/60.0, 3.0*HeartRate/60.0, 4.0*HeartRate/60.0, 5.0*HeartRate/60.0, 6.0*HeartRate/60.0, 7.0*HeartRate/60.0, 8.0*HeartRate/60.0, 9.0*HeartRate/60.0, 10.0*HeartRate/60.0],
                [64.32,  43.08,  21.48,   7.68,   2.64,  1.8,    0.96,   0.84,   0.84,   0.48],
                [-1.375319451, -2.138028334,-2.998475655,-3.548254369,-3.394665395,-3.185225885,-3.131120678,-2.448696941,-2.602285915,-2.441715624]])


#Option to plot waveform to screen
plotv='y'

## Define time points at which you want to plot flows
StartTime=0.0
EndTime=60.0/HeartRate #end of a single breath
dt=0.01 #time step for plotting


##BASELINE VALUES FOR COMPARISON - Not necessary to change unless you change the structure of the geometry(!)
#Definition of geometry
#Generation |Number of vessels at this level | Vessel Radius (mm) | Vessel length (mm) | 
vessels_bl = np.array([(1, 1, 2.0, 100.0,'Uterine'),(2, 1, 2.0, 18.0,'Arcuate'),(3, 50, 0.2, 6.0,'Radial'),(4, 50, 0.2, 9.0,'Anastomose')],
                  dtype=[('generation', 'i4'),('number', 'i4'),('radius', 'f8'),('length', 'f8'),('vessel_type', 'U10')])
#spirals and IVS are defined by  resistance and compliance [R|C]
SA_IVS_bl = np.array([1.6,1e-8,1]);

