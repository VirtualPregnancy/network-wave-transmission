#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import VesselDefinition as params
import FunctionDefinitions as funcs

## Calculate total resistance of the system and compare to baseline (flow decreases by this factor as resistance increases assuming a constant driving pressure)
[TotalResistance,VenousResistance]=funcs.total_resistance(params.vessels,params.SA_IVS)
[BaselineResistance,BaselineVenous]=funcs.total_resistance(params.vessels_bl,params.SA_IVS_bl)
SteadyFlow=params.SteadyFlow*BaselineResistance/TotalResistance

## Calculate characteristic admittance,wave propogation constant, and uterine artery compliance
[CharacteristicAdmittance,WaveProp,UtCompliance]=funcs.characteristic_admittance(params.vessels,params.SA_IVS)
#Calculate effective admittance and reflection coefficient for each vessel
[EffectiveAdmittance,ReflectionCoefficients]=funcs.effective_admittance(params.vessels,params.SA_IVS,CharacteristicAdmittance,WaveProp,VenousResistance)

#Calculate phase and amplitude offset for waveform at insonation site vessel
#Output is [[amplitude],[phase]]
reflect_coeff=np.transpose(np.column_stack((np.absolute(ReflectionCoefficients[0][:]),np.angle(ReflectionCoefficients[0][:]))))
char_admit=np.transpose(np.column_stack((np.absolute(CharacteristicAdmittance[0][:]),np.angle(CharacteristicAdmittance[0][:]))))

#Real and imaginary parts of wave propogation constant
wave_prop_constant=np.transpose(np.column_stack((WaveProp[0][:].real,WaveProp[0][:].imag)))

## Incident flow profile, insonation site flow (incident, reflected, total), insonation site pressure, insonation site velocity
[InsonationSiteVelocity,time]=funcs.timecourse(params.StartTime,params.EndTime,params.dt,reflect_coeff,char_admit,wave_prop_constant,SteadyFlow,UtCompliance,'Uterine')



#Output to screen key properties of the flow velocity waveform    
funcs.flow_velocity_properties(InsonationSiteVelocity)    
#If selected as an option in VesselDefinition.py, plot flow velocity waveform
if(params.plotv=="y"): 
    plt.xlabel('time (seconds)')
    plt.ylabel('velocity (cm/s)') 
    plt.title('Uterine artery flow velocity waveform')      
    plt.plot(time,InsonationSiteVelocity)
    plt.show()


    
    