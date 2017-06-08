#!/usr/bin/python
import numpy as np
import VesselDefinition as params
##Function definitions
def total_resistance(vessels,terminals):
    resistance=np.zeros(np.size(vessels))
    total_resistance=0.0
    
    for i in range(0,np.size(vessels)):
        resistance[i]=81.0*params.mu*vessels['length'][i]/(8.0*np.pi* vessels['radius'][i]**4.0)/vessels['number'][i]
        if(vessels['type'][i]=='Anastomose'):
            anast_index=i
        else:
            total_resistance=total_resistance+resistance[i]

    if(vessels['length'][anast_index]==0.0):
        venous_resistance=0.0
        #Only IVS contribution to resistance (as in Mo et al.)
        parallel_resistance=terminals[0]
    else:        
        venous_resistance=total_resistance/2.0 #Assuming veins are half as resistive as arteries
        #add terminal and anastomosis resistance in parallel
        parallel_resistance=1.0/(1.0/(resistance[anast_index]+venous_resistance)+terminals[2]*1.0/terminals[0])
        
    total_resistance=total_resistance+parallel_resistance
    return [total_resistance,venous_resistance]
    
    
def characteristic_admittance(vessels,terminals):
    char_admit=np.zeros((np.size(vessels),params.NHar),dtype=complex)
    prop_const=np.zeros((np.size(vessels),params.NHar),dtype=complex)
    C=np.zeros(np.size(vessels))
    L=np.zeros(np.size(vessels))
    R=np.zeros(np.size(vessels))
    G=np.zeros(np.size(vessels))
    for i in range(0,np.size(vessels)):
        hbar=params.h*vessels['radius'][i]
        C[i]=3.0*np.pi*vessels['radius'][i]**3/(2.0*hbar*params.E)# Compliance per unit length
        L[i]=9.0*params.rho/(4.0*np.pi*vessels['radius'][i]**2)#per unit length
        R[i]=81.0*params.mu/(8.0*np.pi*vessels['radius'][i]**4) #laminar resistance per unit length
        G[i]=0.0
        for j in range(0,params.NHar):
            omega=(j+1)*2.0*np.pi*params.HeartRate/60.0
            char_admit[i][j]=np.sqrt(G[i]+np.complex(0.0,1.0)*omega*C[i])/np.sqrt(R[i]+np.complex(0.0,1.0)*omega*L[i])*vessels['number'][i] #summed in parallel
            prop_const[i][j]=np.sqrt((G[i]+np.complex(0.0,1.0)*omega*C[i])*(R[i]+np.complex(0.0,1.0)*omega*L[i]))
            
    utcomp=C[0]
    return [char_admit,prop_const,utcomp]
    
def effective_admittance(vessels,terminals,char_admit,prop_const,v_resist):
    eff_admit=np.zeros((np.size(vessels),params.NHar),dtype=complex)
    reflect=np.zeros((np.size(vessels),params.NHar),dtype=complex)
    #First consider the terminal admitance, which is the admittance of the asastomoses and veins (in series) added in parallel to the SA/IVS admittance
    for i in range(0,np.size(vessels)):
        if(vessels['type'][i]=='Anastomose'):
            for j in range(0,params.NHar):
                eff_admit[i][j]=char_admit[i][j]/(1.0+char_admit[i][j]*v_resist) #adding venous resistance in series
            for j in range(0,params.NHar):
                omega=(j+1)*2.0*np.pi*params.HeartRate/60.0
                IVS_admit=terminals[2]*(1.0+np.complex(0.0,1.0)*omega*terminals[0]*terminals[1])/terminals[0]
                eff_admit[i][j]=eff_admit[i][j]+IVS_admit
                
    #Step backward through network and calculate effecitive admittances 
    for i in range(np.size(vessels)-2,-1,-1):
        daughter_admit=eff_admit[i+1][:]
        reflect[i][:]=np.divide(char_admit[i][:]-daughter_admit,char_admit[i][:]+daughter_admit)
        eff_admit[i][:]=char_admit[i][:]*(1.0-reflect[i][:]*np.exp(-2.0*prop_const[i][:]*vessels['length'][i]))\
            /(1.0+reflect[i][:]*np.exp(-2.0*prop_const[i][:]*vessels['length'][i]))
    
    return [eff_admit,reflect]
    
def flow_velocity_properties(velocity):
    print "Velocity vaveform properties "
    print "============================="
    print 'S/D = ' + str(np.max(velocity)/np.min(velocity))
    print 'RI = ' + str((np.max(velocity)-np.min(velocity))/np.max(velocity))
    print 'PI = ' + str((np.max(velocity)-np.min(velocity))*np.size(velocity)/np.sum(velocity))
    
    #check for notch - note that this fails for complicated waveforms but works for most physiologically realistic parameter sets
    point1=1
    point2=1
    checksize=np.size(velocity)-2
    while(velocity[point1+1]>velocity[point1] and point1<checksize): #expect one monotonically increasing then decreasing waveform prior to notch
        point1=point1+1
    while(velocity[point1+1]<velocity[point1] and point1<checksize):
        point1=point1+1
    point2=point1
    while(velocity[point2+1]>velocity[point2] and point2<checksize):
        point2=point2+1
    if(point1 < point2):
        print "Notch present"
        print "Notch height " + str(velocity[point2]-velocity[point1])
        print "Notch ratio "+ str((velocity[point2]-velocity[point1])/(np.max(velocity)-np.min(velocity)))
        
    else:
        print "No notch present"
        
def timecourse(StartTime,EndTime,dt,reflect_coeff,char_admit,wave_prop_constant,SteadyFlow,UtCompliance):
    LengthOfUterine=params.vessels['length'][0] #mm
    UterineRadius=params.vessels['radius'][0]#mm
    NTime=int(np.ceil((EndTime-StartTime)/dt))
    time=np.zeros(NTime)
    IncidentFlow=np.zeros(NTime)
    InsonationIncidentFlow=np.zeros(NTime) #ml/min
    InsonationReflectFlow=np.zeros(NTime) #ml/min
    InsonationSitePressureIn=np.zeros(NTime) #Pa
    InsonationSitePressureRef=np.zeros(NTime) #Pa
    InsonationSitePressure=np.zeros(NTime) #Pa
    InsonationSiteVelocity=np.zeros(NTime) #cm/s
    InsonationSiteTotalFlow=np.zeros(NTime)
    for i in range(0,NTime):
        time[i]=dt*i
        for j in range(0,params.NHar):

            IncidentFlow[i]=IncidentFlow[i]+params.IWavHar[1,j]*np.cos(2*np.pi*time[i]*params.IWavHar[0,j]+params.IWavHar[2,j])
        
            InsonationIncidentFlow[i]=InsonationIncidentFlow[i]+\
            params.IWavHar[1,j]*np.exp(-params.Dist2Inson*wave_prop_constant[0,j])*\
            np.cos(2*np.pi*time[i]*params.IWavHar[0,j]+params.IWavHar[2,j]-wave_prop_constant[1,j]*params.Dist2Inson)
        
            InsonationReflectFlow[i]=InsonationReflectFlow[i]-\
            reflect_coeff[0,j]*params.IWavHar[1,j]*np.exp((-params.Dist2Inson-LengthOfUterine)*wave_prop_constant[0,j])*\
            np.cos(2*np.pi*time[i]*params.IWavHar[0,j]+params.IWavHar[2,j]-wave_prop_constant[1,j]*(params.Dist2Inson+LengthOfUterine)+reflect_coeff[1,j])
        
            #Incident pressure is characteristic impedance times incident flow
            InsonationSitePressureIn[i]=InsonationSitePressureIn[i]+\
            params.IWavHar[1,j]*np.exp(-params.Dist2Inson*wave_prop_constant[0,j])*1000.0/60.0*np.cos(2*np.pi*time[i]*params.IWavHar[0,j]+params.IWavHar[2,j]-wave_prop_constant[1,j]*params.Dist2Inson-char_admit[1,j])/char_admit[0,j]
        
            #Reflected pressure is characteristic impedance times reflcted flow
            InsonationSitePressureRef[i]=InsonationSitePressureRef[i]-\
            reflect_coeff[0,j]*params.IWavHar[1,j]*np.exp((-params.Dist2Inson-LengthOfUterine)*wave_prop_constant[0,j])*1000.0/60.0*\
            np.cos(2*np.pi*time[i]*params.IWavHar[0,j]+params.IWavHar[2,j]-wave_prop_constant[1,j]*(params.Dist2Inson+LengthOfUterine)+reflect_coeff[1,j]-char_admit[1,j])/char_admit[0,j]
        
                
            InsonationSitePressure[i]=InsonationSitePressureIn[i]+InsonationSitePressureRef[i]+80*133   
            InsonationSiteTotalFlow[i]= InsonationIncidentFlow[i]+InsonationReflectFlow[i]+SteadyFlow
            InsonationSiteVelocity[i]=(InsonationSiteTotalFlow[i])*(1000.0/60.0)/(np.pi*UterineRadius*UterineRadius+UtCompliance*(InsonationSitePressure[i]))/10.0
    
    return [InsonationSiteVelocity,time]