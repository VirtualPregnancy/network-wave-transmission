# Uterine Artery Wave Transmission
A library to solve biological wave transmission problems in simple networks of arteries representing the uterine vasculature.

## Project description
This a set of python scripts that allows simulation software allows simulation of wave transmission in the uterine arteries. Users can define properties of different classes of uterine arteries (specifically uterine, radial, arteriovenous anastomoses) and a terminal load representing the placental bed. The scripts will output metrics that represent Doppler ultrasound markers of placental health: The ratio of systolic to diastolic flow velocities (S/D), the resistance index (RI) and the pulsatility index (PI). If the model predicts dicrotic notching, properties of the notch are also output to screen. A plot of the velocity waveform predicted at the ultrasound insonation site is also output.

### Pre-requisites
- The ability to run python scripts from the command line
- The python numpy package python (http://www.numpy.org)
- The python matplotlib package (https://matplotlib.org)


### Running

From the home directory for the project simply type the following into a command line window
```
  python RUNME.py
```

### Expected output

Output to your command line window will summarise predicted velocity waveform properties:
```
Velocity waveform properties 
=============================
S/D = 1.64514486781
RI = 0.392150795003
PI = 0.52725811502
No notch present
```
An (optional) window will also open, showing the velocity waveform at the insonation site.

### Input
The model geometry and parameterisation is defined in *VesselDefinition.py*, and parameter definitions are given within this script. The user can change the physical parameters (blood density and viscosity), elasticity of blood vessels and the size and number of vessels at each generation in the uterine artery structure.

### What does each file actually do (in brief)?
Code structure:

- *VesselDefinition.py:* Defines the properties of the vessels 
- *RUNME.py:* Calls functions required to set up model geometry, conduct simulations and plot outputs
- *FunctionDefinitions.py:* contains all the routines required to calculate admittance, sum through a network structure, apply boundary conditions and export solutions.
