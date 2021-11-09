# mic-vib-matrices
Several acoustic duct transfer matrices useful in microphone vibration

## usage

### straightpipe.m 
```T = straightpipe(kL)``` creates the transfer matrix for an acoustical duct with a uniform area using ```kL=2*pi*f*length/c``` from ```f```, the frequency of the sound propagating through the duct and the total ```length``` of the pipe segment. 
### junction.m
```T = junction(area_ratio)``` creates the transfer matrix for a junction between two different pipes (also called an area discontinuity). If the first pipe area is ```s_1``` and the second pipe area is ```s_2```, area ratio should be ```s_1/s_2```. 

### movingarea.m
```T = movingarea(area_ratio)``` creates the transfer matrix accounting for the volume velocity displaced by the moving area at an area discontinuity. If the first pipe area is ```s_1``` and the second pipe area is ```s_2```, area ratio should be ```s_1/s_2```. 

### realmicrophone.m 
This script produces the symbolic transfer matrices for the back volume (```B```) and the whole microphone (```T```) as described in Walsh 2021 (under review as of 11/08/2021). 
Afterward, it produces the analytical function ```P_a``` which is the pressure related acceleration sensitivity of the microphone (in ```Pa/9m/s^2)```) and compares that to the first order truncated Taylor series approximation of ```P_a``` (which implies P_a is proporional to ```L_1```, ```L_2``` and ```L_3``` if the lines match up) 
