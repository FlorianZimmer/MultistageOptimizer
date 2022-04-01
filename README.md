# Mutlistage optimizer for rockets

Calculate optimal mass for each stage given a deltaV requirement on a multistage rocket.

Finds the lowest total mass for a rocket by optimizing the weight distribution between all stages, without altering total deltaV of the rocket.

Rocket and engines are specified in a json.\
Sample files for each are in the config folder.

Mathematical background, which was also the source for this project can be found here:\
http://www.projectrho.com/public_html/rocket/multistage.php

## Config

Options can be found in config/config.json.

useMultiCore	-> if program should utilize all available cores.

precision 		-> the number specifies the combined steps the program should take in splitting the deltaV;\
higher equals better, but cant go above 255 due to implementing the array that holds these numbers as chars.

i.e. for a two stage Rocket with a precision of 150 the first calculation will be:\
1st Stage 1/150 of total deltaV\
2nd Stage 149/150 of total deltaV
                                                                                     
maxRAM 				-> maximum allowed RAM usage in bytes. Will warn user if defined parameters will cause more RAM usage than specified here.
  
verbose 			-> should only be used for debugging, as it slows down the programm alot and in single-threaded mode, as the output isnt orderd when using multithreading.
  
enginesPath 	-> path to json where engines are defined.
  
rocketPath		-> path to json where rocket is defined.

## Limitations

Upwards of six stages we get into realms of impossible amount of necessary RAM, with higher precisions, as the number of different distributions is calculated with n choose r where n is the precision - 1 and r is the number of stages - 1.\
To mitigate the effect, you can lower the precision.
You can try out what what number of distributions exist here: https://www.calculatorsoup.com/calculators/discretemathematics/combinations.php

Currently the program assumes that the engines always work at their vaccuum efficiency.
I tried using the sl isp of engines for the first stage only. But that always leads to a tiny first stage, because it is so "inefficient" and therefore shouldnt be very big according to the program. Because of that, the first stage engine burns for only a very short time, causing the second stage to also burn near the surface.

## KSP Mods

For using this easily and read out all necessary information in KSP, one should use:
- Kerbal Engineer Redux -> read out total mass and deltaV
- Real Fuels -> dry and wet mass of tank
- Procedural Parts -> easily and very finely adjust the size of tanks
