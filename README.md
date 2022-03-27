# Mutlistage optimizer for rockets

Calculate optimal mass for each stage given a deltaV requirement on a multistage rocket.

Finds the lowest total mass for a rocket by optimizing the weight distribution between all stages, without altering deltaV.

Rocket and engines are specified in a json.\
Sample files for each are in the config folder.

## Config

Options can be found in config/config.json.

useMultiCore	-> if program should utilize all available cores.

precision 		-> the number specifies the combined steps the program should take in splitting the deltaV; higher equals better but cant go above 255 due to implementing the array that holds these numbers as chars.

i.e. for a two stage Rocket with a precision of 150 the first calculation will be:\
1st Stage 1/150 of total deltaV\
2nd Stage 149/150 of total deltaV
                                                                                     
maxRAM 				-> maximum allowed RAM usage in bytes. Will warn user if defined parameters will cause more RAM usage than specified here.
  
verbose 			-> should only be used for debugging, as it slows down the programm alot and in single-threaded mode, as the output isnt orderd when using multithreading.
  
enginesPath 	-> path to json where engines are defined.
  
rocketPath		-> path to json where rocket is defined.
