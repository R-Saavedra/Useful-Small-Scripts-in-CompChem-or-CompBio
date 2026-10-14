In this repository you will find some scripts that I have written over the years that have helped me out to analyze ouput files or generate input files. 

The dist_scan.py generates geometry conformations of a SCAN between 2 points increasing or decreasing its distance, which is pretty useful to run manual SCANs in Gaussian/ORCA if the overall algorithm of these softwares is not what you're looking for. For me it has been particularly useful as they've changed the G09 algorithm in the new G16.

The extract_optimized_coordinates is a simple shell script that prints the geometry of an optimization or the geometries of a optimized SCAN in a Gaussian output (.out/.log).

The manual_counterpoise_input_generator generates Gaussian inputs that correspond to a counter-poise calculation. It is useful to include solvation effects as the current Counterpoise implementation in Gaussian cannot deal with implicit solvents.

Hope you find them useful as well!! Feel free to comment if you have any doubt.
