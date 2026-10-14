In this repository you will find some scripts that I have written over the years that have helped me out to analyze ouput files or generate input files. While I try to build myself the scripts, it is undeniably a time-consuming task and I have realied on the usage of agentic AI to improve or complete some scripts [Marked as "Help from AI"].

The dist_scan.py generates geometry conformations of a SCAN between 2 points increasing or decreasing its distance, which is pretty useful to run manual SCANs in Gaussian/ORCA if the overall algorithm of these softwares is not what you're looking for. For me it has been particularly useful as they've changed the G09 algorithm in the new G16. [Help from AI]

The extract_optimized_coordinates is a simple shell script that prints the geometry of an optimization or the geometries of a optimized SCAN in a Gaussian output (.out/.log).

The manual_counterpoise_input_generator generates Gaussian inputs that correspond to a counter-poise calculation. It is useful to include solvation effects as the current Counterpoise implementation in Gaussian cannot deal with implicit solvents. [Help from AI]

Hope you find them useful as well!! Feel free to comment if you have any doubt.
