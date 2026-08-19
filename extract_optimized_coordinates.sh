#!/bin/bash

## This script extracts the geometries of optimized structures from nosymm scan using gaussian
#; Using nosymm in gaussian can be benefitial though a standard orientation is not generated
#; thus difficulting the job when extracting the coordinates. Here the optimized parameters
#; section is searched and the input geometry of that specific converged structure is obtained.

#! The script takes std-in-1 as the gaussian file and std-in-2 as the total number of atoms in the system

scan_flag=false
output_flag=false
usage() {
cat << EOF
Usage:: bash extract_optimized_coordinates.sh [-s] [-o OUTPUT_NAME.xyz] GAUSSIAN_OUTPUT.out NUMBER_OF_ATOMS
        -s              Select this flag to obtain only optimized coordinates from scan.
			This will generate a directory named '{GAUSSIAN_OUTPUT}_EXTRACT_GEOM/'
			that contains all geometries.

        -h              Displays this message.
EOF
exit 1
}

while getopts "s:o:" opt; do
        case "${opt}" in
                s)
                        scan_flag=true
                        ;;
		o)
			output_name="$OPTARG"	
			output_flag=true
			;;
                *)
                        usage
                        exit 0
                        ;;
        esac
done
shift $((OPTIND - 1))

if [ "$scan_flag" == "true" ]; then
	number_of_geometries=$(grep "Optimized Parameters" $1 | wc -l)
	#; This command assigns to the variable number_of_geometries the quantity of times that 
	#; optimized parameters appear within the gaussian file. Simply using grep to output
	#; all lines where "PATTERN" shows and then counting the number of lines with "wc -l"
	
	mkdir ${1%%.*}_EXTRACT_GEOM
	
	#! A Loop has to be run from 1 to the last geometry of the scan
	#! This corresponds to the number of scan geometries +1 being the original geometry
	i=1
	while [ $i -le $number_of_geometries ]; do
		#?echo $i
		starting_line=$(grep -n "Optimized Parameters" $1 | cut -d: -f1 | awk -v "i=$i" 'NR==i')
		#; This command will $( Search "PATTERN" and show each line and its number within the file
		#; | only use the number of the line | only use the i^th number of the line)
		#?echo $starting_line
		geometry_line=$(head -n $starting_line $1 | grep -n "Input orientation" | cut -d: -f1 | tail -n 1)
		#; Now the guassian file is $(only read until the Optimized Parameters line | search for "Input orientation"
		#; PATTERN and show all lines containing it plus its number within the guassian file | take only the number
		#;  use only the last number) this the corresponding converged geometry is obtained.
		#?echo $geometry_line
		geometry_start=$((geometry_line + 5))
		#?echo $geometry_start
		geometry_end=$((geometry_start + $2 - 1))
		#?echo $geometry_end
		#; From how gaussian output files are written the geometry starts 4 lines (5 lines if counting the Input one) 
		#; after the "Input orientation" PATTERN the last line of the geometry would then correspond to the starting line 
		#; + the number of atoms in the system
		sed -n "${geometry_start},${geometry_end}p" $1 | awk '{ print $2, $4, $5, $6 }' > ${1%%.*}_EXTRACT_GEOM/${1%%.*}-${i}.xyz
		{ echo -e "$2\nExtracted optimized geometry from $1 scan at iteration $i" ; cat ${1%%.*}_EXTRACT_GEOM/${1%%.*}-${i}.xyz; } > tmp.t && mv tmp.t ${1%%.*}_EXTRACT_GEOM/${1%%.*}-${i}.xyz
		#; These 2 command lines 1st: print the gaussian lines that contain the optimized geometry | only prints lines containing
		#; atomic number - x - y - z ;; 2nd: append "the text" to the beggining of the file, that's so it is in xyz format.
		((i++))
	done
else
	#!?[[ $(grep -c "Optimized Parameters" $1) -eq 1 ]] || ( echo "The number of 'Optimized Parameters' is not exactly 1. Aborting..." ; exit 1 )
	#?echo "Passed the opt test"
	starting_line=$(grep -n "Optimized Parameters" $1 | cut -d: -f1)
	if grep -q "Standard orientation" $1 ; then  
		geometry_line=$(tail -n +"$starting_line" $1 | grep -n "Standard orientation" | tail -n 1 | cut -d: -f1)
	else
		geometry_line=$(tail -n +"$starting_line" $1 | grep -n "Input orientation" | tail -n 1 | cut -d: -f1)
	fi
	geometry_start=$((starting_line + geometry_line + 4))
	#?echo "Geometry start is located at: $geometry_start"
	geometry_end=$((geometry_start + $2 - 1))
	#?echo "Geometry end is located at: $geometry_end"
	if [ "$output_flag" == "true" ]; then
		[[ "$output_name" == *.xyz ]] || output_name="${output_name}.xyz"
		sed -n "${geometry_start},${geometry_end}p" $1 | awk '{ print $2, $4, $5, $6 }' > $output_name
		{ echo -e "$2\nFrom $1: $(grep "SCF Done" $1 | tail -n 1 | awk '{print $5}') in Hartree" ; cat $output_name; } > tmp.t && mv tmp.t $output_name 
	else
		sed -n "${geometry_start},${geometry_end}p" $1 | awk '{ print $2, $4, $5, $6 }' > ${1%%.*}.xyz
		{ echo -e "$2\nFrom $1: $(grep "SCF Done" $1 | tail -n 1 | awk '{print $5}') in Hartree" ; cat ${1%%.*}.xyz; } > tmp.t && mv tmp.t ${1%%.*}.xyz
	fi
fi

exit 0

