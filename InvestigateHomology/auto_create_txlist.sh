#!/bin/bash

for f in *	#loops through all the files in the working directy
do
	ANFILE="$f"
	if [[ "$ANFILE" == *".CDS.output"* ]];then
		echo "$ANFILE"
		python3 /home/leann/lib/ortho_finder_butterflies/scripts/maps_getlongestcolumn.py "$ANFILE" "${ANFILE}.LIST"

	fi

#	if [[ "$ANFILE" == *".gtf"* ]];then
#		python3 /home/leann/lib/ortho_finder_butterflies/scripts/gtf_longestCDSinfo.py "$ANFILE" "${ANFILE}.CDS.output"
#	fi
done
