#!/bin/bash

for f in *	#loops through all the files in the working directy
do
	ANFILE="$f"
	echo $ANFILE
	if [[ "$ANFILE" == *".gff"* ]];then
		python3 /home/leann/lib/ortho_finder_butterflies/scripts/gff_longestCDSinfo.py "$ANFILE" "${ANFILE}.CDS.output"

	fi

	if [[ "$ANFILE" == *".gtf"* ]];then
		python3 /home/leann/lib/ortho_finder_butterflies/scripts/gtf_longestCDSinfo.py "$ANFILE" "${ANFILE}.CDS.output"
	fi
done
