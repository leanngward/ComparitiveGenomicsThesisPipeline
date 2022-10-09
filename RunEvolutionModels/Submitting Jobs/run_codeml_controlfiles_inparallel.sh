#!/bin/bash
#$ -pe threads 16

filename=$1
parallel --will-cite /home/leann/software/paml4.9j/bin/codeml {} < $filename
