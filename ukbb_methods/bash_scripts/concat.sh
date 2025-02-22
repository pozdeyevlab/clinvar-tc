#!/bin/bash

# Read in named command line args 
while getopts ":f:o:" opt; do
  case $opt in
    f) file_list="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Loop through txt files and concat them into one file
while IFS= read -r file
do
   cat $file >> $output
done < $file_list