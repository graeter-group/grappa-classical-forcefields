#!/bin/bash

# Replace 'desired_file.txt' with the name of the file you are checking for
desired_file="topol.top"
searchdir='datasets/ss-interaction/'

# Find all directories in the current directory and its subdirectories
# Then use grep to filter out directories that contain the desired file
# Finally, use dirname to extract the directory names from the file paths
find $searchdir -maxdepth 1 -type d \
  | while read dir; do
      if ! test -e "$dir/$desired_file"; then
          echo "$dir"
      fi
    done

