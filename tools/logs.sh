#!/bin/bash

for file in *.o*; do
    mv "$file" "$(basename $file .$(echo $file | rev | cut -d'.' -f1 | rev)).log"
done

rm *.e*

