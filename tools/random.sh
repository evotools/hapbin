#!/bin/bash
#Usage: random.sh [# of individuals] [# of source individuals]

randrange() {
    awk -v loop=$1 -v range=$2 'BEGIN{
    srand()
    do {
        num = int(rand() * range)
        if (!(num in prev)) {
            print 2*num+1
            print 2*num+2
            prev[num] = 1
            count++
        }
    } while (count<loop)
    }' | sort -n;
}

randrange $1 $2 > numbers.$1.txt

