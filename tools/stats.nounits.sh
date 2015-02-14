#!/bin/bash

stats() {
    if [ ! -f  h.$1.$2.$3 ]; then
        echo "h.$1.$2.$3 not found"
    else
        cat h.$1.$2.$3 | grep real | cut -f 2 | tr 'm' ' ' | tr 's' ' ' | 
        awk -v snps="$1" -v threads="$2" '
        { 
            time = $1*60*1000+$2*1000; 
            h = int(time/60/60/1000);
            m = int(time/60/1000-h*60);
            s = $2;
            ps = time/489301; 
            pst = ps*threads;
            pss = ps/snps;
            psst = pss*threads;
        } END {
            printf "%s %s %sh%sm%ss %s %s %s %s %s\n", snps, threads, h, m, s, time, ps, pst, pss, psst
        }'
    fi
}

for i in 50 100 200 300 400 500
do
    for j in 1 2 4 6 8 10 12 24 48ht 48 72 96 144 192 240
    do
        stats $i $j log
    done
    echo "-------------------------------"
done

for i in 1000 2000
do
    for j in 1 2 4 6 8 10 12 24 48ht 48 72 96 144 192 240 480
    do
        stats $i $j log
    done
    echo "-------------------------------"
done

