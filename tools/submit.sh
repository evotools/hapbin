#!/bin/bash

sed "s/500/$1/g" hapbin.500.$2.pbs > hapbin.$1.$2.pbs.tmp
mv hapbin.$1.$2.pbs.tmp hapbin.$1.$2.pbs
cat hapbin.$1.$2.pbs
qsub hapbin.$1.$2.pbs

