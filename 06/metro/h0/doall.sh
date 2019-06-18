#!/bin/bash
for i in `cat T.dat`; do  mkdir t$i && cp input.dat t$i && sed -i -e "s/XXX/$i/g" t$i/input.dat  ; done

for i in `cat T.dat`; do cd t$i  ; ../../../ISING_1D/Monte_Carlo_ISING_1D.exe; cd ..   ; done

