#!/bin/bash

if [[ "$1" == "-h" ]] 
then
    echo "Case study 01 – Harmonic oscillator"
    echo "simul++"
    echo "Test various integrators in NVE"
    echo "For more information see the cs01-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.leg

for integrator in verlet k3m2e k4m2e k4m4 k5m4e k5m5 
do
    cp cs01.config cs01${integrator}.config
    cp cs01.field cs01${integrator}.field
    cat cs01.control | sed 's|integrator [km1-9]*|integrator '"$integrator"'|g' > cs01${integrator}.control
    simul++ cs01${integrator}
    cp cs01${integrator}.cpa results/
    rm -f cs01${integrator}.*
    echo "cs01$integrator.cpa:0:1" >> results/integrators.toplot
    echo "$integrator" >> results/integrators.leg
done

cd results

plot.py $(cat integrators.toplot | oneline ) -l[$(cat integrators.leg | oneline | repl ' ' ',')]
cd ..

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

exit
