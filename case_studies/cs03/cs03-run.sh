#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 03 – Argon in PBC"
    echo "simul++"
    echo "Test of various integrators without SHAKE"
    echo "For more information see the cs03-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.leg

for integrator in verlet k3m2e k4m2e k4m4 k5m5 k5m4e; do
    cp cs03.config cs03${integrator}.config
    cp cs03.field cs03${integrator}.field
    rm results/${integrator}.slope
    rm results/${integrator}.errc
    for tstep in  0.01 0.005 0.0025 0.002 0.001; do
        cat cs03.control | sed 's|integrator [km1-9verlet]*|integrator '"$integrator"'|g' | sed 's|timestep [0-9.]*|timestep '"$tstep"'|g' | sed 's|steps [0-9]*|steps '"$(ev "(10/$tstep)")"'|g' >cs03${integrator}.control
        simul++ cs03${integrator}
        # echo "run" > gdb.cmds
        # gdb -batch -x gdb.cmds --args simul++ "cs03${integrator}"
        if [[ "$?" != "0" ]]; then
            exit
        fi
        cp cs03${integrator}.cpa results/cs03${integrator}_${tstep}.cpa
        cd results
        tail -n +10 cs03${integrator}_${tstep}.cpa >pokus.cpa
        if [[ "$(tabproc "c1" <pokus.cpa | sumetc -e)" == "0" ]]; then
            echo "$tstep 0 0 $(ev "0-log($tstep)") 0" >>${integrator}.slope
        else
            echo $tstep $(eng_slope.sh pokus.cpa 1 7) | tabproc "A" "B" "C" "0-log(abs(A))" "0-log(abs(B))">>${integrator}.slope
        fi
        rm pokus.cpa
        cd ..
    done
    rm -f cs03${integrator}.*
    echo "$integrator.slope:4:5" >> results/integrators.toplot
    echo "$integrator" >> results/integrators.leg
done


cd results

plot.py $(cat integrators.toplot | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$-\ln h$'] -y['$-\ln \frac{\partial E}{\partial t}$'] -p["lower%right"] -e["eng_slope.eps"]

cd ..

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

exit
