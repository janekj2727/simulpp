#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 02 – Dumb-bell"
    echo "simul++"
    echo "Test of SHAKE with various integrators"
    echo "For more information see the cs02-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.toplot2
rm results/integrators.toplot3
rm results/integrators.leg

for integrator in verlet k3m2e k4m2e k4m4 k5m5 k5m4e; do
    cp cs02.config cs02${integrator}.config
    cp cs02.field cs02${integrator}.field
    rm results/${integrator}.slope
    rm results/${integrator}.errc
    for tstep in  0.25 0.2 0.125 0.1 0.08 0.05; do
        cat cs02.control | sed 's|integrator [km1-9]*|integrator '"$integrator"'|g' | sed 's|timestep [0-9.]*|timestep '"$tstep"'|g' | sed 's|steps [0-9]*|steps '"$(ev "(50/$tstep)")"'|g' >cs02${integrator}.control
        simul++ cs02${integrator}
        # echo "run" > gdb.cmds
        # gdb -batch -x gdb.cmds --args simul++ "cs02${integrator}"
        if [[ "$?" != "0" ]]; then # if the previous simulation failed
            exit
        fi
        cp cs02${integrator}.cpa results/cs02${integrator}_${tstep}.cpa
        if [[ "$tstep" == "0.1" ]]; then
            cp cs02${integrator}.cpa results/${integrator}_toplot.cpa
        fi
        cd results
        tail -n15 cs02${integrator}_${tstep}.cpa >pokus.cpa
        if [[ "$(tabproc "c1" <pokus.cpa | sumetc -e)" == "0" ]]; then
            echo "$tstep 0 0 $(ev "log($tstep)") -12" >>${integrator}.slope
        else
            echo $tstep $(eng_slope.sh pokus.cpa 1 6) | tabproc "A" "B" "C" "log(abs(A))" "log(abs(B))">>${integrator}.slope
        fi
        if [[ "$(tabproc "c5" <pokus.cpa | tail -n1)" == "0" ]]; then
            echo "$tstep 0 $(ev "log($tstep)") -16" >>${integrator}.errc
        else
            echo $tstep $(tabproc "c5" < pokus.cpa | tail -n1) | tabproc "A" "B" "log(abs(A))" "log(abs(B))">>${integrator}.errc
        fi
        rm pokus.cpa
        cd ..
    done
    rm -f cs02${integrator}.*
    echo "$integrator.slope:4:5" >> results/integrators.toplot
    echo "$integrator.errc:3:4" >> results/integrators.toplot2
    echo "${integrator}_toplot.cpa:6:1" >> results/integrators.toplot3
    echo "$integrator" >> results/integrators.leg
done


cd results

plot.py $(cat integrators.toplot | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \left|\frac{\partial E_\mathrm{tot}}{\partial t}\right|$'] -p["lower%right"] -e["eng_slope.eps"]
plot.py $(cat integrators.toplot2 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \frac{\left| \Delta l \right|}{l}$'] -p["lower%right"] -e["maxErrC.eps"]
plot.py $(cat integrators.toplot3 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$t$'] -y['$E_\mathrm{tot}$'] -p["lower%right"] -e["Etot_vs_time.eps"]

cd ..

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

exit
