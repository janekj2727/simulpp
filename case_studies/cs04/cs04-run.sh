#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 04 – Nitrogen in PBC"
    echo "simul++"
    echo "Test of various integrators with SHAKE"
    echo "For more information see the cs04-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.leg

for integrator in verlet k3m2e k4m2e k4m4 k5m5 k5m4e; do
    cp cs04.config cs04${integrator}.config
    cp cs04.field cs04${integrator}.field
    rm results/${integrator}.slope
    rm results/${integrator}.errc
    for tstep in 0.01 0.005 0.0025 0.002 0.001; do
        cat cs04.control | sed 's|integrator [km1-9verlet]*|integrator '"$integrator"'|g' | sed 's|timestep [0-9.]*|timestep '"$tstep"'|g' | sed 's|steps [0-9]*|steps '"$(ev "(10/$tstep)")"'|g' >cs04${integrator}.control
        simul++ cs04${integrator}
        # echo "run" > gdb.cmds
        # gdb -batch -x gdb.cmds --args simul++ "cs04${integrator}"
        if [[ "$?" != "0" ]]; then
            exit
        fi
        cp cs04${integrator}.cpa results/cs04${integrator}_${tstep}.cpa
        if [[ "$tstep" == "0.002" ]]; then
            cp cs04${integrator}.cpa results/cs04${integrator}_toplot.cpa
        fi
        cd results
        tail -n +10 cs04${integrator}_${tstep}.cpa >pokus.cpa
        if [[ "$(tabproc "c1" <pokus.cpa | sumetc -e)" == "0" ]]; then
            echo "$tstep 0 0 $(ev "log($tstep)") -12" >>${integrator}.slope
        else
            echo $tstep $(eng_slope.sh pokus.cpa 1 10) | tabproc "A" "B" "C" "log(abs(A))" "log(abs(B))" >>${integrator}.slope
        fi
        if [[ "$(tabproc "c7" <pokus.cpa | tail -n1)" == "0" ]]; then
            echo "$tstep 0 $(ev "0-log($tstep)") 0" >>${integrator}.errc
        else
            echo $tstep $(tabproc "c7" <pokus.cpa | tail -n1) | tabproc "A" "B" "log(abs(A))" "log(abs(B))" >>${integrator}.errc
        fi
        # plot pokus.cpa:10:7
        rm pokus.cpa
        cd ..
    done
    rm -f cs04${integrator}.*
    echo "$integrator.slope:4:5" >>results/integrators.toplot
    echo "$integrator.errc:3:4" >>results/integrators.toplot2
    echo "cs04${integrator}_toplot.cpa:10:1" >>results/integrators.toplot3
    echo "$integrator" >>results/integrators.leg
done

cd results

plot.py $(cat integrators.toplot | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \left|\frac{\partial E_\mathrm{tot}}{\partial t}\right|$'] -p["lower%right"] -e["nitrogen-eng_slope.eps"] -f14
plot.py $(cat integrators.toplot2 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \frac{\left| \Delta l \right|}{l}$'] -p["lower%right"] -e["nitrogen-maxErrC.eps"] -f14
plot.py $(cat integrators.toplot3 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$t$'] -y['$E_\mathrm{tot}$'] -p["lower%right"] -e["nitrogen-Etot_vs_time.eps"]
cd ..

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

exit
