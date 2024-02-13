#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 07 – Nitrogen in PBC for thermostat tests"
    echo "simul++"
    echo "Test of universal integration scheme used for NVT Nosé–Hoover"
    echo "For more information see the cs07-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.toplot2
rm results/integrators.toplot3
rm results/integrators.leg

for integ in verlet 'verlet trvp2' 'k3m2e trvp2' k4m2e k4m4 k5m5 k5m4e; do
    integrator=$(echo $integ | sed 's| |_|')
    cp cs07.config cs07${integrator}.config
    cp cs07.field cs07${integrator}.field
    rm results/${integrator}.slope
    rm results/${integrator}.errc

    for tstep in 0.01 0.005 0.0025 0.002 0.001; do
        cat cs07.control | sed 's|integrator [km1-9verletrp]*|integrator '"$integ"'|g' | sed 's|timestep [0-9.]*|timestep '"$tstep"'|g' | sed 's|steps [0-9]*|steps '"$(ev "(10/$tstep)")"'|g' >cs07${integrator}.control
        simul++ cs07${integrator}
        # echo "run" > gdb.cmds
        # gdb -batch -x gdb.cmds --args simul++ "cs07${integrator}"
        if [[ "$?" != "0" ]]; then
            exit
        fi
        cp cs07${integrator}.cpa results/cs07${integrator}_${tstep}.cpa
        if [[ "$tstep" == "0.002" ]]; then
            cp cs07${integrator}.cpa results/cs07${integrator}_toplot.cpa
        fi
        cd results
        tail -n +10 cs07${integrator}_${tstep}.cpa >pokus.cpa
        if [[ "$(tabproc "c1" <pokus.cpa | sumetc -e)" == "0" ]]; then
            echo "$tstep 0 0 $(ev "log($tstep)") -12" >>${integrator}.slope
        else
            echo $tstep $(eng_slope.sh pokus.cpa 1 13) | tabproc "A" "B" "C" "log(abs(A))" "log(abs(B))" >>${integrator}.slope
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
    rm -f cs07${integrator}.*
    echo "$integrator.slope:4:5" >>results/integrators.toplot
    echo "$integrator.errc:3:4" >>results/integrators.toplot2
    echo "cs07${integrator}_toplot.cpa:13:1" >>results/integrators.toplot3
    echo "$integ" | sed 's| |+|' >>results/integrators.leg
done

cd results

plot.py $(cat integrators.toplot | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \left|\frac{\partial E_\mathrm{tot}}{\partial t}\right|$'] -p["lower%right"] -e["nvtnose-eng_slope.eps"] -f14
plot.py $(cat integrators.toplot2 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \frac{\left| \Delta l \right|}{l}$'] -p["lower%right"] -e["nvtnose-maxErrC.eps"] -f14
plot.py $(cat integrators.toplot3 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$t$'] -y['$E_\mathrm{tot}$'] -p["lower%right"] -e["nvtnose-Etot_vs_time.eps"]
cd ..

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

exit
