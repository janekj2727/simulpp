#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 08 – Nitrogen in PBC for thermostat tests"
    echo "simul++"
    echo "Test of universal integration scheme used for NPT Nosé–Hoover (MTK)"
    echo "For more information see the cs08-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.toplot2
rm results/integrators.toplot3
rm results/integrators.leg

loop_function(){ # defined as a function to be run in parallel
    local integ="$1"
    local integrator=$(echo $integ | sed 's| |_|')
    cp cs08.config cs08${integrator}.config
    cp cs08.field cs08${integrator}.field
    rm results/${integrator}.slope
    rm results/${integrator}.errc
    mkdir results/${integrator}
    local tstep
    for tstep in 0.01 0.005 0.0025 0.002 0.001; do
    # for tstep in 0.002; do
        cat cs08.control | sed 's|integrator [km1-9verletrp]*|integrator '"$integ"'|g' | sed 's|timestep [0-9.]*|timestep '"$tstep"'|g' | sed 's|steps [0-9]*|steps '"$(ev "(10/$tstep)")"'|g' >cs08${integrator}.control
        simul++ cs08${integrator}
        # echo "run" > gdb.cmds
        # gdb -batch -x gdb.cmds --args simul++ "cs08${integrator}"
        if [[ "$?" != "0" ]]; then
            exit
        fi
        cp cs08${integrator}.cpa results/${integrator}/cs08${integrator}_${tstep}.cpa
        if [[ "$tstep" == "0.002" ]]; then
            cp cs08${integrator}.cpa results/cs08${integrator}_toplot.cpa
        fi
        cd results/${integrator}
        # echo "${integrator}: $PWD"
        tail -n +10 cs08${integrator}_${tstep}.cpa >pokus${integrator}.cpa
        if [[ "$(tabproc "c1" <pokus${integrator}.cpa | sumetc -e)" == "0" ]]; then
            echo "$tstep 0 0 $(ev "log($tstep)") -12" >>../${integrator}.slope
        else
            echo $tstep $(eng_slope.sh pokus${integrator}.cpa 1 17) | tabproc "A" "B" "C" "log(abs(A))" "log(abs(B))" >>../${integrator}.slope
        fi
        if [[ "$(tabproc "c9" <pokus${integrator}.cpa | tail -n100 | sumetc -M)" == "0" ]]; then
            echo "$tstep 0 $(ev "0-log($tstep)") -9" >>../${integrator}.errc
        else
            echo $tstep $(tabproc "c9" <pokus${integrator}.cpa | tail -n100 | sumetc -M) | tabproc "A" "B" "log(abs(A))" "log(abs(B))" >>../${integrator}.errc
        fi
        # plot pokus${integrator}.cpa:10:7
        rm pokus${integrator}.cpa
        cd ../..
        # echo "${integrator}: $PWD"
    done
    rm -f cs08${integrator}.*
    # rm -f ${integrator}/*
    echo "$integrator.slope:4:5" >>results/integrators.toplot
    echo "$integrator.errc:3:4" >>results/integrators.toplot2
    echo "cs08${integrator}_toplot.cpa:17:1" >>results/integrators.toplot3
    echo "$integ" | sed 's| |+|' >>results/integrators.leg
}

for integ in 'verlet' 'verlet trvp2' 'k3m2e trvp2' 'k4m2e' 'k4m4' 'k5m4e' 'k5m5'; do
    ( loop_function "$integ"; ) & # run in a subshell in background
done

wait
cd results

plot.py $(cat integrators.toplot | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \left|\frac{\partial E_\mathrm{tot}}{\partial t}\right|$'] -p["lower%right"] -e["nptnose-eng_slope.eps"] -f14
plot.py $(cat integrators.toplot2 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$\log h$'] -y['$\log \frac{\left| \Delta l \right|}{l}$'] -p["lower%right"] -e["nptnose-maxErrC.eps"] -f14
plot.py $(cat integrators.toplot3 | oneline) -l[$(cat integrators.leg | oneline | repl ' ' ',')] -x['$t$'] -y['$E_\mathrm{tot}$'] -p["lower%right"] -e["nptnose-Etot_vs_time.eps"]
cd ..

# rm *.cpa *.mol *.plb *~ *.prt *.configfinal
rm *.mol *.plb *~ *.prt *.configfinal
exit
