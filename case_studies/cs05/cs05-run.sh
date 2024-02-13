#!/bin/bash

if [[ "$1" == "-h" ]]; then
    echo "Case study 05 – Nitrogen in PBC (bigger system)"
    echo "simul++"
    echo "Test of medium-size system with SHAKE"
    echo "For more information see the cs05-readme.md file"
fi

mkdir results
rm results/integrators.toplot
rm results/integrators.leg

integrator_list="verlet k3m2e k4m2e k4m4 k5m5 k5m4e"

for integrator in $integrator_list; do
    cp cs05.config cs05${integrator}.config
    cp cs05.field cs05${integrator}.field
    cat cs05.control | sed 's|integrator [km1-9verlet]*|integrator '"$integrator"'|g'  >cs05${integrator}.control
    simul++ cs05${integrator} & # run parallel with multiple integrators
done
wait 
for integrator in $integrator_list; do
    cp cs05${integrator}.cpa results/cs05${integrator}.cpa
    rm -f cs05${integrator}.*
    echo "$integrator" >>results/integrators.leg
    echo "cs05$integrator.cpa " >> results/integrators.toplot
done
echo " " >> results/integrators.toplot

rm *.cpa *.mol *.plb *~ *.prt *.configfinal

echo ""
echo "--------------------------------------------------------------------"
echo "To show results use (substitute N by column number):"
echo "    - MACSIMUS plot:"
echo '            cd results; plot :0:N $(cat integrators.toplot | oneline); cd ..'
echo "    - PYTHON plot.py (my script):"
echo '            cd results; plot.py  $(cat integrators.toplot | oneline | repl '"'  ' ':0:N ') "'-l[$(cat integrators.leg | oneline | repl '"' ' ',')] ; cd .."
exit
