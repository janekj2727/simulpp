#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Program for executing scripts using mymath class from simul++ package"
    echo "Call by:"
    echo "    $0 [script_name] [arguments_for_script]"
    echo "    where:"
    echo "          script_name is a file with a piece of c++ code, if none is given then stdin is read"
    echo "          arguments_for_script are used for args, argv in the c++ code"
fi

mainfile="$SIMULHOME/mymath/src/mymath_main.cpp"
echo "// This file was created automatically with $0" >"$mainfile"

cat <<EOT >>"$mainfile"
#include "math_utils.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "ARMA.hpp"
#include "NewtonMethod.hpp"
#include "LinearSystem.hpp"
#include "LinearInterpolator.hpp"
#include "Linear3PInterpolator.hpp"
#include "Linear4PInterpolator.hpp"
#include "HermiteCubSplines.hpp"
#include "NaturalCubSplines.hpp"
#include "MacsimusQuadSplines.hpp"
#include "MacsimusHyperbSplines.hpp"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <ctime>
#include <chrono>

int main(int argc, char **argv)
{
EOT

if [[ $# -gt 0 ]]; then
    input_script="$1"
    shift
    while read -r line; do
        echo "$line" >>"$mainfile"
    done <$input_script
else
    while read -r line; do
        echo "$line" >>"$mainfile"
    done
fi

echo 'return 0;' >>"$mainfile"
echo '}' >>"$mainfile"

currdirr=$(pwd)

cd $SIMULHOME/build-debug

make -j4 mymath_script 1>&2

cd "$currdirr"

$SIMULHOME/bin-debug/mymath_script "$@"


