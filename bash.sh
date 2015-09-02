#! /bin/bash

echo "compiling..."

g++ simulation.cpp mt64.c rando3.c -lm
echo "running the program..."
./a.out

