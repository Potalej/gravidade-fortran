#!/bin/bash

# Compila
cmake -B build -G Ninja
ninja -C build
echo '\nCompilado!\n'

# Roda
cd build && ./gravidade preset=\"../presets/exemplo1.txt\" && cd ..
