#!/bin/bash

# Compila
cmake -B build -G Ninja
ninja -C build
echo '\nCompilado!\n'

# Roda
cd build && ./gravidade vi=\"../presets/lemniscata.txt\" && cd ..