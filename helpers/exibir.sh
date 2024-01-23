#!/bin/bash

# Compila
cmake -B build -G Ninja
ninja -C build
echo '\nCompilado!\n'

# Roda
cd build && ./gravidade exibir=\"../data/lemniscata_teste.csv\" && cd ..