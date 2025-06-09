#!/bin/bash

echo '\nGerando os compilaveis com o CMake'
apptainer exec --nv /opt/apptainer/ml-verse_latest.sif cmake -B build

echo '\nCompilando com o Make'
cd build
apptainer exec --nv /opt/apptainer/ml-verse_latest.sif make
cd ..