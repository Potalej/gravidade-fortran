#!/bin/bash

sh ./helpers/omp_threads.sh

cmake -B build -G Ninja
ninja -C build
echo '\nCompilado!\n'