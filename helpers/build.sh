#!/bin/bash

cmake -B build -G Ninja
ninja -C build -t clean
echo '\nCompilado!\n'