#!/bin/bash

cmake -B build -G Ninja
ninja -C build
echo '\nCompilado!\n'