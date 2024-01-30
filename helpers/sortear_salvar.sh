#!/bin/bash

# Compila
sh ./helpers/build.sh

# Roda
./gravidade -sv ./presets/condicionar/exemplo.txt
