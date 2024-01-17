#!/bin/sh -x

docker build -t local:gravidade .

echo "Para acessar o diretorio do gravidade-fortran, use `cd /src`"

docker run -it -v .:/src local:gravidade

