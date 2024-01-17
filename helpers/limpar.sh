#!/bin/sh -x

find . -type d -name "*_*" |\
while read d
do
	cd $d
	rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt core Makefile gravidade
	cd -
done
