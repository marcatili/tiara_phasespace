#!/bin/env bash

prototypefile=$1

Njobs=3
a=56
b=12


for i in $(seq 0 1 $Njobs)
do
	seed1=$(($a+$i))
	seed2=$(($b+$i))	
	#cat main_prototype.mac | sed -e s/"s1"/"$seed1"/g | sed -e s/"s2"/"$seed2"/g | sed s/"NN"/"$i"/g > macro/main_$i.mac
	cat $prototypefile | sed -e s/"s1"/"$seed1"/g | sed -e s/"s2"/"$seed2"/g | sed s/"NN"/"$i"/g > macro/main_$i.mac
done
