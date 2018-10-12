#!/bin/bash

for i in `seq 1001 1200`;
do
#echo $i
./CT.intel 1 $i;
done 
