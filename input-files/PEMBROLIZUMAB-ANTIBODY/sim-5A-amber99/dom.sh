#!/bin/bash

rm dom.txt
touch dom.txt

dom1=$(for i in {1..6601}; do echo -n $i ""; done) 
echo $dom1 >> dom.txt

dom2_1st=$(for i in {6848..10099}; do echo -n $i ""; done)
dom2_2nd=$(for i in {16942..20217}; do echo -n $i ""; done)

dom2=$dom2_1st$dom2_2nd >> dom.txt
echo $dom2 >> dom.txt


dom3=$(for i in {10113..16698}; do echo -n $i ""; done) >> dom.txt
echo $dom3 >> dom.txt
