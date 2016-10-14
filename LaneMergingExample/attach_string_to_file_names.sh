#!/bin/bash
ATTACH=$1

rm ./figures/*Courteous_$ATTACH.png
rm ./OutFiles/MatlabSimfile_Discrete_$ATTACH.txt
rename -v "s/Courteous.png/Courteous_$ATTACH.png/" ./figures/*.png
rename -v "s/MatlabSimfile_Discrete.txt/MatlabSimfile_Discrete_$ATTACH.txt/" ./OutFiles/*.txt