#!/bin/bash 

sources=("modules.f" "grsm.f" "aermod.f" "setup.f" "coset.f" "soset.f" "reset.f" "meset.f" "ouset.f" "inpsum.f" "metext.f" "iblval.f" "siggrid.f" "tempgrid.f"
	 "windgrid.f" "calc1.f" "calc2.f" "prise.f" "arise.f" "prime.f" "sigmas.f" "pitarea.f" "uninam.f" "output.f" "evset.f" "evcalc.f" "evoutput.f"
	 "rline.f" "bline.f" )

for file in "${sources[@]}"; do 
	echo "compiling $file"
	gfortran -c -fbounds-check -Wuninitialized -O2 -static -w $file
done

echo "COMPILATION COMPLETE!"

gfortran -o aermod -static -O2 $(find . -iregex ".*\.o")

echo "LINKING COMPLETE"

rm *.o
rm *.mod

echo "INTERMEDIATE FILES DELETED"

