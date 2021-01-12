#!/bin/bash
rm *.o *.x
gfortran -c -g -O3 stsubs.F90
gfortran -c -g -O3 combinatoric.F90
gfortran -c -g -O3 parse_prim.F90
gfortran -c -g -O3 chem_parser.F90
gfortran chem_parser.o parse_prim.o combinatoric.o stsubs.o -o test_chem_parser.x