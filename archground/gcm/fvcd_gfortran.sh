#! /bin/bash

###
### fvcd_gfortran.sh
###
### ---> This should be used to make VCD fortran stuff directly accessible in python 
### ---> This script is for gfortran, but it is easy to adapt to your other compilers
### ---> A file fvcd.so should be created
###

### PATHS to NETCDF and VCD software: Adapt these to your local settings
NETCDF="/usr"
wherevcd="/srv/workspace/data/vcd/VCD_2.3/vcd"
version="2.3"


#module load gnu/7.2.0 netcdf4/4.5.0-gfortran72

### LOG FILE
\rm -f fvcd.log

### COPY/PREPARE SOURCES
### perform changes that makes f2py not to fail (i.e. remove all '!' commments in codes)
sed s/"\!\!'"/"'"/g $wherevcd/VCD_var.F90 | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.VCD_var.F90
sed s/"\!\!'"/"'"/g $wherevcd/VCD.F90 | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.VCD.F90
sed s/"\!\!'"/"'"/g $wherevcd/julian.F90 | sed s/"\!'"/"'"/g | sed s/"\!"/"\n\!"/g > tmp.julian.F90

### BUILD THROUGH f2py WHAT IS NECESSARY TO CREATE THE PYTHON FUNCTIONS
\rm -f fvcd.pyf
f2py -h fvcd.pyf -m fvcd tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90 > fvcd.log 2>&1

### BUILD
f2py -c fvcd.pyf -m fvcd tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90 --fcompiler=gnu95 \
  -L$NETCDF/lib -lnetcdff -lnetcdf \
  -lm -I$NETCDF/include \
  --f90flags="-fPIC -ffree-form -ffree-line-length-none" \
  --verbose \
  > fvcd.log 2>&1

### CLEAN THE PLACE
\rm -f tmp.VCD_var.F90 tmp.VCD.F90 tmp.julian.F90

