#!/bin/sh

# Example script to process CAPRIA or CAPRIA+S data using oxasl. We assume here
# that: 1) the CAPRIA perfusion images have already been calibrated and sit
# within the current directory; 2) there is a structural image (which could be a
# separate MP-RAGE or a CAPRIA+S structural image) processed with fsl_anat in
# ../../T1.anat; and 3) there is a struc2CAPRIA.mat affine matrix to register
# the two in the current directory.
#
# Tom Okell, July 2025

# Make the output directory
mkdir oxasl_out

# Set up custom options CAPRIA for fabber
echo '--capria' > oxasl_out/fabber_opts.txt
echo '--capriafa1=2.0' >> oxasl_out/fabber_opts.txt
echo '--capriafa2=9.0' >> oxasl_out/fabber_opts.txt
echo '--capriatr=0.0091' >> oxasl_out/fabber_opts.txt

# Call oxasl
oxasl  -i asldata_calib \
    --iaf=diff --ibf=rpt --plds=0.1593,0.4869,0.8145,1.1421,1.4697,1.7972 \
    --bolus=1.8 --casl \
    --fslanat=../../T1.anat \
    --struc2asl=struc2CAPRIA.mat \
    --fixbolus \
    -o oxasl_out \
    --spatial-off \
    --save-all \
    --batsd=2.0 \
    --basil-options=oxasl_out/fabber_opts.txt \
    --output-stddev  --overwrite --debug  |  tee oxasl_out/oxasl_log.txt
