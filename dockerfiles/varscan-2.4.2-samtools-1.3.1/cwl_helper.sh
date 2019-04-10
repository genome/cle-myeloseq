#!/bin/bash

set -o errexit
set -o nounset

if [ $# -lt 3 ]
then
    echo "Usage: $0 [TUMOR_BAM] [NORMAL_BAM] [REFERENCE] [roi_bed?]"
    exit 1
fi

TUMOR_BAM="$1"
NORMAL_BAM="$2"
REFERENCE="$3"
OUTPUT="${HOME}/output"

if [ -z ${4+x} ]
then
    #run without ROI
    java -jar /opt/varscan/VarScan.jar somatic \
        <(/usr/bin/samtools mpileup --no-baq -f "$REFERENCE" "$NORMAL_BAM" "$TUMOR_BAM") \
        $OUTPUT \
        --mpileup 1 \
        --output-vcf
else
    ROI_BED="$4"
    java -jar /opt/varscan/VarScan.jar somatic \
        <(/usr/bin/samtools mpileup --no-baq -l "$ROI_BED" -f "$REFERENCE" "$NORMAL_BAM" "$TUMOR_BAM") \
        $OUTPUT \
        --mpileup 1 \
        --output-vcf
fi
