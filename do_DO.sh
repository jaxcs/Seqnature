#! /bin/bash

# A sample script for creating an individualized reference genome
# and annotation gtf file.

date
echo Building Genomes...

ID=$1
PFX=$2
OUT=$3
if [[ ${#TMPDIR} == 0 ]]; then
    TMPDIR=$OUT
else
    TMPDIR=${TMPDIR}/
fi

if [[ ${#OUT} > 0 ]]; then
    rm -rf $OUT
    mkdir -p $OUT
fi

./build_individualized_genome.py \
    -c ${PFX}${ID}.founder.blocks.L.csv \
    -r NCBIM37_genome \
    -i 20111102-indels-all.annotated.vcf \
    -s 20111102-snps-all.annotated.vcf ${ID} L \
         > ${TMPDIR}${ID}_L.fa 2> ${OUT}${ID}_L_stderr.txt &

./build_individualized_genome.py \
    -c ${PFX}${ID}.founder.blocks.R.csv \
    -r NCBIM37_genome \
    -i 20111102-indels-all.annotated.vcf \
    -s 20111102-snps-all.annotated.vcf ${ID} R \
        > ${TMPDIR}${ID}_R.fa 2> ${OUT}${ID}_R_stderr.txt &

wait

# Bring everything together into one .fa
cat ${TMPDIR}${ID}_{L,R}.fa > ${OUT}${ID}.fa &

date
echo Done building genomes.  Now on to the gtfs...

# Build the DO-specific gtfs 
./adjust_individualized_annotations.py \
    -c 1 -t 3 -s 4 -e 5 -n 5 \
    -o ${TMPDIR}${ID}_L.gtf \
    -x ${PFX}${ID}.founder.blocks.L.csv \
    Mus_musculus.NCBIM37.67.gtf ${ID} L \
        > ${OUT}${ID}_L_adjust_stdout.txt 2> ${OUT}${ID}_L_adjust_stderr.txt &

./adjust_individualized_annotations.py \
    -c 1 -t 3 -s 4 -e 5 -n 5 \
    -o ${TMPDIR}${ID}_R.gtf \
    -x ${PFX}${ID}.founder.blocks.R.csv \
    Mus_musculus.NCBIM37.67.gtf ${ID} R \
        > ${OUT}${ID}_R_adjust_stdout.txt 2> ${OUT}${ID}_R_adjust_stderr.txt &
wait

cat ${TMPDIR}${ID}_{L,R}.gtf > ${OUT}${ID}.gtf
mv ${ID}*offsets_chr* ${OUT}

touch ${OUT}/done.marker
date
echo Done.
