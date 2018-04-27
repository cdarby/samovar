#!/bin/sh

BAMFILE=
VCFFILE=
TMPDIR=./tmp
REFGENO= /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa
RESULTDIR=./results
SPLITREFGENO= /work-zfs/mschatz1/resources/grch38

SAMTOOLS=samtools
CNVNATOR=cnvnator

mkdir -p ${TMPDIR}

if [ ! -f ${REFGENO}.fai ]; then
	${SAMTOOLS} faidx ${REF}
fi

if [ ! -f ${BAMFILE}.bai ]; then
	${SAMTOOLS} index ${BAMFILE}
fi

python pipeline/0.setup/generateVarfile.py --vcf ${VCFFILE} --out ${TMPDIR}/mosaic.varfile --variantspacing 100000 --fai ${REFGENO}.fai

#Mask CNVs
cnvnator -root ${TMPDIR}/out.root -tree ${BAMFILE}
cnvnator -root ${TMPDIR}/out.root -his 100 -d ${SPLITREFGENO}
cnvnator -root ${TMPDIR}/out.root -stat 100
cnvnator -root ${TMPDIR}/out.root -partition 100
cnvnator -root ${TMPDIR}/out.root -call 100 > ${TMPDIR}/cnvnator.tsv
#TODO make py3 compatible
python pipeline/0.setup/cnvnator2bed.py ${TMPDIR}/cnvnator.tsv ${TMPDIR}/cnvnator.bed ${REFGENO}.fai
bedtools complement -g ${REFGENO}.fai -i ${TMPDIR}/cnvnator.bed > ${RESULTDIR}/regions.bed

#No CNV masking
bedtools window -w 1000000 -g ${REFGENO}.fai > ${RESULTDIR}/regions.bed


python pipeline/1.simulate/simulate.py --bam ${BAMFILE} --varfile ${TMPDIR}/mosaic.varfile --simulate --nproc 8 > ${TMPDIR}/mosaic.features.tsv

grep "PASS" ${VCFFILE} | grep "1|0\|0|1" | shuf > ${TMPDIR}/het.vcf
python pipeline/1.simulate/simulate.py --bam ${BAMFILE} --varfile ${TMPDIR}/het.vcf --nproc 8 > ${TMPDIR}/het.features.tsv

grep "PASS" ${VCFFILE} | grep "1|1\|0|0" | shuf > ${TMPDIR}/hom.vcf
python pipeline/1.simulate/simulate.py --bam ${BAMFILE} --varfile ${TMPDIR}/hom.vcf --nproc 8 > ${TMPDIR}/hom.features.tsv

cat ${TMPDIR}/hom.features.tsv ${TMPDIR}/het.features.tsv > ${TMPDIR}/germline.features.tsv


python pipeline/2.train/train.py --germline ${TMPDIR}/germline.features.tsv --mosaic ${TMPDIR}/mosaic.features.tsv --out ${TMPDIR}/clf.pkl


python pipeline/3.classify/classify.py --bam ${BAMFILE} --bed ${RESULTDIR}/regions.bed --nproc 32 --clf ${TMPDIR}/clf.pkl > ${RESULTDIR}/scan.bed

#TODO homopolymer bed file filter

python pipeline/4.filter/computeFilterFeatures.py --bam ${BAMFILE} --bed ${RESULTDIR}/scan.bed
#TODO make python
Rscript pipeline/4.filter/filter.R ${TMPDIR}/filterfeatures.tsv ${RESULTDIR}/filteredsites.tsv #TODO output as bed file
#TODO put in earlier step?
python pipeline/4.filter/linkageFilter.py --bam ${BAMFILE} --bed ${RESULTDIR}/filteredsites.tsv --ref ${REFGENO} --vcfavoid ${TMPDIR}/het.vcf > ${RESULTDIR}/linkage.tsv
#TODO reporting step










