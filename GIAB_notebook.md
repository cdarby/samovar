#GIAB TRIO: CHILD

## 0. General

```
samtools flagstat NA24385_GRCh38_phased_possorted_bam.bam
2109450050 + 0 in total (QC-passed reads + QC-failed reads)
134219446 + 0 secondary
0 + 0 supplementary
365163089 + 0 duplicates
2048778634 + 0 mapped (97.12% : N/A)
1975230604 + 0 paired in sequencing
987615302 + 0 read1
987615302 + 0 read2
1682157276 + 0 properly paired (85.16% : N/A)
1899730002 + 0 with itself and mate mapped
14829186 + 0 singletons (0.75% : N/A)
50407834 + 0 with mate mapped to a different chr
38582413 + 0 with mate mapped to a different chr (mapQ>=5)
 
samtools view ../NA24385_GRCh38_phased_possorted_bam.bam | grep -c "HP"
1412038957
```
1412038957/2048778634 = 69% of short reads are phased by Longranger

samtools calmd

## 1. bamsurgeon 

Generate varfile for bamsurgeon: random uniformly spaced variant sites that are not in the un-filtered VCF from longranger

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/0.setup/generateVarfile.py --fai /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa.fai --out out.varfile --vcf ../NA24385_GRCh38_phased_variants.vcf --variantspacing 90000

wc -l out.varfile
34457 out.varfile
```

Picard tools needs this version of java  
`module load java/1.8.0_45`

Run bamsurgeon on the whole file at once (slurm); code modified to only mutate reads with HP=1 tag  

```
sbatch --job-name=bs --partition=lrgmem --nodes=1  --ntasks-per-node=48 --time=48:0:0 --export=ALL --workdir=/work-zfs/mschatz1/cdarby/giabtrio/NA24385_GRCh38/bamsurgeon --mem=400G --wrap="python /work-zfs/mschatz1/cdarby/bamsurgeon/bin/asddsnv.py -v /work-zfs/mschatz1/cdarby/giabtrio/NA24385_GRCh38/bamsurgeon/out.bed -f /work-zfs/mschatz1/cdarby/giabtrio/NA24385_GRCh38/NA24385_GRCh38_phased_possorted_bam_GRCh38_calmd.bam -r  /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa -o /work-zfs/mschatz1/cdarby/giabtrio/NA24385_GRCh38/bamsurgeon/bamsurgeon.bam --picardjar /work-zfs/mschatz1/resources/picard.jar --minmutreads 1 --ignoresnps --force --tagreads --ignorepileup --aligner mem --seed 0 -p 48 --skipmerge" 

ls addsnv.tmp | grep -v "bai" > BAMFILES.txt  
wc -l BAMFILES.txt  
28221 BAMFILES.txt  
```
Made 28,221 mutations

Process and merge the reads with a variant of bamsurgeon's `replacereads.py` that preserves original file's AS and XS tags.

```
samtools calmd --reference /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa -b MERGED.bam > MERGE_calmd_grch38.bam
python /work-zfs/mschatz1/cdarby/bamsurgeon/bamsurgeon/markreads.py MERGE_calmd_grch38.bam MERGE_calmd_grch38_BStag.bam
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/replacereads.py --bam ../ --replacebam MERGE_calmd_grch38_BStag.bam --outputbam MERGE_ALL.bam --all --progress --keepsecondary --keepsupplementary
```
Sort and index the merged file



## 3. How much training data?

Generate varfile for edited-sites simulator: random uniformly spaced variant sites that are not in the un-filtered VCF from longranger

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/0.setup/generateVarfile.py --fai /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa.fai --out 10K.varfile --vcf ../NA24385_GRCh38_phased_variants.vcf --variantspacing 10000

wc -l 10K.varfile
310438 10K.varfile
```

Simulate mosaicism at the sites

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/1.simulate/simulate.py --varfile 10K.varfile --simulate --bam ../NA24385_GRCh38_phased_possorted_bam.bam --nproc X > mosaic.features.tsv
wc -l mosaic.features.tsv
   251630 mosaic.features.tsv
```

Calculate feature vectors at called het. sites

```
grep "PASS" ../NA24385_GRCh38_phased_variants.vcf | grep "0|1\|1|0" > het.vcf
wc -l het.vcf
2799134 het.vcf
shuf -o het.vcf het.vcf

python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/1.simulate/simulate.py --varfile het.vcf --max 100000 --bam ../NA24385_GRCh38_phased_possorted_bam.bam --nproc X > het.features.tsv

wc -l het.features.tsv
   173005 het.features.tsv
```

Calculate feature vectors at called hom. sites

```
grep "PASS" ../NA24385_GRCh38_phased_variants.vcf | grep "0|0\|1|1" > hom.vcf
wc -l hom.vcf
1785278 hom.vcf
shuf -o hom.vcf hom.vcf

python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/1.simulate/simulate.py --varfile hom.vcf --max 100000 --bam ../NA24385_GRCh38_phased_possorted_bam.bam --nproc X > hom.features.tsv
wc -l hom.features.tsv
   278261 hom.features.tsv
```

At simulated mosaic and germline sites, subsample a fraction of the reads. ALL = 86X (estimated), so 0.23 = 20X; 0.35 = 30X; 0.47=40X; 0.58=50X; 0.7=60X; 0.81=70X

```
python ../../../HCC1954T/pipeline/5.experiments/simulateSubsample.py --nproc X --simulate --subsample 0.23 --varfile 10K.varfile --bam ../NA24385_GRCh38_phased_possorted_bam_GRCh38_calmd.bam > 20X/mosaic.features.tsv

   179798 20X/mosaic.features.tsv
   210566 30X/mosaic.features.tsv
   226678 40X/mosaic.features.tsv
   235884 50X/mosaic.features.tsv
   242518 60X/mosaic.features.tsv
   246605 70X/mosaic.features.tsv

python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/simulateSubsample.py --varfile het.vcf --subsample 0.23 --bam ../NA24385_GRCh38_phased_possorted_bam_GRCh38_calmd.bam --nproc X --max 100000 > 30X/het.features.tsv&

python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/simulateSubsample.py --varfile hom.vcf --subsample 0.23 --bam ../NA24385_GRCh38_phased_possorted_bam_GRCh38_calmd.bam --nproc X --max 100000 > 30X/hom.features.tsv&
```

Give classifier various amounts of training data and calculate AUC of TPR/FPR curve (ranked by score) on cross-validation data (equal number of positive/negative examples in training and validation sets)

```
for COV in `seq 20 10 70`; do for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasizeUnphased.py --mosaic ${COV}X/mosaic.features.tsv --het ${COV}X/het.features.tsv --hom ${COV}X/hom.features.tsv --outpfx ${COV}X/ ; for NUM in `seq 1000 1000 20000`; do Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R ${COV}X/${NUM}.txt | cut -f 2 -d" "  >> ${COV}X/AUC${TRIAL}.txt ; done; done; done

for COV in `seq 20 10 70`; do for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --mosaic ${COV}X/mosaic.features.tsv --het ${COV}X/het.features.tsv --hom ${COV}X/hom.features.tsv --outpfx ${COV}X/P ; for NUM in `seq 1000 1000 20000`; do Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R ${COV}X/P${NUM}.txt | cut -f 2 -d" "  >> ${COV}X/AUCP${TRIAL}.txt ; done; done; done


for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasizeUnphased.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv ; for NUM in `seq 1000 1000 20000`; do Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R ${NUM}.txt | cut -f 2 -d" "  >> AUC${TRIAL}.txt ; done; done

for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv ; for NUM in `seq 1000 1000 20000`; do Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R ${NUM}.txt | cut -f 2 -d" "  >> AUC${TRIAL}.txt ; done; done
```

## 4. Paired 

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/pairedTagBam.py --out paired.bam --vcf ../NA24385_GRCh38_phased_variants.vcf --bam ../bamsurgeon/MERGE_ALL.bam

```
Reads written: 2048811406   
Number of ends where other end in conflict or multiple overlapped variants in conflict: 2470047   
No haplotype assigned: 1824403923  

11% of reads phased by overlapping at least one het. variant (self or mate)

**263145 calls**  
![](testdata_DO_NOT_COMMIT/giab/pairedend/multi_pairedonly.png)

**617538 calls**  
![](testdata_DO_NOT_COMMIT/giab/pairedend/multi_nophasing.png)

**617538 calls, of which 263145 paired**
![](testdata_DO_NOT_COMMIT/giab/pairedend/multi_hybridscore.png)


**103792 calls**
![](testdata_DO_NOT_COMMIT/giab/bamsurgeontrain/multi_repeatfilter.png)


**28733 calls**
![](testdata_DO_NOT_COMMIT/giab/bamsurgeontrain/multi_featurefiltered.png)

## MosaicHunter

Mostly default configurations (min reads 5/5%; depth 25-150; no CNV; no blat)
```
java -jar ~/scratch/MosaicHunter/build/mosaichunter.jar -C genome.properties 
java -jar ~/scratch/MosaicHunter/build/mosaichunter.jar -C trio.properties
```