#HCC1954 PAIR

## 0. General

```
samtools flagstat orig_bams/HCC1954N_calmd.bam
891669864 + 0 in total (QC-passed reads + QC-failed reads)
50455372 + 0 secondary
0 + 0 supplementary
73594732 + 0 duplicates
862920745 + 0 mapped (96.78% : N/A)
841214492 + 0 paired in sequencing
420607246 + 0 read1
420607246 + 0 read2
707643459 + 0 properly paired (84.12% : N/A)
804467151 + 0 with itself and mate mapped
7998222 + 0 singletons (0.95% : N/A)
21011308 + 0 with mate mapped to a different chr
17849315 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view ../orig_bams/HCC1954N_calmd.bam | grep -c "HP"
522328004
```
522328004/862920745 = 60% of short reads are phased by Longranger [Normal]

```
samtools flagstat orig_bams/HCC1954T_calmd.bam
894364684 + 0 in total (QC-passed reads + QC-failed reads)
53290106 + 0 secondary
0 + 0 supplementary
83665156 + 0 duplicates
855013358 + 0 mapped (95.60% : N/A)
841074578 + 0 paired in sequencing
420537289 + 0 read1
420537289 + 0 read2
699787209 + 0 properly paired (83.20% : N/A)
794060443 + 0 with itself and mate mapped
7662809 + 0 singletons (0.91% : N/A)
19581676 + 0 with mate mapped to a different chr
16609727 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view orig_bams/HCC1954T_calmd.bam | grep -c "HP"
485258216
```

485258216/855013358 = 56% of short reads are phased by Longranger [Tumor]


## 1. bamsurgeon

### Normal

Generate varfile for bamsurgeon: random uniformly spaced variant sites that are not in the un-filtered VCF from longranger

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/0.setup/generateVarfile.py --fai /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa.fai --out out.varfile --vcf ../orig_VCFs/HCC1954N_WGS_210_phased_variants.vcf --variantspacing 90000

wc -l out.varfile
34769 out.varfile
```

Picard tools needs this version of java  
`module load java/1.8.0_45`

Run bamsurgeon on the whole file at once (slurm); code modified to only mutate reads with HP=1 tag  

```
sbatch --job-name=bsN --partition=lrgmem --nodes=1  --ntasks-per-node=48 --time=48:0:0 --export=ALL --workdir=/work-zfs/mschatz1/cdarby/HCC1954N/bamsurgeon --mem=400G --wrap="python /work-zfs/mschatz1/cdarby/bamsurgeon/bin/addsnv.py -v /work-zfs/mschatz1/cdarby/HCC1954N/bamsurgeon/out.bed -f /work-zfs/mschatz1/cdarby/HCC1954N/orig_bams/HCC1954N_hg19_calmd.bam -r  /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa -o /work-zfs/mschatz1/cdarby/HCC1954N/bamsurgeon/bamsurgeon.bam --picardjar /work-zfs/mschatz1/resources/picard.jar --minmutreads 1 --ignoresnps --force --tagreads --ignorepileup --aligner mem --seed 0 -p 48 --skipmerge"
```

Made 16,464 mutations  

Process and merge the reads with a variant of bamsurgeon's `replacereads.py` that preserves original file's AS and XS tags.

Sort and index the merged file

Calculate feature vectors at mosaic sites

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/1.simulate/simulate.py --bam ../orig_bams/bamsurgeon_HCC1954N_hg19_calmd.bam --varfile sites.vcf --nproc 4 > bamsurgeon.features.tsv
```

### Tumor

Generate varfile for bamsurgeon: random uniformly spaced variant sites that are not in the un-filtered VCF from longranger

```
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/0.setup/generateVarfile.py --fai /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa.fai --out out.varfile --vcf ../orig_VCFs/HCC1954T_WGS_210_phased_variants.vcf --variantspacing 90000

wc -l out.varfile
34774 out.varfile
```

Picard tools needs this version of java  
`module load java/1.8.0_45`

Run bamsurgeon on the whole file at once (slurm); code modified to only mutate reads with HP=1 tag  

```
sbatch --job-name=bsT --partition=lrgmem --nodes=1  --ntasks-per-node=48 --time=48:0:0 --export=ALL --workdir=/work-zfs/mschatz1/cdarby/HCC1954T/bamsurgeon --mem=400G --wrap="python /work-zfs/mschatz1/cdarby/bamsurgeon/bin/addsnv.py -v /work-zfs/mschatz1/cdarby/HCC1954T/bamsurgeon/out.bed -f /work-zfs/mschatz1/cdarby/HCC1954T/orig_bams/HCC1954T_hg19_calmd.bam -r  /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa -o /work-zfs/mschatz1/cdarby/HCC1954T/bamsurgeon/bamsurgeon.bam --picardjar /work-zfs/mschatz1/resources/picard.jar --minmutreads 1 --ignoresnps --force --tagreads --ignorepileup --aligner mem --seed 0 -p 48 --skipmerge"
```

Made 12,248 mutations

Process and merge the reads with a variant of bamsurgeon's `replacereads.py` that preserves original file's AS and XS tags.

Sort and index the merged file

Calculate feature vectors at mosaic sites


## 2. Performance

Generate varfile for edited-sites simulator: random uniformly spaced variant sites that are not in the un-filtered VCF from longranger

Simulate mosaicism at the sites

```
python ../pipeline/1.simulate/simulate.py --simulate --bam ../orig_bams/HCC1954T_hg19_calmd.bam --varfile out.varfile --nproc 2 > mosaic.features.tsv

wc -l mosaic.features.tsv
25214 mosaic.features.tsv
```

Calculate feature vectors at called het. sites
Calculate feature vectors at called hom. sites

Experiments 1-3: 10K mosaic; 10K germline test/holdout  
Experiment 4: 3K mosaic; 3K germline test/holdout

```
for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --mosaic bamsurgeon/bamsurgeon.features.tsv --testmosaic editedsites/mosaic.features.tsv --het bamsurgeon/het.features.tsv --hom bamsurgeon/hom.features.tsv --outpfx 2eval/tmp ; Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R 2eval/tmpprediction.txt | cut -f 2 -d" "  >> 2eval/{T,N}trainBStestE.txt ; done


for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --testmosaic bamsurgeon/bamsurgeon.features.tsv --mosaic editedsites/mosaic.features.tsv --het bamsurgeon/het.features.tsv --hom bamsurgeon/hom.features.tsv --outpfx 2eval/tmp ; Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R 2eval/tmpprediction.txt | cut -f 2 -d" "  >> 2eval/{T,N}trainEtestBS.txt ; done


for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --mosaic editedsites/mosaic.features.tsv --het bamsurgeon/het.features.tsv --hom bamsurgeon/hom.features.tsv --outpfx 2eval/tmp ; Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R 2eval/tmpprediction.txt | cut -f 2 -d" "  >> 2eval/{T,N}trainEtestE.txt ; done

for TRIAL in `seq 1 11`; do python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/traindatasize.py --mosaic bamsurgeon/bamsurgeon.features.tsv --het bamsurgeon/het.features.tsv --hom bamsurgeon/hom.features.tsv --outpfx 2eval/tmp ; Rscript /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/5.experiments/AUCcl.R 2eval/tmpprediction.txt | cut -f 2 -d" "  >> 2eval/{T,N}trainBStestBS.txt ; done

```

Scan

```
bedtools makewindows -w 1000000 -g /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa.fai > windows_ALL.bed
python ../../HCC1954T/pipeline/3.classify/classify.py --bam ../orig_bams/bamsurgeon_HCC1954N_hg19_calmd.bam --clf clf.pkl --nproc 8 --bed windows_ALL.bed > scan.tsv
```

CNV filter (cnvnator)  
Simple repeat filter  
Extract filtering-features

`python ../../HCC1954T/pipeline/4.filter/computeFilterFeatures.py --bam ../orig_bams/MERGE.bam --bed scan_norepeat_nocnv.tsv --nproc 8 > HCCNBS.features.tsv`

Filter based on features and linkage  

```
join -t $'\t' -j1 variants.txt HCCTE.linkagefiltered.tsv2 | gawk '{print "1\t" $2 "." $8 "\t" $9}' FS='\t' > HCCTE_labeled.txt
join -t $'\t' -j1 sitesT.txt HCCTE.linkagefiltered.tsv1 | gawk '{print "1\t" $2 "." $8 "\t" $9}' FS='\t' >> HCCTE_labeled.txt
cat HCCTE.linkagefiltered.tsv2 | gawk '{print "0\t" $2 "." $8 "\t" $9}' FS='\t' >> HCCTE_labeled.txt
sort -u -k2 -o HCCTE_labeled.txt HCCTE_labeled.txt
sort -nr -k3 -o HCCTE_labeled.txt HCCTE_labeled.txt

join -t $'\t' -j1 sitesN.txt HCCNE.linkagefiltered.tsv1 | gawk '{print "1\t" $2 "." $8 "\t" $9}' FS='\t' > HCCNE_labeled.txt
cat HCCNE.linkagefiltered.tsv1 | gawk '{print "0\t" $2 "." $8 "\t" $9}' FS='\t' >> HCCNE_labeled.txt
sort -u -k2 -o HCCNE_labeled.txt HCCNE_labeled.txt
sort -nr -k3 -o HCCNE_labeled.txt HCCNE_labeled.txt

awk '{print $1 "." $7 "\t" $0}' < scanTE_norepeat_linkage.features.tsv > scanTE_norepeat_linkage.features1.tsv
mv scanTE_norepeat_linkage.features1.tsv scanTE_norepeat_linkage.features.tsv
sort -k1 -o scanTE_norepeat_linkage.features.tsv scanTE_norepeat_linkage.features.tsv
join -t $'\t' -j1 ../variants0idx.txt scanTE_norepeat_linkage.features.tsv | awk '{print "2\t" $0}' > scanTE_norepeat_linkage_labeled.features.tsv
join -t $'\t' -j1 ../sitesT.txt scanTE_norepeat_linkage.features.tsv | awk '{print "1\t" $0}' >> scanTE_norepeat_linkage_labeled.features.tsv
cat scanTE_norepeat_linkage.features.tsv | awk '{print "0\t" $0}' >> scanTE_norepeat_linkage_labeled.features.tsv
sort -u -k2 -o scanTE_norepeat_linkage_labeled.features.tsv scanTE_norepeat_linkage_labeled.features.tsv
```
**4342 calls**
![](testdata_DO_NOT_COMMIT/HCCPair/completefeatures/multi_real_TEfiltered.png)
**5662 calls** 
![](testdata_DO_NOT_COMMIT/HCCPair/completefeatures/multi_bs_TEfiltered.png)


## 3. HapMuC

Run each chromosome individually; needed to restart with segfault and to remove a few windows on chr8 that ran for hours.

```
for NUM in `seq 1 22`; do mkdir ${NUM}; samtools mpileup -B -f /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa -r chr${NUM} ../HCC1954N/orig_bams/HCC1954N_hg19_calmd.bam > ${NUM}/normal.pileup; samtools mpileup -B -f /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa -r chr${NUM} ../HCC1954T/orig_bams/HCC1954T_hg19_calmd.bam > ${NUM}/tumor.pileup; /work-zfs/mschatz1/resources/hapmuc/utils/make_windows.sh ${NUM}/tumor.pileup ${NUM}/normal.pileup ${NUM}; rm ${NUM}/*pileup; done

for NUM in `seq 1 22`; do /work-zfs/mschatz1/resources/hapmuc/bin/hapmuc -a /work-zfs/mschatz1/cdarby/HCC1954T/orig_bams/HCC1954T_hg19_calmd.bam -b /work-zfs/mschatz1/cdarby/HCC1954N/orig_bams/HCC1954N_hg19_calmd.bam -f /work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa -w ${NUM}/windowsSNP -o ${NUM}/${NUM}; done
```

**159979 calls**  
![](testdata_DO_NOT_COMMIT/HCCPair/completefeatures/multi_hapmuc_fisherscore.png)
**36774 calls with a hapmuc score**  
![](testdata_DO_NOT_COMMIT/HCCPair/completefeatures/multi_hapmuc_hmscore.png)


## 4. MosaicHunter 

`awk '{print $1 "\t" $2-1 "\t" $2 "\t" $21}' < pairedMH.tsv > pairedMH.bed`
** 3745 calls **  
![](testdata_DO_NOT_COMMIT/HCCPair/MHunter/multi_MHunter_34.png)