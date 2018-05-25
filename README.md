# samovar
Somatic (mosaic) SNV caller for 10X Genomics data using random forest classification and feature-based filters

Requires: Python 2 or 3 with [pyfaidx](https://github.com/mdshw5/pyfaidx); [scikit-learn](http://scikit-learn.org/); [simplesam](http://simplesam.readthedocs.io/en/latest/); [fisher](https://pypi.org/project/fisher/); and installations of [samtools](http://www.htslib.org/), [bedtools](https://github.com/arq5x/bedtools2/releases). Some steps are compatible with [pypy](https://pypy.org/)  
[MARCC](https://www.marcc.jhu.edu/) The SLURM sbatch parameter --ntasks N should be set to accompany Python script parameter --nproc N  

# Main Pipeline

## 0.setup/generateVarfile.py
**Generates a tab-separated file with fields** `contig \t site \t VAF` **to specify where mosaic-like sites will be simulated**

```
python generateVarfile.py --out out.varfile --vcf sample.vcf --fai genome.fa.fai
```  

`--out` output varfile name [out.varfile]  
`--vcf` sites to avoid in simulation, e.g. the VCF from your sample (takes comma-separated multiple file names)  
`--fai` FASTA file index of the genome from `samtools faidx`  
`--variantspacing` spacing along the genome of simulated sites [80000]  
`--vafspacing	` increment of variant allele fraction [0.05]  

**Notes:**   

*  25,000-30,000 mosaic-like sites are suggested. Running this step with --variantspacing 80000 usually results in enough sites output from the simulation step. Decrease --variantspacing to increase the number of mosaic-like sites.  
*  This step takes only a few minutes, mostly spent reading the VCF file. It is not parallelized and pypy is not necessary.  
*  If you do not want to use miscellaneous contigs and/or sex chromosomes, make a new file for the `--fai` argument that contains only the contigs you want to use for mosaic-like sites.


## 1.simulate/simulate.py
**Simulates mosaic-like training examples at sites specified in varfile where there is at least one haplotype-discordant read**

```
python simulate.py --bam sample.bam --varfile out.varfile --simulate --nproc 4 > mosaic.features.tsv
```

`--simulate` mosaic site simulation mode  
`--bam` bam file to simulate sites from  
`--varfile` varfile from `0.setup/generateVarfile.py` of sites and VAF  
`--nproc` uses `multiprocessing.Pool` (value <= 1 does not import the module and runs in serial)  


**Notes:**   

* Records are printed to stdout, the output order may differ from the input order if multiple processes are used  
* With pypy, 4 processes, and 30,000 sites in the varfile, this step takes 30 minutes. (SLURM shared partition default memory 5G/task is sufficient)  
 
**Computes features at known homozygous/heterozygous sites from VCF that have at least one haplotype-discordant read**  

```
python simulate.py --bam sample.bam --varfile sample.vcf --het --max 50000 --nproc 4 > het.features.tsv
```

`--het` will extract features from sites with genotype `1|0\|0|1` that PASS filter and are SNP [i.e. len(ref) == 1 and len(alt) == 1]  
`--hom` will extract features from sites with genotype `0|0\|1|1` that PASS filter and are SNP [i.e. len(ref) == 1 and len(alt) == 1]  
`--bam` bam file to simulate sites from  
`--varfile` sites to extract features from, e.g. the VCF from your sample  
`--max` maximum number of sites to extract features from  
`--nproc` uses `multiprocessing.Pool` (value <= 1 does not import the module and runs in serial)  
Must catch output from stdout.

**Notes:**   

* Records are printed to stdout, the output order may differ from the input order if multiple processes are used 
* Based on the number of mosaic-like sites, you should set --max to half that number in each of --hom and --het mode.  
* Because more sites must be examined to find enough examples with haplotype-discordant reads, with pypy, 4 processes, and --max 15000, this step takes 30-90 minutes for modes --het and --hom. (SLURM shared partition default memory 5G/task is sufficient)  
 
## 2.train/train.py
**Uses feature vectors from mosaic-like, heterozygous, and homozygous sites to train a random forest classifier**

```
python train.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv
```

`--mosaic` file from `--simulate` in `1.simulate/simulate.py`  
`--het` file from `--het` in `1.simulate/simulate.py`  
`--hom` file from `--simulate` in `1.simulate/simulate.py`   
Equal number of het and hom examples are randomly selected to give equal number of mosaic and germline training examples.
`--out` classifier file in Python's pickle binary format [clf.pkl]
`--mindepth` read coverage required to use a site as a training example [16]  
`--nestimators` [100] `maxleafnodes` [50] scikit-learn random forest hyperparameters  

**Notes:** 
 
* Ensure that the python and scikit-learn versions used to train and dump the classifier are the same as used to load and use the classifier later on.
* This step should just take a few minutes, and is not parallelized or pypy-compatible.

## 3.classifyAndFilter/filter.py
**Scans genome, calculates features for sites that pass all filters. Outputs "vectors" for passing sites which will later be classified and ranked by the random forest**

```
python classifyAndFilter/filter.py --bam sample.bam --nproc 48 > vectors.txt 2> intervalsComplete.txt
```

`--bam` bam file of your sample  
`--nproc` uses `multiprocessing.Pool` (value <= 1 does not import the module and runs in serial)  


**Notes:**  

* Records are printed to stdout, the output order may not be in order along the genome if multiple processes are used 
* Regions are printed to stderr when they are completed.
* To decrease memory usage at possible loss of speed, decrease the "performance parameters" in the script. WINDOW\_SIZE is the number of bases each task scans at a time. READBATCH and SITEBATCH control how often the main data structure is pruned.    
* To adjust filter thresholds, change the "filter parameters" in the script. To consider sites with more (less) evidence of mosaicism, increase (decrease) MIN\_MAF and MIN\_HAPDISCORD\_READS. To consider sites with more (less) phased reads, increase (decrease) MIN\_FRAC\_PHASED. To consider sites with higher (lower) quality phasing, increase (decrease) MIN_FRAC_HAP1, the haplotype imbalance of all phased reads, and decrease (increase) MAX\_HAPDISCORD\_READS\_HAP, the haplotype imbalance of haplotype-discordant reads.  
* If you want to include any contigs besides the autosomes and chrX, add the names to the CONTIGS list. (if the sample is male, the reads on chrX will not have the HP tag from longranger and no calls will be made) 
* Processes have been observed that use >100G of memory on a single interval, and this has not yet been investigated and solved. [MARCC] Use SLURM `--partition=lrgmem --nodes=1 --ntasks-per-node=48 --mem=800G --time=12:0:0` which takes 3-5 hours to scan the whole genome, producing 50,000-100,000 candidate vectors (30-70MB output file).  

## 3.classifyAndFilter/classify.py
**Uses the random forest model to rank the feature vectors, outputting those with a higher-than-threshold probability in an output file with format:** `contig \t position \t position+1 \t Samovar score \t read depth \t Number of haplotype discordant read \t MAF \t minor allele base`

```
python classify.py --clf clf.pkl --vectors vectors.txt > predictions.tsv
```

`--clf` classifier trained from 2.train/train.py  
`--vectors` output file from filter step   

**Notes:**  

* Records are printed to stdout
* Ensure that the python and scikit-learn versions used to train and dump the classifier earlier the same as used to load and use the classifier here.
* This step should just take a few minutes, and is not parallelized or pypy-compatible.  
* An invocation of `eval()` in the python script may report an unexpected EOF in the middle of the file. Replacing this character with a newline solves the issue, which has not yet been investigated and solved.  
* If you do not want to use miscellaneous contigs and/or sex chromosomes, make a new file for the `--fai` argument that contains only the contigs you want to use for mosaic-like sites.  
* Decrease (set to 0) MIN\_CLF\_SCORE to output more (all) of the sites.  
* About 20,000-30,000 sites should remain after this step with a probability threshold of 0.9.

## Use bedtools to apply blacklist region filters
**Filter out simple repeats and CNV**

```
bedtools intersect -v -a predictions.tsv -b repeatsb37.bed | bedtools intersect -v -a stdin -b CNVNATOR.bed > regionfiltered.tsv
```

**Notes:**  

* Our suggested blacklist region is 1,2,3,4-nt repeats at least 4bp long with at least 3 copies of the "unit" with a 2bp buffer, and CNVNATOR calls with a 5bp buffer.  
* This step should only take a few minutes but because there are so many intervals of simple repeats, the memory consumption could be up to 30G unless both files are sorted and the `-sorted` flag is used in `bedtools intersect`.  
* About 5,000-10,000 sites should remain after this step.

## 3.classifyAndFilter/linkageFilter.py
**Calculates (and appends to output line) the minimum Fisher score for association between minor allele reads and mismatches, indels, or alignment endpoints**

```
python linkageFilter.py --bam sample.bam --bed regionfiltered.tsv --ref genome.fa --vcfavoid sample.vcf --nproc 4 > predictionsFinal.tsv
```

`--bam` BAM file of your sample   
`--bed` variant calls  
`--vcfavoid` VCF of sites NOT to consider for possible linkage (i.e. VCF from your sample)  
`--ref` reference genome FASTA  

**Notes:**  

* Change MIN\_FISHER\_PVAL to adjust the threshold below which sites will be rejected.  
* This step takes less than 10 minutes with pypy, 4 processes for 8000 sites.
* About 4,000-8,000 sites should remain after this step with a probability threshold of 0.05.

#Utilities


### 0.setup/cnvnator2bed.py
```
python cnvnator2bed.py [infile] [outfile] [.fai of genome]
```  
Converts cnvnator output into a .bed file with +/- 5 buffer around CNV regions


#Example

### example/EXAMPLE.txt

Reads and variants from [Genome in a Bottle NA24385 GRCh38](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics/) extracted by:

```
samtools view -b [input.bam] chr1:200000000-201000000 > example.bam
samtools index example.bam
vcftools --vcf [input.vcf] -c --recode --chr chr1 --from-bp 200000000 --to-bp 201000000 --remove-indels --remove-filtered-all > example.vcf
```

Can replace "python" with "pymp" and/or add "--nproc N" to commands where permitted  
Without pymp or parallelism, the filter step will take about 6 minutes; the rest will take a few seconds each. 
 
```
python ../1.simulate/simulate.py --bam example.bam --varfile out.varfile --simulate > outs/mosaic.features.tsv
python ../1.simulate/simulate.py --bam example.bam --varfile example.vcf --hom  > outs/hom.features.tsv
python ../1.simulate/simulate.py --bam example.bam --varfile example.vcf --het  > outs/het.features.tsv
python ../2.train/train.py --mosaic outs/mosaic.features.tsv --het outs/het.features.tsv --hom outs/hom.features.tsv --nestimators 1 --maxleafnodes 5 --out outs/clf.pkl
python ../3.classifyAndFilter/filter.py --bam example.bam
python ../3.classifyAndFilter/classify.py --clf outs/clf.pkl --vectors outs/vectors.txt > outs/predictions.txt
python ../3.classifyAndFilter/linkageFilter.py --bam example.bam --bed outs/predictions.txt --vcfavoid example.vcf --ref [GRCh38 reference genome FASTA] > outs/finalPredictions.txt
```

  
#Experiments


### experiments/bamsurgeon2vcf.py
```
python bamsurgeon2vcf.py < [infile]
```
Converts bamsurgeon logfile output into a minimal .vcf file with genotype 0|0 to be used with `simulate.py --hom` 

### experiments/subsampleLinkedBam.py
```
python subsampleLinkedBam.py [infile] [outfile]
```
Subsamples half of the barcodes in the file (BX tag)

### experiments/pairedTagBam.py
```
python pairedTagBam.py --bam sample.bam --out paired.bam --vcf sample.vcf
```
Based on phased VCF, annotates reads in bam file with PH tag of 1 or 2 if it or mate can be phased.  
`--bam` bam file with paired-end reads  
`--vcf` phased VCF (GT 0|1\|1|0)  
`--out` output bam file with PH tags  


### experiments/simulatePaired.py
```
python simulatePaired.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 --subsample 0.5 > mosaic.features.tsv
```

Substitutes PH for HP tag and doesn't produce the ASXS features.  


### experiments/simulateSubsample.py
```
python simulateSubsample.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 --subsample 0.5 > mosaic.features.tsv
```

`--subsample` fraction of depth to retain  


### experiments/trainUnphased.py
```
python train.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv
```

Trains random forest only on features that have no phasing information by subsetting the feature vector   

### experiments/classifyAndFilterUnphased.py
```
python classifyAndFilterUnphased.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8 > scan.tsv
```

Classifies and filters based on model trained with no phasing information   

### experiments/classifyAndFilterPaired.py
```
python classifyAndFilterPaired.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8 > scan.tsv
```

Classifies and filters based on model trained with phasing computed using paired-end reads   


### experiments/comparisonTables.R

Figures
