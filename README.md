# samovar
Somatic (mosaic) SNV caller for 10X Genomics data using random forest classification and feature-based filters

Requires: Python 2 or 3 with [pymp](https://github.com/classner/pymp) (`pip install pymp-pypi`); [scikit-learn](http://scikit-learn.org/); [pysam](http://pysam.readthedocs.io/en/latest/index.html)  
SLURM sbatch parameter --ntasks N should be set to accompany Python script parameter --nproc N
ß
#Main Pipeline

### 0.setup/generateVarfile.py
```
python generateVarfile.py --out out.varfile --vcf sample.vcf --fai genome.fa.fai --variantspacing 100000 --vafspacing 0.05
```  

`--out` output varfile name: columns are contig; position; VAF [out.varfile]  
`--vcf` sites to avoid, e.g. the VCF from your sample (takes comma-separated multiple file names)  
`--fai` fasta file index of the genome  
`--variantspacing` spacing along the genome of simulated sites [100000]  
`--vafspacing	` increment of variant allele fraction [0.05]  


### 1.simulate/simulate.py
```
python simulate.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 > mosaic.features.tsv
```

`--simulate` mosaic site simulation mode  
`--bam` bam file to simulate sites from  
`--varfile` varfile from `0.setup/generateVarfile.py` of sites and VAF  
`--nproc` uses `pymp` for multiprocessing (value <= 1 sets if_=False)  
Must catch output from stdout.

```
python simulate.py --bam sample.bam --varfile sample.vcf --het --max 50000 --nproc 8 > het.features.tsv
```

`--het` will extract features from sites with genotype `1|0\|0|1` that PASS filter and are SNP [i.e. len(ref) == 1 and len(alt) == 1]  
`--hom` will extract features from sites with genotype `0|0\|1|1` that PASS filter and are SNP [i.e. len(ref) == 1 and len(alt) == 1]  
`--bam` bam file to simulate sites from  
`--varfile` sites to extract features from, e.g. the VCF from your sample  
`--max` maximum number of sites to extract features from  
`--nproc` uses `pymp` for multiprocessing (value <= 1 sets if_=False)  
Must catch output from stdout.
 
 
### 2.train/train.py
```
python train.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv
```

`--mosaic` file from `--simulate` in `1.simulate/simulate.py`  
`--het` file from `--het` in `1.simulate/simulate.py`  
`--hom` file from `--simulate` in `1.simulate/simulate.py`   
Equal number of het and hom examples are randomly selected to give equal number of mosaic and germline training examples.
`--out` classifier file [clf.pkl]
`--mindepth` read coverage required to use a site as a training example [16]  
`--nestimators` [100] `maxleafnodes` [50] scikit-learn random forest hyperparameters  

**Ensure that the python and scikit-learn versions used to train and dump the classifier are the same as used to load and use the classifier**  


### 3.classify/classifyAndFilter.py
```
python classifyAndFilter.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8 > scan.tsv
```

`--bam` bam file of your sample  
`--bed` intervals to divide among the parallel processes, e.g. from `bedtools makewindows -g genome.fa.fai -w 1000000` or `bedtools complement -g genome.fa.fai -i bad_regions.bed`   
`--clf` classifier file from `2.train/train.py`  
`--nproc` uses `pymp` for multiprocessing (value <= 1 sets if_=False)  
Must catch output from stdout. Output columns are contig; position (0-indexed); position+1; classifier score  

* **To decrease memory usage at possible loss of speed, decrease the READBATCH and CLFBATCH parameters in the script**  
* **If you have a very large input bed file, to decrease memory usage at possible loss of parallelism efficiency, decrease the INTERVALBATCH parameter in the script**
* **To print more (or all) of the sites and their scores, decrease (set to 0.0) the MIN\_CLF_SCORE parameter in the script**
* **To adjust filter thresholds, change parameters in the script** 


### 4.filter/linkageFilter.py
```
python linkageFilter.py --bam sample.bam --bed scan.features.tsv --ref genome.fa --vcfavoid sample.vcf
```

`--bam` bam file of your sample   
`--bed` variant calls  
`--vcfavoid` VCF of variants NOT to use as linkage (i.e. VCF from your sample)  
`--ref` reference genome .fasta  
Output columns are contig; position with minimum p-value for linkage; minimum p-value; minimum p-value * number of sites tested; 2x2 linkage table [reprints input line]

#Utilities


### 0.setup/cnvnator2bed.py
```
python cnvnator2bed.py [infile] [outfile] [.fai of genome]
```  
Converts cnvnator output into a .bed file with +/- 5 buffer around CNV regions


#Example

### example/EXAMPLE.txt

reads from Genome in a Bottle NA24385 GRCh38 extracted by:
```
samtools view -b [input.bam] chr1:200010000-200036000 > example.bam
samtools index example.bam
```

```
python ../1.simulate/simulate.py --bam example.bam --varfile out.varfile --simulate > outs/mosaic.features.tsv
python ../1.simulate/simulate.py --bam example.bam --varfile example.vcf --hom  > outs/hom.features.tsv
python ../1.simulate/simulate.py --bam example.bam --varfile example.vcf --het  > outs/het.features.tsv
python ../2.train/train.py --mosaic outs/mosaic.features.tsv --het outs/het.features.tsv --hom outs/hom.features.tsv --nestimators 1 --maxleafnodes 5 --out outs/clf.pkl
python ../3.classify/classifyAndFilter.py --bam example.bam --bed window.bed --clf outs/clf.pkl
python ../4.filter/linkageFilter.py --bam example.bam --bed example_scan.tsv --ref /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa --vcfavoid example.vcf
```

#Experiments


### 5.experiments/bamsurgeon2vcf.py
```
python bamsurgeon2vcf.py < [infile]
```
Converts bamsurgeon logfile output into a minimal .vcf file with genotype 0|0 to be used with `simulate.py --hom` 

### 5.experiments/subsampleLinkedBam.py
```
python subsampleLinkedBam.py [infile] [outfile]
```
Subsamples half of the barcodes in the file (BX tag)

### 5.experiments/pairedTagBam.py
```
python pairedTagBam.py --bam sample.bam --out paired.bam --vcf sample.vcf
```
Based on phased VCF, annotates reads in bam file with PH tag of 1 or 2 if it or mate can be phased.  
`--bam` bam file with paired-end reads  
`--vcf` phased VCF (GT 0|1\|1|0)  
`--out` output bam file with PH tags  


### 5.experiments/simulatePaired.py
```
python simulatePaired.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 --subsample 0.5 > mosaic.features.tsv
```

Substitutes PH for HP tag and doesn't produce the ASXS features.  


### 5.experiments/simulateSubsample.py
```
python simulateSubsample.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 --subsample 0.5 > mosaic.features.tsv
```

`--subsample` fraction of depth to retain  


### 5.experiments/trainUnphased.py
```
python train.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv
```

Trains random forest only on features that have no phasing information by subsetting the feature vector   

### 5.experiments/classifyAndFilterUnphased.py
```
python classifyAndFilterUnphased.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8 > scan.tsv
```

Classifies and filters based on model trained with no phasing information   

### 5.experiments/classifyAndFilterPaired.py
```
python classifyAndFilterPaired.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8 > scan.tsv
```

Classifies and filters based on model trained with phasing computed using paired-end reads   


### 5.experiments/comparisonTables.R

Figures
