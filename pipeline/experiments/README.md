### experiments/bamsurgeon2vcf.py
```
python bamsurgeon2vcf.py < [infile]
```
Converts bamsurgeon logfile output into a minimal .vcf file with genotype 0|0 to be used with `simulate.py --hom` to calculate features at these sites  

### experiments/subsampleLinkedBam.py
```
python subsampleLinkedBam.py [infile] [outfile]
```
Subsamples half of the barcodes in the file (BX tag); uses pysam

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
python simulatePaired.py --bam sample.bam --varfile out.varfile --simulate --nproc 8 > mosaic.features.tsv
```

Substitutes PH for HP tag and doesn't produce the ASXS features. Also has hom; het modes.  


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

### experiments/filterUnphased.py
```
python filterUnphased.py --bam sample.bam --nproc 8 > vectors.txt 2> intervalsComplete.txt
```

Filters using no phasing information   

### experiments/filterPaired.py
```
python filterPaired.py --bam sample.bam --nproc 8 > vectors.txt 2> intervalsComplete.txt
```

Filters using phasing computed using paired-end reads   

### experiments/classifyUnphased.py
```
python classifyUnphased.py --clf clf.pkl --vectors vectors.txt > predictions.tsv
```

Output line: chrom, pos, pos+1, score, depth, MAF

### experiments/linkageFilterPaired.py  
```
python linkageFilterPaired.py --bam sample.bam --bed regionfiltered.tsv --ref genome.fa --vcfavoid sample.vcf --nproc 4 > predictionsFinal.tsv
```
