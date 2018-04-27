#SARCOMA STUDY

## 0. General

Samples 0, 10, 1, 2, 3, 6, 9, C1

Sort bam file by name first, otherwise not compatible with bamtofastq

```
~/scratch/bamtofastq --gemcode /scratch/groups/mschatz1/mschatz/sidow/sarcoma/10XSARCOMA0.bam fastqs
~/scratch/longranger-1.3.1/longranger run --sex=m --id=SARCOMA0 --fastqs=fastqs --fastqprefix=bamtofastq --nopreflight --jobmode=slurm --reference=/work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/
samtools view phased_possorted_bam.bam | grep -c "HP"
```

`10XSARCOMA0.bam`

```
891330438 + 0 in total (QC-passed reads + QC-failed reads)
14215374 + 0 secondary
0 + 0 supplementary
6377971 + 0 duplicates
865356029 + 0 mapped (97.09% : N/A)
877115064 + 0 paired in sequencing
438557532 + 0 read1
438557532 + 0 read2
787332988 + 0 properly paired (89.76% : N/A)
846449050 + 0 with itself and mate mapped
4691605 + 0 singletons (0.53% : N/A)
38012988 + 0 with mate mapped to a different chr
30309591 + 0 with mate mapped to a different chr (mapQ>=5)
```

371270109/865356029 42.9% reads phased

`10XSARCOMA10.bam`

```
867495683 + 0 in total (QC-passed reads + QC-failed reads)
13424189 + 0 secondary
0 + 0 supplementary
8080161 + 0 duplicates
842939239 + 0 mapped (97.17% : N/A)
854071494 + 0 paired in sequencing
427035747 + 0 read1
427035747 + 0 read2
767801556 + 0 properly paired (89.90% : N/A)
824154180 + 0 with itself and mate mapped
5360870 + 0 singletons (0.63% : N/A)
32816500 + 0 with mate mapped to a different chr
26547530 + 0 with mate mapped to a different chr (mapQ>=5)
```

501127641/842939239 = 59.5% reads phased

`10XSARCOMA1.bam`

```
877219894 + 0 in total (QC-passed reads + QC-failed reads)
12682666 + 0 secondary
0 + 0 supplementary
8444013 + 0 duplicates
852953731 + 0 mapped (97.23% : N/A)
864537228 + 0 paired in sequencing
432268614 + 0 read1
432268614 + 0 read2
783360690 + 0 properly paired (90.61% : N/A)
836100798 + 0 with itself and mate mapped
4170267 + 0 singletons (0.48% : N/A)
29104260 + 0 with mate mapped to a different chr
22802571 + 0 with mate mapped to a different chr (mapQ>=5)
```

480053123/852953731 = 56.3% reads phased

`10XSARCOMA2.bam`

```
879931172 + 0 in total (QC-passed reads + QC-failed reads)
13529062 + 0 secondary
0 + 0 supplementary
7137760 + 0 duplicates
857219366 + 0 mapped (97.42% : N/A)
866402110 + 0 paired in sequencing
433201055 + 0 read1
433201055 + 0 read2
781684690 + 0 properly paired (90.22% : N/A)
839350980 + 0 with itself and mate mapped
4339324 + 0 singletons (0.50% : N/A)
33636470 + 0 with mate mapped to a different chr
26361711 + 0 with mate mapped to a different chr (mapQ>=5)
```

448121593/857219366 = 52.3% reads phased

`10XSARCOMA3.bam`

```
866768894 + 0 in total (QC-passed reads + QC-failed reads)
12635812 + 0 secondary
0 + 0 supplementary
7643573 + 0 duplicates
843169690 + 0 mapped (97.28% : N/A)
854133082 + 0 paired in sequencing
427066541 + 0 read1
427066541 + 0 read2
771637304 + 0 properly paired (90.34% : N/A)
825618120 + 0 with itself and mate mapped
4915758 + 0 singletons (0.58% : N/A)
29739826 + 0 with mate mapped to a different chr
22309349 + 0 with mate mapped to a different chr (mapQ>=5)
```

414229584/843169690 = 49.1% reads phased

`10XSARCOMA6.bam`

```
873456428 + 0 in total (QC-passed reads + QC-failed reads)
12989166 + 0 secondary
0 + 0 supplementary
7873398 + 0 duplicates
850013840 + 0 mapped (97.32% : N/A)
860467262 + 0 paired in sequencing
430233631 + 0 read1
430233631 + 0 read2
780514540 + 0 properly paired (90.71% : N/A)
832742342 + 0 with itself and mate mapped
4282332 + 0 singletons (0.50% : N/A)
29466838 + 0 with mate mapped to a different chr
22671785 + 0 with mate mapped to a different chr (mapQ>=5)
```

438575936/850013840 = 51.6% reads phased

`10XSARCOMA9.bam`

```
887581382 + 0 in total (QC-passed reads + QC-failed reads)
13497006 + 0 secondary
0 + 0 supplementary
6707799 + 0 duplicates
862999268 + 0 mapped (97.23% : N/A)
874084376 + 0 paired in sequencing
437042188 + 0 read1
437042188 + 0 read2
788651686 + 0 properly paired (90.23% : N/A)
844806626 + 0 with itself and mate mapped
4695636 + 0 singletons (0.54% : N/A)
33349752 + 0 with mate mapped to a different chr
25558301 + 0 with mate mapped to a different chr (mapQ>=5)
```

393211473/862999268 = 45.6% reads phased

`10XSARCOMAC1.bam`

```
853518499 + 0 in total (QC-passed reads + QC-failed reads)
12142641 + 0 secondary
0 + 0 supplementary
7138794 + 0 duplicates
827130048 + 0 mapped (96.91% : N/A)
841375858 + 0 paired in sequencing
420687929 + 0 read1
420687929 + 0 read2
762771160 + 0 properly paired (90.66% : N/A)
809533986 + 0 with itself and mate mapped
5453421 + 0 singletons (0.65% : N/A)
26834738 + 0 with mate mapped to a different chr
19396304 + 0 with mate mapped to a different chr (mapQ>=5)
```

373956625/827130048 = 45.2% reads phased


**Subsample for Ginkgo CNV clustering analysis**  

`samtools view -q 20 -s 0.04 -b 6/SARCOMA6/outs/phased_possorted.bam | bedtools bamtobed -i stdin | cut -f 1-3`

[Ginkgo results](http://qb.cshl.edu/ginkgo?q=results/4bsmN6FxrQtk700znu0s)

**CNVnator**

```
cnvnator -root ${NUM}/out.root -tree ${NUM}/SARCOMA${NUM}/outs/phased_possorted_bam.bam
cnvnator -root ${NUM}/out.root -his 100 -d /work-zfs/mschatz1/resources/grch38/
cnvnator -root ${NUM}/out.root -stat 100
cnvnator -root ${NUM}/out.root -partition 100
cnvnator -root ${NUM}/out.root -call 100 > ${NUM}/cnvnator.tsv'
python /work-zfs/mschatz1/cdarby/HCC1954T/pipeline/0.setup/cnvnator2bed.py ${NUM}/cnvnator.tsv ${NUM}/cnvnator.bed /work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa.fai ; done
```
Base pairs affected by CNVnator +5bp buffer
```
0
2005470598
1
1881236238
10
1864930808
2
1881274616
3
1727553649
6
1692908916
9
1801646122
C1
1370152982
```

**Insert size of mate pair library**  
Control: mean 6766; median 6606

##MosaicHunter

Single: 7 10X libraries / 3 Mate Pair / 7 Illumina (17 possibilities)
Paired: There is separate control of Mate Pair, 10X, and Illumina (17 possibilities, assuming each datatype is paired with its own control)


##HapMuc

Paired: 17 possibilities, assuming each datatype is paired with its own control



##OTHER

`MP_sarcoma_0.combined.sorted.bam`

```
1059331136 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2056262 + 0 supplementary
64522716 + 0 duplicates
914926425 + 0 mapped (86.37% : N/A)
1057274874 + 0 paired in sequencing
528637437 + 0 read1
528637437 + 0 read2
729382518 + 0 properly paired (68.99% : N/A)
773716234 + 0 with itself and mate mapped
139153929 + 0 singletons (13.16% : N/A)
35262860 + 0 with mate mapped to a different chr
28934950 + 0 with mate mapped to a different chr (mapQ>=5)
```

`MP_sarcoma_10.combined.sorted.bam`

```
1045963696 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1868608 + 0 supplementary
73846192 + 0 duplicates
903131298 + 0 mapped (86.34% : N/A)
1044095088 + 0 paired in sequencing
522047544 + 0 read1
522047544 + 0 read2
723796454 + 0 properly paired (69.32% : N/A)
763022398 + 0 with itself and mate mapped
138240292 + 0 singletons (13.24% : N/A)
31219442 + 0 with mate mapped to a different chr
25489855 + 0 with mate mapped to a different chr (mapQ>=5)
```

`MP_sarcoma_9.combined.sorted.bam`

```
1060689350 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2131848 + 0 supplementary
89345001 + 0 duplicates
919439551 + 0 mapped (86.68% : N/A)
1058557502 + 0 paired in sequencing
529278751 + 0 read1
529278751 + 0 read2
742855748 + 0 properly paired (70.18% : N/A)
780479394 + 0 with itself and mate mapped
136828309 + 0 singletons (12.93% : N/A)
28438598 + 0 with mate mapped to a different chr
22213355 + 0 with mate mapped to a different chr (mapQ>=5)
```

`MP_sarcoma_C.combined.sorted.bam`

```
1008691374 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2571318 + 0 supplementary
83896479 + 0 duplicates
865799144 + 0 mapped (85.83% : N/A)
1006120056 + 0 paired in sequencing
503060028 + 0 read1
503060028 + 0 read2
684896654 + 0 properly paired (68.07% : N/A)
724707928 + 0 with itself and mate mapped
138519898 + 0 singletons (13.77% : N/A)
29849134 + 0 with mate mapped to a different chr
23045241 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_0.bam`

```
911541982 + 0 in total (QC-passed reads + QC-failed reads)
5307320 + 0 secondary
0 + 0 supplementary
110488268 + 0 duplicates
907818900 + 0 mapped (99.59% : N/A)
906234662 + 0 paired in sequencing
453117331 + 0 read1
453117331 + 0 read2
885488874 + 0 properly paired (97.71% : N/A)
900236272 + 0 with itself and mate mapped
2275308 + 0 singletons (0.25% : N/A)
9415962 + 0 with mate mapped to a different chr
4124070 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_10.bam`

```
870840749 + 0 in total (QC-passed reads + QC-failed reads)
5127139 + 0 secondary
0 + 0 supplementary
85646514 + 0 duplicates
867453044 + 0 mapped (99.61% : N/A)
865713610 + 0 paired in sequencing
432856805 + 0 read1
432856805 + 0 read2
845856582 + 0 properly paired (97.71% : N/A)
860014160 + 0 with itself and mate mapped
2311745 + 0 singletons (0.27% : N/A)
8726120 + 0 with mate mapped to a different chr
3678561 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_1.bam`

```
894939121 + 0 in total (QC-passed reads + QC-failed reads)
5057145 + 0 secondary
0 + 0 supplementary
99391298 + 0 duplicates
890042702 + 0 mapped (99.45% : N/A)
889881976 + 0 paired in sequencing
444940988 + 0 read1
444940988 + 0 read2
866160716 + 0 properly paired (97.33% : N/A)
881420464 + 0 with itself and mate mapped
3565093 + 0 singletons (0.40% : N/A)
9737848 + 0 with mate mapped to a different chr
4464647 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_2.bam`

```
910886290 + 0 in total (QC-passed reads + QC-failed reads)
5266552 + 0 secondary
0 + 0 supplementary
103655998 + 0 duplicates
907287041 + 0 mapped (99.60% : N/A)
905619738 + 0 paired in sequencing
452809869 + 0 read1
452809869 + 0 read2
884846336 + 0 properly paired (97.71% : N/A)
899731016 + 0 with itself and mate mapped
2289473 + 0 singletons (0.25% : N/A)
9230680 + 0 with mate mapped to a different chr
3964410 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_3.bam`

```
905238060 + 0 in total (QC-passed reads + QC-failed reads)
5481320 + 0 secondary
0 + 0 supplementary
102431303 + 0 duplicates
902386506 + 0 mapped (99.68% : N/A)
899756740 + 0 paired in sequencing
449878370 + 0 read1
449878370 + 0 read2
880244278 + 0 properly paired (97.83% : N/A)
895228914 + 0 with itself and mate mapped
1676272 + 0 singletons (0.19% : N/A)
9219610 + 0 with mate mapped to a different chr
3877255 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_6.bam`

```
870500643 + 0 in total (QC-passed reads + QC-failed reads)
5417973 + 0 secondary
0 + 0 supplementary
86249164 + 0 duplicates
867859429 + 0 mapped (99.70% : N/A)
865082670 + 0 paired in sequencing
432541335 + 0 read1
432541335 + 0 read2
846072432 + 0 properly paired (97.80% : N/A)
860971136 + 0 with itself and mate mapped
1470320 + 0 singletons (0.17% : N/A)
9143468 + 0 with mate mapped to a different chr
3774512 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_9.bam`

```
879089373 + 0 in total (QC-passed reads + QC-failed reads)
5159837 + 0 secondary
0 + 0 supplementary
82253328 + 0 duplicates
876372688 + 0 mapped (99.69% : N/A)
873929536 + 0 paired in sequencing
436964768 + 0 read1
436964768 + 0 read2
855316598 + 0 properly paired (97.87% : N/A)
869653152 + 0 with itself and mate mapped
1559699 + 0 singletons (0.18% : N/A)
8789820 + 0 with mate mapped to a different chr
3689228 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_C2-2.bam`

```
913818342 + 0 in total (QC-passed reads + QC-failed reads)
6103242 + 0 secondary
0 + 0 supplementary
103752342 + 0 duplicates
911199153 + 0 mapped (99.71% : N/A)
907715100 + 0 paired in sequencing
453857550 + 0 read1
453857550 + 0 read2
886412960 + 0 properly paired (97.65% : N/A)
903566420 + 0 with itself and mate mapped
1529491 + 0 singletons (0.17% : N/A)
10558768 + 0 with mate mapped to a different chr
4360841 + 0 with mate mapped to a different chr (mapQ>=5)
```

`Sarcoma_C2-2.chr22.bam`

```
9938795 + 0 in total (QC-passed reads + QC-failed reads)
36079 + 0 secondary
0 + 0 supplementary
914413 + 0 duplicates
9918220 + 0 mapped (99.79% : N/A)
9902716 + 0 paired in sequencing
4952449 + 0 read1
4950267 + 0 read2
9784730 + 0 properly paired (98.81% : N/A)
9861566 + 0 with itself and mate mapped
20575 + 0 singletons (0.21% : N/A)
67238 + 0 with mate mapped to a different chr
33872 + 0 with mate mapped to a different chr (mapQ>=5)
```
