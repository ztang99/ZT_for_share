Updated by Zitian [[10/12/2023]].

---

# [[10-10-2023]]

### ***Record for paths to BAMs & GVCFs***

***Parent Directory: all files listed below is in this parent directory:*** 
`/storage1/fs1/jin810-2/Active/Projects/SMaHT/data`

***File path for each sample***
- *BAM Paths*:
	- HG002:  `./bams_300x/HG002.bam`
	- HG005:  `./bams_300x/HG005.bam`
- *GVCF Paths*:
	- HG002:  `./gvcfs_300x/hg002/hg002_pbgen_allvar.g.vcf`
		- Num lines (excluding headers): 26, 144, 798
	- HG005:  `./gvcfs_300x/hg005/hg005_pbgen_allvar.g.vcf`
		- Num lines (excluding headers): 52, 433, 826

***Joint Calling***
- Joint-called VCF:
	- **Path**: `./gvcfs_300x/combined_25.vcf.gz`
		- Num lines (excluding headers): 7, 302, 096
	- (gunzipped `./gvcfs_300x/combined_25.vcf` is stored in the same directory)
- True mutation set:
	- VCF file: `./gvcfs_300x/filtered_combined_25.vcf`
		- Num lines (excluding headers): 331, 782
	- BED file: `./gvcfs_300x/true_mutation_set.bed`
- False positive set (tms +- 10bp):
	- BED file: `./gvcfs_300x/false_positive_set.bed`

# [[10-04-2023]]

### ***Record for size, running time, etc. for HG002 & 005***

*STEP 1 - bwa mem alignment* (FASTQ -> SAM)

| SampleID | SAM | FASTQ | MEM | Num CPU | TIME |
|---|---|---|---|---|---|
|HG002|2.5T|~700G|~13G|40|~30 hrs|
|HG005|2.5T|~750G|~18G|32|~38 hrs|

*STEP2 - Sort, addRG (SAM -> BAM)*

|SampleID | SAM | BAM | MEM_req | CPU_req | TIME |
|---|---|---|---|---|---|
|HG002|2.5T|752G|250GB|20|157h (6.5days)|
|HG005|2.5T|837G|250GB|40|108h (4.5days)|

*STEP3 - Variant calling (BAM -> GVCF)*

- Comparison between using parabricks vs. using GATK directly (using chr1)

|Software | CPU_req | GPU_req | Mem_req | TIME | Size | Homozygous? |
|---|---|---|---|---|---|---|
|Parabricks|20|2|100GB|11.5 min|87M|No
|GATK|20|2|100GB|741.4 min|249M|Yes

**Conclusion: Parabricks is way faster than GATK, meaning that GATK is not actually using GPUs when doing parallel running.**

|SampleID | BAM | GVCF | CPU_req | GPU_req | MEM_req | TIME |
|---|---|---|---|---|---|---|
|HG002|752G|3.0G|20|2|250GB|248 min (4.13h)|
|HG005|837G|5.3G|20|2|100GB|66 1min (11h)|