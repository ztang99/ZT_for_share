>Created by Zitian.
>Last modified [11-07-2023](11-07-2023).

---
>[!Start]
>As part of your rotation, you will be analyzing germline variants out of whole genome sequencing (WGS) data of one(1) family/trio within a congenital anomalies in kidney and urinary tract (CAKUT) disease cohort.

# Useful Resources

Before we start, make sure to take a look at the [on-boarding slides](https://docs.google.com/presentation/d/1aT26I-gsQVd8RJbeNTUzcj9Tfp-lgLXpVIToHwHc8AQ/edit?usp=sharing).

>[!Note]
>Multiple steps in our pipeline incorporates the usage of Genomic Analysis Toolkit (GATK) variant calling tools.

As a starting point, take a close look at [GATK introductory videos](https://www.youtube.com/playlist?list=PLjiXAZO27elDHGlQwfd06r7coiFtpPkvI) before actually dig into utilizing the pipeline. The videos contain relevant material in understanding how whole exome sequencing (WES) and whole genome sequencing (WGS) data is processed and analyzed.

> For general introduction: # 1, 2, 3, 4, 7, 8 
> For germline variant calling: # 5, 9, 10, 11, 12, 13 
> For pipelining: # 16, 17

# 1. Data paths

See below the paths to the fastq files for the one(1) trios mentioned above.

>[!Note] 
>How does fastq files look like? cd into each directory below and visualize the fastq files.

```bash
#data path
/storage1/fs1/jin810/Active/References/Jeannie_Basta_WES_data_06-08-23/MGI_exome_sequencing_072021/RLTO_12
```

# 2. Before Variant Calling

>[!Note]
>What is a fastq file? What is a BAM file? What is each type of file used for in bioinformatics analyses? How do you visualize a BAM/CRAM file?


## Script to use

>[!Note]
>For each and every script listed below, take a closer look at each of them, make necessary modifications to the output path(s), and try to understand the commands used before actually running the code.


- You can find **ALL** the scripts under this directory: `/storage1/fs1/jin810/Active/testing/ztang/misc/test_RoTATion`.
- You can either **COPY** the scripts above to your own directory or make your own bash scripts.
- If you’re copying over the scripts, make sure you have the correct **permission** (reading/writing) to run them.
- If you’re making your own bash scripts, **make them executable before running the script**: `chmod +x /path/to/script.sh` .
- Run each script using `./path/to/script.sh`
- If you are unsure about any input/output files, you can refer to `/storage1/fs1/jin810/Active/testing/ztang/misc/test_RoTATion/results` this folder (where I stored all output files when I tested this pipeline).

### Fastq → Bam → VCF

>[! Do the following before running ANY of the script]
>- `cd` to your folder where you want to strore/run the script
>- Export LSF docker volumes:
>	```bash
>	export LSF_DOCKER_VOLUMES="/storage1/fs1/bga/Active:/storage1/fs1/bga/Active /scratch1/fs1/jin810:/scratch1/fs1/jin810 /storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active $HOME:$HOME"
>	```
>- Request Docker
>	- The command for requesting a docker in **interactive mode** can be found around the top few lines in each script.
>	- If a job takes more than a couple hours to run, you can submit it instead in **uninteractive mode**
>	```bash
>	#interactive 
>	bsub -Is -G compute-jin810 -q general-interactive -n 4 -R 'rusage[mem=32GB]' -a 'docker(nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1)' /bin/bash 
>	#uninteractive 
>	bsub -G compute-jin810 -q general -n 4 -R 'rusage[mem=32GB]' -a 'docker(nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1)' bash parabricksFsq2Bam.sh
>	```

##### `Fsq2Bam.sh`

```bash
#!/bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa"
KS1="/storage1/fs1/bga/Active/gmsroot/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-19443-e48c595a620a432c93e8dd29e4af64f2/snvs.hq.vcf.gz"
KS2="/storage1/fs1/bga/Active/gmsroot/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-20267-00cb8ff552914c17ad66d86031e10d30/indels.hq.vcf.gz"
KS3="/storage1/fs1/bga/Active/gmsroot/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-20211-26b393cc7ab04120ac68cc2cbd4a15df/indels.hq.vcf.gz"

data_path="/storage1/fs1/jin810/Active/References/Jeannie_Basta_WES_data_06-08-23/MGI_exome_sequencing_072021/RLTO_12"
samplenames=("RLTO_12_Child_C" "RLTO_12_Dad_B" "RLTO_12_Mom_A")

### Change the out_path accordingly ###
out_path="/path/to/results"
[ ! -d $out_path ] && mkdir -p $out_path
echo "Output results to: ${out_path}"
#######################################

for sample_name in "${samplenames[@]}"; do
	fastq_r1=$(find "$data_path/$sample_name" -name "*_R1.fastq.gz" -type f)
	fastq_r2=$(find "$data_path/$sample_name" -name "*_R2.fastq.gz" -type f)
		out_bam="${out_path}/${sample_name}.cram"
		out_recal="${out_path}/${sample_name}_report.txt"
		out_duplicate_metrics="${out_path}/${sample_name}_dup_metrics.txt"
		temp_dir="${out_path}/${sample_name}_temp"

		pbrun fq2bam --x4 --ref "${REF}" \\
		--in-fq "${fastq_r1}" "${fastq_r2}" \\
		--knownSites "${KS1}" \\
		--knownSites "${KS2}" \\
		--knownSites "${KS3}" \\
		--out-bam "${out_bam}" \\
		--out-recal-file "${out_recal}" \\
		--out-duplicate-metrics "${out_duplicate_metrics}" \\
		--bwa-options='-Y' \\
		--read-group-pl ILLUMINA \\
		--read-group-sm "${sample_name}" \\
		--tmp-dir "${temp_dir}" \\
		--num-gpus 1 \\
		--memory-limit 80
done
```

##### `GATK_BQSR_HC.sh`

```bash
#!/bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa"
BED="/storage1/fs1/jin810/Active/yung-chun/database/capture_bedfile/Exome-IDT_V1V2_span50bp.bed"

samplenames=("RLTO_12_Child_C" "RLTO_12_Dad_B" "RLTO_12_Mom_A")

### Change the out_path accordingly ###
out_path="/path/to/results"
[ ! -d $out_path ] && mkdir -p $out_path
#######################################

for sample_name in "${samplenames[@]}"; do
	/gatk/gatk ApplyBQSR -R "${REF}" -I "${out_path}/${sample_name}.cram" --bqsr-recal-file "${out_path}/${sample_name}_report.txt" -O "${out_path}/${sample_name}_applyBQSR.cram" && \\
	/gatk/gatk HaplotypeCaller -I "${out_path}/${sample_name}_applyBQSR.cram" -O "${out_path}/${sample_name}_BP_germline.g.vcf.bgz" -R "${REF}" -L "${BED}" -ERC BP_RESOLUTION
done
```

##### `ReEND.sh`

```bash
#!/bin/bash

### Input/Output paths ###
samplenames=("RLTO_12_Child_C" "RLTO_12_Dad_B" "RLTO_12_Mom_A")

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

[ ! -d $out_path ] && mkdir -p $out_path

source /path/to/function/function_lst_ver1
# source /storage1/fs1/jin810/Active/testing/ztang/code/function/function_lst_ver1

for sample_name in "${samplenames[@]}"; do
	re-END "${out_path}/${sample_name}_BP_germline.g.vcf.bgz"
done
done
```


### Joint call → GVCF

##### `GATK_DBI_Joint_vqsr.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 12 -M 160GB -R 'gpuhost rusage[mem=160GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(broadinstitute/gatk:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"
BED="/storage1/fs1/jin810/Active/yung-chun/database/capture_bedfile/Exome-IDT_V1V2_span50bp.bed"

### Change the out_path accordingly ###
out_path="/path/to/results"
map_path="/path/to/Combine_DBI/DBImapfile.map"
#######################################

DBI_path="${out_path}/Combine_DBI"
tmp_dir="${DBI_path}/tmp"
[ ! -d $out_path ] && mkdir -p $out_path
[ ! -d $tmp_dir ] && mkdir -p $tmp_dir
[ ! -d $DBI_path ] && mkdir -p $DBI_path

# cd "${DBI_path}"

/gatk/gatk GenomicsDBImport --java-options "-Xmx160G -Xms160G" -R "${REF}" -L "${BED}" --merge-input-intervals --genomicsdb-workspace-path "${DBI_path}/joint_list_DBI_work" --batch-size 10 --reader-threads 10 --consolidate --tmp-dir "${tmp_dir}" --sample-name-map "${map_path}"

/gatk/gatk GenotypeGVCFs --java-options "-Xmx160G -Xms160G" -R "${REF}" -V "gendb:///${DBI_path}/joint_list_DBI_work" -O "${DBI_path}/RLTO_12_joint_list_BP_genotype.vcf.gz" --tmp-dir "${tmp_dir}"

/gatk/gatk --java-options "-Xmx160G -Xms160G" VariantRecalibrator -V "${DBI_path}/RLTO_12_joint_list_BP_genotype.vcf.gz" -O "${DBI_path}/RLTO_12_joint_list_joint_vqsr.recal" --tranches-file "${DBI_path}/RLTO_12_joint_list_joint_vqsr.tranches" --resource:omni,known=false,training=true,truth=true,prior=12.0 /storage1/fs1/jin810/Active/known_sites/1000G_omni2.5.hg38.vcf.gz --resource:1000g,known=false,training=true,truth=false,prior=10.0 /storage1/fs1/jin810/Active/known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /storage1/fs1/jin810/Active/known_sites/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /storage1/fs1/jin810/Active/known_sites/hapmap_3.3.hg38.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=12.0 /storage1/fs1/jin810/Active/known_sites/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an MQ -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 --mode BOTH

/gatk/gatk --java-options "-Xmx160G -Xms160G" ApplyVQSR -R "${REF}" -V "${DBI_path}/RLTO_12_joint_list_BP_genotype.vcf.gz" --recal-file "${DBI_path}/RLTO_12_joint_list_joint_vqsr.recal" --tranches-file "${DBI_path}/RLTO_12_joint_list_joint_vqsr.tranches" -O "${DBI_path}/RLTO_12_joint_list_joint_vqsr.vcf.gz" --mode BOTH

/gatk/gatk --java-options "-Xmx160G -Xms160G" SelectVariants -R "${REF}" -V "${DBI_path}/RLTO_12_joint_list_joint_vqsr.vcf.gz" -O "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs.vcf.gz" -select 'vc.isNotFiltered()'
```

##### `BCF_splitmulti_leftnorm_reID.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 12 -M 160GB -R 'gpuhost rusage[mem=160GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(biocontainers/bcftools:v1.5_cv3)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
tmp_dir="${DBI_path}/tmp"
[ ! -d $out_path ] && mkdir -p $out_path
[ ! -d $tmp_dir ] && mkdir -p $tmp_dir

/opt/conda/bin/bcftools norm -m-both -f "${REF}" -O z "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs.vcf.gz" | /opt/conda/bin/bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz"

/opt/conda/bin/bcftools index -t -f "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz"
```


### Quality Control

##### `KinshipSexcheck.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 4 -M 16GB -R 'gpuhost rusage[mem=16GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(sam16711/plink:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
QC_path="${out_path}/Quality_Control"

[ ! -d $out_path ] && mkdir -p $out_path
[ ! -d $QC_path ] && mkdir -p $QC_path

/bin/plink --double-id --geno 0.01 --genome --hwe 0.001 --maf 0.05 --snps-only --allow-extra-chr --vcf-half-call m --vcf "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz" --out "${QC_path}/RLTO_12_kinship"

/bin/plink --double-id --allow-extra-chr --vcf-half-call m --vcf "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz" --out "${QC_path}/RLTO_12_sexcheck1"

/bin/plink -split-x hg38 --make-bed --allow-extra-chr --bfile "${QC_path}/RLTO_12_sexcheck1" --out "${QC_path}/RLTO_12_sexcheck2"

/bin/plink --impute-sex ycount --make-bed --allow-extra-chr --bfile "${QC_path}/RLTO_12_sexcheck2" --out "${QC_path}/RLTO_12_sexcheck3"
```

##### `BamMetrics.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 4 -M 16GB -R 'gpuhost rusage[mem=16GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(sam16711/bam_metrics:v1)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"
BED="/storage1/fs1/jin810/Active/yung-chun/database/capture_bedfile/Exome-IDT_V1V2_span50bp.bed"

samplenames=("RLTO_12_Child_C" "RLTO_12_Dad_B" "RLTO_12_Mom_A")

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

QC_path="${out_path}/Quality_Control"

for sample_name in "${samplenames[@]}"; do
	/opt/bamMetrics -t 4 -r "${REF}" -b "${BED}" -o "${QC_path}/${sample_name}_bamMetrics.txt" "${out_path}/${sample_name}.cram"
done
```

>[! Do the following before proceed to the next script]
>- Copy-paste the following script(s) to your own folder:
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/liftover_hg38tob37_v2.py`
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/mergeBamMetrics.py`
>- Substitute all `/path/to/...` in the following scripts with the actual paths.

##### `MergeBamMetrics.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 4 -M 16GB -R 'gpuhost rusage[mem=16GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(spashleyfu/ubuntu18_vep104:hail_gsutil)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
QC_path="${out_path}/Quality_Control"

export JAVA_HOME="/opt/conda"

/opt/conda/bin/python3 /path/to/liftover_hg38tob37_v2.py "${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz" "${QC_path}/RLTO_12_PCA.vcf.bgz"

/opt/conda/bin/python3 /path/to/mergeBamMetrics.py /path/to/results/Quality_Control/merge_bammetrics_pathlist.txt "${QC_path}/RLTO_12_All_BamMetrics.tsv"
```

##### `PCA_Laser.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -n 4 -M 16GB -R 'gpuhost rusage[mem=16GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(dreammaerd/laser_trace_v2.04:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"
GDB="/storage1/fs1/jin810/Active/yung-chun/database/reference/LASER/HGDP_938.geno"
GRC="/storage1/fs1/jin810/Active/yung-chun/database/reference/LASER/HGDP_938.RefPC.coord"

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
QC_path="${out_path}/Quality_Control"

vcf2geno --inVcf "${QC_path}/RLTO_12_PCA.vcf.bgz" --out "${QC_path}/RLTO_12_PCA_All" 

trace -p "${QC_path}/my_parameterfile" -k 10 -K 20 -s "${QC_path}/RLTO_12_PCA_All.geno" -g "${GDB}" -c "${GRC}" -o "${QC_path}/RLTO_12_tracePCA"
```


<aside> 💡 Do the following before proceed to the next script:

- Copy-paste this script to your own folder: `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/hdgp_PCA.R`
    
- Substitute all `/path/to/...`in the following script with the actual paths. </aside>
    
- `PCA_Laser_R.sh`
    
    ```bash
    #!/bin/bash
    
    ## Request Docker ##
    # bsub -Is -n 4 -M 16GB -R 'gpuhost rusage[mem=16GB] span[hosts=1]' -G compute-jin810 -q general-interactive -a 'docker(dreammaerd/r-jupyter:4.2.1)' /bin/bash
    
    ### Reference paths DO NOT CHANGE ###
    REF="/storage1/fs1/jin810/Active/reference/Homo_sapiens_assembly38.fasta"
    GDB="/storage1/fs1/jin810/Active/yung-chun/database/reference/LASER/HGDP_938.geno"
    GRC="/storage1/fs1/jin810/Active/yung-chun/database/reference/LASER/HGDP_938.RefPC.coord"
    
    ### Change the out_path accordingly ###
    out_path="/path/tp/results"
    #######################################
    
    DBI_path="${out_path}/Combine_DBI"
    QC_path="${out_path}/Quality_Control"
    
    /opt/conda/bin/Rscript /path/to/hdgp_PCA.R -i "${QC_path}/RLTO_12_tracePCA.ProPC.coord" -r "${GRC}"
    ```
    

---

<aside> 💡 How do you interpret your quality control results? Take a look at [PLINK](https://zzz.bwh.harvard.edu/plink/summary.shtml) and [Hail](https://hail.is/docs/0.2/methods/relatedness.html) references for interpretation of kinship statistics, PCA, and many more. Google and anyone with QC experience in lab would be available for you to ask questions!

</aside>

## Monitor your jobs on RIS

Some useful commands for you to operate/monitor your jobs:

```bash
#see all submitted jobs
bjobs
#terminate jobs by job ID
bkill [job_id_1] [job_id_2] ...
```

# 3. De Novo Variant Calling

>[!Note]
>Now you’ve done with all the pre-steps. We’re going to move on to call the denovo variants and annotate them using various tools. But first, what are De Novo variants? Why are they rare and important?


## Script to use

>[!Do the following before proceed to the next script:]
>- Copy the following 3 scripts to your testing folder:
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/hail_denovo_v2.py`
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/VCF_to_VariantTable_with_ParentalGT_v4.py`
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/Test_Hail_annotation_v2.py`
>- You will see that we used pedigree files as inputs for some of the following steps, and here’s a [link to their technical documentation for PED files](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). These files have already been generated (by me), but all the information within the pedigree files should be based on your sexcheck results during quality control.
    
##### `Denovo.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -G compute-jin810 -q general-interactive -n 8 -R 'rusage[mem=32GB]' -a 'docker(spashleyfu/hail_vep_gnomad:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
Date=$(date +%Y%m%d)
jobname="RLTO_12_denovo_"$Date

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
Denovo_path="${out_path}/Denovo"
[ ! -d $Denovo_path ] && mkdir -p $Denovo_path

vcf="${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz"
ped="${Denovo_path}/RLTO_12_ped.fam"
out="${Denovo_path}/"$jobname
mode="ht"

/opt/conda/bin/python /path/to/hail_denovo_v2.py -i $vcf -o $out -p $ped -m $mode
```


# 4. Annotaion

>[!Note]
>Once we obtained the variants, we need to further annotate them before we can manually interpret whether they’re disease causing (pathogenic) or not. In this pipeline, we use CADD, VEP, etc. to annotate/filter the variants gotten to leave only those that are likely pathogenic.

## Script to use

### Annotate CADD

##### `VCF2VariantTable.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -G compute-jin810 -q general-interactive -n 8 -R 'rusage[mem=32GB]' -a 'docker(spashleyfu/hail_vep_gnomad:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
Date=$(date +%Y%m%d)
jobname="RLTO_12_denovo_"$Date

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
Denovo_path="${out_path}/Denovo"

vcf="${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID.vcf.gz"
ped="${Denovo_path}/RLTO_12_VariantTable.ped"

/opt/conda/bin/python /path/to/VCF_to_VariantTable_with_ParentalGT_v4.py $vcf $ped
```

##### `Annotate_cadd_null.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -G compute-jin810 -q general-interactive -n 8 -R 'rusage[mem=32GB]' -a 'docker(spashleyfu/hail_vep_gnomad:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###
pheno="/storage1/fs1/jin810/Active/testing/ztang/misc/test_RoTATion/results/Annotation/RLTO_12_pheno.txt"

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

DBI_path="${out_path}/Combine_DBI"
Denovo_path="${out_path}/Denovo"
Annotate_path="${out_path}/Annotation"
[ ! -d $Anno_path ] && mkdir -p $Anno_path

variant_table="${DBI_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID_variantTable.txt"
mode="step1"

# interactive mode for debugging
/opt/conda/bin/python /path/to/Test_Hail_annotation_v2.py -i $variant_table -o $Annotate_path -m $mode -p $pheno

# un-interactive mode for actual job execution
bsub -G compute-jin810 -q general -n 8 -R 'rusage[mem=32GB]' -oo /path/to/logs/%J.log -a 'docker(spashleyfu/hail_vep_gnomad:latest)' /opt/conda/bin/python /path/to/Test_Hail_annotation_v2.py -i $variant_table -o $Annotate_path -m $mode -p $pheno
```

### Upload to [CADD web](https://cadd.gs.washington.edu/)

```bash
cd /path/to/Annotation

# Output from previous step
file="20230613_GATK_DBImport_joint_CAKUT_WES_Germline_BP_20230612_joint_list_joint_vqsr_gs_sp_ln_reID_variantTable_annotation_step1_cadd_null.tsv"

## Split by size, maxium is 2MB for cadd upload
split -n 2 "$file" -d "${file::-4}"_ --suffix-length=3 --additional-suffix=".tsv"

## Upload the split files to <https://cadd.gs.washington.edu/score>
GRCh38-v1.6_58c333a3e6f7d0c6ef8244aad2c26a53.tsv.gz
GRCh38-v1.6_18c56bc037dc0a21ca7876f182b7d3f1.tsv.gz

## Need to merge two cadd_file together
cd /path/to/Annotation/

#f=20230404_joint_variantTable_annotation_cadd_null_fromcaddweb.tsv
of=20230613_CAKUT_WES_variant_table_annotation_cadd_null_caddweb.tsv

cat GRCh38-v1.6_58c333a3e6f7d0c6ef8244aad2c26a53.tsv > $of
cat GRCh38-v1.6_18c56bc037dc0a21ca7876f182b7d3f1.tsv | tail -n +3 >> $of

/opt/conda/bin/bgzip -@4 $of
# output: 20230613_CAKUT_WES_variant_table_annotation_cadd_null_caddweb.tsv.gz
```

### Annotate CADD web and Variant types

##### `Cadd_web_variant_type.sh`

```bash
#!/bin/bash

## Request Docker ##
# bsub -Is -G compute-jin810 -q general-interactive -n 8 -R 'rusage[mem=32GB]' -a 'docker(spashleyfu/hail_vep_gnomad:latest)' /bin/bash

### Reference paths DO NOT CHANGE ###

### Change the out_path accordingly ###
out_path="/path/to/results"
#######################################

Denovo_path="${out_path}/Denovo"
Annotate_path="${out_path}/Annotation"

## CADD web ##

hail_table="${Annotate_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID_variantTable_annotation_step1.ht" # output from Annotate_cadd_null.sh
cadd_web="/path/to/annotation_caddweb.tsv.gz"

/opt/conda/bin/python /path/to/Test_Hail_annotation_v2.py -i $hail_table -m "step2" -c $cadd_web

## Variant Types ##

step2_out="${Annotate_path}/RLTO_12_joint_list_joint_vqsr_gs_sp_ln_reID_variantTable_annotation_step1_annotation_step2.ht"
denovo_out="${Denovo_path}/RLTO_12_denovo_20230905.ht"

/opt/conda/bin/python /path/to/Test_Hail_annotation_v2.py -i $step2_out -m "step3" -d $denovo_out
```


# 5. Plotreads

>[! Do the following before you proceed to next step]
>- Copy over the following scripts to your own directory
>	- `/storage1/fs1/jin810/Active/testing/ztang/misc/test_RoTATion/Plotread_CMD.sh`
>		- Change the paths in ==line 134 and 149== to your actual path
>	- `/storage1/fs1/jin810/Active/testing/ztang/code/CAKUT/Execute_JOB_bsub_CMD.sh`
>		- Change the `path_to_scripts` below when you call this file to whichever you put in ==line149== you put in the file above.

##### Run Plotreads
```bash
export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /storage1/fs1/bga:/storage1/fs1/bga $HOME:$HOME'

step3file="/path/to/Annotation/20230613_GATK_DBImport_joint_CAKUT_WES_Germline_BP_20230612_joint_list_joint_vqsr_gs_sp_ln_reID_variantTable_annotation_step1_annotation_step2_annotation_step3_001_forplotread.txt"
mode="trio"
Date=$(date +%Y%m%d)
name="CAKUT_WES_Denovo_Plotread_"$Date

bash Plotread_CMD.sh -d $step3file -m $mode -n $name

bash Execute_JOB_bsub_CMD.sh -d "/path/to/CAKUT_WES_Denovo_Plotread_20230703_scripts"
```

# 6. Interpretation

#### Some rules you could use to manually filter out false positives:
- Check if only proband has variants at the specific chromosome position but not parents (because they should be unaffected disease wise).
- See whether alignment directions are random from top to bottom (can further check whether variants have ~50/50 ratio of forward and reverse sequences).
- Check quality of variants: are they full squares/rectangles? Full size squares mean ~99 GQ.
- Check whether there is repetitions around the variant site (TTTTT, or AGTAGTAGTAGT, etc.)
- Check whether the region is C/G heavy (A→red; C→blue; G→green; T→yellow)
- Check on [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat) whether there’s only ONE match with 100% (with total size).

If total number of variants is not too large, let someone else to check it (validate the results).