# Read mapping and variant calling

# 1. Read Mapping

Read mapping is the process of aligning short reads to a reference genome. This step is necessary to determine the origin of each read and to identify variations in the sequenced sample.

## 1.1 Index

```bash
srun --cpus-per-task=10 --mem=16G --time=2:00:00 bwa index -p /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz
```

- `bwa`: This is calling the BWA software.
- `index`: This is the subcommand to tell BWA to create an index of the reference genome.
- `p`: This option allows you to specify a prefix for the names of the index files that will be created. In this case, the prefix is `/scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1`.
- `/scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz`: This is the file path of the reference genome you are indexing. The reference genome should be in FASTA format. In this case, it's a gzipped FASTA file.

## 1.2 Mapping

```bash
module load bwa
module load SAMtools
```

Please see the command to map the trimmed reads against the first 10Mb of one chromosome. This step should take about 1h30. 

```bash
srun --cpus-per-task=40 --mem=80G --time=16:00:00 bwa mem -t 20 -r 'CM051614.1:1-10000000' /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114296.pair1.subsampled.truncated.gz /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114296.pair2.subsampled.truncated.gz | samtools view -F 4 -b - | samtools sort -@20 -O BAM -o ~/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1_sorted.bam -
```

- `bwa mem`: Perform alignment using BWA-MEM algorithm.
    - `t 20`: Use 20 threads.
    - `/scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1`: Path to the reference index.
    - `/scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN190240.pair1.truncated.gz`: Path to the first paired-end FASTQ file.
    - `/scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN190240.pair2.truncated.gz`: Path to the second paired-end FASTQ file.
- `|`: Pipe the output of the previous command to the next command.
- `samtools view`:
    - `F 4`: Skip reads with the 4 flag set (i.e., unmapped reads).
    - `b`: Output in BAM format.
- `samtools sort`:
    - `@20`: Use 20 threads for sorting.
    - `O BAM`: Specify output format as BAM.
    - `o ~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted.bam`: Specify output file.

```bash
srun --cpus-per-task=40 --mem=80G --time=16:00:00 bwa mem -t 20 -r 'CM051614.1:1-10000000' /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114297.pair1.subsampled.truncated.gz /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114297.pair2.subsampled.truncated.gz | samtools view -F 4 -b - | samtools sort -@20 -O BAM -o ~/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1_sorted.bam -
```

```bash
srun --cpus-per-task=40 --mem=80G --time=16:00:00 bwa mem -t 20 -r 'CM051614.1:1-10000000'  /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114393.pair1.subsampled.truncated.gz /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN114393.pair2.subsampled.fastq.gz | samtools view -F 4 -b - | samtools sort -@20 -O BAM -o ~/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1_sorted.bam -
```

```bash
srun --cpus-per-task=40 --mem=80G --time=16:00:00 bwa mem -t 20 -r 'CM051614.1:1-10000000'  /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN115950.pair1.subsampled.truncated.gz /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN115950.pair2.subsampled.truncated.gz | samtools view -F 4 -b - | samtools sort -@20 -O BAM -o ~/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1_sorted.bam -
```

```bash
srun --cpus-per-task=40 --mem=80G --time=16:00:00 bwa mem -t 20 -r 'CM051614.1:1-10000000'  /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1 /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN190240.pair1.subsampled.truncated.gz /scratch/projects/hpc-workshop/clouded_leopard/fastq_trimmed/NN190240.pair2.subsampled.truncated.gz | samtools view -F 4 -b - | samtools sort -@20 -O BAM -o ~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted.bam -
```

Table 1: List of read mapping software

| Software | Description | Algorithm(s) | Advantages |
| --- | --- | --- | --- |
| BWA | Widely used for mapping short reads to a reference genome. | BWA-backtrack, BWA-SW, BWA-MEM | Fast, efficient, and accurate. |
| Minimap2 | Versatile aligner for both DNA and RNA sequences, effective for long reads and spliced alignments. | Minimizer-based chaining and DP-based extension | Efficient for long reads, supports spliced alignment. |
| Bowtie2 | Ultrafast and memory-efficient tool for aligning short reads to long reference sequences. | FM-Index based on the Burrows-Wheeler transform | Extremely fast, requires little memory. |
| STAR | Specialized for aligning spliced RNA-seq reads to a reference genome. | Maximal Mappable Prefix (MMP) | Very fast, suited for RNA-seq, discovers novel splices. |
| HISAT2 | Efficient in handling spliced alignments and is often used for RNA-seq data. | Hierarchical Graph FM index | Efficient for spliced alignments, uses less memory. |

This table summarizes the main software tools used for genome mapping, describing their algorithms, advantages, and example commands. Note that SAMtools is not an aligner but is included here as it is crucial for manipulating the output files of alignment tools.

# 2. Post-processing and Variant Calling with GATK and Picard

The Genome Analysis Toolkit (GATK) and Picard are used for post-alignment processing and variant calling.

Once reads are aligned, the alignment file needs to be post-processed before variant calling. This includes sorting the file, marking duplicates, and recalibrating base quality scores. Variant calling is then performed to identify SNPs, indels, and other forms of variations.

## 2.1 Sorting and Marking Duplicates with GATK

Sorting the BAM file helps order the alignments by their reference genome coordinates, which makes downstream analysis faster and more efficient.

Marking duplicates is a critical step in the pipeline. During the sequencing process, PCR amplification is often used to generate a sufficient amount of DNA. However, this can lead to the creation of identical or near-identical duplicates. These duplicates can distort the variant calling process by making it seem like there are more reads supporting a variant than there actually are.

```bash
module load GATK
module load java/OpenJDK/17.0.7.0.7
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk MarkDuplicates I=/scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1_sorted.bam O=~/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1_sorted_marked.bam M=~/clouded_leopard/bam_files/NN114296_marked_dup_metrics.txt
```

- **`picard MarkDuplicates`** is the Picard command for marking duplicates.
- **`I=sorted.bam`** specifies the input file, which is the sorted BAM file from the previous step.
- **`O=marked.bam`** specifies the output file, which will contain the same alignments but with duplicates marked.
- **`M=marked_dup_metrics.txt`** specifies the file to which metrics on the numbers of duplicates should be written. This can be useful for quality control and checking how much of your data consists of duplicates.

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk MarkDuplicates I=/scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1_sorted.bam O=~/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1_sorted_marked.bam M=~/clouded_leopard/bam_files/NN114297_marked_dup_metrics.txt
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk MarkDuplicates I=/scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1_sorted.bam O=~/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1_sorted_marked.bam M=~/clouded_leopard/bam_files/NN114393_marked_dup_metrics.txt
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk MarkDuplicates I=/scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1_sorted.bam O=~/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1_sorted_marked.bam M=~/clouded_leopard/bam_files/NN115950_marked_dup_metrics.txt
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk MarkDuplicates I=/scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted.bam O=~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted_marked.bam M=~/clouded_leopard/bam_files/NN190240_marked_dup_metrics.txt
```

## 2.2 Create an index for the reference genome and marked BAM file:

```bash
module load SAMtools
```

**DO NOT RUN THE COMMAND BELOW FOR THE TUTORIAL**

```bash
samtools faidx /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz

```

```bash
srun --cpus-per-task=1 --mem=16G --time=4:00:00 samtools index ~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted_marked.bam
```

## 2.3 Base Recalibration with GATK

**IMPORTANT**: This steps is only performed for datasets where a well curated reference panel already exists. **We are not going to perform this step since we do not have a reference SNP panel.**

Base quality scores are estimates of the probability of a base being called incorrectly by the sequencer. However, these scores are not always accurate, and systematic errors can occur due to various factors related to the sequencing process and the sequencing machine itself. BQSR adjusts these scores based on the actual data and known variant sites, improving the accuracy of the variant calls in the next step of the pipeline.

```
gatk BaseRecalibrator -I marked.bam -R reference.fa --known-sites sites.vcf -O recal_data.table
gatk ApplyBQSR -R reference.fa -I marked.bam --bqsr-recal-file recal_data.table -O final.bam

```

BaseRecalibrator

- **`gatk BaseRecalibrator`** runs the BaseRecalibrator tool in GATK.
- **`I marked.bam`** specifies the input file, which should be the BAM file with duplicates marked.
- **`R reference.fa`** specifies the reference genome.
- **`-known-sites sites.vcf`** specifies a VCF file with known variant sites, which are used to distinguish sequencing errors from true variants.
- **`O recal_data.table`** specifies the output file, which is a table of recalibration parameters.

ApplyBQSR

- **`gatk ApplyBQSR`** runs the ApplyBQSR tool in GATK.
- **`R reference.fa`** specifies the reference genome.
- **`I marked.bam`** specifies the input file, which should be the BAM file with duplicates marked.
- **`-bqsr-recal-file recal_data.table`** specifies the recalibration file that was created in the previous step.
- **`O final.bam`** specifies the output file, which is a BAM file with adjusted base quality scores.

## 2.4 Variant Calling with GATK

Finally, variant calling is performed with GATK's HaplotypeCaller:

```bash
module load java/OpenJDK/17.0.7.0.7
module load GATK
```

### 2.4.1 HaplotypeCaller

HaplotypeCaller with gVCF as output

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk HaplotypeCaller -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -I /scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1_sorted_marked.bam -O ~/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1.g.vcf -ERC GVCF
```

1. `gatk`: Calls the Genome Analysis Toolkit (GATK).
2. `HaplotypeCaller`: Specifies the tool within GATK being used.
3. `R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz`: Specifies the reference genome to be used in the analysis. `R` stands for "reference".
4. `I input.bam`: Specifies the input file containing the aligned sequencing data. `I` stands for "input". You will need to replace `input.bam` with the path to your actual input file.
5. `O output.g.vcf`: Specifies the output file where the called variants will be written. `O` stands for "output". In this case, the output is a genomic VCF file (gVCF).
6. `ERC GVCF`: Stands for "Emit Ref Confidence GVCF". This option tells HaplotypeCaller to output a GVCF or genomic VCF file, which is a kind of VCF file that, in addition to the variants, includes reference confidence information (basically a quality score for the reference bases). This is useful especially in the context of joint genotyping across multiple samples, as it maintains information about which regions in the genome were confidently called as reference.

The gVCF format is particularly useful in large cohort studies where you are processing many samples. It allows for a more sophisticated and accurate joint variant calling on all your samples. The idea is that you first call variants in each individual sample using HaplotypeCaller in GVCF mode and then combine these gVCF files and perform joint genotyping on them in a subsequent step.

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk HaplotypeCaller -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -I /scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1_sorted_marked.bam -O ~/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1.g.vcf -ERC GVCF
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk HaplotypeCaller -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -I /scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1_sorted_marked.bam -O ~/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1.g.vcf -ERC GVCF
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk HaplotypeCaller -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -I /scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1_sorted_marked.bam -O ~/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1.g.vcf -ERC GVCF
```

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk HaplotypeCaller -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -I /scratch/projects/hpc-workshop/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1_sorted_marked.bam -O ~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1.g.vcf -ERC GVCF
```

You can specify whether you want **`HaplotypeCaller`** to produce a VCF or a gVCF using the **`-ERC`** (EmitRefConfidence) parameter. To produce a gVCF, you would use the argument **`-ERC GVCF`**.

It is common practice to split the genome into chromosomes to perform SNP calling step. The haplotype calling can take a long time. Please see below an **example** on how to run HaplotypeCaller for each chromosome.

```bash
#!/bin/bash

# Path to the reference genome
REFERENCE="/scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz"

# Path to the input BAM file
INPUT_BAM="input.bam"

# Array of chromosome names (adjust this to match the chromosome names in your reference genome)
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# Loop through each chromosome
for CHROMOSOME in ${CHROMOSOMES[@]}
do
    # Output file for this chromosome
    OUTPUT_VCF="output_${CHROMOSOME}.g.vcf"
    
    # Run GATK HaplotypeCaller for this chromosome
    gatk HaplotypeCaller -R $REFERENCE -I $INPUT_BAM -O $OUTPUT_VCF -ERC GVCF -L $CHROMOSOME
done
```

### 2.4.2 Merge VCF

Before joint genotyping, you'll need to combine the individual gVCFs into a single gVCF.

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk CombineGVCFs -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -V ~/clouded_leopard/bam_files/NN190240_clouded_leopard_CM051614.1.g.vcf -V ~/clouded_leopard/bam_files/NN114297_clouded_leopard_CM051614.1.g.vcf -V ~/clouded_leopard/bam_files/NN114296_clouded_leopard_CM051614.1.g.vcf -V ~/clouded_leopard/bam_files/NN114393_clouded_leopard_CM051614.1.g.vcf -V ~/clouded_leopard/bam_files/NN115950_clouded_leopard_CM051614.1.g.vcf -O /scratch/projects/hpc-workshop/clouded_leopard/bam_files/clouded_leopard_CM051614.1.g.vcf.gz
```

Repeat the `-V` option for each sample gVCF.

- `gatk CombineGVCFs`: This invokes the GATK tool `CombineGVCFs` for combining multiple gVCF files.
- `R`: This specifies the reference genome file that was used during the variant calling process. The reference genome should be in FASTA format.
- `V`: You can use this flag multiple times to specify each of the gVCF files you want to combine. In this example, two gVCF files named `sample1.g.vcf` and `sample2.g.vcf` are being combined.
- `O`: This specifies the output file where the combined gVCF will be written. In this example, the output file is named `combined.g.vcf`.

Instead of using `CombineGVCFs`, you can opt for a more scalable solution using `GenomicsDBImport`  to create a GenomicsDB, which is more efficient for handling large datasets with many samples. This is particularly useful when dealing with a large cohort. The command below is an **example**.

```bash
gatk GenomicsDBImport \
    --genomicsdb-workspace-path my_genomicsdb \
    -V sample1.g.vcf \
    -V sample2.g.vcf \
    --intervals my_interval_list.bed
```

### 2.4.3 GenotypeGVCF

Now, perform joint genotyping on the combined gVCF.

```bash
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk GenotypeGVCFs -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -V /scratch/projects/hpc-workshop/clouded_leopard/bam_files/clouded_leopard_CM051614.1.g.vcf.gz -O ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_final.vcf.gz
```

- `gatk GenotypeGVCFs`: This invokes the GATK tool `GenotypeGVCFs` for joint genotyping.
- `R`: This specifies the reference genome file that was used during the variant calling process. The reference genome should be in FASTA format.
- `V`: Specifies the input gVCF file that contains the combined genomic data from multiple samples or chromosomes. In this example, the combined file is `combined.g.vcf`.
- `O`: Specifies the output file where the final VCF containing the variant calls will be written. In this example, the output file is named `final.vcf`.

This command takes the combined gVCF file, and for each genomic position, it performs joint genotyping across all samples to produce final variant calls. This process accounts for sample-specific information and increases the power and accuracy of variant detection compared to calling variants in each sample separately.

If you have the database you can use the following command: 

```bash
gatk GenotypeGVCFs \
   -R /path/to/reference.fasta \
   -V gendb://my_genomicsdb \
   -O output.vcf
```

After this command is executed, the `final.vcf` file will contain the final variant calls and can be used for downstream analysis such as variant annotation, filtering, and various types of genetic analyses.

# 3. Variant Analysis with BCFtools

BCFtools  are used for further analysis and manipulation of variant call files.

Once variants are called, they can be filtered, annotated, and manipulated for downstream analyses. This can involve separating SNPs and indels, calculating variant statistics, and more.

## 3.1 Filtering variants with BCFtools

```bash
module load bcftools
```

SNPs and indels can be separated with BCFtools:

```
srun --cpus-per-task=10 --mem=16G --time=8:00:00 bcftools view -v snps /scratch/projects/hpc-workshop/clouded_leopard/bam_files/clouded_leopard_CM051614.1_final.vcf.gz > ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_snps.vcf.gz
```

```bash
srun --cpus-per-task=10 --mem=16G --time=8:00:00 bcftools view -v indels /scratch/projects/hpc-workshop/clouded_leopard/vcf/clouded_leopard_CM051614.1_final.vcf.gz > ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_indels.vcf.gz
```

# 4. Hard Filtering of Variants with GATK

After variant calling, the next step is variant filtering, which is done to remove potential false positives. One way to do this is through hard filtering, where variants are filtered out based on certain criteria.

Hard filtering involves setting thresholds for various metrics and then filtering out any variant that doesn't meet these thresholds. While this approach is simple and straightforward, it might not be the best way to filter variants, as it doesn't take into account the context of each variant. Other more sophisticated methods include variant quality score recalibration (VQSR), which uses machine learning to classify variants.

## 4.1 SNPs

To perform hard filtering with GATK, you can use the `VariantFiltration` tool. Here is an example command:

```
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk VariantFiltration -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -V ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_snps.vcf.gz -O ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_snps_filtered.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "my_filter"

```

Here's what each part of the command means:

- `gatk VariantFiltration` runs the VariantFiltration tool in GATK.
- `R reference.fa` specifies the reference genome.
- `V input.vcf` specifies the input VCF file.
- `O filtered.vcf` specifies the output file, which will be a VCF file with failing variants marked as filtered.
- `-filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0"` sets the criteria for filtering. In this case, variants are filtered out if the quality by depth (QD) is less than 2.0, the Fisher strand bias (FS) is more than 60.0, or the mapping quality (MQ) is less than 40.0.
- `-filter-name "my_filter"` assigns a name to this filter, which will be added to the FILTER column for failing variants.

After running this command, variants that do not meet the specified criteria will be marked as filtered in the output VCF file. They are not removed from the file, but their FILTER column will contain the filter name instead of the "PASS" label.

These filtering criteria are just examples, and the actual values and metrics used should be adapted according to the specifics of your experiment and the recommendations of the GATK Best Practices.

## 4.2 Indels

The hard filtering criteria mentioned in the previous example are indeed often used for SNP variants. However, indels (insertions and deletions) have a different error profile and thus typically require different filter expressions.

Here is an example of how you might set up hard filters for indels with GATK's `VariantFiltration`:

```
srun --cpus-per-task=40 --mem=16G --time=8:00:00 gatk VariantFiltration -R /scratch/projects/hpc-workshop/clouded_leopard/reference/mNeoNeb1.pri_genomic.fasta.gz -V ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_indels.vcf.gz -O ~/clouded_leopard/bam_files/clouded_leopard_CM051614.1_indels_filtered.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "indel_filter"

```

In this command:

- `gatk VariantFiltration` runs the VariantFiltration tool in GATK.
- `R reference.fa` specifies the reference genome.
- `V indels.vcf` specifies the input VCF file, which in this case should contain indel calls.
- `O filtered_indels.vcf` specifies the output file.
- `-filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"` sets the criteria for filtering indels. For indels, the filter for Fisher Strand Bias (FS) is often set higher than for SNPs.
- `-filter-name "indel_filter"` gives a name to this filter, which will be added to the FILTER column for failing variants.

Again, these are example filter settings. The actual thresholds will depend on your specific dataset and should be adjusted based on your understanding of the sequencing and alignment quality, as well as the acceptable false positive and false negative rates for your study.

For the remaining tutorials we are going to use a vcf file containing only the filtered snps. 

Table 2: List of variant calling softwares

| Software | Description | Algorithm(s) | Advantages |
| --- | --- | --- | --- |
| GATK HaplotypeCaller | Used for calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. | Local De Novo Assembly and Pair-HMMs | Highly accurate, widely used, good for both SNPs and indels. |
| FreeBayes | Bayesian haplotype-based polymorphism discovery and genotyping. | Bayesian, haplotype-based | Can handle complex variation. Handles both haploid and diploid data. |
| SAMtools mpileup and BCFtools | Used together for calling variants from aligned reads. | Consensus Caller | Simple, fast, and uses less memory compared to others. |
| VarScan | Variant calling and somatic mutation/CNV detection for exome-scale data. | Heuristic/statistic based | Good for somatic mutation detection, handles tumor-normal pairs. |
| SNVer | Statistical tool for calling common and rare variants via log-likelihood ratio test. | Maximum likelihood estimation (MLE) based | Capable of handling pooled samples. |
| DeepVariant | Uses deep learning to call genetic variants from next-generation DNA sequencing data. | Deep learning (convolutional neural networks) | High accuracy, particularly good for tricky variants. |

This table summarizes the main software tools used for variant calling, describing their algorithms, advantages, and example commands. These tools are widely used for identifying genetic variants such as SNPs and indels from aligned sequencing reads. Note that some of the tools like GATK are more resource-intensive but tend to be more accurate, while others like SAMtools are lighter but might not be as comprehensive.

### What should you do next?

After performing SNP calling, it is essential to perform quality control and filtering to ensure that the data is of high quality before conducting population genetic analyses. Below are the steps commonly taken:

1. **Quality Filtering**: You might want to filter out low-quality SNPs and reads. You can use software like VCFtools, BCFtools, or GATK to do so. Typically, you might set a threshold for the quality values (e.g., QD, FS, MQ) and only retain SNPs above those thresholds.
2. **Missing Data Filtering**: It is also important to filter out SNPs with too much missing data, as they can bias your results. You might want to remove SNPs that are not genotyped in a certain percentage of individuals.
3. **Minor Allele Frequency Filtering (MAF)**: SNPs with very low allele frequencies might be due to sequencing errors. Filtering based on Minor Allele Frequency (e.g., MAF > 0.01) can be a good practice.
4. **Linkage Disequilibrium (LD) Pruning**: This is the process of removing SNPs that are in high linkage disequilibrium with each other. Too much LD can bias certain population genetics statistics. You can use tools like PLINK to remove SNPs in high LD.
5. **Filtering Related Individuals**: If your dataset includes related individuals, you might want to remove them or account for them in your analysis, as they can bias the results of association studies. Software like KING can be used to estimate relatedness and help filter related individuals.
6. **Population Stratification**: It is important to account for population structure in your dataset, as it can also bias the results. Software like EIGENSOFT can be used to correct for population stratification.
7. **Heterozygosity Filtering**: It might be a good idea to remove individuals or SNPs with an unusually high or low rate of heterozygosity.
8. **Spatial Filtering**: If you are working with geographically structured populations, it might be beneficial to also consider spatial filtering techniques to remove SNPs that are highly spatially autocorrelated.
9. **Functional Annotation Filtering**: Depending on your study, you may want to retain only SNPs within certain genomic features (e.g., coding regions, regulatory regions).

Once you've filtered your data and ensured its quality, you can move on to various population genetic analyses like calculating fixation indices (F_ST), performing population structure analysis (e.g., PCA, ADMIXTURE), demographic inference, phylogenetic analysis, and genome-wide association studies (GWAS).

# References:

- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18):3094-3100. doi: 10.1093/bioinformatics/bty191.
- McKenna A et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9):1297-1303. doi: 10.1101/gr.107524.110.
- Broad Institute (2023). Picard. GitHub. https://github.com/broadinstitute/picard
- Li H et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16):2078-2079. doi: 10.1093/bioinformatics/btp352.
- Danecek P et al. (2011). The variant call format and VCFtools. Bioinformatics, 27(15):2156-2158. doi: 10.1093/bioinformatics/btr330.
