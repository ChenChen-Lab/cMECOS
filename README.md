# cMECOS
This document outlines the computational pipeline used for processing metagenomic data, reconstructing microbial genomes, and analyzing phage-host interactions. 

```bash

# Core tools
conda install -c bioconda \
  soapfilter=2.2 \
  bowtie2=2.2.5 \
  spades=3.15.3 \
  bwa=0.7.15 \
  samtools=1.9 \
  metabat2=2.12.1 \
  checkm-genome=1.2.0 \
  prokka=1.13.3 \
  mmseqs2=13.45111 \
  blast=2.9.0 \
  mash=2.0 \
  gtdbtk=1.7.0 \
  fasttree=2.1.10

# Specialized tools
conda install -c bioconda \
  athena=1.2 \
  drep=3.4.0 \
  barrnap=0.9 \
  eggnog-mapper=2.0.1 \
  crisprcasfinder=2.0.2 \
  carveme=1.5.1 \
  smetana=1.0

# Python packages
pip install igraph==0.9.8 scikit-learn==0.24.2

```

## 1. Reads Processing and Genome Binning

### 1.1 Quality Control & Host DNA Removal
- **Tool**: SOAPfilter v2.2  
  - Removes low-quality reads and adapter sequences
- **Host Filtering**:  
  - Bowtie2 v2.2.5 (`--very-sensitive` mode)  
  - Reference: Human genome hg38
  
```{bash eval=FALSE,include=FALSE, echo=TURE}

# Quality filtering
SOAPfilter_v2.2 -F CTGTCTCTTATACACATCT -R CTGTCTCTTATACACATCT \
  -q 20 -p 0.1 -M 2 -f input.fq -o clean.fq
# Host DNA removal
bowtie2 --very-sensitive -x hg38_index -U clean.fq \
  --un nohost.fq > host_mapped.sam 2> bowtie2.log  
  
```


### 1.2 Metagenome Assembly
1. **Primary Assembly**:  
   MetaSPAdes v1.1.3 (default parameters)
2. **Assembly Refinement**:  
   - Read mapping with BWA v0.7.15 (`mem` mode)  
   - BAM sorting with Samtools v1.9  
   - Co-barcode optimization using ATHENA v1.2
   
```{bash eval=FALSE,include=FALSE, echo=TURE}

# Primary assembly
metaspades.py -k 21,33,55 -t 32 -o assembly_dir \
  --pe1-1 read1.fq --pe1-2 read2.fq

# Co-barcode refinement
bwa mem -t 16 contigs.fa reads.fq | samtools sort -@4 -o aligned.bam
athena --threads 32 --config config.json
  
```
### 1.3 Genome Binning
- **Binner**: MetaBAT2 v2.12.1  
- **Quality Control**:  
  CheckM v1.2.0 with criteria:  
Completeness >50%
Contamination <5%
Quality Score = Completeness - 5×Contamination >50
```{bash eval=FALSE,include=FALSE, echo=TURE}
# Coverage calculation
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

# MetaBAT2 binning
metabat2 -i contigs.fa -a depth.txt -o bins_dir/bin \
  -m 1500 -t 24 --seed 42

# Quality assessment
checkm lineage_wf -x fa bins_dir/ checkm_out/ \
  --tab_table -f checkm_results.tsv --threads 16
```  

## 2. Species/Strain Clustering

### 2.1 Reference Database
- UHGG database (v2.0): [Download Link](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0)

### 2.2 Clustering Workflow
1. **Species-Level (95% ANI)**:  
 dRep v3.4.0 `-pa 0.95 -sa 0.99 -nc 0.30 -cm larger`
2. **Strain-Level (99% ANI)**:  
 - MASH v2.0 (k=21, sketch=10k)  
 - Agglomerative clustering with scikit-learn (single linkage)

```{bash eval=FALSE,include=FALSE, echo=TURE}

second_cluster_main.py my_MAG.txt Cdb.csv UHGG_28W_with_us/ /disk_RAID5_27T/UHGG2/genomes-all_metadata.tsv  strain_cluster_res 0.99 190
```

### 2.3 Taxonomic Assignment
- **Tool**: GTDB-Tk v1.7.0 (r202 database)
- **Phylogenetics**:  
- 120 marker genes alignment  
- FastTree v2.1.10 for ML tree  
- pplacer for taxonomic placement
```{bash eval=FALSE,include=FALSE, echo=TURE}
# GTDB-Tk workflow
gtdbtk classify_wf --genome_dir genomes/ --out_dir gtdbtk_out \
  --cpus 32 --prefix output

```

## 3. rRNA Gene Analysis
- **Prediction**: barrnap v0.9
- **Database Integration**:  
Silva + Greengene + rrnDB → MMseqs2 (80% coverage, 100% identity)
- **Annotation**: BLAST+ v2.9.0

```{bash }
for i in `find raw -name '*.ffn'`;do x=`echo ${i}|sed -e 's/\.ffn//g' -e 's/.*\///g'`;echo "~/anaconda3/bin/seqkit grep -j 192 -n -r -p  '16S ribosomal RNA$' ${i}|~/anaconda3/bin/seqkit replace -p '.*' -r \"${x}:_:\" |~/anaconda3/bin/seqkit rename|sed '/^>/ s/ $//g'";done|parallel -j 192  > all_full_len_16S.fa


makeblastdb -in combine3.fa -dbtype nucl  -out all_16S_db

blastn -db /disk_RAID6_180T/Workstation/Chenchen2023_2/16S_database_punblic/all_16S_db -num_threads 192 -task megablast -query all_full_len_16S.fa  -max_target_seqs 3 -out all_16s.m6 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus'
awk -F "\t" '!a[$1]++{print}' 1.blast > uniq.blast

```

## 4. Gene Catalog Construction

### 4.1 CDS Prediction
- **Method**: Prokka v1.13.3
- Prodigal v2.6.3 (`-p single -c -m`)


for i in `ls 5006_Mine_MAG/*.fa`;do x=`echo $i|sed -e 's/^5006_Mine_MAG\/pre\.//g' -e 's/\.fa//g'`;echo "prokka --outdir 5006_mine_annote_dir/raw/prokka_${x} --prefix ${x} ${i} --cpus 4";done| ~/anaconda3/bin/parallel -j 50



### 4.2 Protein Clustering
- **Method**: MMseqs2 
- Clustering levels: 100%, 95%, 90% identity
- Minimum coverage: 80% of shorter sequence


mkdir -p mmseq_dir/{DB,Clu,tsv}
mmseqs createdb all_raw.faa mmseq_dir/DB/mm_seq_DB

#100%
mmseqs linclust mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_100 mmseq_dir/tmp_mmseq_100 --kmer-per-seq 80 --cov-mode 1 -c 0.8 --min-seq-id 1 --threads 192
mmseqs createtsv mmseq_dir/DB/mm_seq_DB mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_100 mmseq_dir/tsv/mm_seq_clu_100.tsv
mmseqs createsubdb mmseq_dir/Clu/mm_seq_clu_100 mmseq_dir/DB/mm_seq_DB mmseq_dir/mm_seq_clu_100_rep
mmseqs convert2fasta mmseq_dir/mm_seq_clu_100_rep mmseq_dir/mm_seq_clu_100_rep.faa
#95%
mmseqs linclust mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_95 mmseq_dir/tmp_mmseq_95 --kmer-per-seq 80 --cov-mode 1 -c 0.8 --min-seq-id 0.95 --threads 192
mmseqs createtsv mmseq_dir/DB/mm_seq_DB mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_95 mmseq_dir/tsv/mm_seq_clu_95.tsv
mmseqs createsubdb mmseq_dir/Clu/mm_seq_clu_95 mmseq_dir/DB/mm_seq_DB mmseq_dir/mm_seq_clu_95_rep
mmseqs convert2fasta mmseq_dir/mm_seq_clu_95_rep mmseq_dir/mm_seq_clu_95_rep.faa
#90%
mmseqs linclust mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_90 mmseq_dir/tmp_mmseq_90 --kmer-per-seq 80 --cov-mode 1 -c 0.8 --min-seq-id 0.9 --threads 192
mmseqs createtsv mmseq_dir/DB/mm_seq_DB mmseq_dir/DB/mm_seq_DB mmseq_dir/Clu/mm_seq_clu_90 mmseq_dir/tsv/mm_seq_clu_90.tsv
mmseqs createsubdb mmseq_dir/Clu/mm_seq_clu_90 mmseq_dir/DB/mm_seq_DB mmseq_dir/mm_seq_clu_90_rep
mmseqs convert2fasta mmseq_dir/mm_seq_clu_90_rep mmseq_dir/mm_seq_clu_90_rep.faa


### 4.3 Functional Annotation
- **Method**: eMapper 
emapper.py -m diamond --no_annot --no_file_comments --data_dir /dev/shm --seed_ortholog_evalue 0.00001 --cpu 192 -i novel_90.faa -o ./emapper --override
emapper.py --annotate_hits_table emapper.emapper.seed_orthologs --no_file_comments -o ./anno_ --cpu 192 --data_dir /dev/shm --override

## 5. Phage-Host Network

CRISPR Spacer Analysis
- **Prediction**: CRISPRCasFinder v2.0.2 (evidence levels 3-4)
- **Phage Matching**:  
BLAST+ vs MGV database  
Criteria: ≤1 mismatch/gap over ≥95% spacer length

conda activate CRISPRCasFinder
mkdir CRISPR_FIND_ALL
for i in `find -name '*.fa'`;do x=`echo ${i}|sed -e 's/^\.\/pre\.//g' -e 's/\.fa$//g' -e 's/^\.\///g'`;echo "extract_spacers_and_count_CRISPER.sh ${i} CRISPR_FIND_ALL/${x} 3 1";done|parallel -j 180 


## 6. Metabolic Modeling
- **GEM Reconstruction**: CarveMe v1.5.1
- **Interaction Metrics**:  
SMETANA v1.0 (MIP, MRO, SMETANA scores)
- **Validation**:  
100× subsampling (n=2-10 species)  
Wilcoxon test (FDR <0.05)
```{bash}
# GEM reconstruction
carve genome.faa --gapfill -v -o model.xml

# SMETANA analysis
smetana -d . -o interactions/ \
  --mediadb media.tsv --flavor cobra
```

---

## 7. Visualization  
The R script [`FiguresFinal.R`](analysis/FiguresFinal.R) contains the complete code for generating all figures presented in this study, including taxonomic composition bar plots, phage-host interaction networks, metabolic heatmaps, and clinical association boxplots.
