<h2>Nextflow pipeline for whole genome sequence analysis of Human genome (Chr2) using NanoPore/Illumina Hybrid sequence </h2>

<h3>Software used in the pipeline.</h3>

**FastQC** is widely used and is robust, efficient, and versatile quality control software for a varied range of raw genetic data. It outputs a quality report which can be viewed to give an indication on how good the respective reads are. For more information, please visit the FastQC website https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

**MultiQC** is used create a single report with interactive plots for multiple bioinformatics analyses across many samples. It reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control. For more information, please visit the MultiQC website https://github.com/MultiQC/MultiQC.

**NanoQC** is a plotting tool for long read sequencing data and alignments. For more information, please check its website: https://biocontainers.pro/tools/nanoplot and https://github.com/wdecoster/NanoPlot. <br>
<img width="400" height="400" alt="LengthvsQualityScatterPlot_legacy_kde" src="https://github.com/user-attachments/assets/1bcb8557-c3c7-451b-a7db-7ea96fab3aa5" />

**Hifiasm** 
is a fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads. It produces arguably the best single-sample telomere-to-telomere assemblies combing HiFi, ultralong and Hi-C reads, and it is one of the best haplotype-resolved assemblers for the trio-binning assembly given parental short reads. For more information, please visit https://github.com/chhylp123/hifiasm. <br> 
&ensp; **Completeness:** <br>
  &emsp;  **-Merqury** is a tool for validating genomes assembled using long read sequencing. The algorithm uses k-mer<br> &emsp; &ensp;to evaluate base accuracy and completeness of a genome by comparing the de novo assembled genome <br> &emsp; &ensp;with high accuracy reads that are not used in the assembly. It does not require another reference assembly <br> &emsp; &ensp;for its evaluation. It is able to evaluate the haplotype-specific accuracy,completeness, phase block  <br>
  &emsp; &ensp;continuity and switch errors in trio binned assemblies https://github.com/marbl/merqury. <br>
  <img width="400" height="400" alt="LengthvsQualityScatterPlot_legacy_kde" src="https://github.com/user-attachments/assets/234120a5-6654-428a-b638-a9a43b229940" />

  
  &emsp; **-QUAST** (QUality ASsessment Tool) is a widely used software for evaluating and comparing genome  <br>
  &emsp; &ensp;assemblies in bioinformatics, providing metrics on assembly quality, such as contiguity and completeness.  <br>
  &emsp; &ensp;The tool can be used with or without a reference genome and generates detailed reports, tables, and plots  <br>
  &emsp; &ensp;for analysis https://github.com/ablab/quast.

**Maker2** 
is a genome annotation and data management tool designed for second-generation genome projects. It is a multi-threaded, parallelized application that can process second-generation datasets of virtually any size. It can produce accurate annotations for novel genomes where training-data are limited, of low quality or even non-existent. . For more information, please visit the website https://github.com/Yandell-Lab/maker.

**BUSCO**
 (Benchmarking Universal Single-Copy Orthologs) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness based on evolutionarily informed expectations of gene content from near-universal single-copy orthologs, please visit https://busco.ezlab.org/.


