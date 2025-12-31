#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------
// Default parameters
// ----------------------------
params.illumina_dir       = 'data/illumina'
params.ont_dir            = 'data/NanoPore'
params.ref_dir            = 'data/Reference'
params.outdir             = 'results'

params.threads            = 8
params.racon_rounds       = 2
params.racon_map_preset   = 'map-ont'
params.pilon_rounds       = 1
params.kmer               = 21
params.busco_lineage      = ''

params.targetLongReadCov  = 100
params.genomeSize         = 3000000000
params.targetBases        = params.targetLongReadCov * params.genomeSize

/********************************************************************************
 * Channels
 ********************************************************************************/
Channel.fromFilePairs("${params.illumina_dir}/*_R{1,2}.fastq.gz")
       .set { illumina_pairs }  // set(sample_id, [R1, R2])

Channel.fromPath("${params.ont_dir}/*.fastq.gz")
       .ifEmpty { error "No ONT reads found in ${params.ont_dir}" }
       .set { ont_reads }

Channel.fromPath("${params.ref_dir}/*.fasta")
       .first()
       .set { ref_fasta }

/********************************************************************************
 * 1) Illumina QC (FastQC)
 ********************************************************************************/
process FASTQC {
    tag { sample_id }
    cpus 2
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_output", emit: reports

    script:
    """
    mkdir -p fastqc_output
    fastqc -q -t ${task.cpus} -o fastqc_output ${reads.join(' ')}
    """
}

/********************************************************************************
 * 2) Raw ONT QC (NanoPlot)
 ********************************************************************************/
process NANOPLOT_RAW {
    tag { ont_files.baseName }
    cpus 4
    publishDir "${params.outdir}/nanoplot_raw", mode: 'copy'

    input:
    path ont_files

    output:
    path "nanoplot_raw_output", emit: reports

    script:
    """
    NanoPlot --fastq ${ont_files} \
		--outdir nanoplot_raw_output \
		--threads ${task.cpus} \
        --loglength --plots dot kde --verbose
    """
}

/********************************************************************************
 * 3) Filtlong (combine all ONT reads)
 ********************************************************************************/
process FILTLONG {
    tag "filtlong"
    cpus 4
    publishDir "${params.outdir}/filtlong", mode: 'copy'

    input:
    path ont_files

    output:
    path "long_reads.filtered.fastq.gz", emit: filtered_long

    conda '/home/anaconda3/envs/FILTLONG'

    script:
    """
    filtlong \
        --min_length 1000 \
        --keep_percent 90 \
        --length_weight 0.5 \
        --target_bases ${params.targetBases} \
        ${ont_files} \
        > long_reads.filtered.fastq

    gzip -c long_reads.filtered.fastq > long_reads.filtered.fastq.gz
    """
}

/********************************************************************************
 * 4) NanoPlot filtered ONT reads
 ********************************************************************************/
process NANOPLOT_FILTERED {
    tag "nanoplot_filtered"
    cpus 4
    publishDir "${params.outdir}/nanoplot_filtered", mode: 'copy'

    input:
    path filtered_fastq

    output:
    path "nanoplot_filtered_output", emit: reports

    script:
    """
    NanoPlot --fastq ${filtered_fastq} \
        --outdir nanoplot_filtered_output \
        --threads ${task.cpus} \
        --loglength --plots dot kde --verbose
    """
}

/********************************************************************************
 * 5) Ratatosk hybrid correction
 ********************************************************************************/
process RATATOSK {
    tag { sample_id }
    publishDir "${params.outdir}/ratatosk", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(longreads)

    output:
    path "long_reads.corrected.fastq.gz", emit: corrected_long

    conda '/home/anaconda3/envs/RATATOSK'

    script:
    """
    echo "Running Ratatosk on sample: ${sample_id}"

    Ratatosk correct -v -G -c 16 \
        -l ${longreads} \
        -s ${reads[0]} ${reads[1]} \
        -o long_reads.filtered.corrected

    if [ -f long_reads.filtered.corrected.fastq ]; then
        gzip -c long_reads.filtered.corrected.fastq > long_reads.corrected.fastq.gz
    elif [ -f long_reads.filtered.corrected.fastq.gz ]; then
        mv long_reads.filtered.corrected.fastq.gz long_reads.corrected.fastq.gz
    else
        echo "Ratatosk did not produce expected output files"
        exit 1
    fi
    """
}

/********************************************************************************
 * 6) NanoPlot corrected ONT reads
 ********************************************************************************/
process NANOPLOT_CORRECTED {
    tag "nanoplot_corrected"
    cpus 4
    publishDir "${params.outdir}/nanoplot_corrected", mode: 'copy'

    input:
    path corrected

    output:
    path "nanoplot_corrected_output", emit: reports

    script:
    """
    NanoPlot --fastq ${corrected} --outdir nanoplot_corrected_output --threads ${task.cpus} --loglength --plots dot kde --verbose
    """
}

/********************************************************************************
 * 7) Flye assembly
 ********************************************************************************/
process FLYE {
    tag "flye"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path longcorr

    output:
    path("assembly.fasta"), emit: draft_assembly

    conda '/home/anaconda3/envs/flye_env'

    script:
    """
    mkdir -p flye_out
    flye --nano-raw ${longcorr} --out-dir flye_out --threads ${task.cpus} --genome-size ${params.genomeSize} --asm-coverage 50
    cp flye_out/assembly.fasta assembly.fasta
    """
}

/********************************************************************************
 * 8) Racon polishing
 ********************************************************************************/
process RACON {
    tag "racon"
    publishDir "${params.outdir}/polish", mode: 'copy'

    input:
    path assembly
    path longreads

    output:
    path "assembly.racon2.fasta", emit: racon_polished

    conda '/home/anaconda3/envs/RACON'

    script:
    """
    # Round 1
    minimap2 -x map-ont ${assembly} ${longreads} > round1.paf

    racon -m 8 -x -6 -g -8 -w 500 ${longreads} round1.paf ${assembly} > round1.fasta

    # Round 2
    minimap2 -x map-ont round1.fasta ${longreads} > round2.paf

    racon -m 8 -x -6 -g -8 -w 500 ${longreads} round2.paf round1.fasta > assembly.racon2.fasta
    """
}

/********************************************************************************
 * 9) Pilon polishing
 ********************************************************************************/
process PILON {
    tag { sample_id }
    publishDir "${params.outdir}/polish", mode: 'copy'

    input:
    path assembly
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.pilon.fasta", emit: pilon_polished
    path "${sample_id}.pilon.changes", emit: pilon_changes

    conda '/home/anaconda3/envs/PILON'

    script:
    """
    set -euo pipefail

    bwa index ${assembly}
    samtools faidx ${assembly}

    bwa mem -t ${task.cpus} ${assembly} ${reads[0]} ${reads[1]} |
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam

    set +e
    java -Xmx64G -jar ~/anaconda3/envs/PILON/share/pilon-1.24-0/pilon.jar \\
        --genome ${assembly} \\
        --frags ${sample_id}.sorted.bam \\
        --output ${sample_id}.pilon \\
        --threads ${task.cpus} \\
        --changes
    PILON_EXIT=\$?
    set -e

    if [[ \$PILON_EXIT -ne 0 && \$PILON_EXIT -ne 1 ]]; then
        echo "Pilon failed with exit code \$PILON_EXIT"
        exit \$PILON_EXIT
    fi
    """
}

/********************************************************************************
 * 10) Purge_Dups
 ********************************************************************************/
process PURGE_DUPS {
    tag "purge_dups"
    cpus 1
    publishDir "${params.outdir}/purge_dups", mode: 'copy'

    input:
    path assembly
    path longreads

    output:
    path "assembly.purged.fasta", emit: purged_assembly
    path "purged.bed", emit: purged_bed

    conda '/home/anaconda3/envs/purge_dups'

    script:
    """
    set -euo pipefail

    # ---- Assembly size check ----
    ASM_SIZE=\$(awk '{sum+=length(\$0)} END {print sum}' ${assembly})

    if [ "\$ASM_SIZE" -lt 50000000 ]; then
        echo "Assembly too small for purge_dups (\$ASM_SIZE bp). Skipping."
        cp ${assembly} assembly.purged.fasta
        touch purged.bed
        exit 0
    fi

    # ---- Map ONT reads ----
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${longreads} \
        | samtools sort -@ ${task.cpus} -o aln.sorted.bam
    samtools index aln.sorted.bam

    # ---- Self-alignment ----
    minimap2 -x asm5 -t ${task.cpus} ${assembly} ${assembly} > self.paf

    # ---- Coverage stats ----
    pbcstat aln.sorted.bam || true
    calcuts PB.stat > cutoffs || true

    # ---- Purge duplicates ----
    purge_dups -T cutoffs self.paf > purged.bed || true

    # ---- Extract sequences ----
    if [ -s purged.bed ] && command -v get_seqs >/dev/null 2>&1; then
        get_seqs purged.bed ${assembly} > assembly.purged.fasta || cp ${assembly} assembly.purged.fasta
    else
        cp ${assembly} assembly.purged.fasta
        touch purged.bed
    fi
    """
}

/********************************************************************************
 * 11) RagTag scaffolding
 ********************************************************************************/
process RAGTAG {
    tag "ragtag"
    publishDir "${params.outdir}/ragtag", mode: 'copy'

    input:
    tuple path(asm), path(ref)

    output:
    path "assembly.ragtag.fasta", emit: ragtag_assembly

    conda '/home/anaconda3/envs/ragtag'

    script:
    """
    ragtag.py scaffold ${ref} ${asm} -o ragtag_out -t ${task.cpus}
    cp ragtag_out/ragtag.scaffold.fasta assembly.ragtag.fasta
    """
}

/********************************************************************************
 * 12) TGS-GapCloser
 ********************************************************************************/
process TGS_GAPCLOSER {
    tag "tgs-gapcloser"
    publishDir "${params.outdir}/gapcloser", mode: 'copy'

    input:
    path asm
    path longreads

    output:
    path("assembly.gapclosed.fasta"), emit: gapclosed_assembly

    conda '/home/anaconda3/envs/TGS_GAPCLOSER'

    script:
    """
    echo -e "[general]\\nproject=gapcloser\\nthread=${task.cpus}\\nref=${asm}\\nreads=${longreads}" > gapcloser.cfg
    TGS-GapCloser -i gapcloser.cfg -o gapcloser_out -t ${task.cpus} || true

    if [ -f gapcloser_out/final.fasta ]; then
        cp gapcloser_out/final.fasta assembly.gapclosed.fasta
    else
        cp ${asm} assembly.gapclosed.fasta
    fi
    """
}

/********************************************************************************
 * 13) QUAST
 ********************************************************************************/
process QUAST {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple path(asm), path(ref)

    output:
    path "quast_report", emit: reports

    script:
    """
    if [ -f "${ref}" ]; then
        quast.py ${asm} -R ${ref} -o quast_report -t ${task.cpus}
    else
        quast.py ${asm} -o quast_report -t ${task.cpus}
    fi
    """
}

/********************************************************************************
 * 14) Merqury
 ********************************************************************************/
process MERQURY {
    tag { sample }
    publishDir "${params.outdir}/merqury/${sample}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(r1), path(r2), path(final_assembly)

    output:
    path "${sample}.meryl", emit: meryl_db
    path "merqury_output", emit: reports

    conda '/home/anaconda3/envs/MERQURY'

    script:
    """
    set -euo pipefail

    k=21

    # Build and merge meryl DBs in one logical flow
    meryl k=\$k count output ${sample}.meryl.count ${r1} ${r2}
    meryl union-sum output ${sample}.meryl ${sample}.meryl.count

    # Run Merqury (non-trio)
    mkdir -p merqury_output
    merqury.sh ${sample}.meryl ${final_assembly} merqury_output/MerquryEval || echo "Merqury completed with warnings"
    
    # Clean up intermediate files
    rm -rf ${sample}.meryl.count
    """
}

/********************************************************************************
 * 16) MultiQC - SIMPLE ORIGINAL VERSION
 ********************************************************************************/
process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    multiqc . -n multiqc_report.html || echo "MultiQC completed"
    """
}

/********************************************************************************
 * Workflow
 ********************************************************************************/
workflow {
    // Collect ONT reads
    ont_all = ont_reads.collect()

    // 1. FastQC on Illumina
    fastqc_out = FASTQC(illumina_pairs)

    // 2. Raw NanoPlot ONT
    nanoplot_raw_out = NANOPLOT_RAW(ont_all)

    // 3. Filtlong
    filtlong_out = FILTLONG(ont_all)

    // 4. NanoPlot filtered
    nanoplot_filtered_out = NANOPLOT_FILTERED(filtlong_out)

    // 5. Ratatosk hybrid correction
    ratatosk_input = illumina_pairs.combine(filtlong_out)
    ratatosk_out = RATATOSK(ratatosk_input)

    // 6. NanoPlot corrected ONT
    nanoplot_corrected_out = NANOPLOT_CORRECTED(ratatosk_out)

    // 7. Flye assembly
    flye_out = FLYE(ratatosk_out)

    // 8. Racon polishing
    racon_out = RACON(flye_out, ratatosk_out)

    // 9. Pilon polishing
    pilon_out = PILON(
        assembly = racon_out,
        reads = illumina_pairs
    )

    // 10. Purge_Dups
    purge_out = PURGE_DUPS(
        pilon_out.pilon_polished,
        ratatosk_out.corrected_long
    )

    // 11. RagTag scaffolding
    ragtag_ch = purge_out.purged_assembly
        .ifEmpty { error "No purged assemblies found!" }
        .combine(ref_fasta)

    ragtag_out = RAGTAG(ragtag_ch)

    // 12. TGS-GapCloser
    gapcloser_out = TGS_GAPCLOSER(
        asm = ragtag_out.ragtag_assembly,
        longreads = ratatosk_out.corrected_long
    )

    // 13. QUAST
    quast_ch = gapcloser_out.gapclosed_assembly.combine(ref_fasta)
    quast_out = QUAST(quast_ch)

    // 14. Merqury
    merqury_ch = illumina_pairs
        .map { sample, reads -> tuple(sample, reads[0], reads[1]) }
        .combine(gapcloser_out.gapclosed_assembly)
        .map { sample, r1, r2, asm -> tuple(sample, r1, r2, asm) }

    merqury_out = MERQURY(merqury_ch)

    // 15. MultiQC aggregation
    // Collect all outputs into one channel
    all_reports = Channel.empty()
        .mix(fastqc_out.reports)
        .mix(nanoplot_raw_out.reports)
        .mix(nanoplot_filtered_out.reports)
        .mix(nanoplot_corrected_out.reports)
        .mix(quast_out.reports)
        .mix(merqury_out.reports)

    multiqc_out = MULTIQC(all_reports.collect())
}
