//
// Nanopore de novo assembly with Dragonflye, host removal with Kraken2, and QC with QUAST.
// Assembler is selected via --assemblers; currently only 'dragonflye' is supported for Nanopore.
// Guard with --skip_assembly / --assemblers as done for the Illumina assemblers.
//

include { KRAKEN2_KRAKEN2 } from '../../../modules/nf-core/kraken2/kraken2/main'
include { DRAGONFLYE      } from '../../../modules/nf-core/dragonflye/main'
include { QUAST           } from '../../../modules/nf-core/quast/main'

workflow ASSEMBLY_DRAGONFLYE {

    take:
    ch_fastq      // channel: [ val(meta), path(fastq) ]  -- single-end Nanopore reads
    ch_kraken2_db // channel: path(db)
    ch_fasta      // channel: path(genome.fasta)
    ch_gff        // channel: [ val(meta), path(genome.gff) ] or [ [:], [] ]

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Kraken2 – remove host reads, keep unclassified (non-host) reads for assembly
    //
    ch_assembly_fastq = ch_fastq
    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2 (
            ch_fastq,
            ch_kraken2_db,
            true,  // save_output_fastqs – save classified & unclassified FASTQs
            false  // save_reads_assignment
        )
        ch_multiqc_files  = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect { it[1] }.ifEmpty([]))
        ch_versions       = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

        // Use unclassified (non-host) reads for assembly
        ch_assembly_fastq = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
    }

    //
    // MODULE: Dragonflye – de novo assembly from Nanopore long reads
    //
    DRAGONFLYE (
        ch_assembly_fastq.map { meta, fastq -> [ meta, [], fastq ] } // no short-read polishing
    )
    ch_versions = ch_versions.mix(DRAGONFLYE.out.versions.first())

    //
    // MODULE: QUAST – assembly quality assessment across all samples
    //
    ch_quast_results = channel.empty()
    if (!params.skip_assembly_quast) {
        DRAGONFLYE.out.contigs
            .collect { meta, fa -> fa }
            .map { contigs -> tuple([id: "dragonflye_quast"], contigs) }
            .set { ch_to_quast }

        QUAST (
            ch_to_quast,
            ch_fasta.map { fasta -> [ [:], fasta ] },
            ch_gff
        )
        ch_quast_results = QUAST.out.results
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect { it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(QUAST.out.versions)
    }

    emit:
    assembly_fastq   = ch_assembly_fastq  // channel: [ val(meta), path(fastq) ]
    contigs          = DRAGONFLYE.out.contigs   // channel: [ val(meta), path(*.fa) ]
    quast_results    = ch_quast_results         // channel: [ val(meta), path(dir) ]
    multiqc_files    = ch_multiqc_files         // channel: mixed multiqc inputs
    versions         = ch_versions              // channel: [ versions.yml ]
}
