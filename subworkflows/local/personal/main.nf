//
// snippy - Rapid variant calling from sequence reads
//
//include { initOptions } from '../../../lib/nf/functions'
// import saveFiles needs for local processes - might be better to make it it's own module
include { initOptions; saveFiles } from '../../../lib/nf/functions'

// snippy options
snippy_opts = initOptions(params.containsKey("options") ? params.options : [:], 'snippy')
snippy_opts.is_module = params.wf == 'snippy' ? true : false
snippy_opts.args = [
    params.bwaopt ? "--bwaopt ${params.bwaopt}" : "",
    params.fbopt ? "--fbopt ${params.fbopt}" : "",
    "--mapqual ${params.mapqual}",
    "--basequal ${params.basequal}",
    "--mincov ${params.mincov}",
    "--minfrac ${params.minfrac}",
    "--minqual ${params.minqual}",
    "--maxsoft ${params.maxsoft}",
    params.snippy_opts ? "${params.snippy_opts}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
snippy_opts.subdir = params.run_name
snippy_opts.logs_use_prefix = true

// snippy-core options
MASK = params.mask ? file(params.mask) : []
core_opts = initOptions(params.containsKey("options") ? params.options : [:], 'snippy-core')
core_opts.is_module = false
core_opts.args = [
    "--maxhap ${params.maxhap}",
    "--mask-char ${params.mask_char}",
    params.snippy_core_opts ? "${params.snippy_core_opts}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
core_opts.publish_to_base = [".full.aln.gz", ".samples.txt"]
core_opts.suffix = "core-snp"

include { SNIPPY_RUN as SNIPPY_RUN_MODULE }  from '../../../modules/nf-core/snippy/run/main' addParams( options: snippy_opts )
include { SNIPPY_CORE as SNIPPY_CORE_MODULE }  from '../../../modules/nf-core/snippy/core/main' addParams( options: core_opts )
include { GUBBINS } from '../gubbins/main' addParams( options: [suffix: 'core-snp', ignore: [".aln.gz"]] )
include { IQTREE } from '../iqtree/main' addParams( options: [suffix: 'core-snp', ignore: [".aln.gz"]] )
include { SNPDISTS } from '../../../modules/nf-core/snpdists/main' addParams( options: [suffix: 'core-snp.distance'] )




process REFFIRSTFASTA {

    // Not strictly necessary to be samtools - should run on pure bash, but will call non-existing docker image otherwise
    container 'quay.io/biocontainers/samtools:1.9--h91753b0_8'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta")                  , emit: fasta
    path "*.{log,err}"                           , optional: true, emit: logs
    path ".command.*"                                            , emit: nf_logs
    path "versions.yml"                                          , emit: versions

    script:
    """
    sequence_count=\$(grep -c "^>" "$fasta")

    awk '/^>/{n++} {if (n == '\$sequence_count') print}' ${fasta} > reffirst-${fasta}.fasta

    awk '/^>/{n++; if (n == '\$sequence_count') exit} {if (n < '\$sequence_count') print}' ${fasta} >> reffirst-${fasta}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version) | sed -n \'s/.*version \\([0-9.()]\\+\\)-release.*/\\1/p\')
    END_VERSIONS
    """
}

workflow PERSONAL {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    reference // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run Snippy per-sample
    SNIPPY_RUN_MODULE(reads, reference)
    ch_versions = ch_versions.mix(SNIPPY_RUN_MODULE.out.versions.first())

    // Identify core SNPs
    SNIPPY_RUN_MODULE.out.vcf.collect{meta, vcf -> vcf}.map{ vcf -> [[id:'snippy-core'], vcf]}.set{ ch_merge_vcf }
    SNIPPY_RUN_MODULE.out.aligned_fa.collect{meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'snippy-core'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }
    SNIPPY_CORE_MODULE(ch_snippy_core, reference, MASK)
    ch_versions = ch_versions.mix(SNIPPY_CORE_MODULE.out.versions.first())

    // Assumes last aln from snippy is always Reference. 
    // Flip last to be first entry in alignment because gubbins use first entry as ref for vcf generation.
    REFFIRSTFASTA(SNIPPY_CORE_MODULE.out.clean_full_aln)

    // Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(REFFIRSTFASTA.out.fasta)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)

        // Recombination filtered Per-sample SNP distances
        SNPDISTS(GUBBINS.out.fasta)
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE_MODULE.out.clean_full_aln)
        }
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    // Adding more stuff here should not stale the resume cache

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
