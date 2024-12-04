// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options     = initOptions(params.containsKey("options") ? params.options : [:], 'reffirstfasta')
options.btype = options.btype ?: "comparative"
conda_tools = "bioconda::perl-bioperl=1.7.8"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process REFFIRSTFASTA {
    tag "${meta.id}"
    label "process_single"

    // Not strictly necessary to be samtools - should run on pure bash, but will call non-existing docker image otherwise
    container 'quay.io/biocontainers/perl-bioperl:1.6.924--pl5.22.0_7'
    //containerOptions ''

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("reordered_*")        , emit: fasta
    path "*.{log,err}"                           , optional: true, emit: logs
    path ".command.*"                                            , emit: nf_logs
    path "versions.yml"                                          , emit: versions

    script:
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    fasta_reorder.pl $fasta_name Reference 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version) | grep -oE \'(v[0-9.]+)\')
    END_VERSIONS
    """
}
