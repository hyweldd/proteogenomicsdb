process PYPGATK_COSMICDB {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_single_thread'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/ggrimes/pypgatk_ar:latest' :
        'docker://ghcr.io/ggrimes/pypgatk_ar:latest' }"

    input:
    tuple val(meta), path(cosmic_genes)
    tuple val(meta), path(cosmic_mutations)
    path cosmic_config

    output:
    path '*.fa'           , emit: database
    path  "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = "${prefix}.fa"

    """
    pypgatkc.py cosmic-to-proteindb \\
        --config_file ${cosmic_config} \\
        --input_mutation ${cosmic_mutations} \\
        --input_genes ${cosmic_genes} \\
        --output_db ${name} \\
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgatk: \$(pypgatk --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """

    touch ${prefix}.fa 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgatk: \$(pypgatk --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}