#!/usr/bin/env nextflow

params.output_dir = "_output"

include { R_ALPHA_DIVERSITY } from './microbiome-workflows/subworkflows/r_alpha_diversity_wf.nf'
include { R_BETA_DIVERSITY } from './microbiome-workflows/subworkflows/r_beta_diversity_wf.nf'

process PREPARE_TSES {
    tag "prepare_tses"
    conda params.r_conda_env

    input:
    path table                 // final_table.qza
    path taxonomy             // final_taxonomy.qza
    path rooted_tree          // rooted_tree.qza
    path metadata             // metadata.tsv
    path exclude_samples      // excluded_samples.txt
    path seqs                 // final_seqs.qza
    path prepare_script       // prepare_data.R

    output:
    path "tse_full.rds", emit: tse_full
    path "tse_feces.rds", emit: tse_feces
    path "tse_gut.rds", emit: tse_gut

    script:
    """
    Rscript ${prepare_script} \
        --features ${table} \
        --taxonomy ${taxonomy} \
        --tree ${rooted_tree} \
        --metadata ${metadata} \
        --exclude_samples ${exclude_samples} \
        --seqs ${seqs} \
        --output_dir .
    """
}

workflow {
    main:

    table_ch = Channel.fromPath("${params.output_dir}/merged/final_table.qza")
    taxonomy_ch = Channel.fromPath("${params.output_dir}/merged/final_taxonomy.qza")
    tree_ch = Channel.fromPath("${params.output_dir}/merged/rooted_tree.qza")
    seqs_ch = Channel.fromPath("${params.output_dir}/merged/final_seqs.qza")

    metadata_ch = Channel.fromPath("data/metadata/metadata.tsv")
    exclude_samples_ch = Channel.fromPath("data/metadata/excluded_samples.txt")
    prepare_script_ch = Channel.fromPath("scripts/prepare_data.R")

    prepare_out = PREPARE_TSES(
        table_ch,
        taxonomy_ch,
        tree_ch,
        metadata_ch,
        exclude_samples_ch,
        seqs_ch,
        prepare_script_ch,
    )

    // Build channel of [tse_type, tse_file] for diversity calculations
    tse_data = Channel.empty()
        .mix(
            prepare_out.tse_full.map { tse_file -> ["full", tse_file] },
            prepare_out.tse_gut.map { tse_file -> ["gut", tse_file] },
            prepare_out.tse_feces.map { tse_file -> ["feces", tse_file] }
        )

    // Calculate alpha diversity
    R_ALPHA_DIVERSITY(tse_data, [:])

    // Calculate beta diversity (chain from alpha output)
    alpha_tse_data = R_ALPHA_DIVERSITY.out.tse_with_alpha.map { tse_file ->
        def type = tse_file.name.replace("tse_", "").replace("_with_alpha.rds", "")
        [type, tse_file]
    }

    R_BETA_DIVERSITY(alpha_tse_data, [:])

    publish:
    // TSE objects with diversity calculations
    tse_full_final = R_BETA_DIVERSITY.out.tse_with_beta.filter { it.name.contains("full") }
    tse_feces_final = R_BETA_DIVERSITY.out.tse_with_beta.filter { it.name.contains("feces") }
    tse_gut_final = R_BETA_DIVERSITY.out.tse_with_beta.filter { it.name.contains("gut") }
}

output {
    tse_full_final {
        path "r/tses/final"
    }
    tse_feces_final {
        path "r/tses/final"
    }
    tse_gut_final {
        path "r/tses/final"
    }
}
