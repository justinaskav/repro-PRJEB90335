#!/usr/bin/env nextflow

// Parallel equivalent of R/make_all.R: one task per script in R/figures/.
// Each task gets the repo layout the scripts expect (R/, shared/R/, _output/r/tses/final/)
// and its outputs are copied back to the same places under the launch directory.

process RENDER_FIGURE {
    tag "${script.baseName}"
    conda params.r_conda_env
    cache 'deep'   // shared/R is staged as a directory; hash its contents, not its mtime
    cpus 1
    publishDir ".", mode: 'copy', overwrite: true

    input:
    path script, stageAs: 'R/figures/*'
    path setup, stageAs: 'R/setup.R'
    path helpers, stageAs: 'R/helpers.R'
    path shared_r, stageAs: 'shared/R'
    path tses, stageAs: '_output/r/tses/final/*'

    output:
    path "figures/*.png", emit: figures
    path "_output/r/maaslin3", optional: true, emit: maaslin3

    script:
    """
    Rscript ${script}
    """
}

workflow {
    scripts_ch = Channel.fromPath("R/figures/*.R")
    tses_ch = Channel.fromPath("_output/r/tses/final/tse_*_with_beta.rds").collect()

    RENDER_FIGURE(
        scripts_ch,
        file("R/setup.R"),
        file("R/helpers.R"),
        file("shared/R", type: 'dir'),
        tses_ch,
    )
}
