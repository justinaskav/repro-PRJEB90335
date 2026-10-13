#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WF_PREPARE_HOST_INDEX } from './microbiome-workflows/modules/host_alignment/bowtie2_host_wf.nf'
include { CREATE_SILVA_CLASSIFIER as SILVA_V3V4 } from './microbiome-workflows/modules/reference_databases/create_classifier_wf.nf'
include { AMPLICON_WORKFLOW as RRNA_1 ; AMPLICON_WORKFLOW as RRNA_2 ; AMPLICON_WORKFLOW as RRNA_3 } from './microbiome-workflows/workflows/amplicon_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from './microbiome-workflows/modules/phylogeny/phylogeny_wf.nf'
include { MERGE_FEATURE_TABLES ; MERGE_SEQUENCES ; MERGE_TAXONOMIES } from './microbiome-workflows/modules/data_processing/merge_qiime2_wf.nf'

// Main workflow
workflow {
    main:

    // This run produces only what r_prepare.nf reads: the merged feature table,
    // sequences, taxonomy and tree. Everything else the per-batch workflow can do
    // is switched off here; none of it changes those four artefacts.
    def core_only = [
        skip_qc: true,                    // FastQC/MultiQC reports
        skip_summaries: true,             // QIIME2 .qzv summaries, gram-staining evaluation
        skip_contaminant_taxonomy: true,  // classification of the ASVs decontam removed
        skip_picrust2: true,              // functional prediction
        skip_phylogeny: true,             // per-batch trees; the tree is built once after the merge
        skip_diversity: true,             // QIIME2 core metrics; diversity is computed in r_prepare.nf
    ]

    silva_classifier_v3v4_ch = SILVA_V3V4(
        [id: "silva-${params.silva.version}-${params.silva.target}-V3V4"],
        params.silva.version,
        params.silva.target,
        params.silva.filter_labels,
        params.silva.filter_min_lengths,
        params.silva.derep_mode,
        params.primers.forward,
        params.primers.reverse,
        false,
    ).classifier

    // bowtie2 index built once, shared by the three batches
    host_index_ch = WF_PREPARE_HOST_INDEX(params.host_alignment.host_fasta).index_files

    batch1_results = RRNA_1(
        params.batches[0].name,
        "${projectDir}/../${params.batches[0].metadata_path}",
        params.batches[0].input,
        params.batches[0].trimming,
        params.batches[0].dada2,
        [:],
        [:],
        [classifier: silva_classifier_v3v4_ch],
        params.batches[0].other + core_only,
        [:],
        host_index_ch,
    )

    batch2_results = RRNA_2(
        params.batches[1].name,
        "${projectDir}/../${params.batches[1].metadata_path}",
        params.batches[1].input,
        params.batches[1].trimming,
        params.batches[1].dada2,
        [:],
        [:],
        [classifier: silva_classifier_v3v4_ch],
        params.batches[1].other + core_only,
        [:],
        host_index_ch,
    )

    batch3_results = RRNA_3(
        params.batches[2].name,
        "${projectDir}/../${params.batches[2].metadata_path}",
        params.batches[2].input,
        params.batches[2].trimming,
        params.batches[2].dada2,
        [:],
        [:],
        [classifier: silva_classifier_v3v4_ch],
        params.batches[2].other + core_only,
        [:],
        host_index_ch,
    )

    feature_tables_ch = batch1_results.feature_table.mix(batch2_results.feature_table, batch3_results.feature_table).toSortedList()
    rep_seqs_ch = batch1_results.rep_seqs.mix(batch2_results.rep_seqs, batch3_results.rep_seqs).toSortedList()
    taxonomy_ch = batch1_results.taxonomy.mix(batch2_results.taxonomy, batch3_results.taxonomy).toSortedList()

    // Merge results from all batches
    merged_tables = MERGE_FEATURE_TABLES(feature_tables_ch)
    merged_seqs = MERGE_SEQUENCES(rep_seqs_ch)
    merged_taxa = MERGE_TAXONOMIES(taxonomy_ch)

    // Generate phylogeny from merged sequences
    phylogeny = ALIGN_AND_TREE_QIIME2(merged_seqs.merged_seqs)

    publish:
    merged_table = merged_tables.merged_table
    merged_representative_sequences = merged_seqs.merged_seqs
    merged_taxonomy = merged_taxa.merged_taxonomy
    rooted_tree = phylogeny.rooted_tree_qza
}

output {
    rooted_tree {
        path 'merged'
    }
    merged_table {
        path 'merged'
    }
    merged_representative_sequences {
        path 'merged'
    }
    merged_taxonomy {
        path 'merged'
    }
}
