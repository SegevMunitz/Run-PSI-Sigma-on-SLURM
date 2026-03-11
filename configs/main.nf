nextflow.enable.dsl=2

process PREP {
    tag 'prep'
    output:
    path 'sample_tasks.tsv', emit: sample_tasks
    path "${params.outdir}/.done/prep.done", emit: done

    script:
    """
    set -euo pipefail
    ${params.python} ${params.code_dir}/step1_prep_run.py

    OUTDIR="${params.outdir}"
    SAMPLES_FILE="\${OUTDIR}/samples.txt"
    awk 'NF{print NR"\t"\$0}' "\${SAMPLES_FILE}" > sample_tasks.tsv
    mkdir -p "\${OUTDIR}/.done"
    touch "\${OUTDIR}/.done/prep.done"
    """
}

process STAR_ALIGN {
    tag { "star:${sample_id}" }
    input: tuple val(task_id), val(sample_id)
    output: path "star.${task_id}.done", emit: done

    script:
    """
    set -euo pipefail
    ${params.python} ${params.code_dir}/step2_star_align_core.py --task-id ${task_id}
    touch star.${task_id}.done
    """
}

process MAKE_GROUPS {
    tag 'groups'
    input: val(star_done_markers)
    output: path 'groups.done', emit: done

    script:
    """
    ${params.python} ${params.code_dir}/step3_make_group_files.py
    touch groups.done
    """
}

process FEATURE_COUNTS {
    tag 'featureCounts'
    input: val(star_done_markers)
    output: path 'featurecounts.done', emit: done

    script:
    def pe_flag = params.featurecounts_paired_end ? '--paired-end' : ''
    """
    ${params.python} ${params.code_dir}/step4_featurecounts.py ${pe_flag}
    touch featurecounts.done
    """
}

process PSI_SIGMA {
    tag 'psi_sigma'
    input: path(groups_done)
    path(featurecounts_done)
    output: path 'psisigma.done', emit: done

    script:
    """
    ${params.python} ${params.code_dir}/step5_run_psi_sigma.py
    touch psisigma.done
    """
}

process SALMON_QUANT {
    tag { "salmon:${sample_id}" }
    input: tuple val(task_id), val(sample_id)
    output: path "salmon.${task_id}.done", emit: done

    script:
    """
    set -euo pipefail
    ${params.python} ${params.code_dir}/step6_salmon_quant_core.py --task-id ${task_id}
    touch salmon.${task_id}.done
    """
}

process SALMON_FILTER {
    tag 'salmon_filter'
    input: path(psisigma_done)
    val(salmon_done_markers)
    output: path 'salmon_filter.done', emit: done

    script:
    """
    ${params.python} ${params.code_dir}/step7_psi_sigma_filtered_by_salmon.py
    touch salmon_filter.done
    """
}

workflow {
    prep = PREP()

    sample_tasks_ch = prep.out.sample_tasks
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple((row[0] as Integer), row[1] as String) }

    star = STAR_ALIGN(sample_tasks_ch)
    star_all = star.out.done.collect()

    groups = MAKE_GROUPS(star_all)
    feature_counts = FEATURE_COUNTS(star_all)

    psi_sigma = PSI_SIGMA(groups.out.done, feature_counts.out.done)

    // Prefix star_files with _ to suppress "Parameter not used" warning
    salmon_input = sample_tasks_ch.combine(star_all).map{ task, _star_files -> task }
    salmon = SALMON_QUANT(salmon_input)
    
    salmon_all = salmon.out.done.collect()

    filter = SALMON_FILTER(psi_sigma.out.done, salmon_all)
    filter.out.done.view { result -> "Pipeline completed: ${result}" }
}
