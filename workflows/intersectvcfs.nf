/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_CONCAT } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_ISEC } from '../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_NORM } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT } from '../modules/nf-core/bcftools/sort/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_intersectvcfs_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTERSECTVCFS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Concat Strelka SNVs and Indels
    //

    strelka_vcfs = ch_samplesheet
        .map { meta, strelka_snvs, strelka_snvs_tbi, strelka_indels, strelka_indels_tbi, mutect2, mutect2_tbi ->
            [meta + [caller: 'strelka2'], [strelka_snvs, strelka_indels], [strelka_snvs_tbi, strelka_indels_tbi]]
        }

    BCFTOOLS_CONCAT(strelka_vcfs)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    mutect2_vcfs = ch_samplesheet
        .map { meta, strelka_snvs, strelka_snvs_tbi, strelka_indels, strelka_indels_tbi, mutect2, mutect2_tbi ->
            [meta + [caller: 'mutect2'], mutect2, mutect2_tbi]
        }

    BCFTOOLS_NORM(BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi).mix(mutect2_vcfs))
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    BCFTOOLS_SORT(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    sorted_vcfs = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi)

    ch_intersect_in = sorted_vcfs.map { meta, vcf, tbi ->
                                    [[id:meta.id], vcf, tbi]
                                }.groupTuple(size: 2)

    BCFTOOLS_ISEC(ch_intersect_in)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
