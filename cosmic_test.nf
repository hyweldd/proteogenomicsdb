
include { BCFTOOLS_SORT } from '/exports/eddie/scratch/s2215490/proteogenomicsdb/modules/nf-core/bcftools/sort/main.nf'
 
        params.vcf = '/exports/eddie/scratch/s2215490/paper_results/proteogenomicsdb/bcftools/after_merge.vcf'

workflow {

    vcf = Channel.fromPath(params.vcf)

    BCFTOOLS_SORT (
        vcf.map { [ [id: 'vcf' ], it ] }
    )
}