
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//modules for downloading databases from ENSEMBL
include { PYPGATK_ENSEMBL_DOWNLOAD } from '../../../modules/local/pypgatk/ensembl_downloader/main.nf'
include { CAT_CAT                  } from '../../../modules/nf-core/cat/cat/main.nf'

//modules for optional database generation
include { PYPGATKDNA as PYPGATK_NCRNA       } from '../../../modules/local/pypgatk/dnaseq_to_proteindb/main.nf'
include { PYPGATKDNA as PYPGATK_PSEUDOGENES } from '../../../modules/local/pypgatk/dnaseq_to_proteindb/main.nf'
include { PYPGATKDNA as PYPGATK_ALRORFS     } from '../../../modules/local/pypgatk/dnaseq_to_proteindb/main.nf'

//modules for the generation of the main ENSEMBL database
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main.nf'
include { CRABZ_COMPRESS } from '../../../modules/nf-core/crabz/compress/main.nf'
include { PYPGATK_VCF    } from '../../../modules/local/pypgatk/vcf_to_proteindb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN ENSEMBLDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ENSEMBLDB {

take:

    //inputs from the main workflow 
    ensembl_downloader_config //channel: /path/to/ensembl downloader config
    species_name              //string: ENSEMBL species ID
    ensembl_config            //channel: /path/to/ensembl config
    altorfs_config            //channel: /path/to/altorfs config
    pseudogenes_config        //channel: /path/to/pseudogenes config
    ncrna_config              //channel: /path/to/ncrna config
    skip_proteome   
    skip_ncrna      
    skip_pseudogenes
    skip_altorfs   
    skip_ensembl_vcf

main:

    //creates empty channels for tool versions and the peptide database
    versions_ch     = Channel.empty()
    mixed_databases = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF ENSEMBL DOWNLOAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //creating empty channels used in the ENSEMBL download pathway
    cdna_mixed = Channel.empty()
    total_cdna = Channel.empty()

    //PYPGATK_ENSEMBL uses the ensembl species ID to downloads files from ENSEMBL
    PYPGATK_ENSEMBL_DOWNLOAD (
        ensembl_downloader_config,
        species_name,
        skip_ensembl_vcf
    )
    versions_ch = versions_ch.mix(PYPGATK_ENSEMBL_DOWNLOAD.out.versions).collect()
    cdna_mixed = PYPGATK_ENSEMBL_DOWNLOAD.out.cdna.mix(PYPGATK_ENSEMBL_DOWNLOAD.out.ncrna).collect()
        .map { [ [id: 'Total_cDNA' ], it ] }

    CAT_CAT (
        cdna_mixed
    )
    versions_ch = versions_ch.mix(CAT_CAT.out.versions_cat).collect()
    total_cdna = CAT_CAT.out.file_out.collect()
        .map { meta, it ->
            return [it] }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF ENSEMBL DOWNLOAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF OPTIONAL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//conditional execution based on if skip_proteome is true or false
if (!skip_proteome) {

    //creates an empty channel that will then be populated with the protein files downloaded from ENSEMBL
    reference_proteome = Channel.empty()
    reference_proteome = PYPGATK_ENSEMBL_DOWNLOAD.out.protein.collect()

    //adds the protein files downloaded from ENSEMBL to the mixed database channel
    mixed_databases = mixed_databases.mix(reference_proteome).collect()

}

else {
        //bypass the subworkflow
        log.info "proteome skipped."
    }

//conditional execution based on if skip_ncrna is true or false
if (!skip_ncrna) {

    //PYPGATK_NCRNA takes the total_cdna and the ncrna_config to generate a peptide database
    PYPGATK_NCRNA (
        total_cdna.map { [ [id: 'ncRNA'], it ] },
        ncrna_config
    )
    versions_ch = versions_ch.mix(PYPGATK_NCRNA.out.versions).collect()

    //adds the peptide database generated from ncrna data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_NCRNA.out.database).collect()

}

else {
        //bypass the subworkflow
        log.info "ncrna database skipped."
    }

//conditional execution based on if skip_pseudogenes is true or false
if (!skip_pseudogenes) {

    //PYPGATK_PSEUDOGENES takes the total_cdna and the pseudogenes_config to generate a peptide database
    PYPGATK_PSEUDOGENES (
        total_cdna.map { [ [id: 'pseudogenes'], it ] },
        pseudogenes_config
    )
    versions_ch = versions_ch.mix(PYPGATK_PSEUDOGENES.out.versions).collect()

    //adds the peptide database generated from pseudogenes data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_PSEUDOGENES.out.database).collect()

}

else {
        //bypass the subworkflow
        log.info "pseudogenes database skipped."
    }

//conditional execution based on if skip_altorfs is true or false
if (!skip_altorfs) {

    //creates an empty channel that will then be populated with the cdna files downloaded from ENSEMBL
    cdna_database = Channel.empty()
    cdna_database = PYPGATK_ENSEMBL_DOWNLOAD.out.cdna.collect()
        .map { [ [id: 'Altorfs_database'], it ] }

    //PYPGATK_ALTORFS takes the cdna files and the altorfs_config to generate a peptide database
    PYPGATK_ALRORFS (
        cdna_database,
        altorfs_config
    )
    versions_ch = versions_ch.mix(PYPGATK_ALRORFS.out.versions).collect()

    //adds the peptide database generated from altorfs data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_ALRORFS.out.database).collect()
}

else {
        //bypass the subworkflow
        log.info "altorfs database skipped."
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF OPTIONAL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF ENSEMBL VCF DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//conditional execution based on if skip_ensembl_vcf is true or false
if (!skip_ensembl_vcf) {

    //creating empty channels used in the ENSEMBL vcf database generation pathway
    ensembl_vcf  = Channel.empty()
    reference_ch = Channel.empty()
    fasta_fai_ch = Channel.empty()

    compressed_vcf = Channel.empty()

    total_vcf    = Channel.empty()
    vcf_index_ch = Channel.empty()
    vcf_merge_ch = Channel.empty()
    merged_vcf   = Channel.empty()
    ensembl_gtf  = Channel.empty()

    //populates ensembl_vcf with the vcf file downloaded from ENSEMBL
    ensembl_vcf = PYPGATK_ENSEMBL_DOWNLOAD.out.vcf.flatten()
        .map { [ [id: 'concatenated_vcf' ], it ] }
    
    //populates ensembl_gtf with the gtf file downloaded from ENSEMBL 
    ensembl_gtf = PYPGATK_ENSEMBL_DOWNLOAD.out.gtf.collect()
        .map { [ [id: 'gtf' ], it ] } 

    reference_ch = PYPGATK_ENSEMBL_DOWNLOAD.out.fasta.collect()

    SAMTOOLS_FAIDX (
        reference_ch.map { [ [ id: 'reference' ], it ] },
        [ [], [] ],
        false
    )
    versions_ch = versions_ch.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    fasta_fai_ch = SAMTOOLS_FAIDX.out.fai.collect()

    CRABZ_COMPRESS (
        ensembl_vcf
    )
    versions_ch = versions_ch.mix(CRABZ_COMPRESS.out.versions)
    compressed_vcf = CRABZ_COMPRESS.out.archive

    BCFTOOLS_INDEX (
        compressed_vcf
    )
    versions_ch = versions_ch.mix(BCFTOOLS_INDEX.out.versions_bcftools)
    vcf_index_ch = BCFTOOLS_INDEX.out.tbi
        .map { meta, index ->
            return index }
        .collect()
        .map { index ->
            def meta2 = [ id: 'key' ]
            return [ meta2, index ] }
  
    merged_vcf = compressed_vcf
            .map { meta, vcf ->
                return vcf }
            .collect()
            .map { vcf ->
                def meta2 = [ id: 'key' ]
                return [ meta2, vcf ] }

    vcf_merge_ch = merged_vcf.join(vcf_index_ch).collect()

    BCFTOOLS_MERGE (
        vcf_merge_ch,
        reference_ch.map { [ [ id: 'reference' ], it ] },
        fasta_fai_ch,
        [ [], [] ]
    )
    versions_ch = versions_ch.mix(BCFTOOLS_MERGE.out.versions_bcftools)
    vcf_file_ch = BCFTOOLS_MERGE.out.vcf


    //PYPGATK takes the concatenated vcf file, the gtf file, and the 
    //concatenated cDNA along with the ensembl_config to generate a peptide database
    PYPGATK_VCF (
        vcf_file_ch,
        ensembl_gtf,
        total_cdna.map { [ [id: 'ensembl_vcf'], it ] },
        ensembl_config
    )
    versions_ch = versions_ch.mix(PYPGATK_VCF.out.versions).collect()

    //adds the peptide database generated from the vcf file to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_VCF.out.database).collect()
}

else {
        //bypass the subworkflow
        log.info "ensembl vcf database skipped."
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF MAIN ENSEMBL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


emit:
    
    //emits to the main workflow
    mixed_databases //channel: [ val(meta), [ database ] ]
    versions_ch     //channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
