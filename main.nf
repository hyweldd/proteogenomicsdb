#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteogenomicsdb
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteogenomicsdb
    Website: https://nf-co.re/proteogenomicsdb
    Slack  : https://nfcore.slack.com/channels/proteogenomicsdb
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DATABASE_GENERATION     } from './workflows/database_generation.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.reference = getGenomeAttribute('fasta')
params.annotation = getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PROTEOGENOMICSDB {

take:
    bam_files
    bam_index

main:
    //
    // WORKFLOW: Run pipeline
    //
    DATABASE_GENERATION (
        //RNASEQDB paramiters
        bam_files,                     //channel: [val(meta), path(samplesheet)]
        bam_index,                     //channel: [val(meta), path(samplesheet)]
        params.reference,              //path/to/genome/reference
        params.annotation,             //path/to/genome/annotation
        params.transcripts,            //path/to/transcripts
        params.custom_config,          //path/to/custom_config
        params.dna_config,             //path/to/dna_config              
        params.faidx_get_genome_sizes, //true or false statement 

        //ENSEMBLDB paramaters
        params.ensembl_downloader_config, //path to the ensembl_downloader_config
        params.species_id,                //species id to download from ensembl
        params.ensembl_config,            //path to the ensembl_config
        params.altorfs_config,            //path to the altorfs_config
        params.pseudogenes_config,        //path to the pseudogenes_config
        params.ncrna_config,              //path to the ncrna_config
        
        //COSMICDB paramiters
        params.cosmic_config,        //path/to/cosmic_config
        params.username,             //string: cosmic username           
        params.password,             //string: cosmic password
        params.cosmic_genes_url,     //string: cosmic genes url
        params.cosmic_mutations_url, //string: cosmic mutations url

        //GENCODEDB paramaters
        params.gencode_transcripts_url, //string: gencode transcripts url
        params.gencode_annotations_url, //string: gencode annotations url
        params.gencode_reference_url,   //string: gencode reference url
        params.gnomad_url,              //string: gnomad url
        params.gencode_config,          //path/to/gnomad_config
        
        //CBIOPORTALDB paramaters
        params.cbioportal_url, //string: cbioportal url
        params.grch38_url,     //string: grch38 url
        params.cbio_study,     //string: name of cbioportal study
        params.cbio_config,    //path/to/cbioportal_config
        
        //MERGEDB paramaters
        params.clean_config,        //path/to/clean_config
        params.decoy_config,        //path/to/decoy_config
        params.additional_database, //path/to/additional_database
        
        //skip options 
        params.skip_rnaseqdb,           //boolean: option to skip the rnaseqdb workflow
        params.skip_dnaseq,             //boolean: option to skip RNA-seq canonical database
        params.skip_vcf,                //boolean: option to skip RNA-seq variant database
        params.skip_ensembldb,          //boolean: option to skip the ensembldb workflow
        params.skip_proteome,           //boolean: option to skip ensembl protein sequences
        params.skip_ncrna,              //boolean: option to skip ensembl ncRNA database
        params.skip_pseudogenes,        //boolean: option to skip ensembl pseudogenes database
        params.skip_altorfs,            //boolean: option to skip ensembl altORFs database
        params.skip_ensembl_vcf,        //boolean: option to skip ensembl variant database
        params.skip_cosmicdb,           //boolean: option to skip the cosmicdb workflow
        params.skip_gencodedb,          //boolean: option to skip the gencodedb workflow
        params.skip_cbioportaldb,       //boolean: option to skip the cBioPortaldb workflow
        params.skip_decoy,              //boolean: option to skip decoy
        params.skip_additional_database //boolean: option to skip additional database

    )

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    NFCORE_PROTEOGENOMICSDB (
        PIPELINE_INITIALISATION.out.bam_samplesheet,
        PIPELINE_INITIALISATION.out.bai_samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url     
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
