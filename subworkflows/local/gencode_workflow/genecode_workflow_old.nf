
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WGET   as WGET_TRANSCRIPTS   } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_TRANSCRIPTS } from '../../../modules/nf-core/gunzip/main.nf'
include { WGET as WGET_ANNOTATION      } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_ANNOTATION  } from '../../../modules/nf-core/gunzip/main.nf'

include { WGET as WGET_REFERENCE      } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_REFERENCE  } from '../../../modules/nf-core/gunzip/main.nf'
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main.nf'

include { GSUTIL      } from '../../../modules/local/gsutil/main.nf'
include { CAT_CAT     } from '../../../modules/nf-core/cat/cat/main.nf'
include { PYPGATK_VCF } from '../../../modules/local/pypgatk/vcf_to_proteindb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN GNOMADDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENECODEDB {

take:

    //inputs from the main workflow
    genecode_transcripts_url //string: genecode transcripts url
    genecode_annotation_url  //string: genecode annotations url
    genecode_reference_url
    gnomad_url               //string: gnomad vcf url
    genecode_config          //channel: /path/to/genecode config

main:

    //creates empty channels for tool versions and the peptide database
    versions_ch       = Channel.empty()
    genecode_database = Channel.empty()

    //creates empty channels used in the GENECODEDB workflow
    genecode_transcripts_ch = Channel.empty()
    genecode_annotation_ch  = Channel.empty()
    genecode_reference_ch   = Channel.empty()

    genecode_annotation_gunzipped  = Channel.empty()
    genecode_transcripts_gunzipped = Channel.empty()
    genecode_reference_gunzipped   = Channel.empty()

    fasta_fai_ch = Channel.empty()

    vcf_compressed_ch = Channel.empty()
    vcf_file_ch       = Channel.empty()

    merged_vcf        = Channel.empty()



    //WGET_TRANSCRIPTS downloads the transcript files from genecode 
    WGET_TRANSCRIPTS (
        genecode_transcripts_url.map { [ [id: 'transcripts.fa' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_TRANSCRIPTS.out.versions).collect()
    genecode_transcripts_ch = WGET_TRANSCRIPTS.out.outfile.collect()

    //GUNZIP_TRANSCRIPTS unzipps the transcript files downloaded from genecode
    GUNZIP_TRANSCRIPTS (
        genecode_transcripts_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_TRANSCRIPTS.out.versions_gunzip).collect()
    genecode_transcripts_gunzipped = GUNZIP_TRANSCRIPTS.out.gunzip.collect()

    //WGET_ANNOTATION downloads the annotation files from genecode 
    WGET_ANNOTATION (
        genecode_annotation_url.map { [ [id: 'annotations.gtf' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_ANNOTATION.out.versions).collect()
    genecode_annotation_ch = WGET_ANNOTATION.out.outfile.collect()

    //GUNZIP_ANNOTATION unzipps the annotation files downloaded from genecode
    GUNZIP_ANNOTATION (
        genecode_annotation_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_ANNOTATION.out.versions_gunzip).collect()
    genecode_annotation_gunzipped = GUNZIP_ANNOTATION.out.gunzip.collect()

    WGET_REFERENCE (
        genecode_reference_url.map { [ [id: 'reference.fa' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_REFERENCE.out.versions).collect()
    genecode_reference_ch = WGET_REFERENCE.out.outfile.collect()

    //GUNZIP_ANNOTATION unzipps the annotation files downloaded from genecode
    GUNZIP_REFERENCE (
        genecode_reference_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_REFERENCE.out.versions_gunzip).collect()
    genecode_reference_gunzipped = GUNZIP_REFERENCE.out.gunzip.collect()
    
    //GSUTIL downloads the vcf file from gnomad 
    GSUTIL (
        gnomad_url.map { [ [id: 'gnomad.vcf' ], it ] }
    )
    versions_ch = versions_ch.mix(GSUTIL.out.versions).collect()
    vcf_compressed = GSUTIL.out.vcf.collect()

    CAT_CAT (
        vcf_compressed
        )
    versions_ch = versions_ch.mix(CAT_CAT.out.versions_cat)
    vcf_file_ch = CAT_CAT.out.file_out.collect()

    //PYPGATK takes the unzipped vcf, annotation, and transcript along with the 
    //genecode config to generate a peptide database
    PYPGATK_VCF (
        vcf_file_ch,
        genecode_annotation_gunzipped,
        genecode_transcripts_gunzipped,
        genecode_config  
    )
    versions_ch = versions_ch.mix(PYPGATK_VCF.out.versions).collect()
    
    //adds the peptide database generated from the PYPGATK_VCF to genecode_database
    genecode_database = PYPGATK_VCF.out.database.collect()

emit:
    
    // emits to the main workflow
    genecode_database //channel: [ val(meta), [ database ] ]
    versions_ch       //channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
