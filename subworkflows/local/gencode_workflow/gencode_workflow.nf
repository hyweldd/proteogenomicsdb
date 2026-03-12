
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//download GENCODE transcripts
include { WGET   as WGET_TRANSCRIPTS   } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_TRANSCRIPTS } from '../../../modules/nf-core/gunzip/main.nf'

//download GENCODE annotation
include { WGET as WGET_ANNOTATION      } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_ANNOTATION  } from '../../../modules/nf-core/gunzip/main.nf'

//download GENCODE genome reference
include { WGET as WGET_REFERENCE      } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_REFERENCE  } from '../../../modules/nf-core/gunzip/main.nf'

//generate a index on the GENCODE genome referencce
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main.nf'

//download and merge the GnomAD VCF files
include { GSUTIL                       } from '../../../modules/local/gsutil/main.nf'
include { BCFTOOLS_INDEX               } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_MERGE               } from '../../../modules/nf-core/bcftools/merge/main.nf'

//create the Gencode variant database
include { PYPGATK_VCF                  } from '../../../modules/local/pypgatk/vcf_to_proteindb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN GNOMADDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENCODEDB {

take:

    //inputs from the main workflow
    gencode_transcripts_url //string: gencode transcripts url
    gencode_annotation_url  //string: gencode annotations url
    gencode_reference_url   //string: gencode reference url 
    gnomad_url               //string: gnomad vcf url
    gencode_config          //channel: /path/to/gencode config

main:

    //creates empty channels for tool versions and the peptide database
    versions_ch       = Channel.empty()
    gencode_database = Channel.empty()

    //creates empty channels for storing the GENCODE downloads
    gencode_transcripts_ch = Channel.empty()
    gencode_annotation_ch  = Channel.empty()
    gencode_reference_ch   = Channel.empty()

    //creates empty channels for storing the GENCODE downloads (once extracted)
    gencode_annotation_gunzipped  = Channel.empty()
    gencode_transcripts_gunzipped = Channel.empty()
    gencode_reference_gunzipped   = Channel.empty()

    //empty channel for holding the genome reference index
    fasta_fai_ch = Channel.empty()

    //empty channels for holding the VCF file 
    vcf_compressed_ch = Channel.empty()
    vcf_file_ch       = Channel.empty()
    merged_vcf        = Channel.empty()

    //empty channel for holding the VCF index
    vcf_index_ch = Channel.empty()

    //WGET_TRANSCRIPTS downloads the cDNA transcript files from gencode 
    WGET_TRANSCRIPTS (
        gencode_transcripts_url.map { [ [id: 'transcripts.fa' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_TRANSCRIPTS.out.versions).collect()
    gencode_transcripts_ch = WGET_TRANSCRIPTS.out.outfile.collect()

    //GUNZIP_TRANSCRIPTS unzips the transcript files downloaded from gencode
    GUNZIP_TRANSCRIPTS (
        gencode_transcripts_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_TRANSCRIPTS.out.versions_gunzip).collect()
    gencode_transcripts_gunzipped = GUNZIP_TRANSCRIPTS.out.gunzip.collect()

    //WGET_ANNOTATION downloads the annotation files from gencode 
    WGET_ANNOTATION (
        gencode_annotation_url.map { [ [id: 'annotations.gtf' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_ANNOTATION.out.versions).collect()
    gencode_annotation_ch = WGET_ANNOTATION.out.outfile.collect()

    //GUNZIP_ANNOTATION unzipps the annotation files downloaded from gencode
    GUNZIP_ANNOTATION (
        gencode_annotation_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_ANNOTATION.out.versions_gunzip).collect()
    gencode_annotation_gunzipped = GUNZIP_ANNOTATION.out.gunzip.collect()

    WGET_REFERENCE (
        gencode_reference_url.map { [ [id: 'reference.fa' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_REFERENCE.out.versions).collect()
    gencode_reference_ch = WGET_REFERENCE.out.outfile.collect()

    //GUNZIP_ANNOTATION unzipps the annotation files downloaded from gencode
    GUNZIP_REFERENCE (
        gencode_reference_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_REFERENCE.out.versions_gunzip).collect()
    gencode_reference_gunzipped = GUNZIP_REFERENCE.out.gunzip.collect()

    //generates a fasta index from the gencode genome reference
    SAMTOOLS_FAIDX (
        gencode_reference_gunzipped,
        [ [], [] ],
        false  
    )
    versions_ch = versions_ch.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    fasta_fai_ch = SAMTOOLS_FAIDX.out.fai.collect()
    
    //GSUTIL downloads VCF files from gnomad 
    GSUTIL (
        gnomad_url.map { [ [id: 'gnomad.vcf' ], it ] }
    )
    versions_ch = versions_ch.mix(GSUTIL.out.versions).collect()
    vcf_compressed = GSUTIL.out.vcf
        .map { meta, vcf -> 
            return vcf }
        .flatten()

    //generates a VCF index from the gnomad downloads
    BCFTOOLS_INDEX (
        vcf_compressed.map { [ [ id: 'vcf' ], it ] }
    )
    versions_ch = versions_ch.mix(BCFTOOLS_INDEX.out.versions_bcftools)
    vcf_index_ch = BCFTOOLS_INDEX.out.tbi
        .map { meta, index ->
            return index }
        .collect()
        .map { index ->
            def meta2 = [ id: 'key' ]
            return [ meta2, index ] }
    
        merged_vcf = vcf_compressed.collect()
            .map { vcf ->
                def meta = [ id: 'key' ]
                return [ meta, vcf ] }

    vcf_merge_ch = merged_vcf.join(vcf_index_ch).collect()

    //merges the VCF files togther using the VCF and fasta indexes
    BCFTOOLS_MERGE (
        vcf_merge_ch,
        gencode_reference_gunzipped,
        fasta_fai_ch,
        [ [], [] ]
    )
    versions_ch = versions_ch.mix(BCFTOOLS_MERGE.out.versions_bcftools)
    vcf_file_ch = BCFTOOLS_MERGE.out.vcf

    //PYPGATK takes the merged VCF, annotation, and transcripts along with the 
    //gencode config to generate a variant peptide database
    PYPGATK_VCF (
        vcf_file_ch,
        gencode_annotation_gunzipped,
        gencode_transcripts_gunzipped,
        gencode_config  
    )
    versions_ch = versions_ch.mix(PYPGATK_VCF.out.versions).collect()
    gencode_database = PYPGATK_VCF.out.database.collect()

emit:
    
    // emits to the main workflow
    gencode_database //channel: [ val(meta), [ database ] ]
    versions_ch      //channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
