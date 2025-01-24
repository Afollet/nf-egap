include { samplesheetToList  } from 'plugin/nf-schema'


params.file_path = './assets/schema_input.json'
params.publishDir = './results'

   /**trimmomatic_cmd = ["java", f"-Xmx{RAM_GB}G", "-jar", find_file("trimmomatic-0.39.jar"), 
                            "PE", "-threads", str(CPU_THREADS), "-phred33",
                            ILLUMINA_RAW_F_READS, ILLUMINA_RAW_R_READS, 
                            trimmo_f_pair_path, fwd_unpaired_out,
                            trimmo_r_pair_path, rev_unpaired_out,
                            ILLUMINACLIP, HEADCROP, CROP, SLIDINGWINDOW, MINLEN]
                            **/
process trimmomatic {
  
    publishDir params.publishDir + "/trimmomatic"
    container 'staphb/trimmomatic:0.39'
    label "trimmomatic"

    input:
    val data 

    output:
    tuple val(key), file(trimmedF), file(trimmedR), emit: trimmedF

    script:
    key = data.key
    trimmedF= "${data.illuminaRawFBaseName}_1P.fastq"
    trimmedR= "${data.illuminaRawRBaseName}_2P.fastq"

    """
    trimmomatic \
    PE -phred33 ${data.illuminaRawFR} ${data.illuminaRawFR} \
    ${data.illuminaRawFBaseName}_1P.fastq unpaired_1P.fastq \
    ${data.illuminaRawRBaseName}_2P.fastq unpaired_2P.fastq \
    ${params.ILLUMINACLIP} ${params.HEADCROP} ${params.CROP} ${params.SLIDINGWINDOW} ${params.MINLEN}
    """
}

/**
process filtlong {


    input:
    val data 

    output:
    path("long_read_filtlong.fastq", includeInputs: true)

    script:
    """
    echo '\
    filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    --length_weight 0.5\
    --target_bases  ${target_lr_length} \
    ${lr} > lr_filtlong.fastq '
    """
}**/


workflow {

    mainChannel = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))
    mappedMainChannel = mainChannel.map { 
        // Convert paths to File objects and get base names without extensions
        def ontRawBaseName = new File(it[0]).name.replaceAll(/(\.fastq|\.fq\.gz)$/, '')
        def illuminaRawFBaseName = new File(it[1]).name.replaceAll(/(\.fastq|\.fq\.gz)$/, '')
        def illuminaRawRBaseName = new File(it[2]).name.replaceAll(/(\.fastq|\.fq\.gz)$/, '')
        
        [
            key: it[0],
            ontRawReads: it[1],
            illuminaRawFR: it[2],
            illuminaRawRR: it[3],
            speciesId: it[4],
            organismKingdom: it[5],
            organismKaryote: it[6],
            compleasm1: it[7],
            compleasm2: it[8],
            estSize: it[9],
            refSeq: it[10],
            ontRawBaseName: ontRawBaseName,   // Base name of ONT raw read file
            illuminaRawFBaseName: illuminaRawFBaseName,   // Base name of F-read file
            illuminaRawRBaseName: illuminaRawRBaseName    // Base name of R-read file
        ]
    }
    trimmomaticOutput = trimmomatic(mappedMainChannel)
    //trim long nanpore with porechop?
    //filtlong(mappedMainChannel)
    // View the final channel with appended outputs
}

