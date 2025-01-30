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
  
    publishDir params.publishDir, mode: 'copy'
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

process fastqc {
// Short read quality control
    label "fastqc"
    publishDir "${params.publishDir}/qc", mode: 'copy'
    container 'staphb/fastqc:0.12.1'

    input:
    tuple val(key), path(trimmedF), path(trimmedR)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results
    tuple val(key), env(total_bases), topic: "qc"

    script:

    key = key

    """
    fastqc -q ${trimmedF} ${trimmedR} -o .

    total_bases_fwd=\$(grep "Total Sequences" ${params.publishDir}/qc/shortread/fastqc/${trimmedF}_fastqc.html | awk '{print \$3}')
    total_bases_rev=\$(grep "Total Sequences" ${params.publishDir}/qc/shortread/fastqc/${trimmedF}_fastqc.html | awk '{print \$3}')
    
    total_bases=\$((total_bases_fwd + total_bases_rev))
    """
}

process nanoplot {
// Quality check for nanopore reads and Quality/Length Plots
    label "nanoplot"
    publishDir "${params.outdir}/qc/", mode: 'copy'
    container 'staphb/nanoplot:1.42.0'

    input:
    tuple val(key), file(long_reads)

    output:
    file '*.png'
    file '*.html'
    file '*.txt'
    file("*_NanoStats.txt") 
    env(RAW_ONT_READS), topic: "qc"
    env(RAW_ONT_MEAN_LENGTH), topic: "qc"
    env(RAW_ONT_MEAN_QUAL), topic: "qc"
    env(RAW_ONT_TOTAL_BASES), topic: "qc"

    script:
       /** raw_nanoplot_cmd = [ "NanoPlot", "--fastq", ONT_RAW_READS, "-t", str(CPU_THREADS),
                                "-o", raw_nanoplot_dir, "--plots", "kde", "dot", "--loglength",
                                "--N50", "--title", "Raw ONT Reads: Preliminary Data",
                                "--prefix", "RawONT", "--verbose"] **/
    """

    NanoPlot --fastq  ${long_reads} -t ${params.defaultThreads}  --title ${long_reads} -c darkblue \
    --plots kde dot --loglength --N50 --prefix ${long_reads} --verbose --prefix ${long_reads}_
    nanostats_file="${long_reads}_NanoStats.txt"
    RAW_ONT_READS=\$(grep -oP '(?<=Number of reads: )\\d+' "\$nanostats_file")
    RAW_ONT_MEAN_LENGTH=\$(grep -oP '(?<=Mean read length: )[\\d.]+' "\$nanostats_file")
    RAW_ONT_MEAN_QUAL=\$(grep -oP '(?<=Mean read quality: )[\\d.]+' "\$nanostats_file")
    RAW_ONT_TOTAL_BASES=\$(grep -oP '(?<=Total bases: )\\d+' "\$nanostats_file")
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
            illuminaRawRBaseName: illuminaRawRBaseName, // Base name of R-read file
        ]
    }
    
    //mappedMainChannel.view()

    trimmomaticOutput = trimmomatic(mappedMainChannel)
    fastqc(trimmomaticOutput) 
    
    mappedMainChannel.filter { it.ontRawReads }
        .map { [key: it.key, long_reads: file(it.ontRawReads)] }
        .set { nanoplotInput }

    nanoplotInput.view()

    nanoplotOut = nanoplot(nanoplotInput) 
    //trim long nanpore with porechop?
    //filtlong(mappedMainChannel)
    // View the final channel with appended outputs
}

