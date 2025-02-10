include { samplesheetToList  } from 'plugin/nf-schema'


params.file_path = './assets/schema_input.json'
params.publishDir = './results'

process trimmomatic {
  
    publishDir "$baseDir/results", mode: 'copy'
    container 'staphb/trimmomatic:0.39'
    label "trimmomatic"

    input:
    val data
    file BARCODE_FILE
    tuple val(ILLUMINACLIP), val(HEADCROP), val(CROP), val(SLIDINGWINDOW), val(MINLEN)

    output:
    tuple val(key), file("${illuminaRawFBaseName}_1P.trim.fastq.gz"), file("${illuminaRawRBaseName}_2P.trim.fastq.gz"), emit: trimmedF

    script:
    key = data.key
    illuminaRawFBaseName = file(data.illuminaRawFR).simpleName
    illuminaRawRBaseName = file(data.illuminaRawFR).simpleName

    """
    trimmomatic \
    PE -phred33 ${data.illuminaRawFR} ${data.illuminaRawRR} \
    ${illuminaRawFBaseName}_1P.trim.fastq.gz ${illuminaRawFBaseName}_unpaired_1P.trim.fastq.gz \
    ${illuminaRawRBaseName}_2P.trim.fastq.gz ${illuminaRawRBaseName}_unpaired_2P.trim.fastq.gz \
    ${ILLUMINACLIP} ${HEADCROP} ${CROP} ${SLIDINGWINDOW} ${MINLEN}
    """
}

process bbduk {

    publishDir "$baseDir/results/dedup", mode: 'copy'
    container 'staphb/bbtools:39.13'

    input:
    tuple val(key), path(trimmo_f_pair_path), path(trimmo_r_pair_path)
    path adapters_path

    output:
    tuple val(key), path("${trimmo_f_pair_path.simpleName}.bbduk.fq.gz"), path("${trimmo_r_pair_path.simpleName}.bbduk.fastq.gz")

    script:
    key = key

    """
    bbduk.sh \
        in1=${trimmo_f_pair_path} \
        in2=${trimmo_r_pair_path} \
        out1=${trimmo_f_pair_path.simpleName}.bbduk.fq.gz \
        out2=${trimmo_r_pair_path.simpleName}.bbduk.fastq.gz \
        ref=${adapters_path} \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        qtrim=rl \
        trimq=20
    """
}

process clumpify {

  publishDir "$baseDir/results/dedup", mode: 'copy'
  container 'staphb/bbtools:39.13'

  input:
  tuple val(key), path(bbduk_f_map_path), path(bbduk_r_map_path)

  output:
  tuple val(key), path("${bbduk_f_map_path.simpleName}.dedup.fq.gz"), path("${bbduk_r_map_path.simpleName}.dedup.fq.gz")

  script:
  """
  clumpify.sh in=${bbduk_f_map_path} in2=${bbduk_r_map_path} \
  out=${bbduk_f_map_path.simpleName}.dedup.fq.gz out2=${bbduk_r_map_path.simpleName}.dedup.fq.gz dedupe
  """
}



process fastqc {
// Short read quality control
    label "fastqc"
    publishDir "${params.publishDir}/qc", mode: 'copy'
    container 'staphb/fastqc:0.12.1'

    input:
    tuple val(key), path(forward_read), path(reverse_read)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results
    tuple val(key), env(total_bases), topic: "qc"

    script:

    key = key

    """
    fastqc -q ${forward_read} ${reverse_read} -o .

    total_bases_fwd=\$(grep "Total Sequences" ${baseDir}/results/qc/shortread/fastqc/${forward_read}_fastqc.html | awk '{print \$3}')
    total_bases_rev=\$(grep "Total Sequences" ${baseDir}/results/qc/shortread/fastqc/${reverse_read}_fastqc.html | awk '{print \$3}')
    
    total_bases=\$((total_bases_fwd + total_bases_rev))
    """
}

process nanoplot {
// Quality check for nanopore reads and Quality/Length Plots
    label "nanoplot"
    publishDir "${baseDir}/results/qc/nanoplot", mode: 'copy'
    container 'staphb/nanoplot:1.42.0'

    input:
    tuple val(key), file(long_reads), val(defaultThreads)

    output:
    file '*.png'
    file '*.html'
    file '*.txt'
    file("*_NanoStats.txt") 
    tuple val(key), env(RAW_ONT_READS), env(RAW_ONT_MEAN_LENGTH), env(RAW_ONT_MEAN_QUAL), env(RAW_ONT_TOTAL_BASES), emit: nanoplotStats

    script:
    key = key    

    """
    NanoPlot --fastq  ${long_reads} -t ${defaultThreads}  --title ${long_reads.baseName} -c darkblue \
    --plots kde dot --loglength --N50 --prefix ${long_reads.baseName} --verbose --prefix ${long_reads.baseName}_

    nanostats_file="${long_reads.baseName}_NanoStats.txt"

    RAW_ONT_READS=\$(grep -oP '(?<=Number of reads: )\\s+[\\d.]+' "\$nanostats_file" | tr -d " ") 

    RAW_ONT_MEAN_LENGTH=\$(grep -oP '(?<=Mean read length: )\\s+[\\d.]+' "\$nanostats_file" | tr -d " ") 

    RAW_ONT_MEAN_QUAL=\$(grep -oP '(?<=Mean read quality: )\\s+[\\d.]+' "\$nanostats_file" | tr -d " ") 

    RAW_ONT_TOTAL_BASES=\$(grep -oP '(?<=Total bases: )\\s+[\\d.]+' "\$nanostats_file" | tr -d " ") 
    """
}


process filtlong {

    publishDir "$baseDir/results", mode: 'copy'
    container "staphb/filtlong:0.2.1"

    input:
    tuple val(key), path(clump_f_dedup_path), path(clump_r_dedup_path), path(ONT_RAW_READS), val(defaultThreads)
    

    output:
    tuple val(key), path(clump_f_dedup_path), path(clump_r_dedup_path), path("${ONT_RAW_READS.simpleName}.filtered.fq.gz"), emit: outputForRatatosk
    tuple val(key), path("${ONT_RAW_READS.simpleName}.filtered.fq.gz"), emit: outputForNanoplot

    script:
    key = key

    """
    filtlong --min_length 1000 --min_mean_q 8 --target_bases 500000000 \\
            --trim -1 $clump_f_dedup_path \\
            -2 $clump_r_dedup_path \\
            $ONT_RAW_READS | gzip > ${ONT_RAW_READS.simpleName}.filtered.fq.gz
    """
}

process ratatosk {

    publishDir "$baseDir/results", mode: 'copy'
    container 'ajslee/ratatosk:0.9.0'


    input:
    tuple val(key), path(clump_f_dedup_path), path (clump_r_dedup_path), path(gzipped_filtered_ONT_reads), val(CPU_THREADS)

    output:
    tuple val(key), path("${gzipped_filtered_ONT_reads.simpleName}.corrected.fq.gz")

    script:
    """
    ratatosk correct -G -v -s ${clump_f_dedup_path} -s ${clump_r_dedup_path} \
                      -l ${gzipped_filtered_ONT_reads} \
                      -o ${gzipped_filtered_ONT_reads.simpleName}.corrected.fq.gz \
                      -c ${CPU_THREADS} 
    """
}


workflow {

    mainChannel = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))
    dataChannel = mainChannel.map { 
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
        ]
    }
    
    trimmomaticParams = Channel.of(
        ["${params.ILLUMINACLIP}", "${params.HEADCROP}", "${params.CROP}", "${params.SLIDINGWINDOW}", "${params.MINLEN}"]
        )
    barcodeFileChannel = Channel.fromPath("${params.BARCODE_FILE}")
    
    trimmomaticOutput = trimmomatic(dataChannel, barcodeFileChannel, trimmomaticParams)
    bbdukOutput = bbduk(trimmomaticOutput, file("${params.BARCODE_FILE}"))
    clumpifyOutput = clumpify(bbdukOutput)
    fastqc(clumpifyOutput) 
    
    dataChannel.filter { it.ontRawReads }
        .map { [key: it.key, long_reads: file(it.ontRawReads), defaultThreads: "${params.nanoplotThreads}"] }
        .set { longReadsChannel } 


    flattenedLongReadsChannel = longReadsChannel.map { [it.key, it.long_reads, it.nanoplotThreads]}
    joinedClumpifyAndLongChannel = clumpifyOutput.join(flattenedLongReadsChannel)
    joinedClumpifyAndLongChannel
        .filter { it[0] && it[1] && it[2] }
        .set{ clumpifyAndLongReads }

    filtlongOutput = filtlong(clumpifyAndLongReads)
    correctedLongReads = ratatosk(filtlongOutput.outputForRatatosk.map { it + ["${params.ratatoskThreads}"] })

    filtNanoplotInput =(filtlongOutput.outputForNanoplot.map { it + [params.nanoplotThreads] })
    correctedNanoplotInput = (correctedLongReads.map { it + [params.nanoplotThreads] })
    nanoplot(longReadsChannel.mix(filtNanoplotInput).mix(correctedNanoplotInput))
}

//staphb/masurca:4.1.0