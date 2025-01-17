include { samplesheetToList  } from 'plugin/nf-schema'


params.file_path = './assets/schema_input.json'


workflow {


    ch_input = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    ch_input.view()
} 
