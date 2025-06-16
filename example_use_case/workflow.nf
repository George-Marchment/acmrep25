//Command to execute
//nextflow workflow.nf


process filter_data {
    input:
        path input_file
    output:
        path "adult_patients.csv"

    script:
        """
        awk -F ',' 'NR==1 || \$3 >= 18' $input_file > "adult_patients.csv"
        """
}

process preprocessing {
    label 'python2'
    publishDir 'results/', mode: 'copy'

    input:
        path script
        file data
    output:
        path "*" 

    script:
        """
        python ${script} ${data}
        """
}

process analysis {
    label 'python3'
    

    input:
        path script
        file data
    output:
        path "*.csv" 

    script:
        """
        python3 ${script} ${data}
        """
}

process show_clusters {
    label 'R'
    publishDir 'results/', mode: 'copy'

    input:
        path script
        file data
    output:
        path "*.pdf" 

    script:
        """
        Rscript ${script} ${data}
        """
}



workflow {
    filtered_data = filter_data(file("patients.csv"))
    cleaned_data = preprocessing(Channel.fromPath('src/python2_script.py', followLinks: false), filtered_data)
    points = analysis(Channel.fromPath('src/clustering.py', followLinks: false), cleaned_data)
    show_clusters(Channel.fromPath('src/generate_cluster.R', followLinks: false), points)
}