params.reference = "../data/genome/reference.fa"
params.fastq = "../data/reads"


process indexRef {
    label 'indexRef'

    input:
    file ref

    output:
    tuple val(ref.name), file("${ref.baseName}.*")

    script:
    """
    bwa index ${ref}
    """
    }

process mapping {
    label 'mapping'

    publishDir "results/bams/", mode:'copy'

    input: 
    tuple val(name), file(f1), file(f2)
    tuple val(refName), file(ref)

    output:
    tuple val(name), file("*.bam"), file("*.bai")
    
    script:
    """
    bwa mem -t ${task.cpus} ${refName} ${f1} ${f2} > tmp.sam  
    samtools sort -o ${name}.bam tmp.sam
    samtools index ${name}.bam
    """
    }


process consensus {
    publishDir "results/consensus/", mode: 'copy'
    
    input:
    tuple val(name), file(bam), file(bai)

    output:
    file "${name}.fa"

    script:
    """
    samtools mpileup -d 600000 -A -Q 0 -F 0 ${bam} | ivar consensus -q 20 -t 0 -m 5 -n N -p ${name}
    """
}

process nextcladedata {

    label 'getref'
    
    output:
    path "ncref"
    
    script:
    """
    nextclade dataset get --name 'sars-cov-2' --output-dir 'ncref'
    """
}

process nextclade {
    publishDir "results/", mode: 'copy'

    label 'nextclade'
    
    input:
    path ncref
    path seq
    
    output:
    path "annotations.tsv"
    
    script:
    """
    nextclade run --in-order --input-dataset ${ncref} --output-tsv 'annotations.tsv'  ${seq}
    sed -i 's/\\x0D\$//' annotations.tsv
    """
}

process pangolin {
    publishDir "results/", mode:'copy'

    input:
    file fa

    output:
    file "*.csv"

    script:
    """
    PATH=/opt/conda/envs/pangolin/bin/:\$PATH
    pangolin --analysis-mode usher '${fa}' -t ${task.cpus} --outfile lineage_report.csv
    """
}

workflow {
    referenceFile = file(params.reference)
    data = Channel.fromFilePairs(params.fastq+'/*_{1,2}.fastq.gz', flat:true)

    index = indexRef(referenceFile)
    bam = mapping(data,index)
    consens = consensus(bam)
    allconsens =  consens.collectFile(name:"all_consensus.fasta")
    ncdata=nextcladedata()
    nextclade(ncdata,allconsens)
    pangolin(allconsens)
}

