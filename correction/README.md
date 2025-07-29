# Use-Case Workflow

This Nextflow Workflow is an implementation of the pratical work workflow presented in the [lecture](https://github.com/George-Marchment/acmrep25/blob/main/tutoriel_material/slides.pdf) for ACM REP 2025 tutorial "Computational Reproducibility With Scientific Workflows: Analysing viral genomes with Nextflow".

This workflow is composed of 3 main steps:
1. Mapping the reads onto the reference genome
2. Building the consensus sequence from the mapped reads
3. Identifying the clade of the virus based on the consensus sequence

<img align="center" src="../img/wf.png" width="60%">

To run the analysis (using *Apptainer*), simply run th  e following command:
```
nextflow workflow.nf
```

> To use *Docker* instead of *Apptainer*, simply change "`apptainer`" to "`docker`" in the `nextflow.config` file.