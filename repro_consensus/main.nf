process consensus {

    label 'python3'

    input:
    path script
    path true_value
    path predicted

    output:
    stdout

    script:
    """
    python3 ${script} ${true_value} ${predicted}
    """
}

workflow {
    true_value = Channel.fromPath("results/", type: 'dir')
    predicted_values = Channel.fromPath("results_participants/",  type: 'dir')
    consensus(Channel.fromPath('src/consensus.py', followLinks: false), true_value, predicted_values).view()
}