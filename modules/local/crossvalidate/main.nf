process CROSSVALIDATE_SAMPLES {
    tag "crossvalidate"

    // exec process: no container needed

    input:
    val(input_mapping)
    val(input_pairing)

    output:
    val(valid_samples),   emit: valid_samples
    val(invalid_samples), emit: invalid_samples
    val(valid_pairings),  emit: valid_pairings

    exec:
    // Validate that all samples in pairing have corresponding mapping entries
    def mappingSamples = input_mapping.collect{ it[0] }.unique()
    valid_pairings = input_pairing.findAll{ pair ->
        mappingSamples.contains(pair[0]) && mappingSamples.contains(pair[1])
    }
    def validIds = valid_pairings.flatten().unique()
    valid_samples = input_mapping.findAll{ validIds.contains(it[0]) }
    invalid_samples = input_mapping.findAll{ !validIds.contains(it[0]) }

    if (valid_samples.size() == 0) {
        error "CrossValidateSamples: No valid samples between pairing and mapping files."
    }
}
