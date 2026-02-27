# Tempo 100% Feature Parity Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Achieve 100% feature parity with the original MSKCC Tempo `develop` branch by replacing all `skip_*` flags with the original `--workflows`/`--assayType` control mechanism, adding all 45 missing processes as local modules, and wiring them into the main workflow.

**Architecture:** The original Tempo uses a `--workflows snv,sv,mutsig,...` comma-separated string parsed into `doWF_*` booleans, plus `--assayType exome|genome` to gate WGS-only features. We will replicate this exact control pattern in the nf-core migration, replacing the current `skip_*` flag approach. All 45 missing processes will be added as local modules following the existing module pattern (script + stub blocks, versions.yml, container directives).

**Tech Stack:** Nextflow DSL2, nf-core module conventions, original cmopipeline/* Docker containers, stub-mode testing

**Reference files:**
- Original workflow: `/sessions/sleepy-tender-newton/tempo-dev/dsl2.nf` (322 lines)
- Current nf-core workflow: `workflows/tempo.nf` (1081 lines)
- Process catalog: `/sessions/sleepy-tender-newton/ORIGINAL_PROCESS_CATALOG.md` (1683 lines)
- Original containers: `/sessions/sleepy-tender-newton/tempo-dev/conf/containers.config`

---

## Phase 1: Replace skip_* Flags with workflows/assayType Control

### Task 1.1: Add workflows and assayType params to nextflow.config

**Files:**
- Modify: `nextflow.config` (lines 55-79)

**Step 1: Edit nextflow.config to replace skip params with workflows/assayType**

Replace the `skip_*` block and add `workflows`/`assayType`/`aggregate` params:

```groovy
    // Pipeline control — mirrors original Tempo
    workflows                  = ''         // Comma-separated: snv,sv,mutsig,lohhla,facets,qc,msisensor,germSNV,germSV
    assay_type                 = 'exome'    // 'exome' or 'genome'
    aggregate                  = false      // true or path to aggregate TSV
    svcnv                      = 'hisens'   // 'ascat', 'hisens', or 'purity' (WGS SV CNV source)
    cosmic                     = 'v3'       // Cosmic signature version for tempoSig

    // Keep only these skip flags (not workflow-level)
    skip_qc                    = false
    skip_multiqc               = false
```

Remove all other `skip_*` params (skip_somatic_snv, skip_somatic_sv, skip_germline_snv, skip_germline_sv, skip_facets, skip_msi, skip_lohhla, skip_polysolver, skip_neoantigen, skip_annotation, skip_vcf2maf, skip_somatic_combine, skip_germline_combine, skip_alignment, skip_markduplicates, skip_bqsr).

**Step 2: Run stub test to verify config parses**

Run: `nextflow run . -profile test -stub --outdir /tmp/t1 -with-report false -with-timeline false -with-trace false 2>&1 | head -20`
Expected: No config parse errors

**Step 3: Commit**

```bash
git add nextflow.config
git commit -m "feat: replace skip_* flags with --workflows/--assayType control"
```

---

### Task 1.2: Add doWF_* boolean derivation to workflows/tempo.nf

**Files:**
- Modify: `workflows/tempo.nf` (top of `main:` block, ~line 84)

**Step 1: Add workflow flag parsing at top of main block**

Insert after `ch_versions = Channel.empty()` / `ch_multiqc_files = Channel.empty()`:

```groovy
    // =============================================
    // WORKFLOW CONTROL FLAGS (mirrors original Tempo dsl2.nf)
    // =============================================

    // Parse --workflows string into boolean flags
    def WFs = params.workflows instanceof Boolean ? '' : (params.workflows ?: '')
    def wfList = WFs.split(',').collect{ it.trim().toLowerCase() }.unique().findAll{ it }

    def doWF_SNV        = ['snv', 'mutsig'].any{ it in wfList }
    def doWF_SV         = 'sv' in wfList
    def doWF_manta      = ['snv', 'sv', 'mutsig'].any{ it in wfList }
    def doWF_facets     = ['lohhla', 'facets', 'snv', 'mutsig', 'germsnv'].any{ it in wfList }
    def doWF_loh        = ['lohhla', 'snv', 'mutsig'].any{ it in wfList }
    def doWF_germSNV    = 'germsnv' in wfList
    def doWF_germSV     = 'germsv' in wfList
    def doWF_QC         = 'qc' in wfList
    def doWF_msiSensor  = 'msisensor' in wfList
    def doWF_mutSig     = 'mutsig' in wfList
    def doWF_mdParse    = doWF_manta && doWF_facets && doWF_loh && doWF_SNV && doWF_msiSensor && doWF_mutSig

    // WGS-only: ASCAT/BRASS/HRDetect gated by assay_type == 'genome'
    def isWGS = params.assay_type == 'genome'

    // If SV workflow + WGS + facets-based CNV source, ensure facets runs
    if (doWF_SV && isWGS && ['hisens','purity'].contains(params.svcnv)) {
        doWF_facets = true
    }
```

**Step 2: Replace all `if (!params.skip_*)` blocks with `if (doWF_*)` equivalents**

Mapping:
- `if (!params.skip_somatic_snv)` → `if (doWF_SNV)`
- `if (!params.skip_somatic_sv)` → `if (doWF_SV)`
- `if (!params.skip_somatic_combine)` → stays inside `if (doWF_SNV)` (it's part of SNV)
- `if (!params.skip_facets)` → `if (doWF_facets)`
- `if (!params.skip_msi)` → `if (doWF_msiSensor)`
- `if (!params.skip_polysolver)` → `if (doWF_loh)` (polysolver is part of LOH workflow)
- `if (!params.skip_lohhla)` → stays inside `if (doWF_loh)` block
- `if (!params.skip_annotation && !params.skip_somatic_snv)` → `if (doWF_SNV)` (annotation is inherent to SNV)
- `if (!params.skip_vcf2maf)` → stays inside SNV annotation block
- `if (!params.skip_germline_snv)` → `if (doWF_germSNV)`
- `if (!params.skip_germline_combine)` → stays inside `if (doWF_germSNV)` block
- `if (!params.skip_germline_sv)` → `if (doWF_germSV)`
- `if (!params.skip_qc)` → `if (doWF_QC)`
- `if (!params.skip_multiqc)` → `if (doWF_QC && !params.skip_multiqc)`

**Step 3: Run stub test**

Run:
```bash
PERL5LIB=/sessions/sleepy-tender-newton/perl-lib \
PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/sessions/sleepy-tender-newton/stub-tools" \
NXF_JAVA_HOME=/sessions/sleepy-tender-newton/.local/lib/python3.10/site-packages/jdk4py/java-runtime \
NXF_OFFLINE=true \
/usr/bin/bash /sessions/sleepy-tender-newton/.local/bin/nextflow run . \
  -profile test \
  --workflows 'snv,sv,facets,lohhla,qc,msisensor,germSNV,germSV' \
  --outdir /tmp/t2 \
  -with-report false -with-timeline false -with-trace false \
  -stub
```
Expected: All existing tasks pass (75 tasks, 0 failures)

**Step 4: Commit**

```bash
git add workflows/tempo.nf
git commit -m "feat: wire doWF_* workflow control booleans into tempo.nf"
```

---

### Task 1.3: Update test configs for workflows param

**Files:**
- Modify: `conf/test.config`
- Modify: `conf/test_comprehensive.config`

**Step 1: Update test.config**

Replace all `skip_*` lines with a minimal `workflows` setting:

```groovy
    // Pipeline options
    outdir                = 'results'
    workflows             = ''   // minimal test: alignment only
```

**Step 2: Update test_comprehensive.config**

Replace all `skip_*` lines with the full workflow set:

```groovy
    // ---- Enable ALL calling paths ----
    workflows              = 'snv,sv,mutsig,lohhla,facets,qc,msisensor,germSNV,germSV'

    // Skip neoantigen (no module wired yet)
    // Neoantigen will be enabled when module is added in Phase 2
```

Keep all the reference file paths unchanged.

**Step 3: Run both profiles**

Minimal test:
```bash
nextflow run . -profile test --outdir /tmp/t3 -stub ...
```
Expected: Only alignment tasks run

Comprehensive test:
```bash
nextflow run . -profile test,test_comprehensive --outdir /tmp/t4 -stub ...
```
Expected: All paths enabled, ~75+ tasks pass

**Step 4: Commit**

```bash
git add conf/test.config conf/test_comprehensive.config
git commit -m "feat: update test configs to use --workflows instead of skip_* flags"
```

---

## Phase 2: Add Missing Local Modules (Groups A-F)

Each module follows this pattern (reference: `modules/local/delly/call/main.nf`):
1. Process declaration with tag, label, container
2. `input:` / `output:` blocks with nf-core meta map convention
3. `script:` block with the original command
4. `stub:` block that creates empty outputs + versions.yml
5. `path "versions.yml", emit: versions`

### Task 2.1: Group F — Utility Modules (4 modules)

These have no dependencies on other new modules and are needed earliest.

**Files to create:**
- `modules/local/splitintervals/main.nf` (CreateScatteredIntervals)
- `modules/local/germline/combine_hc_vcf/main.nf` (GermlineCombineHaplotypecallerVcf)
- `modules/local/metadata_parser/main.nf` (MetaDataParser)
- `modules/local/crossvalidate/main.nf` (CrossValidateSamples)

**Step 1: Create modules/local/splitintervals/main.nf**

```groovy
process SPLIT_INTERVALS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.1.0.0' :
        'broadinstitute/gatk:4.1.0.0' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(dict)
    path(intervals)
    val(scatter_count)

    output:
    path("*.interval_list"), emit: interval_lists
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.assay_type == 'genome' ? 'INTERVAL_SUBDIVISION' : 'BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW'
    """
    gatk SplitIntervals \\
        --reference ${fasta} \\
        --intervals ${intervals} \\
        --scatter-count ${scatter_count} \\
        --subdivision-mode ${mode} \\
        --output scattered

    for i in scattered/*.interval_list; do
        BASENAME=\$(basename \$i)
        mv \$i scattered-\$BASENAME
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '2s/.*v//p')
    END_VERSIONS
    """

    stub:
    """
    touch scattered-0001.interval_list scattered-0002.interval_list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.1.0.0
    END_VERSIONS
    """
}
```

**Step 2: Create modules/local/germline/combine_hc_vcf/main.nf**

```groovy
process GERMLINE_COMBINE_HC_VCF {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/bcftools-vt:1.1.1' :
        'cmopipeline/bcftools-vt:1.1.1' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)
    tuple val(meta2), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${prefix}.haplotypecaller.vcf.gz"), path("${prefix}.haplotypecaller.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools concat \\
        --allow-overlaps \\
        ${vcfs} | \\
    bcftools sort | \\
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --check-ref s \\
        --multiallelics -both | \\
    bcftools norm --rm-dup all \\
        --output-type z \\
        --output ${prefix}.haplotypecaller.vcf.gz

    tabix --preset vcf ${prefix}.haplotypecaller.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.haplotypecaller.vcf.gz
    touch ${prefix}.haplotypecaller.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
```

**Step 3: Create modules/local/metadata_parser/main.nf**

```groovy
process METADATA_PARSER {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/metadataparser:0.5.9' :
        'cmopipeline/metadataparser:0.5.9' }"

    input:
    tuple val(meta), path(purity_out), path(maf_file), path(qc_output), path(msi_file), path(mutsig), path(polysolver_file)
    path(coding_bed)

    output:
    tuple val(meta), path("*.sample_data.txt"), emit: metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_metadata_file.py \\
        --sampleID ${prefix} \\
        --tumorID ${meta.tumor_id} \\
        --normalID ${meta.normal_id} \\
        --facetsPurity_out ${purity_out} \\
        --facetsQC ${qc_output} \\
        --MSIsensor_output ${msi_file} \\
        --mutational_signatures_output ${mutsig} \\
        --polysolver_output ${polysolver_file} \\
        --MAF_input ${maf_file} \\
        --coding_baits_BED ${coding_bed}

    mv ${prefix}_metadata.txt ${prefix}.sample_data.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metadataparser: 0.5.9
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sample_data.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metadataparser: 0.5.9
    END_VERSIONS
    """
}
```

**Step 4: Create modules/local/crossvalidate/main.nf**

Note: This is a Nextflow `exec` process (pure Groovy, no container). In the nf-core migration we keep it as an exec process.

```groovy
process CROSSVALIDATE_SAMPLES {
    tag "crossvalidate"

    // exec process: no container needed

    input:
    val(input_mapping)
    val(input_pairing)

    output:
    val(valid_samples),  emit: valid_samples
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
```

**Step 5: Run stub test to verify new modules parse**

```bash
nextflow run . -profile test,test_comprehensive --outdir /tmp/t5 -stub ...
```
Expected: Pipeline parses (new modules not yet wired, so task count unchanged)

**Step 6: Commit**

```bash
git add modules/local/splitintervals/ modules/local/germline/combine_hc_vcf/ modules/local/metadata_parser/ modules/local/crossvalidate/
git commit -m "feat: add Group F utility local modules (splitintervals, combine_hc_vcf, metadata_parser, crossvalidate)"
```

---

### Task 2.2: Group A — SV Pipeline Modules (9 modules)

**Files to create:**
- `modules/local/svaba/somatic/main.nf`
- `modules/local/svaba/germline/main.nf`
- `modules/local/svtools/vcf2bedpe_somatic/main.nf`
- `modules/local/svtools/vcf2bedpe_germline/main.nf`
- `modules/local/iannotatesv/somatic/main.nf`
- `modules/local/iannotatesv/germline/main.nf`
- `modules/local/clustersv/main.nf`
- `modules/local/svcircos/main.nf`
- `modules/local/svclone/main.nf`

Each module follows the same pattern. Key containers from the original:
- SvABA: `cmopipeline/svaba:0.0.1`
- SVVcf2Bedpe: `cmopipeline/svtools:0.0.3`
- AnnotateSVBedpe: `cmopipeline/iannotatesv:0.0.2`
- ClusterSV: `cmopipeline/clustersv:0.0.1`
- SVCircos: `cmopipeline/biocircos:0.0.1`
- SVclone: `cmopipeline/svclone:0.0.1`

**Step 1: Create all 9 module files**

For each module: write the process with nf-core meta map convention, original script block, stub block with empty outputs + versions.yml. Use the source code from `/sessions/sleepy-tender-newton/ORIGINAL_PROCESS_CATALOG.md` for exact commands.

Inputs/outputs must use `tuple val(meta), path(...)` convention instead of the original's `tuple val(idTumor), val(idNormal), ...` pattern. The meta map carries `tumor_id`, `normal_id`, `patient`, etc.

**Step 2: Run parse check (no wiring yet)**

```bash
nextflow run . -profile test,test_comprehensive --outdir /tmp/t6 -stub ...
```

**Step 3: Commit**

```bash
git add modules/local/svaba/ modules/local/svtools/ modules/local/iannotatesv/ modules/local/clustersv/ modules/local/svcircos/ modules/local/svclone/
git commit -m "feat: add Group A SV pipeline local modules (svaba, svtools, iannotatesv, clustersv, svcircos, svclone)"
```

---

### Task 2.3: Group B — WGS-Only Modules (8 modules)

**Files to create:**
- `modules/local/ascat/allelecount/main.nf` (runAscatAlleleCount)
- `modules/local/ascat/run/main.nf` (runAscat)
- `modules/local/brass/generate_bas/main.nf` (generateBasFile)
- `modules/local/brass/input/main.nf` (SomaticRunBRASSInput)
- `modules/local/brass/cover/main.nf` (SomaticRunBRASSCover)
- `modules/local/brass/run/main.nf` (SomaticRunBRASS)
- `modules/local/hrdetect/main.nf` (HRDetect)
- `modules/local/svsignatures/main.nf` (RunSVSignatures)

Key containers:
- ASCAT: `quay.io/wtsicgp/ascatNgs:4.4.0`
- generateBasFile: `quay.io/wtsicgp/pcap-core:5.5.0`
- BRASS: `cmopipeline/brass:0.0.2`
- HRDetect/SVSignatures: `cmopipeline/signaturetoolslib:0.0.1`

**Step 1: Create all 8 module files**

Use catalog source code. All WGS-only modules should include `when: params.assay_type == 'genome'` in addition to the standard `task.ext.when` check.

**Step 2: Commit**

```bash
git add modules/local/ascat/ modules/local/brass/ modules/local/hrdetect/ modules/local/svsignatures/
git commit -m "feat: add Group B WGS-only local modules (ascat, brass, hrdetect, svsignatures)"
```

---

### Task 2.4: Group C — Annotation & Signatures Modules (6 modules)

**Files to create:**
- `modules/local/germline/annotate_maf/main.nf` (GermlineAnnotateMaf)
- `modules/local/somatic/facets_annotation/main.nf` (SomaticFacetsAnnotation)
- `modules/local/germline/facets_annotation/main.nf` (GermlineFacetsAnnotation)
- `modules/local/facets/preview_qc/main.nf` (DoFacetsPreviewQC)
- `modules/local/neoantigen/main.nf` (RunNeoantigen)
- `modules/local/mutsig/main.nf` (RunMutationSignatures)

Key containers:
- GermlineAnnotateMaf: `cmopipeline/vcf2maf:vep88_1.2.7`
- Facets annotation: `cmopipeline/facets-suite-preview-htstools:0.0.1`
- Neoantigen: `cmopipeline/neoantigen:0.3.3`
- MutSig: `cmopipeline/temposig:0.2.3`

**Step 1: Create all 6 module files**

**Step 2: Commit**

```bash
git add modules/local/germline/annotate_maf/ modules/local/somatic/facets_annotation/ modules/local/germline/facets_annotation/ modules/local/facets/preview_qc/ modules/local/neoantigen/ modules/local/mutsig/
git commit -m "feat: add Group C annotation & signatures local modules"
```

---

### Task 2.5: Group D — QC & Reporting Modules (5 modules)

**Files to create:**
- `modules/local/alfred/main.nf` (QcAlfred)
- `modules/local/conpair/all/main.nf` (QcConpairAll)
- `modules/local/multiqc/sample/main.nf` (SampleRunMultiQC)
- `modules/local/multiqc/somatic/main.nf` (SomaticRunMultiQC)
- `modules/local/multiqc/cohort/main.nf` (CohortRunMultiQC)

Key containers:
- Alfred: `cmopipeline/alfred:v0.1.17`
- ConpairAll: `cmopipeline/conpair:v0.3.3`
- MultiQC (all 3): `cmopipeline/multiqc:0.1.3`

**Step 1: Create all 5 module files**

**Step 2: Commit**

```bash
git add modules/local/alfred/ modules/local/conpair/all/ modules/local/multiqc/
git commit -m "feat: add Group D QC & reporting local modules (alfred, conpair_all, multiqc variants)"
```

---

### Task 2.6: Group E — Aggregation Modules (13 modules)

**Files to create:**
- `modules/local/aggregate/somatic_maf/main.nf`
- `modules/local/aggregate/somatic_sv/main.nf`
- `modules/local/aggregate/somatic_facets/main.nf`
- `modules/local/aggregate/somatic_netmhc/main.nf`
- `modules/local/aggregate/somatic_metadata/main.nf`
- `modules/local/aggregate/somatic_lohhla/main.nf`
- `modules/local/aggregate/somatic_hrdetect/main.nf`
- `modules/local/aggregate/somatic_svclone/main.nf`
- `modules/local/aggregate/somatic_svsignatures/main.nf`
- `modules/local/aggregate/germline_maf/main.nf`
- `modules/local/aggregate/germline_sv/main.nf`
- `modules/local/aggregate/qc_bam/main.nf`
- `modules/local/aggregate/qc_conpair/main.nf`

Most aggregate modules use `cmopipeline/bcftools-vt:1.1.1` or `1.2.0` containers and are simple `awk`/`cat` merge scripts.

**Step 1: Create all 13 module files**

**Step 2: Commit**

```bash
git add modules/local/aggregate/
git commit -m "feat: add Group E aggregation local modules (13 aggregate processes)"
```

---

## Phase 3: Wire New Modules into Workflow

### Task 3.1: Wire Group F utilities into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements for Group F modules**

```groovy
include { SPLIT_INTERVALS          } from '../modules/local/splitintervals/main'
include { GERMLINE_COMBINE_HC_VCF  } from '../modules/local/germline/combine_hc_vcf/main'
include { METADATA_PARSER          } from '../modules/local/metadata_parser/main'
```

Note: CrossValidateSamples is not needed in nf-core migration since we use samplesheet-based input validation via nf-schema.

**Step 2: Wire GERMLINE_COMBINE_HC_VCF into the germline SNV block**

Inside `if (doWF_germSNV)`, after GATK4_HAPLOTYPECALLER, when doing scattered calling:

If scattered calling is used, the HC VCFs per interval need combining before going to germline combine. This matches the original flow where HC called per-scatter, then combined.

For the initial migration, keep single-interval HC calling (matching current behavior) and add GERMLINE_COMBINE_HC_VCF as a future enhancement for scatter support.

**Step 3: Wire METADATA_PARSER into the mdParse block**

Inside `if (doWF_mdParse)`, after combining facets + maf + qc + msi + mutsig + polysolver channels:

```groovy
    if (doWF_mdParse) {
        // Combine all required channels by patient/pair
        // ... channel joins ...
        METADATA_PARSER(ch_metadata_input, ch_coding_bed)
    }
```

**Step 4: Run stub test**

**Step 5: Commit**

---

### Task 3.2: Wire Group C annotation modules into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements**

```groovy
include { SOMATIC_FACETS_ANNOTATION  } from '../modules/local/somatic/facets_annotation/main'
include { GERMLINE_FACETS_ANNOTATION } from '../modules/local/germline/facets_annotation/main'
include { GERMLINE_ANNOTATE_MAF      } from '../modules/local/germline/annotate_maf/main'
include { FACETS_PREVIEW_QC          } from '../modules/local/facets/preview_qc/main'
include { NEOANTIGEN                 } from '../modules/local/neoantigen/main'
include { MUTATION_SIGNATURES        } from '../modules/local/mutsig/main'
```

**Step 2: Wire each into appropriate workflow block**

- `SOMATIC_FACETS_ANNOTATION`: Inside `if (doWF_SNV)` after VCF2MAF, joining facets hisens RData + MAF
- `GERMLINE_FACETS_ANNOTATION`: Inside `if (doWF_germSNV)` after germline VCF2MAF, joining facets
- `GERMLINE_ANNOTATE_MAF`: Inside `if (doWF_germSNV)` after GERMLINE_COMBINE_CHANNEL
- `FACETS_PREVIEW_QC`: Inside `if (doWF_facets)` after FACETS
- `NEOANTIGEN`: Inside `if (doWF_loh && doWF_SNV)` after POLYSOLVER + somatic MAF
- `MUTATION_SIGNATURES`: Inside `if (doWF_mutSig)` after somatic MAF

**Step 3: Run stub test**

**Step 4: Commit**

---

### Task 3.3: Wire Group A SV modules into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements for all 9 SV modules**

**Step 2: Wire somatic SV modules**

Inside `if (doWF_SV)`, after SOMATIC_MERGE_SV:

1. `SVABA_SOMATIC` — call alongside DELLY, joining with tumor-normal BAMs
2. Wire SVABA output into SOMATIC_MERGE_SV (extend merge to include svaba+delly+manta)
3. `SV_VCF2BEDPE_SOMATIC` — after SOMATIC_MERGE_SV, convert merged VCF to BEDPE
4. `ANNOTATE_SV_BEDPE_SOMATIC` — after VCF2BEDPE, with repeat_masker/blacklist references
5. `CLUSTER_SV` — after ANNOTATE_SV_BEDPE_SOMATIC
6. `SV_CIRCOS` — after CLUSTER_SV, joining with FACETS CNV

**Step 3: Wire germline SV modules**

Inside `if (doWF_germSV)`, after GERMLINE_MERGE_SV:

1. `SVABA_GERMLINE` — alongside germline Delly/Manta
2. Wire into GERMLINE_MERGE_SV
3. `SV_VCF2BEDPE_GERMLINE` — after merge
4. `ANNOTATE_SV_BEDPE_GERMLINE` — after VCF2BEDPE

**Step 4: Wire SVclone**

Inside `if (doWF_SV && doWF_SNV && isWGS)`, after SV annotation + somatic MAF + FACETS:

```groovy
    SVCLONE(ch_svclone_input, svclone_wrapper_script)
```

**Step 5: Run stub test**

**Step 6: Commit**

---

### Task 3.4: Wire Group B WGS-only modules into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements for ASCAT, BRASS, HRDetect, SVSignatures**

**Step 2: Wire ASCAT**

Inside `if (doWF_SV && isWGS && params.svcnv == 'ascat')`:

```groovy
    ASCAT_ALLELECOUNT(ch_ascat_allelecount_input, fasta, fai, snp_gc_corrections)
    ASCAT_RUN(ch_ascat_run_input, fasta, fai, snp_gc_corrections)
```

**Step 3: Wire BRASS**

Inside `if (doWF_SV && isWGS)`:

```groovy
    GENERATE_BAS_FILE(ch_bas_input, fasta, fai)
    BRASS_INPUT(ch_brass_input_input, fasta, fai, brass_ref_dir, vagrent_ref_dir)
    BRASS_COVER(ch_brass_cover_input, fasta, fai, brass_ref_dir, vagrent_ref_dir)
    BRASS_RUN(ch_brass_run_input, fasta, fai, brass_ref_dir, vagrent_ref_dir)
```

Wire BRASS output into SV merge (extends the merge to include brass+delly+manta+svaba).

**Step 4: Wire HRDetect + SVSignatures**

Inside `if (doWF_SV && doWF_SNV && isWGS)`:

```groovy
    HRDETECT(ch_hrdetect_input, hrdetect_script)
    SV_SIGNATURES(ch_svsig_input, sv_signature_script)
```

**Step 5: Run stub test**

**Step 6: Commit**

---

### Task 3.5: Wire Group D QC modules into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements for QC modules**

**Step 2: Wire QcAlfred**

Inside `if (doWF_QC)`, alongside existing PICARD_COLLECTHSMETRICS and QUALIMAP_BAMQC:

```groovy
    ALFRED(ch_recal_bam_bai, fasta)
```

**Step 3: Wire QcConpairAll**

Inside `if (doWF_QC)`, after CONPAIR_PILEUP tumor-normal pairing:

```groovy
    CONPAIR_ALL(ch_conpair_all_input, [fasta, fai, dict])
```

**Step 4: Wire Sample/Somatic/Cohort MultiQC**

Replace the existing nf-core MULTIQC call with the three Tempo-specific MultiQC processes:

- `SAMPLE_MULTIQC` — per-sample QC (alfred + fastp + qualimap + hsmetrics)
- `SOMATIC_MULTIQC` — per-pair QC (conpair + qualimap T/N + facets QC)
- `COHORT_MULTIQC` — cohort-level (all of the above aggregated)

Keep the nf-core MULTIQC as a fallback when Tempo-specific MultiQC configs are not available.

**Step 5: Run stub test**

**Step 6: Commit**

---

### Task 3.6: Wire Group E aggregation modules into tempo.nf

**Files:**
- Modify: `workflows/tempo.nf`

**Step 1: Add include statements for all 13 aggregate modules**

**Step 2: Wire aggregation block**

At the end of the workflow, add an `if (params.aggregate)` block:

```groovy
    if (params.aggregate) {
        // Collect outputs from all workflows and aggregate
        if (doWF_SNV)      { AGGREGATE_SOMATIC_MAF(ch_all_somatic_mafs) }
        if (doWF_SV)       { AGGREGATE_SOMATIC_SV(ch_all_somatic_bedpes) }
        if (doWF_facets)   { AGGREGATE_SOMATIC_FACETS(ch_all_facets_outputs) }
        if (doWF_loh)      { AGGREGATE_SOMATIC_LOHHLA(ch_all_lohhla_outputs) }
        if (doWF_germSNV)  { AGGREGATE_GERMLINE_MAF(ch_all_germline_mafs) }
        if (doWF_germSV)   { AGGREGATE_GERMLINE_SV(ch_all_germline_bedpes) }
        if (doWF_QC) {
            AGGREGATE_QC_BAM(ch_all_alfred_outputs)
            AGGREGATE_QC_CONPAIR(ch_all_conpair_outputs)
        }
        if (doWF_mdParse)  { AGGREGATE_SOMATIC_METADATA(ch_all_metadata_outputs) }
        if (doWF_SNV && doWF_loh) { AGGREGATE_SOMATIC_NETMHC(ch_all_neoantigen_outputs) }
        // WGS-only aggregations
        if (doWF_SV && doWF_SNV && isWGS) {
            AGGREGATE_SOMATIC_HRDETECT(ch_all_hrdetect_outputs)
            AGGREGATE_SOMATIC_SVCLONE(ch_all_svclone_outputs)
            AGGREGATE_SOMATIC_SVSIGNATURES(ch_all_svsig_outputs)
        }
    }
```

**Step 3: Run stub test with --aggregate**

**Step 4: Commit**

---

## Phase 4: Update Test Data & Final Verification

### Task 4.1: Add missing test reference files

**Files:**
- Create: `test-data/genome/` (various empty/stub reference files)
- Modify: `conf/test_comprehensive.config`

Some new modules need additional reference files that are not yet in `test-data/`. Create stub versions:

```bash
# WGS-specific references (BRASS, ASCAT)
touch test-data/reference/brass/HiDepth.bed.gz
touch test-data/reference/brass/brass_np.groups.gz
# ... etc

# SV blacklist references
touch test-data/genome/sv_blacklist.bed
touch test-data/genome/sv_blacklist.bedpe
touch test-data/genome/sv_blacklist_foldback.bedpe
touch test-data/genome/sv_blacklist_te.bedpe

# Neoantigen references
touch test-data/reference/neoantigen/neoantigen_cdna.fa
touch test-data/reference/neoantigen/neoantigen_cds.fa
```

Add corresponding params to `test_comprehensive.config`.

**Commit**

---

### Task 4.2: Full stub-mode verification

**Step 1: Run comprehensive stub test**

```bash
PERL5LIB=/sessions/sleepy-tender-newton/perl-lib \
PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/sessions/sleepy-tender-newton/stub-tools" \
NXF_JAVA_HOME=/sessions/sleepy-tender-newton/.local/lib/python3.10/site-packages/jdk4py/java-runtime \
NXF_OFFLINE=true \
/usr/bin/bash /sessions/sleepy-tender-newton/.local/bin/nextflow run . \
  -profile test,test_comprehensive \
  --workflows 'snv,sv,mutsig,lohhla,facets,qc,msisensor,germSNV,germSV' \
  --outdir /tmp/results_final \
  -with-report false -with-timeline false -with-trace false \
  -stub
```

Expected: ~120+ tasks (75 existing + ~45 new), 0 failures

**Step 2: Verify all outputs**

```bash
find /tmp/results_final -type f | wc -l
```

Expected: Significantly more output files than before

**Step 3: Run minimal test (alignment only)**

```bash
nextflow run . -profile test --outdir /tmp/results_minimal --workflows '' -stub ...
```

Expected: Only alignment tasks run

**Step 4: Run SNV-only test**

```bash
nextflow run . -profile test,test_comprehensive --outdir /tmp/results_snv --workflows 'snv' -stub ...
```

Expected: Alignment + somatic SNV tasks only

**Step 5: Commit**

```bash
git add -A
git commit -m "feat: complete 100% feature parity with original Tempo pipeline"
```

---

## Phase 5: Process-by-Process Verification Checklist

Run a final audit comparing original `containers.config` processes vs nf-core modules:

| Original Process | nf-core Module | Status |
|---|---|---|
| AlignReads | FASTP + BWAMEM2_MEM | EXISTS |
| MergeBamsAndMarkDuplicates | GATK4_MARKDUPLICATES | EXISTS |
| RunBQSR | GATK4_BASERECALIBRATOR + GATK4_APPLYBQSR | EXISTS |
| CreateScatteredIntervals | SPLIT_INTERVALS | NEW |
| RunMutect2 | GATK4_MUTECT2 | EXISTS |
| SomaticCombineMutect2Vcf | STRELKA2_COMBINE_SOMATIC | EXISTS |
| SomaticRunStrelka2 | STRELKA_SOMATIC | EXISTS |
| SomaticCombineChannel | SOMATIC_COMBINE_CHANNEL | EXISTS |
| SomaticAnnotateMaf | ENSEMBLVEP_VEP + VCF2MAF | EXISTS |
| SomaticFacetsAnnotation | SOMATIC_FACETS_ANNOTATION | NEW |
| SomaticDellyCall | DELLY_CALL_SOMATIC | EXISTS |
| DellyCombine | DELLY_COMBINE | EXISTS |
| SomaticRunManta | MANTA_SOMATIC | EXISTS |
| SomaticMergeSVs | SOMATIC_MERGE_SV | EXISTS |
| SomaticRunSvABA | SVABA_SOMATIC | NEW |
| SomaticSVVcf2Bedpe | SV_VCF2BEDPE_SOMATIC | NEW |
| SomaticAnnotateSVBedpe | ANNOTATE_SV_BEDPE_SOMATIC | NEW |
| SomaticRunClusterSV | CLUSTER_SV | NEW |
| SomaticRunSVCircos | SV_CIRCOS | NEW |
| SomaticRunSVclone | SVCLONE | NEW |
| DoFacets | FACETS | EXISTS |
| DoFacetsPreviewQC | FACETS_PREVIEW_QC | NEW |
| RunMsiSensor | MSISENSORPRO_SCAN + MSISENSORPRO_MSISOMATIC | EXISTS |
| RunPolysolver | POLYSOLVER | EXISTS |
| RunLOHHLA | LOHHLA | EXISTS |
| RunNeoantigen | NEOANTIGEN | NEW |
| RunMutationSignatures | MUTATION_SIGNATURES | NEW |
| MetaDataParser | METADATA_PARSER | NEW |
| runAscatAlleleCount | ASCAT_ALLELECOUNT | NEW |
| runAscat | ASCAT_RUN | NEW |
| generateBasFile | GENERATE_BAS_FILE | NEW |
| SomaticRunBRASSInput | BRASS_INPUT | NEW |
| SomaticRunBRASSCover | BRASS_COVER | NEW |
| runBRASS | BRASS_RUN | NEW |
| HRDetect | HRDETECT | NEW |
| RunSVSignatures | SV_SIGNATURES | NEW |
| GermlineRunHaplotypecaller | GATK4_HAPLOTYPECALLER | EXISTS |
| GermlineCombineHaplotypecallerVcf | GERMLINE_COMBINE_HC_VCF | NEW |
| GermlineRunStrelka2 | STRELKA2_GERMLINE | EXISTS |
| GermlineCombineChannel | GERMLINE_COMBINE_CHANNEL | EXISTS |
| GermlineAnnotateMaf | GERMLINE_ANNOTATE_MAF | NEW |
| GermlineFacetsAnnotation | GERMLINE_FACETS_ANNOTATION | NEW |
| GermlineDellyCall | DELLY_CALL_GERMLINE | EXISTS |
| GermlineRunManta | MANTA_GERMLINE | EXISTS |
| GermlineMergeSVs | GERMLINE_MERGE_SV | EXISTS |
| GermlineRunSvABA | SVABA_GERMLINE | NEW |
| GermlineSVVcf2Bedpe | SV_VCF2BEDPE_GERMLINE | NEW |
| GermlineAnnotateSVBedpe | ANNOTATE_SV_BEDPE_GERMLINE | NEW |
| QcPileup | CONPAIR_PILEUP | EXISTS |
| QcConpair | CONPAIR_CONCORDANCE | EXISTS |
| QcConpairAll | CONPAIR_ALL | NEW |
| QcAlfred | ALFRED | NEW |
| QcCollectHsMetrics | PICARD_COLLECTHSMETRICS | EXISTS |
| QcQualimap | QUALIMAP_BAMQC | EXISTS |
| SampleRunMultiQC | SAMPLE_MULTIQC | NEW |
| SomaticRunMultiQC | SOMATIC_MULTIQC | NEW |
| CohortRunMultiQC | COHORT_MULTIQC | NEW |
| SomaticAggregateMaf | AGGREGATE_SOMATIC_MAF | NEW |
| SomaticAggregateNetMHC | AGGREGATE_SOMATIC_NETMHC | NEW |
| SomaticAggregateFacets | AGGREGATE_SOMATIC_FACETS | NEW |
| SomaticAggregateSv | AGGREGATE_SOMATIC_SV | NEW |
| SomaticAggregateMetadata | AGGREGATE_SOMATIC_METADATA | NEW |
| SomaticAggregateLOHHLA | AGGREGATE_SOMATIC_LOHHLA | NEW |
| SomaticAggregateHRDetect | AGGREGATE_SOMATIC_HRDETECT | NEW |
| SomaticAggregateSVclone | AGGREGATE_SOMATIC_SVCLONE | NEW |
| SomaticAggregateSvSignatures | AGGREGATE_SOMATIC_SVSIGNATURES | NEW |
| GermlineAggregateMaf | AGGREGATE_GERMLINE_MAF | NEW |
| GermlineAggregateSv | AGGREGATE_GERMLINE_SV | NEW |
| QcBamAggregate | AGGREGATE_QC_BAM | NEW |
| QcConpairAggregate | AGGREGATE_QC_CONPAIR | NEW |

**Total: 67 processes → 67 nf-core equivalents = 100% parity**

---

## Execution Notes

- **All new modules use original Docker containers** — no software changes
- **Stub tests are the verification mechanism** — no Docker needed in dev environment
- **Phase 1 is the highest priority** — it restructures the control plane
- **Phase 2 tasks (A-F) are independent** — can be parallelized with subagents
- **Phase 3 tasks depend on Phase 2** — must wire after modules exist
- **WGS-only modules (Group B)** will only run when `--assay_type genome` is set
- **Aggregation (Group E)** only runs when `--aggregate` is set
