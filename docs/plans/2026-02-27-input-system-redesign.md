# Input System Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add TARGET column to samplesheet, remove sex column, add --bamMapping for pre-aligned BAM input, and add aggregate-only mode via --aggregate TSV with PATH column.

**Architecture:** Three independent input modes sharing the same downstream channel shape: (1) FASTQ samplesheet → alignment → recalibrated BAMs, (2) --bamMapping TSV → skip alignment, feed BAMs directly, (3) --aggregate TSV with PATH → skip all per-sample processing, run only cohort aggregation. The `target` field flows through the meta map for per-sample interval/bait selection (not yet wired to process inputs — that's a future task).

**Tech Stack:** Nextflow DSL2, nf-schema plugin, Python 3 (samplesheet validator)

---

### Task 1: Update samplesheet schema — remove sex, add target

**Files:**

- Modify: `assets/schema_input.json`

**Step 1: Edit schema_input.json**

Replace the entire file with this updated schema:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/mskcc/tempo/main/assets/schema_input.json",
  "title": "mskcc/tempo pipeline - params.input schema",
  "description": "Schema for the file provided with params.input",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "patient": {
        "type": "string",
        "pattern": "^\\S+$",
        "meta": ["patient"],
        "errorMessage": "Patient ID must be provided and cannot contain spaces"
      },
      "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "meta": ["sample"],
        "errorMessage": "Sample name must be provided and cannot contain spaces"
      },
      "status": {
        "type": "integer",
        "enum": [0, 1],
        "meta": ["status"],
        "errorMessage": "Status must be 0 (normal) or 1 (tumor)",
        "default": 0
      },
      "target": {
        "type": "string",
        "pattern": "^\\S+$",
        "meta": ["target"],
        "errorMessage": "Target/bait set must be provided (e.g., agilent, idt, wgs)"
      },
      "lane": {
        "type": "string",
        "meta": ["lane"],
        "default": "L001",
        "errorMessage": "Lane identifier (e.g., L001, L002)"
      },
      "fastq_1": {
        "type": "string",
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "format": "file-path",
        "exists": true,
        "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
      },
      "fastq_2": {
        "type": "string",
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "format": "file-path",
        "exists": true,
        "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
      }
    },
    "required": ["patient", "sample", "status", "target", "fastq_1", "fastq_2"]
  }
}
```

**Step 2: Verify the JSON is valid**

Run: `python3 -c "import json; json.load(open('assets/schema_input.json'))"`
Expected: No output (valid JSON)

**Step 3: Commit**

```bash
git add assets/schema_input.json
git commit -m "feat(input): remove sex column, add required target column to samplesheet schema"
```

---

### Task 2: Update samplesheet parsing subworkflow

**Files:**

- Modify: `subworkflows/local/utils_nfcore_tempo_pipeline/main.nf`

**Step 1: Update the samplesheet channel mapping**

In the `.map` closure (around line 71-91), update the meta map construction:

Replace:

```nextflow
def meta = [
    id:      meta_raw.sample,
    patient: meta_raw.patient,
    sample:  meta_raw.sample,
    status:  meta_raw.status instanceof Integer ? meta_raw.status : meta_raw.status.toInteger(),
    sex:     meta_raw.sex ?: 'NA',
    lane:    meta_raw.lane ?: 'L001'
]
```

With:

```nextflow
def meta = [
    id:      meta_raw.sample,
    patient: meta_raw.patient,
    sample:  meta_raw.sample,
    status:  meta_raw.status instanceof Integer ? meta_raw.status : meta_raw.status.toInteger(),
    target:  meta_raw.target,
    lane:    meta_raw.lane ?: 'L001'
]
```

**Step 2: Commit**

```bash
git add subworkflows/local/utils_nfcore_tempo_pipeline/main.nf
git commit -m "feat(input): parse target column instead of sex in samplesheet subworkflow"
```

---

### Task 3: Update workflow to remove sex references and carry target through meta

**Files:**

- Modify: `workflows/tempo.nf`

**Step 1: Update multi-lane merge grouping (line ~231)**

Replace:

```nextflow
def new_meta = meta.subMap('patient', 'sample', 'status', 'sex') + [id: meta.sample]
```

With:

```nextflow
def new_meta = meta.subMap('patient', 'sample', 'status', 'target') + [id: meta.sample]
```

**Step 2: Update tumor-normal pairing (line ~343-351)**

Replace the pair_meta construction:

```nextflow
def pair_meta = [
    id:        "${tumor_meta.sample}__${normal_meta.sample}",
    patient:   patient,
    tumor_id:  tumor_meta.sample,
    normal_id: normal_meta.sample,
    status:    1,
    sex:       tumor_meta.sex ?: 'NA'
]
```

With:

```nextflow
def pair_meta = [
    id:        "${tumor_meta.sample}__${normal_meta.sample}",
    patient:   patient,
    tumor_id:  tumor_meta.sample,
    normal_id: normal_meta.sample,
    status:    1,
    target:    tumor_meta.target
]
```

**Step 3: Commit**

```bash
git add workflows/tempo.nf
git commit -m "feat(input): replace sex with target in workflow meta maps"
```

---

### Task 4: Add --bamMapping parameter and parsing

**Files:**

- Modify: `nextflow.config` (add param)
- Modify: `subworkflows/local/utils_nfcore_tempo_pipeline/main.nf` (add BAM parsing)

**Step 1: Add params.bamMapping to nextflow.config**

After `input = null` (line 17), add:

```nextflow
bamMapping                 = null    // TSV: PATIENT, SAMPLE, STATUS, TARGET, BAM, BAI
```

**Step 2: Add BAM samplesheet parsing to PIPELINE_INITIALISATION**

After the existing `ch_samplesheet` creation (line ~92) and before `emit:`, add a new channel for BAM inputs. The full logic:

```nextflow
    //
    // Create channel from BAM mapping (if provided)
    //
    if (params.bamMapping) {
        Channel
            .fromPath(params.bamMapping, checkIfExists: true)
            .splitCsv(sep: '\t', header: true)
            .map { row ->
                if (!row.PATIENT || !row.SAMPLE || !row.TARGET || !row.BAM || !row.BAI) {
                    error("bamMapping TSV must have columns: PATIENT, SAMPLE, STATUS, TARGET, BAM, BAI. Found: ${row.keySet()}")
                }
                def meta = [
                    id:      row.SAMPLE,
                    patient: row.PATIENT,
                    sample:  row.SAMPLE,
                    status:  row.STATUS ? row.STATUS.toInteger() : 0,
                    target:  row.TARGET
                ]
                def bam = file(row.BAM, checkIfExists: true)
                def bai = file(row.BAI, checkIfExists: true)
                if (!bam.name.endsWith('.bam')) {
                    error("BAM file must end with .bam: ${row.BAM}")
                }
                if (!bai.name.endsWith('.bai') && !bai.name.endsWith('.bam.bai')) {
                    error("BAI file must end with .bai: ${row.BAI}")
                }
                return [ meta, bam, bai ]
            }
            .set { ch_bam_input }
    } else {
        ch_bam_input = Channel.empty()
    }
```

**Step 3: Add ch_bam_input to emit block**

Update the emit block:

```nextflow
    emit:
    samplesheet = ch_samplesheet
    bam_input   = ch_bam_input
```

**Step 4: Commit**

```bash
git add nextflow.config subworkflows/local/utils_nfcore_tempo_pipeline/main.nf
git commit -m "feat(input): add --bamMapping parameter and BAM input channel parsing"
```

---

### Task 5: Wire BAM input into workflow — skip alignment

**Files:**

- Modify: `workflows/tempo.nf` (accept BAM channel, skip alignment for BAMs)
- Modify: `main.nf` (pass BAM channel to TEMPO workflow)

**Step 1: Add ch_bam_input to TEMPO workflow take block**

In `workflows/tempo.nf`, after `ch_pon_tbi` in the `take:` section (line ~133), add:

```nextflow
    ch_bam_input
```

**Step 2: Add BAM input path that skips alignment**

After the existing `ch_recal_bam_bai` creation (line ~324, after SAMTOOLS_INDEX_RECAL join), add the BAM input merge logic:

```nextflow
    // ===============
    // BAM INPUT (skip alignment for pre-aligned BAMs)
    // ===============
    // BAM inputs bypass FASTQC→FASTP→BWAMEM2→SORT→MERGE→MARKDUP→BQSR
    // and feed directly into the recalibrated BAM channel
    ch_bam_input
        .map { meta, bam, bai -> [ meta, bam, bai ] }
        .set { ch_bam_direct }

    // Combine FASTQ-path recalibrated BAMs with direct BAM inputs
    ch_recal_bam_bai
        .mix(ch_bam_direct)
        .set { ch_recal_bam_bai }
```

**Step 3: Update main.nf to pass BAM channel**

In `main.nf`, update the MSKCC_TEMPO workflow:

Add after `take: samplesheet`:

```nextflow
    bam_input   // channel: BAM inputs from --bamMapping
```

Update the TEMPO call to pass the BAM channel:

```nextflow
    TEMPO (
        samplesheet,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwa_index,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_known_indels,
        ch_known_indels_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_intervals,
        ch_pon,
        ch_pon_tbi,
        bam_input
    )
```

Update the entry workflow to pass BAM channel:

```nextflow
    MSKCC_TEMPO (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.bam_input
    )
```

**Step 4: Handle input validation — require at least one of --input or --bamMapping**

In `subworkflows/local/utils_nfcore_tempo_pipeline/main.nf`, update the mandatory check:

Replace:

```nextflow
    if (!params.input) {
        error("Please provide an input samplesheet with --input")
    }
```

With:

```nextflow
    if (!params.input && !params.bamMapping && !(params.aggregate instanceof String && params.aggregate != 'true')) {
        error("Please provide an input samplesheet with --input, a BAM mapping with --bamMapping, or an aggregate TSV with --aggregate")
    }
```

Also make the FASTQ samplesheet channel conditional:

```nextflow
    if (params.input) {
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { row ->
                // ... existing parsing logic ...
            }
            .set { ch_samplesheet }
    } else {
        ch_samplesheet = Channel.empty()
    }
```

**Step 5: Commit**

```bash
git add workflows/tempo.nf main.nf subworkflows/local/utils_nfcore_tempo_pipeline/main.nf
git commit -m "feat(input): wire BAM input to skip alignment and feed directly into variant calling"
```

---

### Task 6: Add aggregate-only mode via --aggregate TSV with PATH

**Files:**

- Modify: `workflows/tempo.nf` (add aggregate-only input path)

**Step 1: Add aggregate-only file resolution**

In the existing aggregation block (line ~1538), update the `aggregateIsFile` branch. When the TSV has a `PATH` column and no FASTQ/BAM inputs are provided, this is aggregate-only mode. The per-pair results are resolved from the PATH column using glob patterns matching original Tempo's directory structure.

Replace the aggregate block opening (lines 1538-1549) with:

```nextflow
    if (params.aggregate) {
        def aggregateIsFile = params.aggregate instanceof String && params.aggregate != 'true' && file(params.aggregate).exists()
        if (aggregateIsFile) {
            Channel.fromPath(params.aggregate)
                .splitCsv(sep: '\t', header: true)
                .map { row -> [ row.COHORT, row.TUMOR_ID, row.NORMAL_ID, row.PATH ?: '' ] }
                .set { ch_aggregate_raw }

            // Separate: rows with PATH (aggregate-only) vs without PATH (just cohort grouping)
            ch_aggregate_raw
                .map { cohort, tid, nid, path -> [ cohort, tid, nid ] }
                .set { ch_aggregate_map }

            // If PATH is present and no FASTQ/BAM inputs, resolve files from previous output
            ch_aggregate_raw
                .filter { cohort, tid, nid, path -> path }
                .set { ch_aggregate_with_path }
        } else {
            ch_tumor_normal_pair
                .map { meta, tbam, tbai, nbam, nbai -> [ "default_cohort", meta.tumor_id, meta.normal_id ] }
                .set { ch_aggregate_map }
            ch_aggregate_with_path = Channel.empty()
        }
```

**Step 2: Add aggregate-only file resolution for each output type**

After the `ch_aggregate_with_path` creation, add file resolution logic that mirrors the original Tempo's `inputAggregate.multiMap`. For each aggregation module, if `ch_aggregate_with_path` has data, resolve files from PATH; otherwise use pipeline-computed outputs as before.

This is complex — each aggregation module's input needs a fallback channel. The pattern for each is:

```nextflow
        // Example for somatic MAF:
        if (doWF_SNV && doWF_facets) {
            // From pipeline
            def ch_maf_from_pipeline = SOMATIC_FACETS_ANNOTATION.out.final_maf
                .map { meta, maf -> [ meta.tumor_id, meta.normal_id, maf ] }
            // From aggregate PATH
            def ch_maf_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def maf_files = file("${path}/somatic/${tid}__${nid}/*/*.final.maf")
                    def maf = maf_files instanceof List ? maf_files[0] : maf_files
                    [ tid, nid, maf ]
                }
            // Combine both sources
            ch_maf_from_pipeline.mix(ch_maf_from_path)
                .combine(ch_aggregate_map, by: [0,1])  // join on [tumor_id, normal_id]
                .groupTuple(by: 2)  // group by cohort
                .map { tids, nids, cohort, mafs -> [ cohort, mafs ] }
                .set { ch_agg_somatic_maf }
            ch_agg_somatic_maf.map { it[0] }.set { ch_agg_maf_cohort }
            ch_agg_somatic_maf.map { it[1] }.set { ch_agg_maf_files }
            AGGREGATE_SOMATIC_MAF ( ch_agg_maf_cohort, ch_agg_maf_files )
        }
```

**NOTE:** This pattern must be repeated for ALL aggregate modules (somatic MAF, somatic SV, FACETS, neoantigen, metadata, LOHHLA, HRDetect, SVClone, SV signatures, germline MAF, germline SV, QC BAM, QC Conpair). The full implementation should follow the glob patterns from original Tempo's `extractCohort` / `inputAggregate.multiMap` block documented in the context exploration.

**Important:** For aggregate-only mode (PATH present, no FASTQ/BAM), the workflow control flags still apply — only enabled workflows get aggregated. The aggregate-only mode must also handle the case where `ch_tumor_normal_pair` is empty (no FASTQ or BAM inputs).

**Step 3: Commit**

```bash
git add workflows/tempo.nf
git commit -m "feat(input): add aggregate-only mode with PATH column for pre-computed outputs"
```

---

### Task 7: Update test samplesheets

**Files:**

- Modify: `tests/csv/fastq_tumor_normal_pair.csv`
- Create: `tests/csv/bam_tumor_normal_pair.tsv` (test BAM mapping)
- Create: `tests/csv/aggregate_cohort.tsv` (test aggregate TSV)

**Step 1: Update FASTQ test samplesheet**

Replace `tests/csv/fastq_tumor_normal_pair.csv`:

```csv
patient,sample,status,target,lane,fastq_1,fastq_2
patient_1234,1234N,0,agilent,L001,test-data/testdata/tiny/normal/tiny_n_L001_R1_xxx.fastq.gz,test-data/testdata/tiny/normal/tiny_n_L001_R2_xxx.fastq.gz
patient_1234,1234N,0,agilent,L002,test-data/testdata/tiny/normal/tiny_n_L002_R1_xxx.fastq.gz,test-data/testdata/tiny/normal/tiny_n_L002_R2_xxx.fastq.gz
patient_1234,1234T,1,agilent,L001,test-data/testdata/tiny/tumor/tiny_t_L001_R1_xxx.fastq.gz,test-data/testdata/tiny/tumor/tiny_t_L001_R2_xxx.fastq.gz
patient_1234,1234T,1,agilent,L002,test-data/testdata/tiny/tumor/tiny_t_L002_R1_xxx.fastq.gz,test-data/testdata/tiny/tumor/tiny_t_L002_R2_xxx.fastq.gz
```

**Step 2: Create BAM test mapping**

Create `tests/csv/bam_tumor_normal_pair.tsv`:

```tsv
PATIENT	SAMPLE	STATUS	TARGET	BAM	BAI
patient_1234	1234N	0	agilent	test-data/testdata/tiny/normal/tiny_n.bam	test-data/testdata/tiny/normal/tiny_n.bam.bai
patient_1234	1234T	1	agilent	test-data/testdata/tiny/tumor/tiny_t.bam	test-data/testdata/tiny/tumor/tiny_t.bam.bai
```

**Step 3: Create aggregate test TSV**

Create `tests/csv/aggregate_cohort.tsv`:

```tsv
TUMOR_ID	NORMAL_ID	COHORT	PATH
1234T	1234N	test_cohort	/path/to/previous/tempo/output
```

**Step 4: Update test profiles**

Check `conf/test.config` and `conf/test_comprehensive.config` for any `sex` references in hardcoded samplesheet paths or inline params. Update as needed.

**Step 5: Commit**

```bash
git add tests/csv/ conf/test*.config
git commit -m "test: update test samplesheets with target column, add BAM and aggregate test inputs"
```

---

### Task 8: Update check_samplesheet.py validator

**Files:**

- Modify: `bin/check_samplesheet.py`

**Step 1: Add target validation to RowChecker**

The nf-schema plugin handles primary validation via `schema_input.json`, but `check_samplesheet.py` provides secondary validation. Update it to expect the new columns. Note: This validator is optional since nf-schema handles the primary validation, but keeping it consistent is good practice.

No `sex` references exist in the current `check_samplesheet.py` — it only validates `sample`, `fastq_1`, `fastq_2`. The schema handles `patient`, `status`, `target`. No changes needed unless we want to add target cross-validation (tumor-normal pairs must share same target). That cross-validation is better done in the workflow Nextflow code.

**Step 2: Skip this task if no changes needed**

If `check_samplesheet.py` doesn't reference `sex` or need `target` validation, skip.

**Step 3: Commit (if changes made)**

```bash
git add bin/check_samplesheet.py
git commit -m "refactor: update samplesheet validator for new column structure"
```

---

### Task 9: Run stub tests to verify all changes

**Files:** None (test execution only)

**Step 1: Run exome stub test**

```bash
NXF_JAVA_HOME=... nextflow run . -profile test,test_comprehensive --outdir /tmp/stub_test -stub
```

Expected: All processes succeed with stub outputs.

**Step 2: Run WGS stub test**

```bash
NXF_JAVA_HOME=... nextflow run . -profile test,test_comprehensive,test_wgs --outdir /tmp/stub_test_wgs -stub
```

**Step 3: Run aggregation stub test**

```bash
NXF_JAVA_HOME=... nextflow run . -profile test,test_comprehensive --outdir /tmp/stub_test_agg --aggregate true -stub
```

**Step 4: Run BAM input stub test (if test BAMs exist)**

```bash
NXF_JAVA_HOME=... nextflow run . --bamMapping tests/csv/bam_tumor_normal_pair.tsv -profile test,test_comprehensive --outdir /tmp/stub_test_bam -stub
```

Expected: Alignment processes are skipped, variant calling processes execute.

**Step 5: Fix any failures and re-run**

Iterate until all 4 test modes pass.

**Step 6: Commit any fixes**

```bash
git add -A
git commit -m "fix: resolve stub test failures from input system changes"
```

---

### Task 10: Verify multi-lane handling

**Files:** None (verification only)

**Step 1: Verify the test samplesheet has >1 lane per sample**

The test CSV already has L001 and L002 for both tumor and normal. This exercises the multi-lane merge path.

**Step 2: Confirm in stub test output that SAMTOOLS_MERGE runs**

Check stub test logs for `SAMTOOLS_MERGE` being invoked. Both samples should trigger the `multiple` branch since they each have 2 lanes.

**Step 3: Verify that adding a 3rd lane row works**

Create a temporary test CSV with 3 lanes for one sample and run stub test. Confirm it merges correctly.

---

## Execution Notes

- Tasks 1-3 are safe to parallelize (schema, subworkflow, workflow — they don't depend on each other until Task 5 wires them together)
- Task 4 depends on Tasks 1-3
- Task 5 depends on Task 4
- Task 6 depends on Task 5
- Tasks 7-8 can be done in parallel
- Task 9 depends on all previous tasks
- Task 10 depends on Task 9

## Key Gotchas

1. **nf-schema `samplesheetToList` ordering**: The non-meta fields are returned in schema property order. After removing `sex` and adding `target`, the order is: `fastq_1` (index 1), `fastq_2` (index 2). Verify the indices in the `.map` closure.

2. **Channel reassignment**: In Nextflow DSL2, `ch_recal_bam_bai` cannot be reassigned after `.set{}`. The BAM mix must use `.mix()` before `.set{}` or use a different channel name. Use `ch_all_recal_bam_bai` if needed.

3. **Empty channels**: When `--input` is not provided (BAM-only or aggregate-only mode), `ch_samplesheet` must be `Channel.empty()` not null. All downstream processes that depend on FASTQ alignment outputs need to handle empty channels gracefully.

4. **Aggregate-only PATH resolution**: The `file()` function with glob patterns may return empty lists if the PATH doesn't exist or files aren't found. Add error handling for missing files.
