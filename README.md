# Structural variant evaluation workflow

Evaluate structural variant (SV) calls against a truth set using [truvari](https://github.com/ACEnglish/truvari) and [sveval](https://github.com/jmonlong/sveval).

1. [Inputs](#inputs)
2. [Outputs](#outputs)
3. [Test locally](#test-locally)
4. [Evaluate against GIAB on GRCh37](#evaluate-against-giab-on-grch37)
4. [Evaluate against GIAB on GRCh38](#evaluate-against-giab-on-grch38)

## Inputs

- `CALLS_VCF`: VCF with SV calls. Must be sorted and bgzipped.
- `TRUTH_VCF`: VCF with SV truthset. Must be sorted and bgzipped.
- `REFERENCE_FASTA`: Genome reference FASTA file
- `CONF_REGIONS_BED`: Optional. BED file defining confident regions to evaluate on
- `SIMPREP_BED`: Optional. Simple repeat annotation to help match SVs in repeats
- `SAMPLE`: Optional. Sample name. Should match between calls and truthset VCF. If empty, it will use the first sample in CALLS_VCF.

## Outputs

- `truvariSummary`: JSON summary from Truvari
- `truvariTarball`: tarball with results from Truvari
- `svevalSummary`: TSV summary from sveval for each quality threshold (first rows correspond to no threshold, i.e all variants considered)
- `svevalTarball`: tarball with results from sveval (inc. a summary PDF)

## Test locally

Quick test with dummy simulated data from [`test`](test).

```sh
miniwdl run --as-me -i test/test.inputs.json workflow.wdl
miniwdl run --as-me -i test/test.minimal.inputs.json workflow.wdl
```

## Evaluate against GIAB on GRCh37

Download the data:

```sh
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed

wget https://raw.githubusercontent.com/jmonlong/sveval/master/docs/simpleRepeat_GRCh37.bed.gz
```

For testing purpose, download a VCF with SV calls:

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/vgsv2019/vcfs/giab5-vg-HG002.vcf.gz
```

Run the workflow:

```sh
miniwdl run --as-me -i giab37.inputs.json workflow.wdl
miniwdl run --as-me -i giab37.minimal.inputs.json workflow.wdl
```

## Evaluate against GIAB on GRCh38

Disclaimer: this uses a lifted version of the GIAB truthset. 
Lifting from GRCh37 to GRCh38 seemed to be unambiguous for the vast majority of SVs (~99%), but of course a few might be lifted incorrectly or be missing.

The lifted GIAB truthset is available in the [`giab-inputs`](giab-inputs) directory. 
Download it locally with:

```sh
wget https://raw.githubusercontent.com/jmonlong/sv_evaluation_wf/master/giab-inputs/giab6_hg38-truth-baseline.vcf.gz

wget https://raw.githubusercontent.com/jmonlong/sv_evaluation_wf/master/giab-inputs/HG002_SVs_Tier1_v0.6.lifted.bed
```

Download the genome and repeat annotation:

```sh
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

wget https://raw.githubusercontent.com/jmonlong/sveval/master/docs/simpleRepeat_GRCh38.bed.gz
```

Assuming the SV calls are in a `calls.vcf.gz` file, run the workflow:

```sh
miniwdl run --as-me -i giab38.inputs.json workflow.wdl
```
