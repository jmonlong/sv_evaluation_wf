## Test dataset

Not supposed to be realistic! It's just to make sure the pipeline/tools run.

See `test.inputs.json` for all the simulated input files.

### Simulation

The dataset was prepared with the following commands.

#### Simulate a small chromosome and some calls/truth SVs.

```sh
python3 sim-test-dataset.py
```

#### Sort them

```sh
bcftools sort calls.vcf | bgzip > calls.sorted.vcf.gz
bcftools sort truth.vcf | bgzip > truth.sorted.vcf.gz
```

