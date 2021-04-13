import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")
adapter=config["ref"]["adapter"]
#print(units)


##### target rules #####

rule all:
    input:
        #expand("logs/trimmomatic/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r1.html", unit=units.itertuples()),
        expand("qc/fastqc_pretrim/{unit.sample}-{unit.unit}_r2.html", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r1.html", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/{unit.sample}-{unit.unit}_r2.html", unit=units.itertuples()),
        "qc/multiqc_report_pretrim.html",
        "qc/multiqc_report_posttrim.html",
        directory(config["ref"]["index"]),
        #expand("star/{unit.sample}-{unit.unit}/Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        #expand("results/{unit.sample}-{unit.unit}.featureCounts", unit=units.itertuples())
        "qc/multiqc_report_all.html",
        "results/featureCounts/all.featureCounts",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples())

##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/fpkm.smk"
