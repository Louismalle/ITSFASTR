import pandas as pd
configfile: "../../../config/config.yaml" #Set config file.

Samplesheet = pd.read_csv(config["Samplesheet"], sep=' ') #Read sample sheet in a dataframe.

localrules: all_mapping

# Select output path depending on run mode.
RefPath = config["RefPath"]
metaPath = ""
print("Running in normal mode!")

if config["Trimmer"] in ["bbduk", "cutadapt"]:
    ProjDirName = config["ProjName"] + "/trimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.
    MapIn = config["TmpDir"] + "/" + config["ProjName"] + "/trimmed" + "/1_trimming/"

elif config["Trimmer"] == "none":
    ProjDirName = config["ProjName"] + "/untrimmed" + metaPath
    tmp_dir = config["TmpDir"] + "/" + ProjDirName + metaPath  # Set TEMPDIR.
    MapIn = config["SamplePath"] + "/"

else:
    print("\nMissing trimming option in config file!\n")
    exit()

rule all_mapping:
    input:
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam",
                sample = Samplesheet["sample_name"]),
        expand(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}",
                    sample = Samplesheet["sample_name"])

if config["Mode"] == "ONT":  # Choose the sequencer used for data.
    rule NanoPlot:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}")
        threads: config["ThreadNr"]
        conda: "../../envs/map_ONT_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}.log"
        shell:
            """
            NanoPlot --bam {input} \
                     -o {output} \
                     --raw \
                     --alength \
                     -t {threads} \
                     --huge \
                     2>> {log}
            """
else:
    rule qualimap:
        input:
            config["OutPath"] + "/" + ProjDirName + "/1_mapping/{sample}.bam"
        output:
            dir=directory(config["OutPath"] + "/" + ProjDirName + "/1_mapping_quality/{sample}")
        threads: config["ThreadNr"]
        resources:
                mem_mb=60000, #round(60/(math.sqrt(12))),
                time=120
        conda: "../../envs/preprocessing_env.yaml"
        benchmark:
            config["OutPath"] + "/benchmark/" + ProjDirName + "/1_mapping_quality/{sample}.tsv"
        log:
            config["OutPath"] + "/logs/" + ProjDirName + "/1_mapping_quality/{sample}.log"
        shell:
            """
            unset DISPLAY; \
            qualimap bamqc -bam {input} \
            --java-mem-size={resources.mem_mb}M \
            -nt {threads} \
            -outdir {output.dir} \
            `# -outfile {output} #> gives errors` \
            &> {log}
            """

