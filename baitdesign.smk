import os

# Load the configuration file
configfile: "config.yaml"

# Extract parameters from the config
puumala_segments = config["puumala_segments"]
orthohanta_segments = config["orthohanta_segments"]
blacklist_files = config["blacklist_files"]
probe_design = config["probe_design"]
analysis_params = config["analysis"]
output_dirs = config["output_dirs"]
max_cores = config["max_cores"]
conda_envs = config["conda_envs"]
logs_dir = config["logs_dir"]

# Parameters for optional behavior and filtering
design_puumala_ends = probe_design.get("design_puumala_ends", True)
use_lsh_minhash = probe_design.get("use_lsh_minhash", False)
lsh_minhash_threshold = probe_design.get("lsh_minhash_threshold", 0.6)

segments = ["L", "M", "S"]

os.makedirs(logs_dir, exist_ok=True)

rule all:
    input:
        os.path.join(output_dirs["probes"], "final_probe_set.fasta"),
        os.path.join(output_dirs["analysis"], "coverage_analysis.txt")

#####################################
# Helper Functions
#####################################

def design_command_params(pl, ps, m, l, output, input_file):
    cmd = [
        "design.py",
        "-pl", str(pl),
        "-ps", str(ps),
        "-m", str(m),
        "-l", str(l),
        "--skip-set-cover"
    ]
    # Include blacklist files
    for blk in blacklist_files:
        cmd.extend(["--avoid-genomes", blk])
    # Add LSH MinHash if enabled
    if use_lsh_minhash:
        cmd.extend(["--filter-with-lsh-minhash", str(lsh_minhash_threshold)])
    # Add output and input
    cmd.extend(["-o", output, input_file])
    return " ".join(cmd)

#####################################
# Step 1: Design
#####################################

rule extract_puumala_segments:
    input:
        genome=lambda w: puumala_segments[w.segment]
    output:
        ends=os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"),
        middles=os.path.join(output_dirs["segments"], "puumala_{segment}_middles.fasta")
    params:
        end_size=probe_design["puumala"]["segment_end_size"],
        overlap=probe_design["puumala"]["overlap"],
        start_middle=lambda w: probe_design["puumala"]["segment_end_size"] + 1 - probe_design["puumala"]["overlap"],
        end_middle=lambda w: -(probe_design["puumala"]["segment_end_size"] + 1 - probe_design["puumala"]["overlap"]),
        output_dir=output_dirs["segments"],
        temp_dir=lambda w: os.path.join(output_dirs["segments"], f"temp_{w.segment}")
    log:
        os.path.join(logs_dir, "extract_puumala_segments_{segment}.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        mkdir -p {params.temp_dir}
        seqkit subseq --region 1:{params.end_size} --update-faidx {input.genome} > {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta 2>> {log}
        seqkit subseq --region -{params.end_size}:-1 --update-faidx {input.genome} > {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta 2>> {log}
        cat {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta > {output.ends} 2>> {log}

        seqkit subseq --region {params.start_middle}:{params.end_middle} --update-faidx {input.genome} > {output.middles} 2>> {log}
        rm -r {params.temp_dir}
        """

rule design_probes_puumala_ends:
    input:
        ends=os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta")
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_ends_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["puumala"]["segment_end_stride"],
        m=probe_design["puumala"]["mismatches"],
        l=probe_design["puumala"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: (
            f"touch {output.probes}" if not design_puumala_ends 
            else design_command_params(
                probe_design["probe_length"],
                probe_design["puumala"]["segment_end_stride"],
                probe_design["puumala"]["mismatches"],
                probe_design["puumala"]["lcf_thres"],
                output.probes,
                input.ends
            )
        )
    log:
        os.path.join(logs_dir, "design_probes_puumala_ends_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

rule design_probes_puumala:
    input:
        genome=lambda w: puumala_segments[w.segment]
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["puumala"]["probe_stride"],
        m=probe_design["puumala"]["mismatches"],
        l=probe_design["puumala"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: design_command_params(
            probe_design["probe_length"],
            probe_design["puumala"]["probe_stride"],
            probe_design["puumala"]["mismatches"],
            probe_design["puumala"]["lcf_thres"],
            output.probes,
            input.genome
        )
    log:
        os.path.join(logs_dir, "design_probes_puumala_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

rule design_probes_orthohanta:
    input:
        genome=lambda w: orthohanta_segments[w.segment]
    output:
        probes=os.path.join(output_dirs["probes"], "orthohanta_{segment}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["orthohanta"]["probe_stride"],
        m=probe_design["orthohanta"]["mismatches"],
        l=probe_design["orthohanta"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: design_command_params(
            probe_design["probe_length"],
            probe_design["orthohanta"]["probe_stride"],
            probe_design["orthohanta"]["mismatches"],
            probe_design["orthohanta"]["lcf_thres"],
            output.probes,
            input.genome
        )
    log:
        os.path.join(logs_dir, "design_probes_orthohanta_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

#####################################
# Step 2: Combine-Filter
#####################################

def conditional_puumala_ends():
    if design_puumala_ends:
        return expand(os.path.join(output_dirs["probes"], "puumala_{segment}_ends_probes.fasta"), segment=segments)
    else:
        return []

rule combine_candidate_probes:
    input:
        puumala_ends=conditional_puumala_ends(),
        puumala=expand(os.path.join(output_dirs["probes"], "puumala_{segment}_probes.fasta"), segment=segments),
        orthohanta=expand(os.path.join(output_dirs["probes"], "orthohanta_{segment}_probes.fasta"), segment=segments)
    output:
        combined=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta")
    params:
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "combine_candidate_probes.log")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.puumala_ends} {input.puumala} {input.orthohanta} > {output.combined} 2> {log}
        """

rule filter_probes:
    input:
        candidate_probes=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta"),
        puumala_ends=expand(os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"), segment=segments),
        puumala_middles=expand(os.path.join(output_dirs["segments"], "puumala_{segment}_middles.fasta"), segment=segments),
        orthohanta=expand("{genome}", genome=orthohanta_segments.values())
    output:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        m=max(probe_design["puumala"]["mismatches"], probe_design["orthohanta"]["mismatches"]),
        l=min(probe_design["puumala"]["lcf_thres"], probe_design["orthohanta"]["lcf_thres"]),
        c=1.0,
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "filter_probes.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        design.py \
            {input.puumala_ends} \
            {input.puumala_middles} \
            {input.orthohanta} \
            -pl {params.pl} \
            -m {params.m} \
            -l {params.l} \
            -c {params.c} \
            --filter-from-fasta {input.candidate_probes} \
            -o {output.filtered_probes} > {log} 2>&1
        """

rule final_probe_set:
    input:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes.fasta")
    output:
        final=os.path.join(output_dirs["probes"], "final_probe_set.fasta")
    params:
        output_dir=output_dirs["probes"]
    shell:
        """
        mkdir -p {params.output_dir}
        cp {input.filtered_probes} {output.final}
        """

#####################################
# Step 3: Validate
#####################################

rule analyze_probe_coverage:
    input:
        probes=os.path.join(output_dirs["probes"], "final_probe_set.fasta"),
        puumala_genomes=expand("{genome}", genome=puumala_segments.values()),
        orthohanta_genomes=expand("{genome}", genome=orthohanta_segments.values()),
        puumala_flanks=expand(os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"), segment=segments)
    output:
        coverage=os.path.join(output_dirs["analysis"], "coverage_analysis.txt")
    params:
        m_puumala=analysis_params["mismatches_puumala"],
        l_puumala=analysis_params["lcf_thres_puumala"],
        m_orthohanta=analysis_params["mismatches_orthohanta"],
        l_orthohanta=analysis_params["lcf_thres_orthohanta"],
        output_dir=output_dirs["analysis"]
    log:
        os.path.join(logs_dir, "analyze_probe_coverage.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        echo "=== Coverage Analysis Report - $(date) ===" >> {log}
        
        echo "Analyzing coverage for Puumala virus full genomes..." >> {log}
        analyze_probe_coverage.py \
            -d {input.puumala_genomes} \
            -f {input.probes} \
            -m {params.m_puumala} \
            -l {params.l_puumala} \
            --print-analysis > {params.output_dir}/puumala_full_coverage.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus full genomes." >> {log}
        
        echo "\nAnalyzing coverage for Puumala virus flanks..." >> {log}
        analyze_probe_coverage.py \
            -d {input.puumala_flanks} \
            -f {input.probes} \
            -m {params.m_puumala} \
            -l {params.l_puumala} \
            --print-analysis > {params.output_dir}/puumala_flanks_coverage.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus flanks." >> {log}
        
        echo "\nAnalyzing coverage for other Orthohantaviruses..." >> {log}
        analyze_probe_coverage.py \
            -d {input.orthohanta_genomes} \
            -f {input.probes} \
            -m {params.m_orthohanta} \
            -l {params.l_orthohanta} \
            --print-analysis > {params.output_dir}/orthohanta_coverage.txt 2>> {log}
        echo "Completed coverage analysis for other Orthohantaviruses." >> {log}
        
        echo "\nCombining all coverage reports into a single file..." >> {log}
        cat {params.output_dir}/puumala_full_coverage.txt \
            {params.output_dir}/puumala_flanks_coverage.txt \
            {params.output_dir}/orthohanta_coverage.txt > {output.coverage} 2>> {log}
        echo "Combined coverage analysis report created at {output.coverage}" >> {log}
        echo "=== Coverage Analysis Completed ===" >> {log}
        """