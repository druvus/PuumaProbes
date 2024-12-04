# Snakefile

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

# Create lists of segments
segments = ["L", "M", "S"]

# Ensure log directory exists
os.makedirs(logs_dir, exist_ok=True)

# Define the final outputs of the pipeline
rule all:
    input:
        os.path.join(output_dirs["probes"], "final_probe_set.fasta"),
        os.path.join(output_dirs["analysis"], "coverage_analysis.txt")

# Rule to extract segment ends for Puumala virus
rule extract_puumala_segment_ends:
    input:
        genome=lambda wildcards: puumala_segments[wildcards.segment]
    output:
        ends=os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta")
    params:
        end_size=probe_design["puumala"]["segment_end_size"],
        output_dir=output_dirs["segments"],
        temp_dir=os.path.join(output_dirs["segments"], "temp_{segment}")
    log:
        os.path.join(logs_dir, "extract_puumala_segment_ends_{segment}.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        mkdir -p {params.temp_dir}
        # Extract first {params.end_size} bp
        seqkit subseq --region 1:{params.end_size} --update-faidx {input.genome} > {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta 2>> {log}
        
        # Extract last {params.end_size} bp
        seqkit subseq --region -{params.end_size}:-1 --update-faidx {input.genome} > {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta 2>> {log}
        
        # Combine the ends
        cat {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta > {output.ends} 2>> {log}
        
        # Clean up temporary files
        rm -r {params.temp_dir}
        """

# Rule to design probes for Puumala segment ends with higher density
rule design_probes_puumala_ends:
    input:
        ends=os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"),
        blacklist=blacklist_files  # Accept multiple blacklist files
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_ends_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["puumala"]["segment_end_stride"],
        m=probe_design["puumala"]["mismatches"],
        l=probe_design["puumala"]["lcf_thres"],
        c=probe_design["puumala"]["coverage"],
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "design_probes_puumala_ends_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        design.py \
            -pl {params.pl} \
            -ps {params.ps} \
            -m {params.m} \
            -l {params.l} \
            -c {params.c} \
            --avoid-genomes {input.blacklist} \
            -o {output.probes} \
            {input.ends} > {log} 2>&1
        """

# Rule to design probes for full Puumala segments
rule design_probes_puumala_full:
    input:
        genome=lambda wildcards: puumala_segments[wildcards.segment],
        blacklist=blacklist_files  # Accept multiple blacklist files
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_full_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["puumala"]["probe_stride"],
        m=probe_design["puumala"]["mismatches"],
        l=probe_design["puumala"]["lcf_thres"],
        c=probe_design["puumala"]["coverage"],
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "design_probes_puumala_full_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        design.py \
            -pl {params.pl} \
            -ps {params.ps} \
            -m {params.m} \
            -l {params.l} \
            -c {params.c} \
            --avoid-genomes {input.blacklist} \
            -o {output.probes} \
            {input.genome} > {log} 2>&1
        """

# Rule to design probes for other orthohantavirus segments
rule design_probes_orthohanta:
    input:
        genome=lambda wildcards: orthohanta_segments[wildcards.segment],
        blacklist=blacklist_files  # Accept multiple blacklist files
    output:
        probes=os.path.join(output_dirs["probes"], "orthohanta_{segment}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        ps=probe_design["orthohanta"]["probe_stride"],
        m=probe_design["orthohanta"]["mismatches"],
        l=probe_design["orthohanta"]["lcf_thres"],
        c=probe_design["orthohanta"]["coverage"],
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "design_probes_orthohanta_{segment}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        design.py \
            -pl {params.pl} \
            -ps {params.ps} \
            -m {params.m} \
            -l {params.l} \
            -c {params.c} \
            --avoid-genomes {input.blacklist} \
            -o {output.probes} \
            {input.genome} > {log} 2>&1
        """

# Rule to combine all probes into a single FASTA file
rule combine_probes:
    input:
        puumala_ends=expand(os.path.join(output_dirs["probes"], "puumala_{segment}_ends_probes.fasta"), segment=segments),
        puumala_full=expand(os.path.join(output_dirs["probes"], "puumala_{segment}_full_probes.fasta"), segment=segments),
        orthohanta=expand(os.path.join(output_dirs["probes"], "orthohanta_{segment}_probes.fasta"), segment=segments)
    output:
        combined=os.path.join(output_dirs["probes"], "combined_probes.fasta")
    params:
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "combine_probes.log")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.puumala_ends} {input.puumala_full} {input.orthohanta} > {output.combined} 2> {log}
        """

# Rule to deduplicate probes and log counts before and after deduplication
rule deduplicate_probes:
    input:
        combined=os.path.join(output_dirs["probes"], "combined_probes.fasta")
    output:
        deduplicated=os.path.join(output_dirs["probes"], "deduplicated_probes.fasta")
    params:
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "deduplicate_probes.log")
    conda:
        os.path.join(conda_envs, "seqkit.yaml")
    shell:
        """
        # Ensure the output directory exists
        mkdir -p {params.output_dir}

        # Start Deduplication Report
        echo "=== Deduplication Report ===" >> {log}

        # Count probes before deduplication
        echo "Number of probes before deduplication:" >> {log}
        BEFORE=$(seqkit stats {input.combined} | awk 'NR==2 {{print $4}}')
        echo $BEFORE >> {log}

        # Total length before deduplication
        TOTAL_BEFORE=$(seqkit stats {input.combined} | awk 'NR==2 {{print $5}}')
        echo "Total length of probes before deduplication: $TOTAL_BEFORE" >> {log}

        # Perform deduplication using seqkit
        seqkit rmdup {input.combined} -o {output.deduplicated} >> {log} 2>&1

        # Count probes after deduplication
        echo "Number of probes after deduplication:" >> {log}
        AFTER=$(seqkit stats {output.deduplicated} | awk 'NR==2 {{print $4}}')
        echo $AFTER >> {log}

        # Total length after deduplication
        TOTAL_AFTER=$(seqkit stats {output.deduplicated} | awk 'NR==2 {{print $5}}')
        echo "Total length of probes after deduplication: $TOTAL_AFTER" >> {log}

        # Calculate duplicates removed
        DUPLICATES_REMOVED=$((BEFORE - AFTER))
        PERCENT_REMOVED=$(echo "scale=2; ($DUPLICATES_REMOVED/$BEFORE)*100" | bc)
        echo "Number of duplicate probes removed: $DUPLICATES_REMOVED" >> {log}
        echo "Percentage of duplicate probes removed: $PERCENT_REMOVED%" >> {log}

        # End Deduplication Report
        echo "============================" >> {log}
        """

# Rule to finalize the probe set
rule final_probe_set:
    input:
        deduplicated=os.path.join(output_dirs["probes"], "deduplicated_probes.fasta")
    output:
        final=os.path.join(output_dirs["probes"], "final_probe_set.fasta")
    params:
        output_dir=output_dirs["probes"]
    shell:
        """
        mkdir -p {params.output_dir}
        cp {input.deduplicated} {output.final}
        """

# Rule to analyze probe coverage, including coverage for Puumala flanks
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
        # Ensure the output directory exists
        mkdir -p {params.output_dir}
        
        # Start Coverage Analysis Report with timestamp
        echo "=== Coverage Analysis Report - $(date) ===" >> {log}
        
        # Analyze coverage for Puumala virus full genomes
        echo "Analyzing coverage for Puumala virus full genomes..." >> {log}
        analyze_probe_coverage.py \
            -d {input.puumala_genomes} \
            -f {input.probes} \
            -m {params.m_puumala} \
            -l {params.l_puumala} \
            --print-analysis > {params.output_dir}/puumala_full_coverage.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus full genomes." >> {log}
        
        # Analyze coverage for Puumala virus flanks
        echo "\nAnalyzing coverage for Puumala virus flanks..." >> {log}
        analyze_probe_coverage.py \
            -d {input.puumala_flanks} \
            -f {input.probes} \
            -m {params.m_puumala} \
            -l {params.l_puumala} \
            --print-analysis > {params.output_dir}/puumala_flanks_coverage.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus flanks." >> {log}
        
        # Analyze coverage for other Orthohantaviruses
        echo "\nAnalyzing coverage for other Orthohantaviruses..." >> {log}
        analyze_probe_coverage.py \
            -d {input.orthohanta_genomes} \
            -f {input.probes} \
            -m {params.m_orthohanta} \
            -l {params.l_orthohanta} \
            --print-analysis > {params.output_dir}/orthohanta_coverage.txt 2>> {log}
        echo "Completed coverage analysis for other Orthohantaviruses." >> {log}
        
        # Combine all coverage reports into a single file
        echo "\nCombining all coverage reports into a single file..." >> {log}
        cat {params.output_dir}/puumala_full_coverage.txt \
            {params.output_dir}/puumala_flanks_coverage.txt \
            {params.output_dir}/orthohanta_coverage.txt > {output.coverage} 2>> {log}
        echo "Combined coverage analysis report created at {output.coverage}" >> {log}
        
        # End Coverage Analysis Report
        echo "=== Coverage Analysis Completed ===" >> {log}
        """