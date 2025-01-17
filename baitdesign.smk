import os
import glob

configfile: "config.yaml"

################################################
# 1) Load Configuration
################################################

puumala_segments = config["puumala_segments"]
orthohanta_segments = config["orthohanta_segments"]
blacklist_files = config["blacklist_files"]
probe_design = config["probe_design"]
analysis_params = config["analysis"]
output_dirs = config["output_dirs"]
conda_envs = config["conda_envs"]
logs_dir = config["logs_dir"]

design_puumala_ends = probe_design.get("design_puumala_ends", True)
use_lsh_minhash = probe_design.get("use_lsh_minhash", False)
lsh_minhash_threshold = probe_design.get("lsh_minhash_threshold", 0.6)

segments = ["L", "M", "S"]

puumala_mismatches = probe_design["puumala"]["mismatches"]
puumala_strides = probe_design["puumala"]["probe_stride"]
orthohanta_mismatches = probe_design["orthohanta"]["mismatches"]
orthohanta_strides = probe_design["orthohanta"]["probe_stride"]

os.makedirs(logs_dir, exist_ok=True)

################################################
# 2) Utility Functions
################################################

def design_command_params(pl, ps, m, l, output, input_file):
    """
    Constructs the command for design.py based on the given parameters.
    """
    cmd = [
        "design.py",
        "-pl", str(pl),
        "-ps", str(ps),
        "-m", str(m),
        "-l", str(l),
        "--skip-set-cover"
    ]
    for blk in blacklist_files:
        cmd.extend(["--avoid-genomes", blk])
    if use_lsh_minhash:
        cmd.extend(["--filter-with-lsh-minhash", str(lsh_minhash_threshold)])
    cmd.extend(["-o", output, input_file])
    return " ".join(cmd)

def conditional_puumala_ends_probes():
    """
    Returns expanded Puumala ends probe FASTA files if design_puumala_ends is True,
    otherwise returns an empty list.
    """
    if design_puumala_ends:
        return expand(
            os.path.join(output_dirs["probes"], "puumala_{segment}_ends_ps{ps}_m{m}_probes.fasta"),
            segment=segments,
            ps=puumala_strides,
            m=puumala_mismatches
        )
    else:
        return []

def get_unique_mismatches():
    """Get unique mismatch values from both Puumala and Orthohanta configs."""
    puumala_mismatches = set(probe_design["puumala"]["mismatches"])
    orthohanta_mismatches = set(probe_design["orthohanta"]["mismatches"])
    return sorted(list(puumala_mismatches.union(orthohanta_mismatches)))
    
################################################
# 3) Expansions for Designing Probes
################################################

all_puumala_probes = expand(
    os.path.join(output_dirs["probes"], "puumala_{segment}_ps{ps}_m{m}_probes.fasta"),
    segment=segments,
    ps=puumala_strides,
    m=puumala_mismatches
)

all_orthohanta_probes = expand(
    os.path.join(output_dirs["probes"], "orthohanta_{segment}_ps{ps}_m{m}_probes.fasta"),
    segment=segments,
    ps=orthohanta_strides,
    m=orthohanta_mismatches
)

# Collect all design FASTAs (Puumala, Orthohanta, optional Puumala ends)
all_design_probes = (
    all_puumala_probes +
    all_orthohanta_probes +
    conditional_puumala_ends_probes()
)

################################################
# 4) Identify Genome References for Mapping
################################################

genome_folder = config["genome_folder"]
genome_files = glob.glob(os.path.join(genome_folder, "*", "assembly", "*.fa.gz"))

# Build a list of dicts for host, assembly, and file path
genomes = [
    {
        "host": os.path.basename(os.path.dirname(os.path.dirname(filepath))),
        "assembly": os.path.splitext(os.path.basename(filepath))[0],
        "filepath": filepath
    }
    for filepath in genome_files
]

# Unique IDs like "Aedes_aegypti_GCA_002204515.1"
genome_ids = [
    f"{g['host']}_{g['assembly'].replace('.', '_')}" 
    for g in genomes
]


################################################
# 5) Construct expansions for mapping each design
################################################

def design_file_stem(path):
    """
    Given an absolute path like:
        results/probes/puumala_L_ps60_m2_probes.fasta
    Return:
        puumala_L_ps60_m2_probes
    """
    base = os.path.basename(path)       # puumala_L_ps60_m2_probes.fasta
    stem, _ = os.path.splitext(base)    # puumala_L_ps60_m2_probes
    return stem

# For each design FASTA, map to each genome -> a PAF
paf_files = [
    os.path.join(
        output_dirs["mapping"],
        f"{design_file_stem(design_file)}__{genome_id}.paf"  # Note the double underscore
    )
    for design_file in all_design_probes
    for genome_id in genome_ids
]

# Filtered matching oligos for each PAF
matching_oligos_txt = [
    paf.replace(".paf", "_matching_oligos.txt")
    for paf in paf_files
]

# print("Number of genome files found:", len(genome_files))
# print("Genome IDs:", genome_ids)
# print("Number of PAF files to generate:", len(paf_files))
# print("Number of matching oligos files to generate:", len(matching_oligos_txt))

# Add at the beginning of the file after the config loading
def debug_print_design_files():
    print("\nDesign files:")
    for df in all_design_probes:
        print(f"  {design_file_stem(df)}")
    
    print("\nSample PAF files to be generated:")
    for i, paf in enumerate(paf_files[:5]):  # Print first 5 for sample
        print(f"  {paf}")

#debug_print_design_files()


################################################
# 6) rule all
################################################

rule all:
    """
    The 'all' rule ensures:
      1) All design probes are generated
      2) They are mapped individually to each genome
      3) Matching oligos are filtered per design+genome
      4) Candidate probes are combined & filtered for each mismatch value
      5) Final sets of probes are produced for each mismatch value
      6) Coverage analysis is completed for each mismatch value
    """
    input:
        # 1) All design FASTAs
        all_design_probes,
        # 2 & 3) Mapping outputs (PAF + matching oligos per design)
        paf_files,
        matching_oligos_txt,
        # 4) Combined candidate probes
        os.path.join(output_dirs["probes"], "all_candidate_probes.fasta"),
        # 5) Final sets for each mismatch value
        expand(
            os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta"),
            m=get_unique_mismatches()
        ),
        # 6) Coverage analysis for each mismatch value
        expand(
            os.path.join(output_dirs["analysis"], "coverage_analysis_{m}.txt"),
            m=get_unique_mismatches()
        )

################################################
# Step 1: Extract segments
################################################

rule extract_puumala_segments:
    """
    Extracts the ends and middles from Puumala segments for optional ends-probe design.
    """
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
        seqkit subseq --region 1:{params.end_size} --update-faidx {input.genome} \
          > {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta 2>> {log}
        seqkit subseq --region -{params.end_size}:-1 --update-faidx {input.genome} \
          > {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta 2>> {log}

        cat {params.temp_dir}/puumala_{wildcards.segment}_5prime.fasta \
            {params.temp_dir}/puumala_{wildcards.segment}_3prime.fasta \
          > {output.ends} 2>> {log}

        seqkit subseq --region {params.start_middle}:{params.end_middle} --update-faidx {input.genome} \
          > {output.middles} 2>> {log}

        rm -r {params.temp_dir}
        """

################################################
# Step 2: Design Probes
################################################

rule design_probes_puumala_ends:
    """
    Designs probes specifically for the Puumala segment ends, if enabled.
    Uses same stride and mismatch parameters as other Puumala designs.
    """
    input:
        ends=os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta")
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_ends_ps{ps}_m{m}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        l=probe_design["puumala"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: (
            f"touch {output.probes}" if not design_puumala_ends
            else design_command_params(
                probe_design["probe_length"],
                int(wildcards.ps),
                int(wildcards.m),
                probe_design["puumala"]["lcf_thres"],
                output.probes,
                input.ends
            )
        )
    log:
        os.path.join(logs_dir, "design_probes_puumala_ends_{segment}_ps{ps}_m{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """


rule design_probes_puumala:
    """
    Designs probes for full Puumala segments with various stride/mismatch combos.
    """
    input:
        genome=lambda w: puumala_segments[w.segment]
    output:
        probes=os.path.join(output_dirs["probes"], "puumala_{segment}_ps{ps}_m{m}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        l=probe_design["puumala"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: design_command_params(
            probe_design["probe_length"],
            wildcards.ps,
            wildcards.m,
            probe_design["puumala"]["lcf_thres"],
            output.probes,
            input.genome
        )
    log:
        os.path.join(logs_dir, "design_probes_puumala_{segment}_ps{ps}_m{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

rule design_probes_orthohanta:
    """
    Designs probes for Orthohantavirus segments with various stride/mismatch combos.
    """
    input:
        genome=lambda w: orthohanta_segments[w.segment]
    output:
        probes=os.path.join(output_dirs["probes"], "orthohanta_{segment}_ps{ps}_m{m}_probes.fasta")
    params:
        pl=probe_design["probe_length"],
        l=probe_design["orthohanta"]["lcf_thres"],
        output_dir=output_dirs["probes"],
        command=lambda wildcards, input, output: design_command_params(
            probe_design["probe_length"],
            wildcards.ps,
            wildcards.m,
            probe_design["orthohanta"]["lcf_thres"],
            output.probes,
            input.genome
        )
    log:
        os.path.join(logs_dir, "design_probes_orthohanta_{segment}_ps{ps}_m{m}.log")
    conda:
        os.path.join(conda_envs, "catch.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        {params.command} > {log} 2>&1
        """

################################################
# Step 3: Map Each Design FASTA to Each Genome
################################################
rule map_probes_to_genome:
    """
    Maps each designed probe FASTA to each genome reference.
    """
    input:
        probes=lambda wc: next(
            path for path in all_design_probes
            if design_file_stem(path) == wc.design_file_stem
        ),
        genome=lambda wc: next(
            g["filepath"] for g in genomes
            if f"{g['host']}_{g['assembly'].replace('.', '_')}" == wc.genome_id
        )
    output:
        paf=os.path.join(output_dirs["mapping"], "{design_file_stem}__{genome_id}.paf")
    params:
        preset=config["minimap2"]["preset"],
        k=config["minimap2"]["k"],
        w=config["minimap2"]["w"],
        threads=config["minimap2"]["threads"],
        output_dir=output_dirs["mapping"]  # Add output_dir to params
    log:
        os.path.join(logs_dir, "map_probes_to_genome_{design_file_stem}__{genome_id}.log")
    conda:
        os.path.join(conda_envs, "minimap2.yaml")
    shell:
        """
        mkdir -p {params.output_dir}
        minimap2 -x {params.preset} -k {params.k} -w {params.w} -t {params.threads} \
          {input.genome} {input.probes} \
          > {output.paf} 2> {log}
        """


rule filter_paf:
    """
    Applies awk to each PAF file to select records meeting coverage/identity thresholds.
    """
    input:
        paf=os.path.join(output_dirs["mapping"], "{design_file_stem}__{genome_id}.paf")  # Note the double underscore
    output:
        matching_oligos=os.path.join(
            output_dirs["mapping"],
            "{design_file_stem}__{genome_id}_matching_oligos.txt"  # Note the double underscore
        )
    params:
        awk_condition=f"'$1 >= {config['minimap2']['min_coverage']} && $10/$11 >= {config['minimap2']['min_identity']}'"
    log:
        os.path.join(logs_dir, "filter_paf_{design_file_stem}__{genome_id}.log")  # Note the double underscore
    shell:
        """
        awk {params.awk_condition} {input.paf} | cut -f1 | sort -u \
          > {output.matching_oligos} 2> {log}
        """

################################################
# Step 4: Combine + Filter
################################################

rule combine_candidate_probes:
    """
    AFTER mapping each design individually, we combine all the design FASTAs into one file.
    """
    input:
        probes=all_design_probes,
        # Add explicit dependency on mapping results
        paf=paf_files,
        matching_oligos=matching_oligos_txt
    output:
        combined=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta")
    params:
        output_dir=output_dirs["probes"]
    log:
        os.path.join(logs_dir, "combine_candidate_probes.log")
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input.probes} > {output.combined} 2> {log}
        """



rule filter_probes:
    """
    Filters the combined candidate probes, using references (Puumala ends, middles, orthohanta).
    """
    input:
        candidate_probes=os.path.join(output_dirs["probes"], "all_candidate_probes.fasta"),
        puumala_ends=expand(
            os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"),
            segment=segments
        ),
        puumala_middles=expand(
            os.path.join(output_dirs["segments"], "puumala_{segment}_middles.fasta"),
            segment=segments
        ),
        orthohanta=expand("{genome}", genome=orthohanta_segments.values())
    output:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes_{m}.fasta")
    params:
        pl=probe_design["probe_length"],
        l=min(probe_design["puumala"]["lcf_thres"], probe_design["orthohanta"]["lcf_thres"]),
        c=1.0,
        output_dir=output_dirs["probes"]
    wildcard_constraints:
        m="|".join(map(str, get_unique_mismatches()))
    log:
        os.path.join(logs_dir, "filter_probes_{m}.log")
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
            -m {wildcards.m} \
            -l {params.l} \
            -c {params.c} \
            --filter-from-fasta {input.candidate_probes} \
            -o {output.filtered_probes} > {log} 2>&1
        """

rule final_probe_set:
    """
    Creates final probe sets for each mismatch value.
    """
    input:
        filtered_probes=os.path.join(output_dirs["probes"], "filtered_probes_{m}.fasta")
    output:
        final=os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta")
    params:
        output_dir=output_dirs["probes"]
    wildcard_constraints:
        m="|".join(map(str, get_unique_mismatches()))
    shell:
        """
        mkdir -p {params.output_dir}
        cp {input.filtered_probes} {output.final}
        """

################################################
# Step 5: Validate (Coverage Analysis)
################################################

rule analyze_probe_coverage:
    """
    Uses analyze_probe_coverage.py to perform coverage analysis for each mismatch value.
    """
    input:
        probes=os.path.join(output_dirs["probes"], "final_probe_set_{m}.fasta"),
        puumala_genomes=expand("{genome}", genome=puumala_segments.values()),
        orthohanta_genomes=expand("{genome}", genome=orthohanta_segments.values()),
        puumala_flanks=expand(
            os.path.join(output_dirs["segments"], "puumala_{segment}_ends.fasta"),
            segment=segments
        )
    output:
        coverage=os.path.join(output_dirs["analysis"], "coverage_analysis_{m}.txt")
    params:
        l_puumala=analysis_params["lcf_thres_puumala"],
        l_orthohanta=analysis_params["lcf_thres_orthohanta"],
        output_dir=output_dirs["analysis"]
    wildcard_constraints:
        m="|".join(map(str, get_unique_mismatches()))
    log:
        os.path.join(logs_dir, "analyze_probe_coverage_{m}.log")
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
            -m {wildcards.m} \
            -l {params.l_puumala} \
            --print-analysis \
          > {params.output_dir}/puumala_full_coverage_{wildcards.m}.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus full genomes." >> {log}
        
        echo "\nAnalyzing coverage for Puumala virus flanks..." >> {log}
        analyze_probe_coverage.py \
            -d {input.puumala_flanks} \
            -f {input.probes} \
            -m {wildcards.m} \
            -l {params.l_puumala} \
            --print-analysis \
          > {params.output_dir}/puumala_flanks_coverage_{wildcards.m}.txt 2>> {log}
        echo "Completed coverage analysis for Puumala virus flanks." >> {log}
        
        echo "\nAnalyzing coverage for other Orthohantaviruses..." >> {log}
        analyze_probe_coverage.py \
            -d {input.orthohanta_genomes} \
            -f {input.probes} \
            -m {wildcards.m} \
            -l {params.l_orthohanta} \
            --print-analysis \
          > {params.output_dir}/orthohanta_coverage_{wildcards.m}.txt 2>> {log}
        echo "Completed coverage analysis for other Orthohantaviruses." >> {log}
        
        echo "\nCombining all coverage reports into a single file..." >> {log}
        cat {params.output_dir}/puumala_full_coverage_{wildcards.m}.txt \
            {params.output_dir}/puumala_flanks_coverage_{wildcards.m}.txt \
            {params.output_dir}/orthohanta_coverage_{wildcards.m}.txt \
          > {output.coverage} 2>> {log}
        echo "Combined coverage analysis report created at {output.coverage}" >> {log}
        echo "=== Coverage Analysis Completed ===" >> {log}
        """
