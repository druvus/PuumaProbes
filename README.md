
# PuumaProbes

**PuumaProbes** is a streamlined Snakemake pipeline designed for the efficient and accurate design of highly specific probes targeting Puumala virus segments. With an enhanced focus on segment flanks, PuumaProbes ensures comprehensive coverage and optimal performance for diagnostic and research applications.


## Introduction

PuumaProbes is developed to facilitate the design of probes for Puumala virus, a member of the hantavirus family. By leveraging the power of Snakemake and Conda for workflow management and environment control, PuumaProbes offers a reproducible and efficient solution for researchers and diagnosticians aiming to develop robust diagnostic tools.

## Features

- **Automated Probe Design:** Streamlines the process of designing probes for multiple virus segments.
- **Enhanced Focus on Segment Flanks:** Ensures higher probe density at the ends of segments for improved coverage.
- **Blacklist Incorporation:** Excludes unwanted sequences to enhance probe specificity.
- **Conda Environment Management:** Ensures reproducibility by managing dependencies through Conda.
- **Comprehensive Coverage Analysis:** Evaluates probe coverage across target genomes.
- **Detailed Logging:** Facilitates easy troubleshooting and monitoring of the pipeline.


## Requirements

- **Conda:** For environment management. [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- **Snakemake:** Workflow management system.
  ```bash
  conda install -c bioconda snakemake
  ```

## Installation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/druvus/PuumaProbes.git
   cd PuumaProbes
   ```

2. **Set Up Conda Environments:**
   
   PuumaProbes utilizes Conda environments for managing dependencies. Ensure that the `envs/` directory contains the necessary environment YAML files (`catch.yaml` and `seqkit.yaml`).

   ```bash
   # The environments will be automatically created by Snakemake.
   ```

## Configuration

PuumaProbes uses a `config.yaml` file to manage parameters and file paths. Below is an example configuration:

```yaml
# config.yaml

# Input files for Puumala virus segments
puumala_segments:
  L: "data/puumala_L.fasta"
  M: "data/puumala_M.fasta"
  S: "data/puumala_S.fasta"

# Input files for other orthohantavirus segments
orthohanta_segments:
  L: "data/orthohanta_L.fasta"
  M: "data/orthohanta_M.fasta"
  S: "data/orthohanta_S.fasta"

# Blacklist FASTA files
blacklist_files:
  - "data/blacklist_1.fasta"
  - "data/blacklist_2.fasta"
  # Add more blacklist files as needed

# Probe design parameters
probe_design:
  probe_length: 120
  puumala:
    mismatches: 2
    lcf_thres: 100
    coverage: 1.0
    probe_stride_full: 60
    probe_stride_ends: 30
    segment_end_size: 500  # Size of the segment ends (flanks) in bp
  orthohanta:
    mismatches: 5
    lcf_thres: 90
    coverage: 1.0
    probe_stride: 60

# Analysis parameters
analysis:
  mismatches_puumala: 2
  lcf_thres_puumala: 100
  mismatches_orthohanta: 5
  lcf_thres_orthohanta: 90

# Output directories
output_dirs:
  probes: "results/probes"
  segments: "results/segments"
  blacklist: "results/blacklist"
  analysis: "results/analysis"

# Conda environments directory
conda_envs: "envs"

# Logging directory
logs_dir: "logs"

# Other parameters
max_cores: 10  # Adjust based on your system
```

**Ensure that:**

- All paths under `puumala_segments` and `orthohanta_segments` correctly point to your respective FASTA files.
- Blacklist FASTA files are correctly listed under `blacklist_files`.
- Output directories are correctly specified and accessible.
- Probe design and analysis parameters match your experimental requirements.

## Usage

PuumaProbes is managed via Snakemake, which orchestrates the workflow based on the `Snakefile` and `config.yaml`.

### Running the Entire Pipeline

Execute the following command from the project root directory:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 10 --use-conda
```

- `--snakefile`: Specifies the path to the `Snakefile`.
- `--configfile`: Specifies the path to the `config.yaml`.
- `--cores`: Number of CPU cores to utilize (adjust based on your system).
- `--use-conda`: Enables Conda environment management as specified in the `Snakefile`.

### Dry Run

To perform a dry run and see the jobs that Snakemake would execute without actually running them:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 10 --use-conda -n
```

### Running Specific Rules

If you wish to run only specific rules, such as the final `analyze_probe_coverage` rule, specify its output as the target:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 10 results/analysis/coverage_analysis.txt --use-conda
```

## Example

Here's an example of designing probes for the **orthohantavirus L segment**:

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 10 --use-conda results/probes/orthohanta_L_probes.fasta
```

This command will execute all necessary rules leading up to the creation of `orthohanta_L_probes.fasta`.

## Troubleshooting

### Common Issues

1. **Conda Environment Errors:**
   - **Solution:** Ensure that Conda is properly installed and that the `envs/` directory contains the necessary YAML files. Use `--use-conda` flag when running Snakemake.

2. **Missing Input Files:**
   - **Solution:** Verify that all input FASTA files are correctly placed in the `data/` directory and that their paths are correctly specified in `config.yaml`.

3. **Permission Issues:**
   - **Solution:** Ensure you have the necessary read/write permissions for all directories and files involved in the pipeline.

4. **Tool-Specific Errors:**
   - **Solution:** Check the corresponding log files in the `logs/` directory for detailed error messages. Address issues as per the error logs.

### Checking Logs

Each rule writes logs to the `logs/` directory. For example, to view the log for designing probes for the orthohantavirus L segment:

```bash
cat logs/design_probes_orthohanta_L.log
```

## Contributing

Contributions are welcome! If you encounter bugs, have feature requests, or wish to contribute code, please follow these steps:

1. **Fork the Repository**
2. **Create a New Branch**
   ```bash
   git checkout -b feature/YourFeatureName
   ```
3. **Commit Your Changes**
   ```bash
   git commit -m "Add your message here"
   ```
4. **Push to the Branch**
   ```bash
   git push origin feature/YourFeatureName
   ```
5. **Open a Pull Request**

Please ensure that your contributions adhere to the project's coding standards and include relevant tests where applicable.

## License

This project is licensed under the [MIT License](LICENSE).



