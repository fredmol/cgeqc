# cgeqc

> **REPOSITORY MOVED**: This repository has been moved to Bitbucket. Please use the new repository at https://bitbucket.org/genomicepidemiology/cgeqc for the latest updates and releases.

A quality control and sequence trimming tool for the CGELabs suite. 
This tool provides QC assessment and trimming functionality for genomic sequencing data from raw Nanopore data, 
with specialized support for bacterial, viral, and metagenomic data analysis.

## Features

- Quality assessment of sequencing data with data type-specific thresholds
- Trimming of reads based on customizable parameters
- Generation of PDF QC reports with visualizations
- Support for bacterial, viral, and metagenomic sequencing data

## Installation

### Using Conda

The following Conda channels are required:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels genomicepidemiology
```

Mamba can be used for faster installation:

```bash
conda install -c conda-forge mamba
mamba install -c genomicepidemiology cgeqc
```

Alternatively, install directly with conda:

```bash
conda install -c genomicepidemiology cgeqc
```

### Using pip

```bash
pip install cgeqc
```

### From source

```bash
git clone https://bitbucket.org/genomicepidemiology/cgeqc.git
cd cgeqc
pip install .
```

## Dependencies

The tool requires KMA (k-mer alignment) for trimming operations. Make sure KMA is installed and available in your PATH:

```bash
conda install -c bioconda kma
```

## Usage

### Basic usage

```bash
cgeqc -i <input_fastq> -o <output_directory> -n <sample_name>
```

### Specifying pipeline type

```bash
# For bacterial data (default)
cgeqc -i <input_fastq> -o <output_directory> -n <sample_name> --pipeline bacterial

# For viral data
cgeqc -i <input_fastq> -o <output_directory> -n <sample_name> --pipeline viral

# For metagenomic data
cgeqc -i <input_fastq> -o <output_directory> -n <sample_name> --pipeline metagenomic
```

### Customizing trimming parameters

```bash
cgeqc -i <input_fastq> -o <output_directory> -n <sample_name> \
  --min_length 500 \
  --min_phred 25 \
  --min_average_quality 15
```

### Full parameter list

```
usage: cgeqc [-h] -i INPUT [-o OUTPUT] -n NAME
             [--pipeline {bacterial,viral,metagenomic}]
             [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
             [--min_phred MIN_PHRED]
             [--min_internal_phred MIN_INTERNAL_PHRED]
             [--min_average_quality MIN_AVERAGE_QUALITY]
             [--trim_5_prime TRIM_5_PRIME] [--trim_3_prime TRIM_3_PRIME]
             [--version]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input FASTQ file
  -o OUTPUT, --output OUTPUT
                        Output directory for trimmed file and QC report
  -n NAME, --name NAME  Sample name identifier
  --pipeline {bacterial,viral,metagenomic}
                        Pipeline type to use appropriate quality thresholds for reporting
  --min_length MIN_LENGTH
                        Minimum read length for trimming
  --max_length MAX_LENGTH
                        Maximum read length for trimming
  --min_phred MIN_PHRED
                        Minimum phred score for trimming
  --min_internal_phred MIN_INTERNAL_PHRED
                        Minimum internal phred score for trimming
  --min_average_quality MIN_AVERAGE_QUALITY
                        Minimum average quality for trimming
  --trim_5_prime TRIM_5_PRIME
                        Number of bases to trim from the 5' end
  --trim_3_prime TRIM_3_PRIME
                        Number of bases to trim from the 3' end
  --version             show program's version number and exit
```

## Output

The tool produces:

1. A trimmed FASTQ file: `<sample_name>.fq`
2. A JSON file with detailed QC metrics: `<sample_name>.json`
3. A comprehensive PDF QC report: `<sample_name>_qc_report.pdf`

The PDF report includes:
- Read quality distribution
- Read length distribution
- Overall quality assessment
- Data type-specific metrics and thresholds
- Recommendations based on data quality

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.


## Contact

For bug reports or feature requests, please use the issue tracker on Bitbucket:
https://bitbucket.org/genomicepidemiology/cgeqc/issues

For other inquiries, contact:
Frederik Duus Møller (freddu@food.dtu.dk)