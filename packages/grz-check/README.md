## Usage

```sh
Checks integrity of sequencing files (FASTQ, BAM).

Use --fastq-paired for paired-end FASTQ, --fastq-single for single-end FASTQ, --bam for BAM files, or --raw for only calculating checksums of any file. These flags can be used multiple times.

By default, the tool will exit immediately after the first error is found. Use --continue-on-error to check all files regardless of errors.

Usage: grz-check [OPTIONS] --output <OUTPUT>

Options:
      --show-progress <SHOW_PROGRESS>
          Flag to show progress bars during processing

          [possible values: true, false]

      --fastq-paired <FQ1_PATH> <FQ2_PATH> <FQ1_READ_LEN> <FQ2_READ_LEN>
          A paired-end FASTQ sample. Provide FQ1, FQ2, FQ1 read length, and FQ2 read length. Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check

      --fastq-single <FQ_PATH> <READ_LEN>
          A single-end FASTQ sample. Provide the file path and read length. Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check

      --bam <BAM_PATH>
          A single BAM file to validate

      --raw <FILE_PATH>
          A file for which to only calculate the SHA256 checksum, skipping all other validation

      --output <OUTPUT>
          Path to write the output JSONL report

      --continue-on-error
          Continue processing all files even if an error is found

      --threads <THREADS>
          Number of threads to use for processing

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Example

```sh
# --fastq-paired R1 R2 read_length_R1 read_length_R2
# --fastq-single R1 read_length_R1
grz-check --show-progress --output report.jsonl --fastq-paired path/to/sample__R1.fastq.gz path/to/sample_R2.fastq.gz 150 150 --fastq-single path/to/sample.fastq.gz 151
```
