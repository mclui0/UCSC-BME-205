### Overview
This Python script performs a permutation test to determine the statistical significance of the overlap between two sets of genomic ranges provided in BED format. The program takes two BED files (`SetA.bed` and `SetB.bed`), a FASTA index file (`genome.fa.fai`), and optionally the number of permutations to perform, and outputs the number of observed overlapping bases and the associated p-value.

### Command Line Usage
```bash
python Firstname_Lastname_TierX.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]
```
- **path/to/SetA.bed**: Path to the first set of genomic ranges in BED format.
- **path/to/SetB.bed**: Path to the second set of genomic ranges in BED format.
- **path/to/genome.fa.fai**: Path to the FASTA index file describing the chromosome sizes.
- **num_permutations** *(optional)*: Number of permutations to perform in the test (default is 10,000).

#### Example:
```bash
python Firstname_Lastname_TierX.py SetA.bed SetB.bed genome.fa.fai 10000
```

### Output
The program outputs the observed number of overlapping bases between the two sets and the associated p-value:
```
Number of overlapping bases observed: <observed_overlap>, p value: <p_value>
```

### Functions

- **`load_bed(file_path)`**:
  - Loads a BED file and returns genomic ranges as a dictionary keyed by chromosome.
  - **Input**: Path to BED file.
  - **Output**: Dictionary of genomic ranges.

- **`load_fai(file_path)`**:
  - Loads a FASTA index file (`.fai`) and returns chromosome sizes as a dictionary.
  - **Input**: Path to the `.fai` file.
  - **Output**: Dictionary of chromosome sizes.

- **`merge_ranges(ranges)`**:
  - Merges overlapping genomic ranges for each chromosome.
  - **Input**: Dictionary of genomic ranges.
  - **Output**: Dictionary of merged genomic ranges.

- **`count_overlapping_bases(setA, setB)`**:
  - Counts the number of overlapping bases between two sets of genomic ranges.
  - **Input**: Two dictionaries of genomic ranges (merged).
  - **Output**: Integer count of overlapping bases.

- **`randomize_bed(set_ranges, chrom_sizes)`**:
  - Randomly redistributes genomic ranges within chromosome boundaries.
  - **Input**: Dictionary of merged genomic ranges and chromosome sizes.
  - **Output**: Dictionary of randomized genomic ranges.

- **`permutation_test(setA, setB, chrom_sizes, num_permutations)`**:
  - Performs a permutation test to determine the significance of observed overlap.
  - **Input**: Two genomic range sets, chromosome sizes, and the number of permutations.
  - **Output**: Observed overlap and p-value.

- **`self_test()`**:
  - Runs a small self-test to validate the permutation test logic.
  - **Output**: Prints the observed overlap and p-value for test data.

### Self-Test Usage
You can run the self-test using the command:
```bash
python Firstname_Lastname_TierX.py self_test
```

### Dependencies
- **Python 3.12.4**
- **Numpy 2.1.1**

Make sure you are using the specified environment:
```bash
conda create -c conda-forge -n bme205-assn3 "python==3.12.4" "numpy==2.1.1"
```