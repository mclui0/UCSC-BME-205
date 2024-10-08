# Gene Expression Normalization and Log2 Fold Change Calculation Script

## Overview
This script processes gene expression data from control and treatment CSV files, normalizes the gene expression values for each gene, and computes the log2 fold change to determine differential expression. It supports large-scale analysis across multiple CSV files, normalizing gene expression by the total expression for all genes in each file.

## Workflow
1. Input and Directory Handling
    - The script accepts two command-line arguments specifying the directories containing the control and treatment gene expression CSV files.
    - It reads and processes each .csv file from the specified directories.

2. Gene Expression Data Normalization
    - The script reads all the files in each directory and extracts gene expression data. For each file, it calculates the total expression for all genes in that file.
    - It then normalizes each gene's expression by dividing it by the total expression of all genes within the same file.
    - The mean and median normalized expression values for each gene are computed across all files in the directory (control or treatment).

3. Log2 Fold Change Calculation
    - After processing the control and treatment data, the script calculates the log2 fold change for each gene. A positive log2 fold change indicates an increase in gene expression in the treatment group compared to the control.
    - Special handling is included for cases where the control or treatment group has a mean expression of zero to avoid division errors:
        - If both the control and treatment means are zero, the fold change is zero.
        - If only the control mean is zero, the fold change is set to infinity (indicating up-regulation).
        - If only the treatment mean is zero, the fold change is set to negative infinity (indicating down-regulation).

4. Output and Sorting
    - The final results are sorted by the log2 fold change (from the lowest to highest value).
    - The output is written to a tab-delimited file (gene_expression_summary.txt), which includes:
        - Gene ID
        - Mean and median normalized expression values for both control and treatment groups
        - The log2 fold change

5. Validation (Optional)
    - The script includes an optional validation function to compare the calculated mean, median, and log2 fold change with the values in the output file, ensuring accuracy.

## Usage
The script requires two directories as input:

    python script.py <control_directory> <treatment_directory>

The output file (gene_expression_summary.txt) will be generated in the current working directory, summarizing the gene expression results.

## Notes
The normalization is done on a per-file basis, where each gene's expression is divided by the total expression across all genes within that file.
Edge cases (such as zero expression) are handled to avoid mathematical errors in fold change calculations.