import csv
import os
import sys
import math


# Define directories
# control_dir = 'control_files'
# treatment_dir = 'treatment_files'
try: 
    control_dir = sys.argv[1]
    treatment_dir = sys.argv[2]
except: 
    print("Please provide the control and treatment directories as arguments.")
    sys.exit(1)

# Function to read and normalize gene expression data from CSV files
def read_and_normalize_data(directory):
    data = {}
    
    # First pass to calculate the total expression for all genes in each file
    for file in os.listdir(directory):
        if file.endswith('.csv'):
            print(f"Processing file: {file} ({directory})")  # Debugging statement
            
            # Variable to track total expression for all genes in the current file
            total_expression_in_file = 0
            file_data = []  # Temporary storage for gene data before normalizing
            
            with open(os.path.join(directory, file), mode='r') as f:
                reader = csv.reader(f)
                next(reader)  # Skip header row

                for row in reader:
                    if not row:  # Skip empty rows
                        continue
                    gene_id = row[0]  # Gene ID is in the first column
                    expression = float(row[1])  # Expression level is in the second column
                    
                    # Add to the total expression across all genes in the file
                    total_expression_in_file += expression
                    file_data.append((gene_id, expression))
                    
            # Now normalize each gene expression using the total expression from this file
            for gene_id, expression in file_data:
                normalized_expression = expression / total_expression_in_file
                
                # Add to data dictionary, aggregating by gene_id across files
                if gene_id not in data:
                    data[gene_id] = {'sum': 0, 'expressions': []}
                
                data[gene_id]['sum'] += normalized_expression
                data[gene_id]['expressions'].append(normalized_expression)

    # After processing all files, calculate the mean and median normalized expressions
    normalized_data = {}
    for gene_id, info in data.items():
        # Calculate mean
        mean = sum(info['expressions']) / len(info['expressions'])
        
        # Calculate median
        sorted_expressions = sorted(info['expressions'])
        n = len(sorted_expressions)
        median = (sorted_expressions[n // 2] if n % 2 != 0 
                  else (sorted_expressions[n // 2 - 1] + sorted_expressions[n // 2]) / 2)
        
        normalized_data[gene_id] = {
            'mean': mean,
            'median': median
        }
    
    return normalized_data

# Function to check calculated values
def validate_output(file_path, control_data, treatment_data):
    with open(file_path, mode='r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  # Skip header
        
        for row in reader:
            gene_id = row[0]
            mean_control = float(row[1])
            median_control = float(row[2])
            mean_treatment = float(row[3])
            median_treatment = float(row[4])
            log_fold_change = float(row[5])

            # Perform validations
            if gene_id in control_data and gene_id in treatment_data:
                # Validate means and medians
                expected_mean_control = control_data[gene_id]['mean']
                expected_mean_treatment = treatment_data[gene_id]['mean']
                
                if not (math.isclose(mean_control, expected_mean_control) and 
                        math.isclose(mean_treatment, expected_mean_treatment)):
                    print(f"Mismatch in means for {gene_id}: expected control {expected_mean_control}, got {mean_control}.")
                    print(f"Expected treatment {expected_mean_treatment}, got {mean_treatment}.")

                # Validate log fold change
                expected_log_fold_change = math.log2(expected_mean_treatment / expected_mean_control) if expected_mean_control > 0 else float('inf')
                if not math.isclose(log_fold_change, expected_log_fold_change):
                    print(f"Mismatch in log fold change for {gene_id}: expected {expected_log_fold_change}, got {log_fold_change}.")
            else:
                print(f"Gene ID {gene_id} not found in control or treatment data.")

# insert main statement
def main():
    # Read and normalize control and treatment data
    control_data = read_and_normalize_data(control_dir)
    treatment_data = read_and_normalize_data(treatment_dir)

    # Prepare the summary results
    summary = []

    # Calculate log2 Fold Change
    for gene_id in control_data.keys():
        if gene_id in treatment_data:
            mean_control = control_data[gene_id]['mean']
            median_control = control_data[gene_id]['median']
            mean_treatment = treatment_data[gene_id]['mean']
            median_treatment = treatment_data[gene_id]['median']
            
            # Log2 Fold Change calculation
            # log_fold_change = math.log2(mean_treatment / mean_control) if mean_control > 0 else float('inf')

            # Handle edge cases for log2 Fold Change calculation (0, inf, -inf) 
            if mean_treatment == 0 and mean_control == 0:
                log_fold_change = 0  # No change when both means are 0
            elif mean_control == 0:
                log_fold_change = float('inf')  # Infinite up-regulation if control is 0 but treatment is not
            elif mean_treatment == 0:
                log_fold_change = float('-inf')  # Infinite down-regulation if treatment is 0 but control is not
            else:
                log_fold_change = math.log2(mean_treatment / mean_control)  # Regular log2 fold change

            summary.append((gene_id, mean_control, median_control, mean_treatment, median_treatment, log_fold_change))
        else:
            print(f"Warning: {gene_id} not found in treatment data")  # Debugging statement

    # Sort the summary by log2 Fold Change (lowest first)
    summary.sort(key=lambda x: x[5])

    # Write the output to a tab-delimited file
    with open('output_Tier1.txt', mode='w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write header
        writer.writerow(['#gene', 
                        '(mean normalized control expression)', 
                        '(median normalized control expression)', 
                        '(mean normalized treatment expression)', 
                        '(median normalized treatment expression)', 
                        '(logFoldChange)'])
        
        # Write data
        for row in summary:
            writer.writerow(row)

    # Final check on the summary data
    print(f"Total genes processed: {len(summary)}")  # Debugging statement
    
    # # Load control and treatment data for validation
    # control_data = read_and_normalize_data('control_files')
    # treatment_data = read_and_normalize_data('treatment_files')
    
    # # Validate the output file
    # validate_output('gene_expression_summary.txt', control_data, treatment_data)
    
    # # Final check on the validation
    # print("Validation complete.")  # Debugging statement
    
if __name__ == '__main__':
    main()