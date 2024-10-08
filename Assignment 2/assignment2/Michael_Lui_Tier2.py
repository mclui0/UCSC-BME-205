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
            'median': median,
            'expressions': info['expressions']  # Store raw normalized expressions for statistical tests
        }
    
    return normalized_data

# Mann-Whitney U Test
def mann_whitney_u_test(control_expressions, treatment_expressions):
    combined = sorted([(expr, 'control') for expr in control_expressions] + 
                      [(expr, 'treatment') for expr in treatment_expressions], key=lambda x: x[0])
    
    rank_sum_control = sum(rank + 1 for rank, (expr, label) in enumerate(combined) if label == 'control')
    n1, n2 = len(control_expressions), len(treatment_expressions)
    u_stat = rank_sum_control - (n1 * (n1 + 1)) / 2  # U-statistic for control group
    u_prime = n1 * n2 - u_stat  # U' (alternative U statistic)
    u = min(u_stat, u_prime)
    
    mean_u = (n1 * n2) / 2
    std_u = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    
    z = (u - mean_u) / std_u  # Z-score
    # Approximate two-tailed p-value for large sample sizes
    p_value = 2 * (1 - (0.5 * (1 + math.erf(abs(z) / math.sqrt(2)))))
    
    return p_value

# Main function to calculate fold changes, apply Mann-Whitney U test, and output results
def main():
    # Read and normalize control and treatment data
    control_data = read_and_normalize_data(control_dir)
    treatment_data = read_and_normalize_data(treatment_dir)

    # Prepare the summary results
    summary = []

    # Calculate log2 Fold Change and p-value
    for gene_id in control_data.keys():
        if gene_id in treatment_data:
            mean_control = control_data[gene_id]['mean']
            median_control = control_data[gene_id]['median']
            mean_treatment = treatment_data[gene_id]['mean']
            median_treatment = treatment_data[gene_id]['median']
            
            # Log2 Fold Change calculation
            if mean_treatment == 0 and mean_control == 0:
                log_fold_change = 0  # No change when both means are 0
            elif mean_control == 0:
                log_fold_change = float('inf')  # Infinite up-regulation if control is 0 but treatment is not
            elif mean_treatment == 0:
                log_fold_change = float('-inf')  # Infinite down-regulation if treatment is 0 but control is not
            else:
                log_fold_change = math.log2(mean_treatment / mean_control)  # Regular log2 fold change

            # Perform Mann-Whitney U test
            p_value = mann_whitney_u_test(control_data[gene_id]['expressions'], treatment_data[gene_id]['expressions'])

            # Append result for this gene
            summary.append((gene_id, mean_control, median_control, mean_treatment, median_treatment, log_fold_change, p_value))
        else:
            print(f"Warning: {gene_id} not found in treatment data")  # Debugging statement

    # Sort the summary by p-value (smallest to largest)
    summary.sort(key=lambda x: x[6])

    # Write the output to a tab-delimited file
    with open('output_Tier2.txt', mode='w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Write header
        writer.writerow(['#gene', 
                        '(mean normalized control expression)', 
                        '(median normalized control expression)', 
                        '(mean normalized treatment expression)', 
                        '(median normalized treatment expression)', 
                        '(logFoldChange)', 
                        '(p-value)'])
        
        # Write data
        for row in summary:
            writer.writerow(row)

    # Final check on the summary data
    print(f"Total genes processed: {len(summary)}")  # Debugging statement

if __name__ == '__main__':
    main()
