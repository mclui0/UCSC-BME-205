import numpy as np
import sys
from collections import defaultdict

def load_bed(file_path):
    ranges = defaultdict(list)
    with open(file_path, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end = line.strip().split()[:3]
            ranges[chrom].append((int(start), int(end)))
    return ranges

def load_fai(file_path):
    chrom_sizes = {}
    with open(file_path, 'r') as fai_file:
        for line in fai_file:
            chrom, size = line.strip().split()[:2]
            chrom_sizes[chrom] = int(size)
    return chrom_sizes

def merge_ranges(ranges):
    merged = defaultdict(list)
    for chrom, chrom_ranges in ranges.items():
        sorted_ranges = sorted(chrom_ranges)
        current_start, current_end = sorted_ranges[0]

        for start, end in sorted_ranges[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged[chrom].append((current_start, current_end))
                current_start, current_end = start, end
        merged[chrom].append((current_start, current_end))
    return merged

def count_overlapping_bases(setA, setB):
    overlap_count = 0
    for chrom in setA:
        if chrom in setB:
            setA_ranges = sorted(setA[chrom])
            setB_ranges = sorted(setB[chrom])
            a_idx, b_idx = 0, 0
            while a_idx < len(setA_ranges) and b_idx < len(setB_ranges):
                a_start, a_end = setA_ranges[a_idx]
                b_start, b_end = setB_ranges[b_idx]
                if a_end <= b_start:
                    a_idx += 1
                elif b_end <= a_start:
                    b_idx += 1
                else:
                    overlap_start = max(a_start, b_start)
                    overlap_end = min(a_end, b_end)
                    overlap_count += max(0, overlap_end - overlap_start)
                    if a_end < b_end:
                        a_idx += 1
                    else:
                        b_idx += 1
    return overlap_count

def randomize_bed(set_ranges, chrom_sizes):
    randomized = defaultdict(list)
    for chrom, ranges in set_ranges.items():
        max_pos = chrom_sizes.get(chrom, 0)
        for start, end in ranges:
            new_start = np.random.randint(0, max_pos - (end - start) + 1)
            randomized[chrom].append((new_start, new_start + (end - start)))
    return randomized

def permutation_test(setA, setB, chrom_sizes, num_permutations=10000):
    merged_setA = merge_ranges(setA)
    merged_setB = merge_ranges(setB)
    observed_overlap = count_overlapping_bases(merged_setA, merged_setB)

    random_overlaps = np.zeros(num_permutations, dtype=int)
    
    for i in range(num_permutations):
        randomized_setA = randomize_bed(merged_setA, chrom_sizes)
        random_overlaps[i] = count_overlapping_bases(randomized_setA, merged_setB)
    
    p_value = (np.sum(random_overlaps >= observed_overlap) + 1) / (num_permutations + 1)
    
    return observed_overlap, p_value

def main():
    if len(sys.argv) < 4:
        print("Usage: python submission.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]")
        sys.exit(1)

    setA_file = sys.argv[1]
    setB_file = sys.argv[2]
    fai_file = sys.argv[3]
    num_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000

    setA = load_bed(setA_file)
    setB = load_bed(setB_file)
    chrom_sizes = load_fai(fai_file)

    observed_overlap, p_value = permutation_test(setA, setB, chrom_sizes, num_permutations)

    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

def self_test():
    # Small test set for validation
    setA = {"chr1": [(10, 15), (13, 18), (20, 25)]}
    setB = {"chr1": [(16, 22)]}
    chrom_sizes = {"chr1": 1000}

    observed_overlap, p_value = permutation_test(setA, setB, chrom_sizes, num_permutations=100)
    print(f"Test - Observed overlap: {observed_overlap}, P-value: {p_value:.4f}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'self_test':
        self_test()
    else:
        main()
