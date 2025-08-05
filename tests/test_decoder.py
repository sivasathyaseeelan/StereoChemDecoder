import statistics
import csv
import os
from StereoChemDecoder import StereoChemDecoder

def read_smiles_from_csv(file_path):
    smiles_list = []
    with open(file_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header row ("SMILES")
        for row in reader:
            if row:  # Ensure row is not empty
                smiles_list.append(row[0])  # First column is SMILES
    return smiles_list

def main():
    # Define paths to CSV files in the tests subdirectory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    smiles_csv_path = os.path.join(base_dir, "smiles_dataset.csv")
    test_smiles_csv_path = os.path.join(base_dir, "test_smiles_dataset.csv")
    
    # Load SMILES datasets from CSV files
    smiles_list = read_smiles_from_csv(smiles_csv_path)
    test_smiles_list = read_smiles_from_csv(test_smiles_csv_path)
    
    # Ensure unique SMILES (preserving original logic, though CSVs contain unique entries)
    smiles_list = list(dict.fromkeys(smiles_list))[:1000]

    # Collect data for statistical analysis
    chiral_counts = []
    ez_counts = []
    valid_molecules = 0
    results = []
    invalid_smiles = []

    decoder = StereoChemDecoder()
    for smiles in smiles_list:
        try:
            chiral_centers = decoder.detect_chirality(smiles)
            potential_ez = decoder.detect_ez(smiles)
            # Use 'or' condition as per original code
            if chiral_centers is not None or potential_ez is not None:
                valid_molecules += 1
                chiral_count = len(chiral_centers) if chiral_centers is not None else 0
                ez_count = len(potential_ez) if potential_ez is not None else 0
                chiral_counts.append(chiral_count)
                ez_counts.append(ez_count)
                results.append((smiles, chiral_count, ez_count))
            else:
                invalid_smiles.append(smiles)
        except ValueError as e:
            print(e)
            invalid_smiles.append(smiles)

    # Statistical analysis
    print("\nStatistical Analysis of Stereochemistry for", valid_molecules, "Molecules:")
    print("--------------------------------------------------")
    print(f"Invalid SMILES strings: {len(invalid_smiles)}")
    if invalid_smiles:
        print("Examples of invalid SMILES:", invalid_smiles[:5])

    # Chiral centers statistics
    chiral_mean = statistics.mean(chiral_counts) if chiral_counts else 0
    chiral_std = statistics.stdev(chiral_counts) if len(chiral_counts) > 1 else 0
    chiral_min = min(chiral_counts) if chiral_counts else 0
    chiral_max = max(chiral_counts) if chiral_counts else 0
    chiral_nonzero = sum(1 for c in chiral_counts if c > 0)

    print("\nChiral Centers:")
    print(f"  Mean number of chiral centers: {chiral_mean:.2f}")
    print(f"  Standard deviation: {chiral_std:.2f}")
    print(f"  Range: {chiral_min} to {chiral_max}")
    print(f"  Molecules with at least one chiral center: {chiral_nonzero} ({chiral_nonzero/valid_molecules*100:.1f}%)")

    # E/Z double bonds statistics
    ez_mean = statistics.mean(ez_counts) if ez_counts else 0
    ez_std = statistics.stdev(ez_counts) if len(ez_counts) > 1 else 0
    ez_min = min(ez_counts) if ez_counts else 0
    ez_max = max(ez_counts) if ez_counts else 0
    ez_nonzero = sum(1 for e in ez_counts if e > 0)

    print("\nE/Z Double Bonds:")
    print(f"  Mean number of E/Z double bonds: {ez_mean:.2f}")
    print(f"  Standard deviation: {ez_std:.2f}")
    print(f"  Range: {ez_min} to {ez_max}")
    print(f"  Molecules with at least one E/Z double bond: {ez_nonzero} ({ez_nonzero/valid_molecules*100:.1f}%)")

    # Additional statistics
    both_chiral_and_ez = sum(1 for _, chiral_count, ez_count in results if chiral_count > 0 and ez_count > 0)
    no_chiral_no_ez = sum(1 for _, chiral_count, ez_count in results if chiral_count == 0 and ez_count == 0)

    print("\nAdditional Statistics:")
    print(f"  Molecules with both chiral centers and E/Z double bonds: {both_chiral_and_ez} ({both_chiral_and_ez/valid_molecules*100:.1f}%)")
    print(f"  Molecules with no chiral centers and no E/Z double bonds: {no_chiral_no_ez} ({no_chiral_no_ez/valid_molecules*100:.1f}%)")

    # Test specific SMILES strings
    print("\nTesting specific SMILES strings:")
    for test_smiles in test_smiles_list:
        chiral_centers = decoder.detect_chirality(test_smiles)
        ez_bonds = decoder.detect_ez(test_smiles)
        print(f"Testing {test_smiles}")
        print(f"  Chirality       : {chiral_centers}")
        print(f"  E/Z double bond : {ez_bonds}")

if __name__ == "__main__":
    main()
