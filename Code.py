# Read the GFF (annotation) file
annotation_file = "annotation.gff"

# Dictionary to store the counts for each chromosome
chromosome_counts = {}

# Open the GFF file
with open(annotation_file, 'r') as file:
    # Iterate over the lines in the file
    for line in file:
        # Skip comment lines
        if line.startswith('#'):
            continue

        # Split the line by tab delimiter
        parts = line.strip().split('\t')

        # Get the chromosome and feature type
        chromosome = parts[0]
        feature_type = parts[2]

        # Count the feature type for each chromosome
        if chromosome not in chromosome_counts:
            chromosome_counts[chromosome] = {"mRNA": 0, "CDS": 0, "UTR": 0}

        if feature_type == "mRNA":
            chromosome_counts[chromosome]["mRNA"] += 1
        elif feature_type == "CDS":
            chromosome_counts[chromosome]["CDS"] += 1
        elif feature_type == "UTR":
            chromosome_counts[chromosome]["UTR"] += 1

# Print the counts for each chromosome
for chromosome, counts in chromosome_counts.items():
    print("Chromosome:", chromosome)
    print("mRNA Count:", counts["mRNA"])
    print("CDS Count:", counts["CDS"])
    print("UTR Count:", counts["UTR"])
    print()
import csv
# Export the counts as a CSV file
output_file = "counts.csv"

with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["Chromosome", "mRNA Count", "CDS Count", "UTR Count"])

    # Write the data rows
    for chromosome, counts in chromosome_counts.items():
        writer.writerow([chromosome, counts["mRNA"], counts["CDS"], counts["UTR"]])

print("Counts exported to", output_file)


import csv
from Bio import SeqIO

# Read the genome file (FASTA)
genome_file = "genome.fasta"

# Read the annotation file (GFF)
annotation_file = "annotation.gff"

# Dictionary to store the mRNA sequences
mRNA_sequences = {}

# Read the genome file
genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Open the GFF file
with open(annotation_file, 'r') as file:
    # Iterate over the lines in the file
    for line in file:
        # Skip comment lines
        if line.startswith('#'):
            continue

        # Split the line by tab delimiter
        parts = line.strip().split('\t')

        # Get the chromosome, feature type, and start/end positions
        chromosome = parts[0]
        feature_type = parts[2]
        start = int(parts[3]) - 1  # Convert to 0-based index
        end = int(parts[4])  # End position is exclusive

        # Process mRNA features
        if feature_type == "mRNA":
            # Check if the chromosome is present in the genome file
            if chromosome in genome_sequences:
                # Get the full sequence of the chromosome
                chromosome_sequence = genome_sequences[chromosome].seq

                # Extract the mRNA sequence
                mRNA_sequence = chromosome_sequence[start:end].upper()

                # Store the mRNA sequence in the dictionary
                if chromosome not in mRNA_sequences:
                    mRNA_sequences[chromosome] = []
                mRNA_sequences[chromosome].append(mRNA_sequence)

# Export the mRNA sequences as a CSV file
output_file = "mRNA_sequences.csv"

with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["Chromosome", "mRNA Sequence"])

    # Write the data rows
    for chromosome, sequences in mRNA_sequences.items():
        for i, sequence in enumerate(sequences):
            writer.writerow([chromosome, sequence])

print("mRNA sequences exported to", output_file)






# Dictionary to store the mRNA sequences
mRNA_sequences = {}

# Read the genome file
genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Open the GFF file
with open(annotation_file, 'r') as file:
    # Iterate over the lines in the file
    for line in file:
        # Skip comment lines
        if line.startswith('#'):
            continue

        # Split the line by tab delimiter
        parts = line.strip().split('\t')

        # Get the chromosome, feature type, and start/end positions
        chromosome = parts[0]
        feature_type = parts[2]
        start = int(parts[3]) - 1  # Convert to 0-based index
        end = int(parts[4])  # End position is exclusive

        # Process mRNA features
        if feature_type == "CDS":
            # Check if the chromosome is present in the genome file
            if chromosome in genome_sequences:
                # Get the full sequence of the chromosome
                chromosome_sequence = genome_sequences[chromosome].seq

                # Extract the mRNA sequence
                mRNA_sequence = chromosome_sequence[start:end].upper()

                # Store the mRNA sequence in the dictionary
                if chromosome not in mRNA_sequences:
                    mRNA_sequences[chromosome] = []
                mRNA_sequences[chromosome].append(mRNA_sequence)

# Export the mRNA sequences as a CSV file
output_file = "CDS_sequences.csv"

with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["Chromosome", "CDS Sequence"])

    # Write the data rows
    for chromosome, sequences in mRNA_sequences.items():
        for i, sequence in enumerate(sequences):
            writer.writerow([chromosome, sequence])

print("CDS sequences exported to", output_file)



