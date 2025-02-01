# Unique Coding Sequences Project

## Introduction
This project focuses on identifying unique coding sequences from a target genome by comparing it with 10 reference genomes. The process involves annotation, sequence conversion, and local BLAST analysis to identify sequences that are unique to the target genome.

## Features
- Annotates the **target genome** and **10 reference genomes**.
- Converts all annotated genomes into **nucleotide sequences**.
- Adds reference genomes and performs **local BLAST searches**.
- Compares the sequences and identifies **unique coding sequences**.


## Repository Structure
```
Unique_Coding_Sequences_Project/
│-- data/ # Raw and processed data files
│-- scripts/ # Python scripts for analysis
│-- results/ # Output and reports
│-- README.md # Project documentation
│-- requirements.txt # Dependencies
│-- Unique_coding_sequences.ipynb # Jupyter Notebook for analysis
```

## Requirements
To run this project, install the following dependencies:
```bash
pip install biopython pandas numpy
```

## Cloning the Repository
Clone the repository from GitHub:
```bash
git clone https://github.com/Mr-Sarmad/Unique_Coding_Sequences_Project.git
cd Unique_Coding_Sequences_Project
```

## Running the Project
1. Open the Jupyter Notebook:
```bash
jupyter notebook Unique_coding_sequences.ipynb
```
2. Follow the steps in the notebook to process data and analyze unique genes.

### Steps for Identifying Unique Coding Sequences
## 1. Annotating Genomes
First, annotate the target genome to identify all coding regions.
Next, annotate 10 reference genomes to perform comparisons.
## 2. Converting Genomes into Nucleotide Sequences
After annotation, convert the identified coding sequences of the target genome and reference genomes into nucleotide sequences.
## 3. Adding Reference Genomes and Running Local BLAST
Add the reference genomes and create a local BLAST database for sequence comparison.
Run BLAST locally to align the target genome sequences against the reference genomes.
## 4. Comparing Sequences and Identifying Unique Sequences
Compare the sequences using BLAST results to identify coding sequences that are unique to the target genome.

## BLAST Setup
BLAST (Basic Local Alignment Search Tool) is used for sequence alignment.

### 1. Install BLAST
Download and install BLAST from the NCBI website:
- [NCBI BLAST+ Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- Extract the files and add the BLAST directory to your system PATH.

### 2. Creating a Local BLAST Database
To create a BLAST database from a FASTA file:
```bash
makeblastdb -in input_sequences.fasta -dbtype nucl -out my_blast_db
```
- `-in input_sequences.fasta` : Input FASTA file.
- `-dbtype nucl` : Specifies nucleotide database.
- `-out my_blast_db` : Output database name.

### 3. Running a Local BLAST Search
To run BLAST locally against the created database:
```bash
blastn -query query_sequences.fasta -db my_blast_db -out results.txt -outfmt 6
```
- `-query query_sequences.fasta` : Query sequences.
- `-db my_blast_db` : Database to search against.
- `-out results.txt` : Output file.
- `-outfmt 6` : Tabular output format.

## Contributing
If you want to contribute:
1. Fork the repository.
2. Create a new branch: `git checkout -b feature-name`
3. Commit your changes: `git commit -m 'Add new feature'`
4. Push to the branch: `git push origin feature-name`
5. Open a Pull Request.

## License
This project is licensed under the MIT License.

## Contact
For queries, contact [Sarmad Jutt](mailto:sarmadjutt136@gmail.com).

