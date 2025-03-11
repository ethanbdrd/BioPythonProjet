# Genomic Analysis of COVID-19 Sequences

## Project Description
This project retrieves, analyzes, and compares genetic sequences of different coronaviruses using data from the NCBI database. It extracts genetic information, analyzes gene compositions, aligns protein sequences, and identifies similarities between human, bat, and pangolin coronaviruses.

## Features
- Fetches genomic sequences from the NCBI database.
- Extracts and analyzes gene-related information.
- Identifies and isolates the Spike (S), Membrane (M), or Nucleocapsid (N) protein sequences.
- Performs multiple sequence alignments using MAFFT.
- Detects sequence variations between species.
- Computes sequence conservation rates between human and other species.
- Provides a custom sequence alignment algorithm.

## Requirements
To run this project, you need:
- Python 3.x
- Biopython (`pip install biopython`)
- MAFFT (for sequence alignment)

## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo-name.git
   cd your-repo-name
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
Run the main script to analyze genomic sequences:
   ```bash
   python main.py
   ```

## License
This project is licensed under the MIT License.
