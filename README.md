# OSC Genome Hi-C Scaffolding Pipeline

Pipeline for chromosome-scale scaffolding of the OSC genome assembly using Hi-C chromatin interaction data.

Part of the **Handler et al., 2025** publication:

**The Drosophila OSC Genome: A Resource for Studies of Transposon and piRNA Biology**

## Overview

This repository contains the workflow for scaffolding the purged OSC genome assembly into chromosome-scale scaffolds using Hi-C sequencing data. 

## Repository Structure

```
├── script-files/         # Core scaffolding scripts
└── scaffold-genome.sh    # Main submission script for Hi-C scaffolding
```

## Pipeline Components

### Hi-C Data Processing
Scripts for processing raw Hi-C sequencing data, including read mapping and contact matrix generation.

### Scaffolding
Tools for ordering and orienting contigs based on Hi-C contact frequencies to generate chromosome-scale scaffolds.

### Quality Assessment
Validation scripts to assess scaffolding accuracy and identify potential misjoins.

## Requirements

- Apptainer
  
## Usage

Run the main scaffolding pipeline using:

```bash
bash scaffold-genome.sh
```

Adjust parameters in the script files based on your Hi-C data quality and assembly characteristics.

## Output

The pipeline produces:
- Chromosome-scale scaffolded assembly
- Hi-C contact maps for quality assessment
- AGP files describing scaffold structure

## Related Resources

### Main Publication Repository
https://github.com/BrenneckeLab/Handler_2025-OSC-genome

### UCSC Genome Browser Hub
https://genome-euro.ucsc.edu/s/Brennecke%2DLab/OSC_r1.01_Handler_et.al._2025

## Citation

Please find the proper citation in https://github.com/BrenneckeLab/Handler_2025-OSC-genome

## Contact

For questions or additional information, please contact:
dominik.handler@imba.oeaw.ac.at

## License

MIT License
