# Mutation Disequilibrium Analysis

> **Warning**
> We are still completing the checklist below. Do not attempt to use this repo until the entire TODO section has been deleted!

These are controlled via the `mdeqasis` command line tool.

# TODO

- [ ] All of the data should be uploaded to Zenodo
- [ ] write a script that copies all the data and results from zenodo into the container if the data and results are not already present
- [ ] connect the script to the dockerfile via `devcontainer.json` file
- [ ] Provide a link to our other instructions on using VS Code to open a directory when you have Docker installed and instantiate a container.

## Running the commands

We provide a Dockerfile for constructing the computational environment necessary to repeat these analyses.

## Sampling of microbial 16S sequences

The following subcommands were run:

`filter-alignments` -> produces sqlitedb of 16S sequences after removing non nucleotide characters (including gaps)

`microbial-fit-gn` -> fits GN model to all microbial alignments

`microbial-gn-stats` -> extracts statistics for sampling the four seed alignments

## Sampling of the ape sequences

The following subcommands were run:

`ape-align-cds` -> aligns CDS sequences for which there is an alignment of intronic sequence. Produces sqlitedb.

`ape-match-cds-intron` -> Filters intron and CDS alignments based on number of aligned positions and sequence complexity. Produces a separate sqlitedb for CDS and introns.

## Sampling of microbial sequences

`filter-alignments` -> Converts nexus formatted aligned sequences into a sqlitedb, dropping aligned columns with gaps.

`microbial-fit-gn` -> Fits GN model to the collection of microbial 16S alignments.

`microbial-gn-stats` -> Extracts statistics from model fits.

`micro-select-alignments` -> Selects alignments from microbial analysis based on numerical attributes of fit.

`microbial-fg-gsn-synthetic` -> generates synthetic alignments under GSN from selected seed alignments.

`microbial-nabla` -> Generates nabla estimates from bootstrap result.

`microbial-aeop-locations` -> Generates a file of pseudo-locations to support aeop analyes.

`microbial-teop-synthetic` -> Generates synthetic alignments for the teop analyses.

## Sampling of Drosophila sequences

`dros-filter-alignments` -> Samples from alignments 3rd codon positions and other hard-coded sampling conditions.

## Sampling of Fxy sequences

The sequences from the 3 rodent genomes were aligned. The evaluation of alignment quality and the subsequent selection of alignments for further analysis was done using `nks/sample_fxy_introns.ipynb`. Selected alignments were written to `data/raw/fxy/introns/selected`.

`fxy-filter-positions` -> Produces the alignments for analysis, filtering positions to exclude gaps etc..

## Notebooks

Assume directory containing latex is sister to current directory

Assumes all data and result files in place

Generates all figs and tables for the manuscript plus supplementary material
