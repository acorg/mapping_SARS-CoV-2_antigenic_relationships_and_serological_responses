
# Mapping SARS-CoV-2 antigenic relationships and serological responses

This repository contains all the data and code for the paper "Mapping SARS-CoV-2 antigenic relationships and serological responses".

## Viewing interactive figures
Interactive versions of each of the figures in the `docs` directory and are available [here](https://acorg.github.io/mapping_SARS-CoV-2_antigenic_relationships_and_serological_responses).

## Running the code
### Dependencies
Code package dependencies are listed in the `DEPENDENCIES` file.

#### titertools (v0.0.0.9002)
Several functions we use for a likelihood-based approach to dealing with censored titer data are included in the "titertools" package available at github.com/shwilks/titertools. Here we use version the package version 0.0.0.9002.

#### Racmacs (v1.2.3)
For performing antigenic cartography we use the "Racmacs" package v1.2.3, available from github.com/shwilks/Racmacs.

#### ablandscapes (v1.2.3)
For plotting antibody landscapes we use the "ablandscapes" package v1.1.2, available from github.com/shwilks/ablandscapes.

To reprocess the raw data and generate the data in `data/generated_data` run the file `code/generate_data.R`. Note that some of the map bootstrapping in particular may take a considerable amount of time to run.

To regenerate all figures run the file `code/generate_figures.R`.

To regenerate all interactive figures run the file `code/generate_docs.R`.

## Directory structure
| Directory | Description |
| ----------- | ----------- |
| `code/figures` | Code for figure generation |
| `code/docs` | Code for interactive figure generation in `docs/` |
| `code/data_generation` | Code for data processing |
| `data/raw_data` | Raw measurement data files |
| `data/generated_data` | Data generated after processing (for example data after model fits) |
| `data/maps` | Antigenic map data |
| `docs` | Interactive figure documents |
| `figures` | Generated figures |
| `functions` | Functions shared across code |

