# Mutation disequilibrium analysis

This repository contains a Dockerfile and a `.devcontainer` configuration for an environment with Python 3.11, Clang, and Zsh with plugins plus the mdeq package and necessary dependencies to reproduce the results in that work.

We advise you to use VS Code for this, as it provides a consistent mechanism for starting and working within a Docker container, including the ability to run Jupyter notebooks within the container. See these (somewhat out-of-date) [instructions](https://github.com/cogent3/SMBE2024Workshop/wiki/Configuring-your-environment) for getting VS Code on the different platforms.

When you have VS Code working, install the official VS Code extensions for Docker and Dev Containers (both authored by Microsoft). Then open the directory containing this file with VS Code. It should trigger a dialog to start the docker container.

## Features of the Docker environment

- Debian base image.
- Python 3.11 and Clang installed.
- uv and a Python virtual environment.
- Zsh with oh-my-zsh for an enhanced terminal experience.
- Zsh plugins: zsh-autosuggestions and zsh-syntax-highlighting for a more interactive terminal experience.
- accupy python package for accurate floating point arithmetic.
- Eigen C++ library for linear algebra (required for accupy).
- "mdeq==2022.6.30"
- "cogent3[extra]==2025.7.10a5"

## Usage

When the docker container has started within VS Code, you will find the virtual environment already active in the terminal. So the entering `mdeq` in the terminal should display the list of sub-commands available.

Entering `mdeqasis` will show the data prep commands etc.

### Getting the data and results

In the terminal, change into the directory corresponding to the one containing this document

```
cd MutationDiseqAnalysis
```

Then execute the script to download and setup both the data sets and the results.
```
python setup_data_results.py
```
