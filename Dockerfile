FROM condaforge/mambaforge:24.9.2-0

# Set the working directory
WORKDIR /work

# Create the conda environment using Mamba
COPY environment.yml .
ENV CONDA_ENV_NAME=busco_phylogenomics
RUN mamba env create -f environment.yml \
    && mamba clean --all --yes
ENV PATH="${CONDA_DIR}/envs/${CONDA_ENV_NAME}/bin:${PATH}"

CMD ["/bin/bash"]
