FROM nextflow/nextflow:25.04.6

# Install system dependencies as root
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install miniforge, install mamba, and set the PATH
ENV CONDA_DIR=/opt/conda
ARG MINIFORGE_VERSION=25.3.1-0
ARG MINIFORGE_URL=https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-Linux-x86_64.sh
RUN wget --quiet ${MINIFORGE_URL} -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniforge.sh \
    && ${CONDA_DIR}/bin/conda install -n base -c conda-forge mamba \
    && ${CONDA_DIR}/bin/conda clean -afy

ENV PATH="${CONDA_DIR}/bin:$PATH"

# Set the working directory
WORKDIR /pipeline

# Copy environment file and set ownership to the non-root 'nextflow' user
COPY --chown=nextflow:nextflow environment.yml .

# Create the conda environment using Mamba
ENV CONDA_ENV_NAME=busco_phylogenomics
RUN mamba env create -f environment.yml \
    && mamba clean --all --yes
ENV PATH="$CONDA_DIR/envs/$CONDA_ENV_NAME/bin:$PATH"

# Switch to a non-root user
USER nextflow

# Copy the rest of the pipeline code, setting ownership
COPY --chown=nextflow:nextflow . .

# Keep Nextflow as the entrypoint
ENTRYPOINT ["nextflow"]
CMD ["run", "main.nf"]
