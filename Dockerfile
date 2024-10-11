FROM mambaorg/micromamba:1.5.8-lunar

LABEL image.author.name "Alessandro Poletti"
LABEL image.author.email "poletti.alessandro0@gmail.com"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba create -n differential-analysis

RUN micromamba config set extract_threads 1 && \
    micromamba install -y -n differential-analysis -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/differential-analysis/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER ./install_dependencies.R /tmp/install_dependencies.R
RUN Rscript /tmp/install_dependencies.R

RUN /opt/conda/envs/differential-analysis/share/homer/configureHomer.pl -install hg38

