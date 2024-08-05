FROM mambaorg/micromamba:1.5.8-lunar

LABEL image.author.name "Alessandro Poletti"
LABEL image.author.email "poletti.alessandro0@gmail.com"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER ./install-dependencies.R /tmp/insstall-dependencies.R

COPY --chown=$MAMBA_USER:$MAMBA_USER ../diffbind_analysis.R /usr/bin/diffbind_analysis.R
COPY --chown=$MAMBA_USER:$MAMBA_USER ../homer2igv.py /usr/bin/homer2igv.py

RUN micromamba create -n differential-analysis

RUN micromamba config set extract_threads 1 && \
    micromamba install -y -n differential-analysis -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/differential-analysis/bin:$PATH

RUN Rscript /tmp/install-dependencies.R

RUN /opt/conda/envs/differential-analysis/share/homer/configureHomer.pl -install hg38                                                                                                                                  482ms  Mon 05 Aug 2024 10:58:07 AM CEST

