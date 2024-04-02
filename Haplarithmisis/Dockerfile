ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        apt \
        bash \
        curl \
        libcurl4-openssl-dev \
        libext2fs2 \
        libfribidi-dev \
        libpam-modules-bin \
        libpam-runtime \
        libss2 \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("BiocManager", "1.30.20")' \
    && Rscript -e 'remotes::install_version("devtools", "2.4.5")' \
    && Rscript -e 'remotes::install_version("vcfR", "1.15.0")' \
    && Rscript -e 'remotes::install_version("textshaping", "0.3.6")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' 

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocGenerics", \
        "BiocVersion", \
        "SNPlocs.Hsapiens.dbSNP144.GRCh38", \
        "VariantAnnotation", \
        "snpStats", \
        "limma", \
        "signal", \
        "plotrix", \
        "MASS", \
        "data.table", \
        "QDNAseq", \
        "Biobase" \
    ), ask = FALSE, force = TRUE)' 
