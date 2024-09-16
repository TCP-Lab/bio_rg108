FROM rocker/tidyverse:4.4

RUN Rscript --vanilla -e "install.packages(c('BiocManager', 'devtools'), repos='https://cloud.r-project.org')" && \
    Rscript --vanilla -e "BiocManager::install('DESeq2')" && \
    Rscript --vanilla -e "devtools::install_github('TCP-Lab/DeNumerator', ref = 'v0.1.2', upgrade='never')"

COPY . .
