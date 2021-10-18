FROM ubuntu:18.04

# This container allows you to run rare variant aggregation tests using MONSTER and SMMAT; for more information run this container with the help command line option.
ENV TZ=Europe/Berlin
ENV PERL_MM_USE_DEFAULT=1
ENV PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
ENV PATH="/usr/local/bin:$PATH"
ENV LC_ALL=C

RUN apt update
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y software-properties-common

# Install R and other libraries
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN DEBIAN_FRONTEND="noninteractive" apt install -y software-properties-common build-essential autoconf libtool bc man git curl wget make moreutils libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev unzip python libgsl-dev r-base libcurl4-openssl-dev axel
	
# Install HTSlib
WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 \
    && tar xjf htslib-1.13.tar.bz2 \
    && rm htslib-1.13.tar.bz2 \
    && cd htslib-1.13 \
    && ./configure && make && make install


# Install BCFTools
WORKDIR /opt
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 \
    && tar xjf bcftools-1.13.tar.bz2 \
    && rm bcftools-1.13.tar.bz2 \
	&& cd bcftools-1.13 \
	&& ./configure && make && make install

RUN Rscript --vanilla -e "install.packages(c('reshape2', 'R.utils', 'parallel', 'Hmisc', 'argparser', 'data.table', 'BiocManager', 'doMC'), repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('SeqArray', 'SeqVarTools'))"

# Install devtools for R
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev \
    && Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"

# Install GMMAT
RUN Rscript -e "library(devtools); remotes::install_github('hanchenphd/GMMAT@406345258de1ae7b24d76933f142e631475856b4')"

# Perl
RUN perl -MCPAN -e 'foreach (@ARGV) { CPAN::Shell->rematein("notest", "install", $_) }' Module::Build DBI Try::Tiny JSON Data::Dumper File::Basename Getopt::Long Data::Types File::Path

# liftOver (UCSC.tools)
WORKDIR /usr/local/bin
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/liftOver \
    && chmod +x liftOver

# MONSTER
WORKDIR /opt
RUN wget http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/MONSTER_v1.3.tar.gz \
    && tar -zxf MONSTER_v1.3.tar.gz \
    && rm MONSTER_v1.3.tar.gz \
    && cd MONSTER \
    && make \
    && chmod +x /opt/MONSTER/MONSTER \
    && ln -s /opt/MONSTER/MONSTER /usr/local/bin/MONSTER

# BEDtools
WORKDIR /opt
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz \
    && tar -zxf bedtools-2.29.1.tar.gz \
    && rm bedtools-2.29.1.tar.gz \
    && cd bedtools2 \
    && make \
    && chmod +x /opt/bedtools2/bin/*
ENV PATH="/opt/bedtools2/bin:$PATH"

# Transpose
WORKDIR /opt
RUN	wget https://sourceforge.net/projects/transpose/files/transpose/transpose-2.0/2.0/transpose-2.0.zip \
    && unzip transpose-2.0.zip \
    && rm transpose-2.0.zip	 \
    && cd transpose-2.0/src \
    && gcc transpose.c -o /usr/local/bin/transpose2 \
    && rm -rf /opt/transpose-2.0


# Burden Testing
WORKDIR /opt
RUN git clone --depth 1 --branch v1.5.3 https://github.com/hmgu-itg/burden_testing
ENV PATH="/opt/burden_testing:/opt/burden_testing/testing:$PATH"

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION

LABEL maintainer1="Arthur Gilly <arthur.gilly@helmholtz-muenchen.de>"
LABEL maintainer2="Young-Chan Park <young-chan.park@helmholtz-muenchen.de>"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.vcs-ref=$VCS_REF
LABEL version=$VERSION