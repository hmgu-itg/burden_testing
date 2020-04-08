FROM ubuntu:18.04
# This container allows you to run rare variant aggregation tests using MONSTER and SMMAT; for more information run this container with the help command line option.
LABEL Author Arthur Gilly, Andrei Barysenka, Daniel Suveges
LABEL Version v1.5
ENV TZ=Europe/Berlin
ENV PERL_MM_USE_DEFAULT=1
ENV PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
ENV PATH=$PATH:/usr/local/bin:/usr/local/bin/burden_testing:/usr/local/bin/burden_testing/testing:/usr/local/bin/.vep/htslib:/usr/local/bin/ensembl-vep:/usr/local/bin/MONSTER:/usr/local/bin/bedtools2/bin/:/usr/local/bin/UCSC.tools
ENV LC_ALL=C
ENV PERL5LIB=$PERL5LIB:/usr/local/bin/.vep
RUN apt update
RUN DEBIAN_FRONTEND="noninteractive" apt install -y software-properties-common build-essential autoconf libtool bc man git curl wget make moreutils libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev unzip python libgsl-dev r-base libcurl4-openssl-dev emacs axel
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && tar -xvjf htslib-1.10.2.tar.bz2 && cd htslib-1.10.2 && make tabix && make bgzip && cp bgzip tabix /usr/bin
RUN Rscript --vanilla -e "install.packages(c(\"reshape2\", \"parallel\", \"Hmisc\", \"argparser\", \"data.table\", \"BiocManager\"),repos = \"http://cran.us.r-project.org\");BiocManager::install(c(\"SeqArray\", \"SeqVarTools\"));install.packages(\"GMMAT\", repos = \"http://cran.us.r-project.org\")"
RUN perl -MCPAN -e 'foreach (@ARGV) { CPAN::Shell->rematein("notest", "install", $_) }' Module::Build DBI Try::Tiny JSON Data::Dumper File::Basename Getopt::Long Data::Types File::Path
RUN cd /usr/local/bin
RUN git clone https://github.com/hmgu-itg/burden_testing
RUN export PERL5LIB=$PERL5LIB:/usr/local/bin/.vep
RUN git clone https://github.com/Ensembl/ensembl-vep.git
RUN cd ensembl-vep
RUN git checkout release/98
RUN mkdir -p /usr/local/bin/.vep && cd /usr/local/bin/.vep && axel -q ftp://ftp.ebi.ac.uk/ensemblorg/pub/release-98/variation/indexed_vep_cache/homo_sapiens_vep_98_GRCh38.tar.gz && echo Untarring... && tar -xzf homo_sapiens_vep_98_GRCh38.tar.gz && rm homo_sapiens_vep_98_GRCh38.tar.gz && cd -
RUN sed 's/ensembl\.org/ebi\.ac\.uk\/ensemblorg/g' INSTALL.pl | sponge INSTALL.pl
RUN perl INSTALL.pl -a ac -n --ASSEMBLY GRCh38 -s homo_sapiens -c /usr/local/bin/.vep -d /usr/local/bin/.vep
RUN cd /usr/local/bin
RUN mkdir UCSC.tools
RUN cd UCSC.tools
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/liftOver
RUN cd /usr/local/bin
RUN wget http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/MONSTER_v1.3.tar.gz
RUN tar -zxf MONSTER_v1.3.tar.gz
RUN rm MONSTER_v1.3.tar.gz
RUN cd MONSTER
RUN make
RUN chmod o+rx /usr/local/bin/MONSTER
RUN cd /usr/local/bin
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
RUN tar -zxf bedtools-2.29.1.tar.gz
RUN rm bedtools-2.29.1.tar.gz
RUN cd bedtools2
RUN make
RUN cd /usr/local/bin
RUN wget https://sourceforge.net/projects/transpose/files/transpose/transpose-2.0/2.0/transpose-2.0.zip
RUN unzip transpose-2.0.zip
RUN rm transpose-2.0.zip
RUN cd transpose-2.0/src
RUN gcc transpose.c -o transpose2
RUN mv transpose2 /usr/local/bin
RUN cd /usr/local/bin
RUN rm -rf transpose-2.0/
RUN CREATIONDATE=`date`
CMD echo "This container was created: $CREATIONDATE"
