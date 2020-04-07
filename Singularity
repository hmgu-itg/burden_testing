Bootstrap: library
From: ubuntu:18.04

%environment
	PERL_MM_USE_DEFAULT=1
	export PERL_MM_USE_DEFAULT
	
	PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"
	export PERL_EXTUTILS_AUTOINSTALL
	
	PATH=$PATH:/usr/local/bin:/usr/local/bin/burden_testing:/usr/local/bin/burden_testing/testing:/usr/local/bin/.vep/htslib:/usr/local/bin/ensembl-vep:/usr/local/bin/MONSTER:/usr/local/bin/bedtools2/bin/:/usr/local/bin/BigWig-Tools:/usr/local/bin/UCSC.tools
	export PATH
	
	LC_ALL=C
	export LC_ALL
	
	PERL5LIB=$PERL5LIB:/usr/local/bin/.vep
	export PERL5LIB

%post
	apt-get install -y software-properties-common
	apt-get install -y build-essential autoconf libtool
	apt-get install -y bc
	#apt-get install -y emacs
	apt-get install -y man
	apt-get install -y git
	apt-get install -y curl
	apt-get install -y wget
	apt-get install -y make
	add-apt-repository universe
	apt-get install -y moreutils
	apt-get install -y tabix
	apt-get install -y libbz2-dev
	apt-get install -y zlib1g-dev
	apt-get install -y libncurses5-dev 
	apt-get install -y libncursesw5-dev
	apt-get install -y liblzma-dev
	apt-get install -y unzip
	apt-get install -y python
	apt-get install -y rsync
	apt-get install -y libgsl-dev
	apt-get install -y r-base
	#apt-get install -y libcurl4-openssl-dev
	
	Rscript --vanilla -e "install.packages(\"data.table\",repos = \"http://cran.us.r-project.org\")"	
	#Rscript --vanilla -e "install.packages(\"remotes\",repos = \"http://cran.us.r-project.org\")"	
	#Rscript --vanilla -e "install.packages(\"BiocManager\",repos = \"http://cran.us.r-project.org\")"	
	#Rscript --vanilla -e "install.packages(\"RCurl\",repos = \"http://cran.us.r-project.org\")"	
	#Rscript --vanilla -e "BiocManager::install(\"SeqArray\")"	
	#Rscript --vanilla -e "BiocManager::install(\"SeqVarTools\")"	
	#Rscript --vanilla -e "remotes::install_github(\"hanchenphd/GMMAT\")"
	
	cpan install Module::Build
	cpan install DBI
	cpan install Try::Tiny
	cpan install JSON
	cpan install Data::Dumper
	#cpan install File::Copy::Recursive
	#cpan install Path::Tiny
	#cpan install Test::File::ShareDir::Dist
	#cpan install DateTime::Locale
	#cpan install DateTime
	cpan install File::Basename
	cpan install Getopt::Long
	cpan install Data::Types
	cpan install File::Path

	cd /usr/local/bin
	git clone https://github.com/hmgu-itg/burden_testing
	
	export PERL5LIB=$PERL5LIB:/usr/local/bin/.vep
	git clone https://github.com/Ensembl/ensembl-vep.git
	cd ensembl-vep
	git checkout release/98
	sed 's/ensembl\.org/ebi\.ac\.uk\/ensemblorg/g' INSTALL.pl | sponge INSTALL.pl
	perl INSTALL.pl -a ac -n --ASSEMBLY GRCh38 -s homo_sapiens -c /usr/local/bin/.vep -d /usr/local/bin/.vep
	
	cd /usr/local/bin
	git clone https://github.com/Carldeboer/BigWig-Tools.git
	cd BigWig-Tools
	git checkout bb83ba1bc28e11e7949884abc26c3f14523dbd0f
	
	cd /usr/local/bin
	mkdir UCSC.tools
	cd UCSC.tools
	rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./
	
	cd /usr/local/bin
	wget http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/MONSTER_v1.3.tar.gz
	tar -zxf MONSTER_v1.3.tar.gz
	rm MONSTER_v1.3.tar.gz
	cd MONSTER
	make
	chmod o+rx /usr/local/bin/MONSTER
	
	cd /usr/local/bin
	wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
	tar -zxf bedtools-2.29.1.tar.gz
	rm bedtools-2.29.1.tar.gz
	cd bedtools2
	make
	
	cd /usr/local/bin
	wget https://sourceforge.net/projects/transpose/files/transpose/transpose-2.0/2.0/transpose-2.0.zip
	unzip transpose-2.0.zip
	rm transpose-2.0.zip	
	cd transpose-2.0/src
	gcc transpose.c -o transpose2
	mv transpose2 /usr/local/bin
	cd /usr/local/bin
	rm -rf transpose-2.0/
	
	CREATIONDATE=`date`

%runscript
	echo "This container was created: $CREATIONDATE"

%labels
	Author HMGU ITG
	Version v1.4.1

%help
	This is a container designed to run burden testing pipeline; for more information run this container with the help command line option


