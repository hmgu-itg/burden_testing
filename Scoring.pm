package Scoring;

## This class was created to add various scores to a selected set of variants.

## Usage:

# Initialize:
#
# use Scoring;
# my $Scoring = Scoring->new($parameters);

# Adding scores:
# $hash = $Scoring->AddScore($hash);

# ASSUMING EIGEN SCORES ARE 37 BASED, CADD 38 BASED

# TODO: ASSUMING SCORE FILES HAVE CHROMOSOME NAMES LIKE 1,2,3 ...

use strict;
use warnings;
use Data::Dumper;

sub new {
    my ( $class, $parameters ) = @_;
    my $self = {};
    bless( $self, $class );
    
    # Storing all the important data:
    $self->{"score"} = $parameters->{"score"};
    $self->{"shift"} = $parameters->{"shift"};
    $self->{"cutoff"} = $parameters->{"cutoff"};
    $self->{"floor"} = $parameters->{"floor"};
    $self->{"build"} = $parameters->{"build"};
    $self->{"verbose"} = $parameters->{"verbose"};
    $self->{"scriptDir"} = $parameters->{"scriptDir"};
    $self->{"tempdir"} = $parameters->{"tempdir"};
    $self->{"smmat"} = $parameters->{"smmat"}; # so that we know if we have to prepend chromosome name with "chr" when lifting over

    foreach my $k (keys(%{$parameters->{"lof_cons"}})){
	$self->{"lof_cons"}->{$k}=$parameters->{"lof_cons"}->{$k};
    }
    
    # Storing paths:
    $self->{"EigenPath"} = $parameters->{"EigenPath"};
    $self->{"caddPath"} = $parameters->{"caddPath"};
    
#    print Dumper $self;
    return $self;
}

sub backticks_bash {
  open(my $pipe, '-|', '/bin/bash', '-c', $_[0])
     or return;
  local $/ = wantarray ? $/ : undef;
  <$pipe>
}

sub AddScore {
    my $self = shift;
    my $hash = shift;

# no need to liftover if we're using CADD scores (they are b38 based)    
    if ($self->{"score"} =~ /Eigen/i || $self->{"score"} eq "Mixed"){
	$hash = _liftover($self, $hash) unless $self->{"build"} eq "37";
    }

    $hash = _get_Eigen_Score($self, $hash) if $self->{"score"} =~ /Eigen/i;
    $hash = _get_CADD($self, $hash) if $self->{"score"} eq "CADD";

    # Returning mixed Eigen/CADD scores
    $hash = _get_mixed($self, $hash) if $self->{"score"} eq "Mixed";

    # Processing scores based on the flooring method and the cutoff:
    $hash = _process_score($self, $hash);
    
    return $hash;
}

sub _get_mixed {
    # REQUIRES CONSEQUENCE !
    my $self = $_[0];
    my %hash = %{$_[1]};
    my @to_delete=();
    
    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){

        # Extracting consequence:
        my $consequence = $hash{$var}{'consequence'};
        (my $chr37 = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;
        (my $chr38 = $hash{$var}{GRCh38}[0] ) =~ s/chr//i;

        # Initialize score:
        $hash{$var}{"score"} = "NA";

	if ($consequence eq "NA"){
	    printf("[Warning] In Scoring->_get_mixed(): variant %s has not been annotated with a consequence; removing it\n",$var);
            push @to_delete,$var;
	    next;
	}
	
        # If the consequence is not severe, we use Eigen, otherwise CADD:
	# Eigen-raw       5
	# Eigen-phred     7
	# Eigen-PC-raw    6
	# Eigen-PC-phred  8
	if (! exists($self->{"lof_cons"}->{$consequence})){
        #if ( $consequence =~ /intron|intergenic|regulatory|non_coding|upstream|downstream/i ){
            my $EigenFile = $self->{"EigenPath"};
	    $EigenFile=~s/\%/$chr37/i;
	    if (! -e $EigenFile){
		print "[Warning] file with Eigen scores for chromosome $chr37 does not exist. Variant $var will be removed.";
		push @to_delete,$var;
		next;
	    }
	    
            my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f 1-4,7 | grep %s", $EigenFile, $chr37, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2], $hash{$var}{alleles}[1]);
            printf "[Info] %s has benign consequences (%s), so EigenPhred scores are used.\n", $var, $consequence;
            print "$tabix_query\n" if $self->{"verbose"};
	    my $lines=backticks_bash($tabix_query);

            foreach my $line (split("\n", $lines)){
                chomp $line;
                # chr     pos     a1      a2      EigenPhred
                my @array = split("\t", $line);

                # Testing alleles:
                if (($array[2] eq $hash{$var}{alleles}[0] and $array[3] eq $hash{$var}{alleles}[1]) or
                    ($array[2] eq $hash{$var}{alleles}[1] and $array[3] eq $hash{$var}{alleles}[0])) {

                    $hash{$var}{"score"} = $array[4]; # Phred scaled Eigen score.
                }
            }

            $hash{$var}{"flag"} = "Non-coding variant, EigenPhred score added.";
        }
        else {
            printf "[Info] %s has severe consequence (%s) so CADD-Phred scores are used.\n", $var, $consequence if $self->{"verbose"};
	    my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f 3,4,6", $self->{"caddPath"}, $chr38, $hash{$var}{GRCh38}[2], $hash{$var}{GRCh38}[2]);
            print "[Info] $tabix_query\n" if $self->{"verbose"};
	    my $lines=backticks_bash($tabix_query);

            # Looping through all ouput lines:
            foreach my $line (split("\n", $lines)){
		# Line: C	G    25.0 : REF ALT PHRED
		chomp $line;
		my ($ref, $alt, $PHRED) = split("\t", $line);

		# Testing alleles:
		if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
		    ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {
                    $hash{$var}{"score"} = $PHRED;
		}
            }
        }
    }

    for my $v (@to_delete){
	print "[Warning]: Variant $v has NA score; removing it";
	delete($hash{$v});
    }
    
    print "[Info] Mixed EigenPhred/CADDPhred scores have been added to variants.\n" if $self->{"verbose"};
    my $fflag=scalar(keys(%hash))==0;
    print("[Warning] No variants remain after adding mixed weights. \n") if($fflag);
    
    return \%hash
}

# CADD are b38 based
# using PHRED CADD score
sub _get_CADD {
    my $self = $_[0];
    my %hash = %{$_[1]};
    my @to_delete=();
    
    # Looping throuh all variants and returning the CADD PHRED score:
    foreach my $var (keys %hash){
        (my $chr = $hash{$var}{GRCh38}[0] ) =~ s/chr//i;
        my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f 3,4,6", $self->{"caddPath"}, $chr, $hash{$var}{GRCh38}[2], $hash{$var}{GRCh38}[2]); # CADD file is 1-based
        print "[Info] $tabix_query" if $self->{"verbose"};
	my $lines=backticks_bash($tabix_query);
        $hash{$var}{"score"} = "NA";

        foreach my $line (split("\n", $lines)){
            # Line: C	G    25.0 : REF ALT phred
            chomp $line;
            my ($ref, $alt, $phred) = split("\t", $line);

            # Testing alleles:
            if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
                ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {
                    $hash{$var}{"score"} = $phred;
            }
        }

	push @to_delete,$var if ($hash{$var}{"score"} eq "NA");
    }
    
    printf "[Info] %s scores have been added to variants.\n", $self->{"score"} if $self->{"verbose"};
    for my $v (@to_delete){
	print "[Warning]: Variant $v has no CADD score; removing it";
	delete($hash{$v});
    }
    
    my $fflag=scalar(keys(%hash))==0;
    print("[Warning] No variants remaining after CADD scoring\n") if($fflag);
    return \%hash
}

# Eigen are b37 based
sub _get_Eigen_Score {
    my $self = $_[0];
    my %hash = %{$_[1]};
    my @to_delete=();
    
    print  "[Info] Adding Eigen scores\n" if $self->{"verbose"};

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){
        my $EigenFile = $self->{"EigenPath"};
	print "EIGENFILE: $EigenFile" if $self->{"verbose"};
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;
	$EigenFile=~s/\%/$chr/i;
	if (! -e $EigenFile){
	    print "[Warning] File with Eigen scores for chromosome $chr does not exist. Variant $var will be removed.";
	    push @to_delete,$var;
	    next;
	}

        # Two tabix queries will be submitted regardless of the output...
        my $tabix_query = sprintf("tabix %s %s:%s-%s | grep %s ", $EigenFile, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2], $hash{$var}{alleles}[1]); # Eigen file is 1-based
        print "$tabix_query" if $self->{"verbose"};
	my $lines=backticks_bash($tabix_query);

        $hash{$var}{"score"} = "NA"; # Initialize Eigen score.

        foreach my $line (split("\n", $lines)){
            chomp $line;
	    
            # chr     pos     a1      a2      Eigen   EigenPC    EigenPhred    EigenPCPhred
            my @array = split("\t", $line);

            # Testing alleles:
            if (($array[2] eq $hash{$var}{alleles}[0] and $array[3] eq $hash{$var}{alleles}[1]) or
                ($array[2] eq $hash{$var}{alleles}[1] and $array[3] eq $hash{$var}{alleles}[0])) {

                $hash{$var}{"score"} = $array[4] || "NA" if $self->{"score"} eq "Eigen"; # Un-scaled raw Eigen score.
                $hash{$var}{"score"} = $array[6] || "NA" if $self->{"score"} eq "EigenPhred"; # Phred scaled Eigen score.
                $hash{$var}{"score"} = $array[5] || "NA" if $self->{"score"} eq "EigenPC"; # Un-scaled Eigen PC score. 
                $hash{$var}{"score"} = $array[7] || "NA" if $self->{"score"} eq "EigenPCPhred"; # Phred scaled Eigen PC score.
            }
        }

        # Reporting if no score has been found for the variant:
        if ( $hash{$var}{"score"} eq "NA") {
            printf ( "[Warning] %s score was not found for variant %s in the non-coding set; removing it\n\n", $self->{"score"}, $var);
	    push @to_delete,$var;
        }
    }

    for my $v (@to_delete){
	delete $hash{$v};
    }
    
    printf "[Info] Eigen scores have been added to variants (Number of variants: %s).\n\n", scalar keys %hash if $self->{"verbose"};
    my $fflag=scalar(keys(%hash))==0;
    print("[Warning] No variants remaining after Eigen score annotation\n") if($fflag);
    return \%hash
}

# liftOver variant coordinates from 38 to 37
sub _liftover {
    my $self = $_[0];
    my %hash = %{$_[1]};
    my $tempFileName = $self->{"tempdir"}."/temp_GRCH38.bed";
    my $tempFileName37 = $self->{"tempdir"}."/temp_GRCH37.bed";
    my $tempFileNameU = $self->{"tempdir"}."/temp_unmapped.bed";

    my $prefix="";
    $prefix="chr" if (defined$self->{"smmat"});
    
    # Checking if the file exists, in which case we delete it:
    `rm $tempFileName` if -e $tempFileName; 
    `rm $tempFileName37` if -e $tempFileName37; 
    `rm $tempFileNameU` if -e $tempFileNameU; 
    
    # Saving the GRCh38 coordinates:
    open( my $tempbed, ">", $tempFileName) or die "[Error] Temporary bedfile could not be opened for writing: $tempFileName\n";
    foreach my $variant (keys %hash){
        printf $tempbed "%s%s\t%s\t%s\t%s\n", $prefix,$hash{$variant}{"GRCh38"}[0], $hash{$variant}{"GRCh38"}[1], $hash{$variant}{"GRCh38"}[2], $variant;
    }

    # Liftover query:
    my $liftover_query = sprintf("liftOver %s %s/hg38ToHg19.over.chain %s %s  2> /dev/null", $tempFileName, $self->{"scriptDir"},$tempFileName37,$tempFileNameU);

    print "Liftover query: ".$liftover_query if $self->{"verbose"};
    
    # Calling liftover:
     backticks_bash($liftover_query);

    # Reading mapped file:
    open(my $lifted, "< $tempFileName37") or die "[Error] After liftover run, the mapped file could not be opened.\n";
    
    my $liftedVarNo  = 0;
    
    # Adding GRCh37
    while (my $line = <$lifted>) {
        chomp $line;
        my ($chr, $start, $end, $SNPID) = split("\t", $line);
        $hash{$SNPID}{"GRCh37"} = [$chr, $start, $end]; # still 0-based
        $liftedVarNo ++;
    }
    close($lifted);
    
    # Removing variants where the liftover failed:
    foreach my $variant (keys %hash){
        unless (exists $hash{$variant}{"GRCh37"}){
            print "[Info] Deleting $variant: could not be mapped to GRCH37.\n";
            delete $hash{$variant}     
        }
    }

    print  "[Info] Number of variants successfully lifted over: $liftedVarNo\n\n" if $self->{"verbose"};
    my $flag = $liftedVarNo == 0 ? 0 : 1;
    return (\%hash);
}

sub _process_score {
    my $self = $_[0];
    my %hash = %{$_[1]}; # hash with data
    
    printf "[Info] Processing %s scores... \n", $self->{"score"} if $self->{"verbose"};

    # Looping through all variants and modifying scores:
    foreach my $var ( keys %hash){

        # Shifting first
        $hash{$var}{"score"} = $hash{$var}{"score"} + $self->{"shift"} if $self->{"shift"};
	
	if(defined($self->{"cutoff"}) && $hash{$var}{"score"} < $self->{"cutoff"}){
	    if (defined($self->{"floor"})){
		printf "[Warning] score of %s is set to %s, because %s score is below the cutoff threshold: %s < %s!\n",$var, $self->{"floor"}, $self->{"score"}, $hash{$var}{"score"}, $self->{"cutoff"};
		$hash{$var}{"score"} = $self->{"floor"};
	    }
	    else{
		printf "[Warning] deleting %s, because %s score is below the cutoff threshold and no floor value is provided: %s < %s!\n", $var, $self->{"score"}, $hash{$var}{"score"}, $self->{"cutoff"};
		delete $hash{$var};		
	    }
	}
    }
    
    print "Done.\n" if $self->{"verbose"};
    return \%hash;
}

1;
