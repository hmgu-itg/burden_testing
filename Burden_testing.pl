#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use JSON;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Data::Types;
use File::Path qw(make_path);

# Version information:
our $version = "v6.0 Last modified: 31.Mar.2020";

# Get script directory:
our $scriptDir = dirname(__FILE__);

use lib dirname(__FILE__);

$\="\n";

##-----------------------------------------------------------------------------------------------------------
#                                   ASSUMING EIGEN SCORES ARE b37 BASED, CADD SCORES ARE b38 BASED (v. 1.5,  https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz)
#
##-----------------------------------------------------------------------------------------------------------

# Loading custom modules:
use GENCODE;
use Scoring;

# Status report:
print "[Info] Script version: $version";
#printf "[Info] Run date: %s", DateTime->now->strftime("%Y. %b %d %H:%M");
printf "[Info] The script was called with the following parameters:\n%s\n", join(" ", $0, @ARGV);

# some default parameter values
my $parameters = {
    "build"   => "38",
    "extend"  => 0, # Basepairs with which the GENCODE features will be extended.
    "MAF"     => 0.05, # The upper threshold of the minor allel frequency.
    "MAC"     => 0, # The lower threshold of the minor allele count.
    "missingthreshold" => 0.01, # The upper threshold of the missingness. (below which the the missingness will be imputed)
    "score"   => "NA", # Applied score to weight variants.
    "cutoff"  => 0, # Score threshold below which the variants will be removed.
    "floor"   => 0, # All the scores below threshold will be set to this value.
    "shift"   => 0, # A value with which the scores of the variants will be shifted.
    "chr_prefix" => "chr", # Chromosome prefix in VCF files or list
    "scriptDir"=>${scriptDir}, 
    "lof" => undef,	
    "smmat" => undef,	
    "vcf" => undef	
};

# LoF VEP consequences
$parameters->{"lof_cons"} = {
    "transcript_ablation"      => 10,
    "splice_acceptor_variant"  => 9,
    "splice_donor_variant"     => 8,
    "stop_gained"              => 7,
    "frameshift_variant"       => 6,
    "stop_lost"                => 5,
    "start_lost"               => 4,
    "transcript_amplification" => 3,
    "inframe_insertion"        => 2,
    "inframe_deletion"         => 1
};

# Command line options without default values:
my ($inputFile,$outputDir,$outputFile, $help);

# Parsing command line options:
GetOptions(
    # Input/Output:
    'input=s' => \$inputFile,
    'output=s' => \$outputFile,
    'output-dir=s' => \$outputDir,

    # Selecting region source:
    'GENCODE=s' => \$parameters->{"input"}->{"GENCODE"},
    'GTEx=s' => \$parameters->{"input"}->{"GTEx"},
    'overlap=s' => \$parameters->{"input"}->{"overlap"},

    # Extending regions with a defined length:
    'extend=s' => \$parameters->{"extend"},

    # Skipping minor transcripts by APPRIS:
    'SkipMinor' => \$parameters->{"minor"},

    # Variant features:
    'MAF=s' => \$parameters->{"MAF"},
    'MAC=s' => \$parameters->{"MAC"},

    # Verbose output:
    'verbose' => \$parameters->{"verbose"},

    # bgizipped, indexed, 1-based input list for SMMAT
    'smmat=s' => \$parameters->{"smmat"},

    # specifying config file:
    'config=s' => \$parameters->{"configName"},

    # Which score do we need:
    'score=s' => \$parameters->{"score"},
    
    'maxVars=s' => \$parameters->{"maxVars"},

    # Do we need LoF variants only:
    'lof' => \$parameters->{"lof"},
    'loftee' => \$parameters->{"loftee"}, # Filters only high and low confident loss of function variants
    'lofteeHC' => \$parameters->{"lofteeHC"}, # Filters for only high confident loss of function variants

    # Accepting missingness filter:
    'missingness=s' => \$parameters->{"missingthreshold"},

    # Changing scores that will be used for weights:
    'shift=s'  => \$parameters->{"shift"},  # The value used to shift eigen scores:
    'cutoff=s' => \$parameters->{"cutoff"}, # hard threshold that will be applied on scores:
    'floor=s'  => \$parameters->{"floor"},  # How do we want to floor Eigen values:

    'chromosome-prefix=s'  => \$parameters->{"chr_prefix"},  # Chromosome prefix in VCF file(s)
    'vcf=s' => \$parameters->{"vcf"},
    'help|h' => \$help
    );

$parameters->{"maxVars"}=1000 unless $parameters->{"maxVars"};

&usage && die "[Error] No output directory specified. Exiting." unless $outputDir;
if (! -d $outputDir){
    die "[Error] Could not create output directory ($outputDir). Exiting." unless make_path($outputDir)==1;
}

$parameters->{"tempdir"}=$outputDir."/prepare_regions_tempfiles";
if (! -d $parameters->{"tempdir"}){
    die "[Error] Could not create temp directory (".$parameters->{"tempdir"}."). Exiting." unless make_path($parameters->{"tempdir"})==1;
}

&usage && die "[Error] No config file specified. Exiting." unless $parameters->{"configName"};
&usage && die "[Error] The specified config file does not exist. Exiting." unless -e $parameters->{"configName"};
&usage && die "[Error] Gene list input file has to be specified with the --input option. Exiting." unless $inputFile;
&usage && die "[Error] The specified input gene list does not exist. Exiting." unless -e $inputFile;
&usage && die "[Error] Output file has to be specified with the --output option. Exiting." unless $outputFile;
&usage && die "[Error] VCF files or a SMMAT input list have to be specified with the --vcf or --smmat option. Exiting." unless(defined($parameters->{"vcf"}) || defined($parameters->{"smmat"}));
&usage && die "[Error] both VCF files and a SMMAT input list are specified. Exiting." if(defined($parameters->{"vcf"}) && defined($parameters->{"smmat"}));

if (! defined($parameters->{"smmat"})){
    &usage && die "[Error] No VCF files exist." unless &checkVCFs($parameters->{"vcf"});
}

$parameters = &readConfigFile($parameters);

# If the score option is not empty, we have to check if it's a valid score, and the
# required files are exists. If any problem found, the score parameter will be set to its
# default value -> no score.
$parameters = &check_scores($parameters) if $parameters->{"score"} ne "NA";

# Processing the requested features:
$parameters = &parseGENCODE($parameters) if $parameters->{"input"}->{"GENCODE"};
$parameters = &parseRegulation($parameters) if $parameters->{"input"}->{"GTEx"} || $parameters->{"input"}->{"overlap"};

our $verbose = $parameters->{"verbose"};

# Report submitted parameters:
&print_parameters($parameters);

# Initializing helper objects:
print "[Info] Initializing GENCODE data" if ($verbose);
my $GENCODE_data = GENCODE->new($parameters);
print "[Info] Initializing score data data" if ($verbose);
my $AddScore = Scoring->new($parameters);

# Open files:
open (my $INPUT, "<", $inputFile) or die "[Error] Input file ($inputFile) could not be opened. Exiting.";

my ($SNPinfo,$SNPfile,$genotypeFile);

if (defined($parameters->{"vcf"})){
    open ($SNPfile, ">", $outputDir."/".$outputFile."_variant_file.txt") or die "[Error] Output file could not be opened.";
    open ($genotypeFile, ">", $outputDir."/".$outputFile."_genotype_file.txt") or die "[Error] Output genotype file could not be opened.";
}

if (defined($parameters->{"vcf"})){
    open ($SNPinfo, ">", $outputDir."/".$outputFile."_SNPinfo_file.txt") or die "[Error] Output SNPinfo file could not be opened.";
}
else{
    open ($SNPinfo, ">", $outputDir."/".$outputFile.".txt") or die "[Error] Output group file could not be opened.";
}

# Processing the input file gene by gene:
my $gene_count = 0;
while ( my $ID = <$INPUT> ){
    next if $ID=~/^\s*$/;
    
    chomp $ID;

    # genomic coordinates of the current gene
    my ($chr, $start, $end, $stable_ID, $name, $CollapsedBed);

    ($chr, $start, $end, $stable_ID, $name) = $GENCODE_data->GetCoordinates($ID);
    if ($start eq "NA") {
	print "[Warning] Gene $ID was not found in the GENCODE data; skipping [NO_GENE]";
	print "";
	next;
    }

    print "[Info] Queried gene: $name (Ensembl ID: $stable_ID), Genomic location: $chr:$start-$end (Input: $ID)";

    # all lines from the Linked_features file that are associated with the current gene
    my $bedlines = &BedToolsQuery($chr, $start, $end, $stable_ID, $parameters->{"Linked_features"});
    # remove some lines we're not interested in
    $CollapsedBed = &FilterLines($bedlines, $stable_ID, $parameters);

    # This should never be a problem, but still be tested:
    unless ( $CollapsedBed ){
	print "[WARNING] Gene $name did not yield any regions. Skipped. [NO_REGION].";
	print "";
	next;
    }

    # CollapsedBed is 0-based
    
    # Once we have the genomic regions, the overlapping variants have to be extracted:
    my $variants = &getVariants($CollapsedBed, $parameters);

    # Filtering variants based on the provided parameters:
    my ($hash, $genotypes) = &processVar($variants, $parameters,$stable_ID);
    # Gene will be skipped if there are no suitable variations left:
    if (scalar keys %{$hash} < 2){
	print "[Warning] Gene $ID is skipped as not enough variants left to test [NOT_ENOUGH_VAR].";
	print "";
	undef $hash;
	undef $genotypes;
	next;
    }

    # The gene will be skipped if there are too many variants:
    if (scalar keys %{$hash} > $parameters->{"maxVars"}){
	print "[Warning] Gene $ID is skipped as more than ".$parameters->{"maxVars"}." variants are in the set [TOO_MANY_VAR].";
	print "";
	undef $hash;
	undef $genotypes;
	next;
    }

    # If we want we can add scores:
    $hash = $AddScore->AddScore($hash) if $parameters->{"score"} ne "NA";
    checkScores($hash) if $parameters->{"score"} ne "NA";

    if ($parameters->{"score"} ne "NA" && scalar keys %{$hash} < 1){
	print "[Warning] Gene $ID is skipped as no variants remaining post-scoring [NO_VAR_REMAIN].";
	print "";
    }

    # We don't save anything unless there at least two variants:
    if (!defined($parameters->{"smmat"})){
	next unless scalar keys %{$hash} > 1;
    }

    # Once we have the scores we have to print out the SNP file:
    my $flag=0;
    $flag=1 if $parameters->{"score"} ne "NA";

    if (defined($parameters->{"smmat"})){
	&print_group_file($hash, $ID, $SNPinfo, $parameters->{"build"},$flag);
    }
    else{# outut for MONSTER from VCFs
	&print_SNPlist($hash, $ID, $SNPfile,$flag);
	&print_SNP_info($hash, $ID, $SNPinfo, $gene_count, $parameters->{"build"},$flag);
	&print_genotypes($genotypes, $genotypeFile, $parameters, $gene_count);
    }
    
    $gene_count ++; # So the header will only be printed once.
}

close $INPUT;
if (defined($parameters->{"vcf"})){
    close $SNPfile;
    close $genotypeFile;
}
close $SNPinfo;

sub check_parameters {
    my $parameters = shift;

    # checking build:
    die "[Error] Only GRCh37 and GRCh38 genome builds are supported (--build 37 or 38)." unless $parameters->{"build"} eq "37" || $parameters->{"build"} eq "38";

    # Checking MAF (a float below 1):
    my $maf=$parameters->{"MAF"};
    die "[Error] $maf must be a float number < 1" unless defined($maf) && is_float($maf) && $maf<1;

    # Checking MAC (an integer):
    my $mac=$parameters->{"MAC"};
    die "[Error] $mac must be an integer number > 0" unless defined($mac) && is_int($mac) && $mac>0;

    # Checking missingness (a float below 1):
    my $miss=$parameters->{"missingthreshold"};
    die "[Error] $miss must be a float number < 1" unless defined($miss) && is_float($miss) && $miss<1;
}

sub check_scores {
    my $params = $_[0];

    # Check if the specified score is a valid, supported score:
    my %acceptedScores = ("CADD" => 1,
			  #                      "Eigen" => 1,
			  #                      "EigenPC" => 1,
			  "EigenPhred" => 1);
#                      "EigenPCPhred" => 1,
#                      "Linsight" => 1,
#                      "Mixed" => 1);

    # Let's report if the specified score is not supported:
    unless (exists $acceptedScores{$params->{"score"}}){
        print "[Info] Supported weighting options: ", join(", ", keys %acceptedScores), "";
        printf "[Warning] The specified score (%s) is currently not supported.", $params->{"score"};
        print "[Warning] No weighting will be applied.";
        $params->{"score"} = "NA";
        return $params;
    }
    # Now, let's check if the requirements of the specified weighting methods are satisified.
    if ($params->{"score"} eq "CADD" ){
        if( ! exists $params->{"caddPath"}){
            print "[Warning] The config file has not entry for 'caddPath', pointing to the genome-wide CADD scores.";
            print "[Warning] CADD scores as weight cannot be used. No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
        elsif (! -e $params->{"caddPath"}){
            printf "[Warning] The specified genome-wide CADD scores are not available %s.", $params->{"caddPath"};
            print "[Warning] CADD scores as weight cannot be used. No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
    }
    elsif ($params->{"score"} =~ /eigen/i){
        if (! exists $params->{"EigenPath"} ) {
            print "[Warning] The config file has not entry for 'EigenPath', pointing to the genome-wide Eigen scores.";
            print "[Warning] Eigen scores as weights cannot be used. No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
    }

    return $params;
}

sub readConfigFile {
    my $params = $_[0];
    open(my $CONF, "<", $params->{"configName"}) or die "[Error] In readConfigFile: config file could not be oppended. Exiting.";

    while ( my $line = <$CONF>) {
        chomp $line;
        next unless $line;
        next if $line =~ /^#/;

        my ($key, $value) = $line =~ /^(\S+)=(\S+)/;
        $params->{$key} = $value if $key && $value;
    }
    return $params;
}

sub parseGENCODE {
    my %AcceptedFeatures = ( "gene" => 1, "exon" => 1, "transcript" => 1, "CDS" => 1, "UTR" => 1);
    my $parameters = $_[0];

    foreach my $feature (split(",", $parameters->{"input"}->{"GENCODE"})){

        # Checking if the provided feature is exists or not:
        if (exists $AcceptedFeatures{$feature}) {
            $parameters->{"GENCODE"}->{$feature} = 1;
        }
        else {
            printf STDERR "[Error] The provided GENCODE feature name ($feature) is not supported: %s. Use these: %s", $feature, join(", ", keys %AcceptedFeatures);
        }
    }
    return $parameters
}

sub parseRegulation {
    # Accepted features: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg
    my %AcceptedFeatures = ( "promoter" => 1, "CTCF" => 1, "enhancer" => 1, "promoterFlank" => 1, "openChrom" => 1, "TF_bind" => 1, "allreg" => 1 );
    my $parameters = $_[0];
    my @reg_class = ("GTEx", "overlap");

    foreach my $class (@reg_class){
        next unless exists $parameters->{"input"}->{$class};
        my %hash = ();
        # Looping through all the submitted features:
        foreach my $feature (split(",",$parameters->{"input"}->{$class})){
            unless (exists $AcceptedFeatures{$feature}){
                printf STDERR "[Error] The provided regulatory feature name ($feature) is not supported: %s. Use these: %s", $feature, join(", ", keys %AcceptedFeatures);
                next;
            }

            $hash{'TF_binding_site'}   = 1 if "TF_bind" eq $feature;
            $hash{'CTCF_binding_site'} = 1 if "CTCF" eq $feature;
            $hash{'enhancer'}          = 1 if "enhancer" eq $feature;
            $hash{'Open chromatin'}    = 1 if "openChrom" eq $feature;
            $hash{'promoter'}          = 1 if "promoter" eq $feature;
            $hash{'allreg'}            = 1 if "allreg" eq $feature;
            $hash{'promoter_flanking_region'} = 1 if "promoterFlank" eq $feature;

        }
        $parameters->{$class} = \%hash;
    }

    return $parameters;
}

sub print_parameters {
    local $\="\n";
    $parameters = $_[0];
    printf "\n[Info] Current genome build: GRCh%s", $parameters->{"build"};
    print "\n[Info] Selected features:";
    printf "\tThe following GENCODE features are considered: %s", join (", ", keys %{$parameters->{"GENCODE"}}) if exists $parameters->{"GENCODE"};
    printf "\tGENCODE features are extended by: %sbp", $parameters->{"extend"} if exists $parameters->{"GENCODE"};
    printf "\tThe following GTEx linked regulatory features are considered: %s", join (", ", keys %{$parameters->{"GTEx"}}) if exists $parameters->{"GTEx"};
    printf "\tThe following overlapping regulatory features are considered: %s", join (", ", keys %{$parameters->{"overlap"}}) if exists $parameters->{"overlap"};

    print "\n[Info] Variant filters:";
    printf "\tUpper MAF threshold: %s", $parameters->{"MAF"};
    printf "\tLower MAC threshold: %s", $parameters->{"MAC"};
    printf "\tUpper missingness threshold: %s", $parameters->{"missingthreshold"};
    print "\tOnly severe variants will be included in the test." if $parameters->{"lof"};
    print "\tOnly loss-of-function variants will be included in the test (loftee HC and LC)." if $parameters->{"loftee"};
    print "\tOnly high-confidence loss-of-function variants will be included in the test (loftee HC)." if $parameters->{"lofteeHC"};

    print "\n[Info] Variant weighting:";
    printf "\tWeighting method: %s\n", $parameters->{"score"};
    if ($parameters->{"score"} ne "NA"){
        printf "\tScore lower cutoff: %s", $parameters->{"cutoff"};
        printf "\tScore floor applied: %s", $parameters->{"floor"};
        printf "\tScore shifted by: %s\n", $parameters->{"shift"};
    }
}

# get relevant lines from Linked_features file
sub BedToolsQuery {
    my ($chr, $start, $end, $stable_ID, $geneBedFile) = @_;
    my $queryString = sprintf("intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\t%s\") -b %s -sorted 2>/dev/null | cut -f9-",$chr, $start-1, $end, $stable_ID, $geneBedFile); # start and end are from GENCODE (1-based), bed is 0-based
    
    print "[Info] IntersectBed query string: $queryString" if $verbose;
    my $query =Scoring::backticks_bash($queryString);
    #print "QUERY: $query" if $verbose;
    return $query;
}

sub formatLines {
    my %hash = %{$_[0]};
    my $ext = $_[1] // 0;
    
    # 0-based
    return sprintf("%s\t%s\t%s\t%s", $hash{"chr"}, $hash{"start"} - $ext, $hash{"end"} + $ext, encode_json(\%hash))
}

sub getVariants {
    # merged is 0-based
    my ($merged, $parameters) = @_;
    my $inputvcfpattern=$parameters->{"vcf"};
    my $prefix=$parameters->{"chr_prefix"};
    my $smmat=$parameters->{"smmat"};
    my $distance = 0;
    my $variants;
    my $source;
    
    # Finding out which chromosome are we on:
    my ($chr) = $merged =~ /^(.+?)\t/;  # WATCHOUT !

    if (defined($smmat)){
	print  "[Info] Extracting variants from the list:" if $verbose;
	$source=$smmat;
	# SMMAT lists have no chromosome prefix
	$prefix="";
    }
    else{
	print  "[Info] Extracting variants from vcf file(s):" if $verbose;
	(my $vcf = $inputvcfpattern ) =~ s/\%/$chr/g;
	
	if (! -e $vcf){
	    print "[Warning] VCF file for chromosome $chr does not exist";
	    return "";
	}
	$source=$vcf;
    }
    
    my $tabix_query = sprintf("tabix %s ", $source);

    # looping through all lines:
    foreach my $line (split("\n", $merged)){
	my ($chr, $start, $end) = split("\t", $line);
	$distance += $end - $start;
	$tabix_query .= sprintf(" %s%s:%s-%s", $prefix, $chr, $start+1, $end); # VCFs and lists are 1-based
    }

    print  "$tabix_query" if $verbose;
    $variants = Scoring::backticks_bash($tabix_query);

    print sprintf("[Info] Total covered genomic regions: %s bp", $distance) if $verbose;
    
    return $variants;
}

sub FilterLines {
    my ($lines, $stable_ID, $parameters) = @_;

    my @output_lines = ();
    my %hash = ();

    # The following hashes read from the parameter set:
    my %GENCODE = exists $parameters->{"GENCODE"} ? %{$parameters->{"GENCODE"}} : ();
    my %GTEx = exists $parameters->{"GTEx"} ? %{$parameters->{"GTEx"}} : ();
    my %overlap = exists $parameters->{"overlap"} ? %{$parameters->{"overlap"}} : ();

    foreach my $line (split (/\n/, $lines)){
	#print "LINE: $line";
        # Decoding json line:
        my %annot_hash = %{decode_json($line)};

        # Filtering for lines which belong to our target gene:
        next unless $annot_hash{"gene_ID"} eq $stable_ID;

        my $source = $annot_hash{"source"};
        my $class =  exists $annot_hash{"class"} ? $annot_hash{"class"} :  "?";
        # print "\n$source $class";

        if ($source eq "GENCODE" && exists $GENCODE{$class}) {

            # If the user has specified, minor transcripts will be excluded:
            next if exists $GENCODE{"minor"} && $annot_hash{"appris"} eq "Minor";

            push (@output_lines, formatLines(\%annot_hash, $parameters->{'extend'}));
            $hash{GENCODE}{$class} ++;
        }
        elsif (($source eq "GTEx" && exists $GTEx{$class})
               || ($source eq "GTEx" && exists $GTEx{'allreg'})){

            push (@output_lines, formatLines(\%annot_hash));
            $hash{GTEx}{$class} ++;
        }
        elsif (($source eq "overlap" && exists $overlap{$class})
               || ($source eq "overlap" && exists $overlap{'allreg'})){

            push (@output_lines, formatLines(\%annot_hash));
            $hash{overlap}{$class} ++;

        }
    }

    if ($verbose) {
        if (exists $hash{GENCODE}) {
            my $report = '';
            foreach my $feature (keys %{$hash{GENCODE}}){
                $report .= "$feature: $hash{GENCODE}{$feature}, ";
            }
            print  "[Info] From the GENCODE dataset the following features were extracted: $report" if $verbose;
        }
        if (exists $hash{GTEx}) {
            my $report = '';
            foreach my $feature (keys %{$hash{GTEx}}){
                $report .= "$feature: $hash{GTEx}{$feature}, ";
            }
            print  "[Info] The following regulatory features were linked to the gene: $report" if $verbose;
        }
        if (exists $hash{overlap}) {
            my $report = '';
            foreach my $feature (keys %{$hash{overlap}}){
                $report .= "$feature: $hash{overlap}{$feature}, ";
            }
            print  "[Info] The following overlapping regulatory features were extracted: $report"  if $verbose;
        }
    }
    
    # Report:
    # print "[Info] Selected lines:\n", join("\n", @output_lines),"" if $verbose;

    # Saving temporary bedfile:
    my $tmpName=$parameters->{"tempdir"}."/filtered_regions.bed";
    open (my $tempbed, "> $tmpName");
    foreach my $line (@output_lines){
        print $tempbed $line;
    }
    close $tempbed;
    
    `sort -k1,1n -k2,2n $tmpName | sponge $tmpName`;

    # Collapsing overlapping features:
    my $queryString = "mergeBed -i $tmpName";
    my $merged = Scoring::backticks_bash($queryString);
    return $merged;
}

sub checkScores{
    my %hash = %{$_[0]};
    for my $v (keys(%hash)){
	die "[Error]: score for $v is not defined" unless (defined($hash{$v}{"score"}));
    }
}

sub getVariantType{
    my ($ref,$alt)=@_;

    return "SNP" if length($ref)==1 && length($alt)==1;
    return "DEL" if length($ref)>length($alt);
    return "INS" if length($ref)<length($alt);
    
    return "NA";    
}

# Using VEP
sub getConsequences{
    my $variants   = $_[0];
    my $parameters   = $_[1];
    my $stable_ID=$_[2];

    my $vepin;
    my $vepout;
    
    my %cons=();

    local $\="\n";
    local $,="\t";
    
    my $fname1=$parameters->{"tempdir"}."/vep_input.txt";

    open ($vepin, ">", $fname1) or die "[Error] Input file for VEP could not be opened.";
    
    my @total_vars = split("\n", $variants);
    my $count=0;
    while (@total_vars){
	my $variant = shift @total_vars;

	# line: CHROM	POS	ID	REF	ALT	...
	my ($chr, $pos, $id, $ref, $alt, @therest) = split(/\t/, $variant);
	
        (my $c = $chr ) =~ s/chr//i;
	
	my $varID=$c."_".$pos."_".$ref."_".$alt;

	# skip multiallelics
	if ($alt=~/,/){
	    $cons{$varID}="NA";
	    next;
	}

	my $vtype=getVariantType($ref,$alt);
	if ($vtype eq "SNP"){
	    print $vepin $c,$pos,$pos,$ref."/".$alt,"+",$varID;
	    $count++;
	}
	elsif($vtype eq "DEL"){
	    if (length($alt)==1){
		my $r=substr($ref,1,length($ref)-1);
		print $vepin $c,$pos+1,$pos+length($ref)-1,$r."/-","+",$varID;
	    }
	    else{
		my $L1=length($alt);
		my $r=substr($ref,$L1,length($ref)-1);
		print $vepin $c,$pos+$L1,$pos+length($ref)-1,$r."/-","+",$varID;
	    }
	    $count++;
	}
	elsif($vtype eq "INS"){
	    if (length($ref)==1){
		my $a=substr($alt,1,length($alt)-1);
		print $vepin $c,$pos+1,$pos,"-/".$a,"+",$varID;
	    }
	    else{
		my $L1=length($ref);
		my $a=substr($alt,$L1,length($alt)-$L1);
		print $vepin $c,$pos+$L1,$pos+$L1-1,"-/".$a,"+",$varID;	    
	    }
	    $count++;
	}
	else{
	    print "[Error] could not determine variant type of $variant";
	    $cons{$varID}="NA";
	    next;
	}
    }
    close($vepin);

    return undef unless($count>0);
    
    my $queryString="vep -i ".$fname1." --dir /usr/local/bin/.vep --dir_cache /usr/local/bin/.vep -o STDOUT --offline --no_stats | grep -v \"^#\" | awk -v g=".$stable_ID." 'BEGIN{FS=\"\\t\";}\$4==g{print \$0;}' | cut -f 1,7";
    print $queryString if $verbose;
    
    my $query =Scoring::backticks_bash($queryString);
    my %max_severity;
    foreach my $line (split (/\n/, $query)){
	print "getConsequences output: ".$line if $verbose;
	my @a=split(/\t/,$line);
	my $ID=$a[0];
	my $effs=$a[1];
	$max_severity{$ID}=0 unless exists($max_severity{$ID});
	
	foreach my $e (split(/,/,$effs)){
	    if (! exists($parameters->{"lof_cons"}->{$e})){ # low severity
		    $cons{$ID}=$e if ($max_severity{$ID}==0);
	    }
	    else{
		if ($parameters->{"lof_cons"}->{$e} > $max_severity{$ID}){
		    $max_severity{$ID}=$parameters->{"lof_cons"}->{$e};
		    $cons{$ID}=$e;
		}
	    }
	}

    }

    return \%cons;
}

sub processVar {
    my $variants   = $_[0]; # List of all overlapping variants
    my $parameters = $_[1]; # All the submitted parameters to filter variants.
    my $stable_ID=$_[2];
    
    my %hash = (); # SNP info container.
    my %genotypeContainer = (); # genotype information container.

    # Which build are we using?
    my $build = "GRCh".$parameters->{"build"};

    print "[Info] Processing variants:" if $verbose;

    my $cons;
    # --------------------------------------------------------
    # we only need consequences if --lof option is provided
    if (defined($parameters->{"lof"})){
	$cons=getConsequences($variants,$parameters,$stable_ID);
    }
    # --------------------------------------------------------
    
    my @total_vars = split("\n", $variants);
    printf "[Info] Total number of overlapping variants: %s\n", scalar(@total_vars) if $verbose;

    if ($parameters->{"smmat"}){
	while (@total_vars){
	    my $consequence="NA"; # default
	    # Removing one variant at a time:
	    my $variant = shift @total_vars;

	    # line: CHROM	POS    	ID REF	ALT
	    my ($chr, $pos, $id,$a1, $a2) = split(/\t/, $variant);
	    my $SNPID = sprintf("%s_%s_%s_%s", $chr, $pos, $a1, $a2);
	    # We don't consider indels if weights are used:
	    if (( length($a2) > 1 || length($a1) > 1 ) && $parameters->{"score"} ne "NA"){
		print  "[Warning] $SNPID will be omitted because it's an indel and we use scores for weighting! ($a1/$a2).";
		next;
	    }
	    
	    # We don't consider multialleleic sites
	    if ( $a2 =~ /,/){
		print  "[Warning] $SNPID will be omitted because it's multiallelic! ($a2).";
		next;
	    }

	    # --------------------------------------------------------
	    if (defined($parameters->{"lof"})){
		(my $chr2 = $chr ) =~ s/chr//i;
		my $varID=$chr2."_".$pos."_".$a1."_".$a2;
		if (defined($cons)){
		    if (exists($cons->{$varID})){
			$consequence=$cons->{$varID};
			print "CONSEQUENCE: ".$varID." ".$consequence if $verbose;
		    }
		}
	    }

	    # If loss of function variants are required, we skip all those variants that are not LoF:
	    if ( $parameters->{"lof"} && ! exists $parameters->{"lof_cons"}->{$consequence} ) {
		printf "[Warning] $SNPID will be omitted because its consequence (%s) is not lof\n", $consequence;
		next;
	    }

	    # --------------------------------------------------------

	    # TODO: check if we need that for SMMAT
	    # Generating variant name (Sometimes the long allele names cause problems):
	    #my $short_a1 = length $a1 > 5 ? substr($a1,0,4) : $a1;
	    #my $short_a2 = length $a2 > 5 ? substr($a2,0,4) : $a2;
	    
	    #$SNPID = sprintf("%s_%s_%s_%s", $chr, $pos, $short_a1, $short_a2);
	    
	    $hash{$SNPID}{"alleles"} = [$a1, $a2];
	    $hash{$SNPID}{$build} = [$chr, $pos - 1, $pos]; # 0-based
	    $hash{$SNPID}{"consequence"} = $consequence;
	}
    }
    else{	
	while (@total_vars){
	    my $consequence="NA"; # default
	    # Removing one variant at a time:
	    my $variant = shift @total_vars;

	    # line: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EGAN00001033155
	    my ($chr, $pos, $id, $a1, $a2, $qual, $filter, $info, $format, @genotypes) = split(/\t/, $variant);

	    # --------------------------------------------------------
	    if (defined($parameters->{"lof"})){
		(my $chr2 = $chr ) =~ s/chr//i;
		my $varID=$chr2."_".$pos."_".$a1."_".$a2;
		if (defined($cons)){
		    if (exists($cons->{$varID})){
			$consequence=$cons->{$varID};
			print "CONSEQUENCE: ".$varID." ".$consequence if $verbose;
		    }
		}
	    }
	    # --------------------------------------------------------

	    # Generating variant name (Sometimes the long allele names cause problems):
	    my $short_a1 = length $a1 > 5 ? substr($a1,0,4) : $a1;
	    my $short_a2 = length $a2 > 5 ? substr($a2,0,4) : $a2;

	    my $SNPID = sprintf("%s_%s_%s_%s", $chr, $pos, $short_a1, $short_a2);

	    # Parsing info field for relevant information:
	    (my $ac )= $info =~ /AC=(.+?)[;\b]/;
	    (my $an )= $info =~ /AN=(.+?)[;\b]/;
	    
	    # If AN and AC values are not found we skip the variant:
	    if (! $an){
		print "[Warning] $SNPID has no AN; skipping";
		next;
	    }

	    if (! $ac){
		print "[Warning] $SNPID has no AC; skipping";
		next;
	    }

	    # --------------------------------------------------------
	    # We add NA if consequence was not found:
	    #(my $consequence ) = $info =~ /consequence=(.+?)[;\b]/;
	    #$consequence = "NA" unless $consequence;
	    # --------------------------------------------------------

	    # We don't consider multialleleic sites
	    if ( $a2 =~ /,/){
		print  "[Warning] $SNPID will be omitted because it's multiallelic! ($a2).";
		next;
	    }

	    print "[Warning] $variant has no AC" if $ac eq "NA";
	    print "[Warning] $variant has no AN" if $an eq "NA";

	    # Calculating values for filtering:
	    my $missingness = (scalar(@genotypes)*2 - $an)/(scalar(@genotypes)*2);
	    my $MAF = $ac/$an;

	    # This flag shows if the non-ref allele is the major one:
	    my $genotypeFlip = 0;
	    if ( $MAF > 0.5 ){
		print "[Info] MAF of $SNPID is $MAF is greater then 0.5, genotype is flipped.";
		$MAF = 1 - $MAF;
		$genotypeFlip = 1;
	    }

	    my $MAC = $ac;
	    $MAC = $an - $ac if $MAC > $an / 2;

	    # We don't consider indels if weights are used:
	    if (( length($a2) > 1 || length($a1) > 1 ) && $parameters->{"score"} ne "NA"){
		print  "[Warning] $SNPID will be omitted because it's an indel and we use scores for weighting! ($a1/$a2).";
		next;
	    }
	    # Filter out variants because of high missingness:
	    if ( $missingness > $parameters->{"missingthreshold"} ){
		print  "[Warning] $SNPID will be omitted because of high missingness ($missingness).";
		next;
	    }
	    # Filter out variant because of high MAF (regardless of the applied weight):
	    if ( $MAF > $parameters->{'MAF'} ){
		printf "[Warning] $SNPID will be omitted because of high MAF (%.3f, cutoff: %s)\n", $MAF, $parameters->{'MAF'};
		next;
	    }
	    # Filter out variant because of high MAF:
	    if ( $ac < $parameters->{'MAC'} ){
		printf "[Warning] $SNPID will be omitted because of low minor allele count ($ac, cutoff: %s)\n", $parameters->{'MAC'};
		next;
	    }
	    # If loss of function variants are required, we skip all those variants that are not LoF:
	    if ( $parameters->{"lof"} && ! exists $parameters->{"lof_cons"}->{$consequence} ) {
		printf "[Warning] $SNPID will be omitted because its consequence (%s) is not lof\n", $consequence;
		next;
	    }

	    # TODO: fix loftee part
	    # If loftee or lofteeHC are enabled, the script exits if no LoF_conf tag is present in the info field.
	    if (( $parameters->{"lofteeHC"} or $parameters->{"loftee"}) and $info !~ /LoF_conf/ ){
		die "[Error] Based on the provided parameters, variant selection based on the loftee prediction was requested.\n\tHowever the provided vcf file does not contain the obligatory LoF_conf tag.Exiting.";
	    }

	    # If loftee is enabled, only low and high-confidence loss-of-function variants will be selected:
	    if ( $parameters->{"loftee"} && $info =~ /LoF_conf\=-/ ) {
		print "[Warning] $SNPID will be omitted because it is not a high- or low-confidence loss of function variant";
		next;
	    }

	    # If lofteeHC is enabled, only high-confidence loss-of-function variants will be selected:
	    if ( $parameters->{"lofteeHC"} && $info !~ /LoF_conf\=HC/ ) {
		print "[Warning] $SNPID will be omitted because it is not a high- confidence loss of function variant";
		next;
	    }

	    # If everything went fine, initializing variant by adding missingness to the hash:
	    $hash{$SNPID}{"missingness"} = $missingness;

	    # Storing variant data for all variants:
	    # TODO: if genotypeFlip == 1, then ac, an and MAF refer to a1, not a2
	    $hash{$SNPID}{"alleles"} = [$a1, $a2];

	    # TODO: for indels end should be different
	    $hash{$SNPID}{$build} = [$chr, $pos - 1, $pos]; # 0-based
	    $hash{$SNPID}{"frequencies"} = [$ac, $an, $MAF];
	    $hash{$SNPID}{"consequence"} = $consequence;
	    $hash{$SNPID}{"rsID"} = $id;

	    # Parsing genotypes:
	    foreach my $gt (@genotypes){
		my @fields = split(":", $gt);
		if    ($fields[0] eq "0/0") { push(@{$genotypeContainer{$SNPID}}, $genotypeFlip == 1 ? 2 : 0) }
		elsif ($fields[0] eq "1/1") { push(@{$genotypeContainer{$SNPID}}, $genotypeFlip == 1 ? 0 : 2) }
		elsif ($fields[0] eq "1/0" or $fields[0] eq "0/1") { push(@{$genotypeContainer{$SNPID}}, 1) }
		else {push(@{$genotypeContainer{$SNPID}}, "-9")}
	    }
	}
    }

    printf  ( "[Info] Number of filtered variants: %s\n", scalar(keys %hash));
    return (\%hash, \%genotypeContainer);
}

sub print_SNPlist {
    local $\=undef;
    my %hash = %{$_[0]};
    my $gene_name = $_[1];
    my $outputhandle = $_[2];
    my $flag = $_[3];

    # Get the list of variants:
    my @snpIDs = keys %hash;

    # The line with the snps will be written anyways:
    print $outputhandle "$gene_name\t$flag\t", join("\t", @snpIDs),"\n";

    # We return if the flag is 0, so there is no scores given:
    return 1 if $flag == 0;

    # If scores are given, we add an extra line:
    print $outputhandle "$gene_name\t0\t";
    my @scores = ();
    foreach my $snpid (keys %hash){
	die "[Error]: variant $snpid has no defined score" unless defined($hash{$snpid}{"score"});
        push (@scores, $hash{$snpid}{"score"});
    }
    print $outputhandle join("\t", @scores),"\n";

    return 1;
}

sub print_SNP_info {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];
    my $outfilehandle = $_[2];
    my $gene_counter = $_[3];
    my $build = $_[4];
    my $flag=$_[5];
    
    # Extract variant names:
    my @variants = keys %hash;

    # Assembling header for the first gene:
    if ( $gene_counter == 0){
        my $header = "Gene_name\tSNPID(GRCh${build})\trsID\tAlleles\tconsequences\tMAF\tAC\tmissigness";
        $header .= "\tScore" if ($flag==1);
        print $outfilehandle $header;
    }

    # printing out info for each gene and variant:
    foreach my $variant (values %hash){
        my $rsID = $variant->{'rsID'};
	die "[Error]: variant $rsID has no defined score" unless($flag==0 || defined($variant->{"score"}));
        # generating all the required fields:
        my $alleleString = join("/", @{$variant->{'alleles'}});
        my $MAF = sprintf "%.4f", $variant->{'frequencies'}->[2];
        my $missingness = sprintf "%.4f", $variant->{'missingness'};
        my $weight = $variant->{'score'} if ($flag==1);
        my $consequence = $variant->{'consequence'};
        my $AC = $variant->{'frequencies'}->[0];
        my $SNPID = $variant->{'GRCh'.$build}->[0].":".$variant->{'GRCh'.$build}->[2];
        my $line = join ("\t", $gene_name, $SNPID, $rsID, $alleleString, $consequence, $MAF, $AC, $missingness);
        $line .= "\t$weight" if ($flag==1);
        print $outfilehandle $line;
    }
}

# if no scores provided, then output score=1
sub print_group_file {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];
    my $outfilehandle = $_[2];
    my $build = $_[3];
    my $flag=$_[4]; # if we  have scores
    
    # Extract variant names:
    my @variants = keys %hash;

    # printing out info for each gene and variant:
    foreach my $variant (values %hash){
	my $weight=1;
	
        my $rsID = $variant->{'rsID'};
	die "[Error]: variant $rsID has no defined score" unless($flag==0 || defined($variant->{"score"}));
	my $chr=$variant->{'GRCh'.$build}[0];
	$chr  =~ s/chr//i; # after liftover in AddScore we will have "chr" prefix

	my $pos=$variant->{'GRCh'.$build}[2];
	my $ref=$variant->{'alleles'}->[0];
	my $alt=$variant->{'alleles'}->[1];
        $weight = $variant->{'score'} if ($flag==1);
#        my $consequence = $variant->{'consequence'};
        my $line = join ("\t", $gene_name, $chr, $pos, $ref,$alt, $weight);
        print $outfilehandle $line;
    }
}

sub print_genotypes {
    my %genotype      = %{$_[0]};
    my $outputhandler = $_[1];
    my $parameters    = $_[2];
    my $gene_counter  = $_[3];

    my $vcfChrFile = &getVCF($parameters->{"vcf"});

    # Assembling header for the first gene:
    if ( $gene_counter == 0){
        # Get list of sample IDs:
        my $samples = `zgrep -m1 "#CHROM"  $vcfChrFile | cut -f10-`;
	chomp($samples);
        print $outputhandler "0\t$samples";
    }

    # Saving data:
    for my $var (keys %genotype){
        print $outputhandler "$var\t", join("\t", @{$genotype{$var}});
    }
}

###

sub usage {
    print "Usage::";
    print("      Required:");
    print("          --input <input file>");
#    print("          --working-dir <working directory containing Linked_features.bed.gz and gencode.basic.annotation.tsv.gz>");
    print("          --output-dir <output directory where output and temporary files will be created>");
    print("          --output <output filename prefix>");
    print("          --vcf <input VCF(s); either --vcf or --smmat is required>");
    print("          --smmat <5 column tab-delimited input list of variants, for SMMAT; either --vcf or --smmat is required>");
    print("          --config <config file>");
    print("      Optional:");    
#    print("          --build <genome build; default: 38>");
    print("          --GENCODE <comma separated list of GENCODE features (gene, exon, transcript, CDS or UTR)>");
    print("          --GTEx <comma separated list of GTEx features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)>");
    print("          --overlap <comma separated list of overlap features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)>");
    print("          --extend <by how many basepairs the GENCODE features should be extended.>");
    print("          --MAF <MAF upper threshold; default: 0.05>");
    print("          --MAC <MAC lower threshold; default: 0>");
    print("          --maxVars <max number of variants in a gene; default: 1000>");
    print("          --SkipMinor <skip minor transcripts by APPRIS>");
    print("          --verbose <increase verbosity>");
    print("          --score <which score to use to weight variants; one of: CADD, EigenPhred >");
    print("          --lof <only select high impact variants: transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification, inframe_insertion, inframe_deletion>");
#    print("          --loftee <only select high and low confident loss of function variants>");
#    print("          --lofteeHC <only select high confident loss of function variants>");
    print("          --missingness <missingness upper threshold; default: 0.01>");
    print("          --shift <shift scores by this value; default: 0>");
    print("          --cutoff <score threshold below which the variants will be removed; default: 0>");
#    print("          --no-filtering <skip MAF, missingness and MAC variant filtering>");
    print("          --floor <scores below this threshold will be set to this value; default: 0>");
    print("          --chromosome-prefix <chromosome prefix in VCF files; default: \"chr\">");
    print("          --help <this help>");
}

sub checkVCFs{
    my $inputvcfpattern=$_[0];
    print $inputvcfpattern;
    for (my $i=1;$i<=22;$i++){
	my $fname=$inputvcfpattern;
	$fname=~s/\%/$i/;
	return 1 if -e $fname;
    }

    return undef;
}

sub getVCF{
    my $inputvcfpattern=$_[0];

    for (my $i=1;$i<=22;$i++){
	my $fname=$inputvcfpattern;
	$fname=~s/\%/$i/;
	return $fname if -e $fname;
    }

    return undef;
}
