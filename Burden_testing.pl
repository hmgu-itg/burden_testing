#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use JSON;
use DateTime;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Data::Types;

# For debugging:
use Devel::Size qw(total_size);

# Version information:
our $version = "v5.1 Last modified: 05.Nov.2019";

# Get script directory:
our $scriptDir = dirname(__FILE__);

use lib dirname(__FILE__);

$\="\n";

##-----------------------------------------------------------------------------------------------------------
#                                   ASSUMING ALL SCORES ARE 37 BASED
# TODO: --config <config.txt>
##-----------------------------------------------------------------------------------------------------------


# Loading custom modules:
use GENCODE;
use Scoring;

# Status report:
print "[Info] Script version: $version";
printf "[Info] Run date: %s", DateTime->now->strftime("%Y. %b %d %H:%M");
printf "[Info] The script was called with the following parameters:\n%s\n", join(" ", $0, @ARGV);

# In the new version, all the parameters will be stored in a hash reference:
my $parameters = {
    "scriptDir" => $scriptDir,
    "build"   => "38", # !!! HARDCODED BUILD !!!
    "extend"  => 0, # Basepairs with which the GENCODE features will be extended.
    "MAF"     => 0.05, # The upper threshold of the minor allel frequency.
    "MAC"     => 0, # The lower threshold of the minor allele count.
    "missingthreshold" => 0.01, # The upper threshold of the missingness. (below which the the missingness will be imputed)
#    "configFileName" => "$scriptDir/config.txt", # The config file from which the source files will be read.
    "score"   => "NA", # Applied score to weight variants.
    "cutoff"  => 0, # Score threshold below which the variants will be removed.
    "floor"   => 0, # All the scores below this threshold will be set to this value.
    "shift"   => 0, # A value with which the scores of the variants will be shifted.
"chr_prefix" => "chr" # Chromosome prefix for VCF files	
};

# This is the list of those consequences that will be retained upon switching on --lof
$parameters->{"lof_cons"} = {
    "transcript_ablation"      => 1,
    "splice_acceptor_variant"  => 1,
    "splice_donor_variant"     => 1,
    "stop_gained"              => 1,
    "frameshift_variant"       => 1,
    "stop_lost"                => 1,
    "start_lost"               => 1,
    "transcript_amplification" => 1,
    "inframe_insertion"        => 1,
    "inframe_deletion"         => 1,
    "splice_region_variant"    => 1,
};

# Command line options without default values:
my ($inputFile, $outputFile, $help);

# Parsing command line options:
GetOptions(
    #'build=s' => \$parameters->{"build"},

    # Input/Output:
    'input=s' => \$inputFile,
    'output=s' => \$outputFile,

    # Selecting region source:
    'GENCODE=s' => \$parameters->{"input"}->{"GENCODE"},
    'GTEx=s' => \$parameters->{"input"}->{"GTEx"},
    'overlap=s' => \$parameters->{"input"}->{"overlap"},

    # Extending regions with a defined length:
    'extend=s' => \$parameters->{"extend"},

    # GENCODE file
    'gencode-file=s' => \$parameters->{"gencode_file"},
    
    # Skipping minor transcripts by APPRIS:
    'SkipMinor' => \$parameters->{"minor"},

    # Variant features:
    'MAF=s' => \$parameters->{"MAF"},
    'MAC=s' => \$parameters->{"MAC"},

    # Verbose output:
    'verbose' => \$parameters->{"verbose"},

    # specifying config file:
    'configFile=s' => \$parameters->{"configFileName"},

    # Which score do we need:
    'score=s' => \$parameters->{"score"},
    
    'maxVars=s' => \$parameters->{"maxVars"},

    # Do we need only loss of function:
    'lof' => \$parameters->{"lof"},
    'loftee' => \$parameters->{"loftee"}, # Filters only high and low confident loss of function variants
    'lofteeHC' => \$parameters->{"lofteeHC"}, # Filters for only high confident loss of function variants

    # Accepting missingness filter:
    'missingness=s' => \$parameters->{"missingthreshold"},

    # Changing scores that will be used for weights:
    'shift=s'  => \$parameters->{"shift"},  # The value used to shift eigen scores:
    'cutoff=s' => \$parameters->{"cutoff"}, # hard threshold that will be applied on scores:
    'floor=s'  => \$parameters->{"floor"},  # How do we want to floor Eigen values:

    'chromosome-prefix=s'  => \$parameters->{"chr_prefix"},  # Chromosome prefix in VCF files

    # vcf File:
    'vcfFile=s' => \$parameters->{"vcfFile"},

    # TEMP DIR
    'tempdir=s' => \$parameters->{"tempdir"},
    
    # Asking for help:
    'help|h' => \$help
    );

print "SCORE: ".$parameters->{"score"};

if ($help || !defined($inputFile) || !defined($outputFile) || !defined($parameters->{"vcfFile"}) || !defined($parameters->{"gencode_file"}) || !defined($parameters->{"tempdir"})){
    &usage();
    exit(1);
}

$parameters->{"maxVars"}=1000 unless $parameters->{"maxVars"}

# Exit unless the absolute necessary input files are exists or specified:
die "[Error] Gene list input file has to be specified with the --input option. Exiting." unless $inputFile;
die "[Error] Output file has to be specified with the --output option. Exiting." unless $outputFile;
die "[Error] VCF files have to be specified with the --vcfFile option. Exiting." unless $parameters->{"vcfFile"};
die "[Error] No VCF files exist." unless &checkVCFs($parameters->{"vcfFile"});
die "[Error] The specified input gene list does not exist. Exiting." unless -e $inputFile;
die "[Error] No config file specified. Exiting." unless $parameters->{"configFileName"};
die "[Error] The specified config file does not exist. Exiting." unless -e $parameters->{"configFileName"};
die "[Error] No temp dir specified. Exiting." unless $parameters->{"tempdir"};
die "[Error] The specified temp dir does not exist. Exiting." unless -d $parameters->{"tempdir"};

# Check stuffs:
#&check_parameters($parameters);

# Open config file:
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
my $GENCODE_data = GENCODE->new($parameters);
my $AddScore = Scoring->new($parameters);

# Open files:
open (my $INPUT, "<", $inputFile) or die "[Error] Input file ($inputFile) could not be opened. Exiting.";
open (my $SNPfile, ">", $outputFile."_variant_file.txt") or die "[Error] Output file could not be opened.";
open (my $genotypeFile, ">", $outputFile."_genotype_file.txt") or die "[Error] Output genotype file could not be opened.";
open (my $SNPinfo, ">", $outputFile."_SNPinfo_file.txt") or die "[Error] Output SNPinfo file could not be opened.";

# Processing the input file gene by gene:
# looping through all the genes in the list:
my $gene_count = 0;
while ( my $ID = <$INPUT> ){
    chomp $ID;

    # At first we have to extract all the associated genomic coordinates for the gene.
    my ($chr, $start, $end, $stable_ID, $name, $CollapsedBed);

    # If the input is not a region a few extra steps will be taken:
    unless ($ID =~ /(\d+)_(\d+)-(\d+)/i){
        ($chr, $start, $end, $stable_ID, $name) = $GENCODE_data->GetCoordinates($ID);
	#print "DATA FROM GENCODE FOR $ID: $chr, $start, $end, $stable_ID, $name" if $verbose;
        # Skipping genes that were not found in the GENCODE dataset.
        if ($start eq "NA") {
            print "[Warning] Gene $ID was not found in the GENCODE data. Is it a valid gene name? This gene will be skipped! [NO_GENE]";
            next;
        }

        print "[Info] Queried gene: $name (Ensembl ID: $stable_ID), Genomic location: $chr:$start-$end (Input: $ID)";

        my $bedlines = &BedToolsQuery($chr, $start, $end, $stable_ID, $parameters->{"Linked_features"});
        $CollapsedBed = &FilterLines($bedlines, $stable_ID, $parameters);

        # This should never be a problem, but still be tested:
        unless ( $CollapsedBed ){
            print "[Error] Gene $name did not yield any regions. Skipped. [NO_REGION].";
            next;
        }
    }

    # If the submitted input is a genomic region, we have to do something else:
    else {
        ($chr, $start, $end) = $ID =~ /(\d+)_(\d+)-(\d+)/i;
        $CollapsedBed = join("\t", $chr, $start, $end);
        $name = $ID;
        printf "\n[Info] Queried region: %s:%s-%s", $chr, $start, $end;
    }

    # Once we have the genomic regions, the overlapping variants have to be extracted:
    my $variants = &getVariants($CollapsedBed, $parameters->{"vcfFile"},$parameters->{"chr_prefix"});

    # Filtering variants based on the provided parameters:
    my ($hash, $genotypes) = &processVar($variants, $parameters);
    # printf STDERR "%s\t%s\t%s", $gene_count, total_size($hash)/1024, total_size($AddScore)/1024; # Debug line.
    # Gene will be skipped if there are no suitable variations left:
    if (scalar keys %{$hash} < 2){
        print "[Warning] Gene $ID is skipped as not enough variants left to test [NOT_ENOUGH_VAR].";
        undef $hash;
        undef $genotypes;
        next;
    }

    # TODO max number of variants per gene as input parameter
    # The gene will be skipped if there are too many variants:
    if (scalar keys %{$hash} > $parameters->{"maxVars"}){
        print "[Warning] Gene $ID is skipped as more than 1000 variants are in the set [TOO_MANY_VAR].";
        undef $hash;
        undef $genotypes;
        next;
    }

    # If we want we can add scores:
    $hash = $AddScore->AddScore($hash) if $parameters->{"score"} ne "NA";

    if ($parameters->{"score"} ne "NA" && scalar keys %{$hash} <1){
		print "[Warning] Gene $ID is skipped as no variants remaining post-scoring [NO_VAR_REMAIN].";
	}

    # We don't save anything unless there at least two variants:
    next unless scalar keys %{$hash} > 1;

    # Once we have the scores we have to print out the SNP file:
    &print_SNPlist($hash, $ID, $SNPfile);

    &print_SNP_info($hash, $ID, $SNPinfo, $gene_count, $parameters->{"build"});
    &print_genotypes($genotypes, $genotypeFile, $parameters, $gene_count);

    $gene_count ++; # So the header will only be printed once.

    # Some diagnostic functionality:
    #print Dumper $hash;  # Look at the snp container
    #print; # Memory usage of the main variables

}

close $inputFile;
close $SNPfile;
close $genotypeFile;
close $SNPinfo;

sub check_parameters {
    my $parameters = shift;

    # checking build:
    die "[Error] Only GRCh37 and GRCh38 genome builds are supported (--build 37 or 38)." unless $parameters->{"build"} eq "37" || $parameters->{"build"} eq "38";

    # Checking MAF (a float below 1):
    my $maf=$parameters->{"MAF"};
    die "[Error] $maf must be a float number < 1" unless defined($maf) && is_float($maf) && $maf<1;

    # Checking MAC (a number):
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
                      "Eigen" => 1,
                      "EigenPC" => 1,
                      "EigenPhred" => 1,
                      "EigenPCPhred" => 1,
                      "Linsight" => 1,
                      "Mixed" => 1);

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
            print "[Warning] Eigen scores as weight cannot be used. No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
        # elsif (! -e $params->{"EigenPath"}){
        #     printf "[Warning] The specified genome-wide Eigen scores are not available %s.", $params->{"EigenPath"};
        #     print "[Warning] Eigen scores as weight cannot be used. No weights will be applied.";
        #     $params->{"score"} = "NA";
        #     return $params;
        # }

    }
    elsif ( $params->{"score"} eq "Linsight"){
        if (! exists $params->{"Linsight"}) {
            print "[Warning] To use linsight scores the config file has to have 'Linsight' entry.";
            print "[Warning] No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
        elsif ( ! -e $params->{"Linsight"}){
            print "[Warning] Linsight scores as weight cannot be used. No weights will be applied.";
            $params->{"score"} = "NA";
            return $params;
        }
    }
    return $params;
}


sub readConfigFile {
    my $params = $_[0];

#    printf "[Info] Config file: %s", $params->{"configFileName"};

    # Reading file:
    open(my $CONF, "<", $params->{"configFileName"}) or die "[Error] Config file could not be oppended. Exiting.";

    # Reading file:
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
            printf STDERR "[Error] The provided GENCODE feature name is not supported: %s. Use these: %s", $feature, join(", ", keys %AcceptedFeatures);
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
                printf STDERR "[Error] The provided regulatory feature name is not supported: %s. Use these: %s", $feature, join(", ", keys %AcceptedFeatures);
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

sub BedToolsQuery {
    my ($chr, $start, $end, $stable_ID, $geneBedFile) = @_;
    my $queryString = sprintf("intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\t%s\") -b %s -sorted | cut -f9-",$chr, $start, $end, $stable_ID, $geneBedFile);
    
    print "[Info] IntersectBed query string: $queryString" if $verbose;
    my $query =Scoring::backticks_bash($queryString);
    #print("RETURN FROM BEDTOOLSQUERY");
    return $query;
}

sub formatLines {
    my %hash = %{$_[0]};
    my $ext = $_[1] // 0;
    return sprintf("%s\t%s\t%s\t%s", $hash{"chr"}, $hash{"start"} - $ext, $hash{"end"} + $ext, encode_json(\%hash))
}

sub getVariants {
    my ($merged, $inputvcfpattern,$prefix) = @_;
    my $distance = 0;

    # Finding out which chromosome are we on:
    my ($chr) = $merged =~ /^(.+?)\t/;  # WATCHOUT !

    print("GETVARIANTS CHR=$chr") if $verbose;
    
    # Print info:
    print  "\n[Info] Extracting variants from vcf files:" if $verbose;
    (my $vcfFile = $inputvcfpattern ) =~ s/\%/$chr/g;
    return "" unless -e $vcfFile;
    my $bcftoos_query = sprintf("tabix %s ", $vcfFile);

    # looping through all lines:
    foreach my $line (split("\n", $merged)){
        my ($chr, $start, $end) = split("\t", $line);
        $distance += $end - $start;
        $bcftoos_query .= sprintf(" %s%s:%s-%s", $prefix, $chr, $start, $end);
    }

    print  "$bcftools_query" if $verbose;
    my $variants = Scoring::backticks_bash($bcftools_query);

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

        # Decoding json line:
        my %annot_hash = %{decode_json($line)};

        # Filtering for lines which belong to our target gene:
        next unless $annot_hash{"gene_ID"} eq $stable_ID;

        my $source = $annot_hash{"source"};
        my $class =  exists $annot_hash{"class"} ? $annot_hash{"class"} :  "?";
        # print "\n$source $class";

        if ($source eq "GENCODE" and exists $GENCODE{$class}) {

            # If the user has specified, minor transcripts will be excluded:
            next if exists $GENCODE{"minor"} && $annot_hash{"appris"} == "Minor";

            push (@output_lines, formatLines(\%annot_hash, $parameters->{'extend'}));
            $hash{GENCODE}{$class} ++;
        }
        elsif (($source eq "GTEx" and exists $GTEx{$class})
               or ($source eq "GTEx" and exists $GTEx{'allreg'})){

            push (@output_lines, formatLines(\%annot_hash));
            $hash{GTEx}{$class} ++;
        }
        elsif (($source eq "overlap" and exists $overlap{$class})
               or ($source eq "overlap" and exists $overlap{'allreg'})){

            push (@output_lines, formatLines(\%annot_hash));
            $hash{overlap}{$class} ++;

        }
        # It is easy to expand this design.
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
    print "[Info] Selected lines:\n", join("\n", @output_lines),"" if $verbose;

    # Saving temporary bedfile:
    open (my $tempbed, "> filtered_regions.bed");
    foreach my $line (@output_lines){
        print $tempbed $line;
    }
    `sort -k1,1 -k2,2n filtered_regions.bed | sponge filtered_regions.bed`;
    close $tempbed;

    # Collapsing overlapping features:
    my $queryString = "mergeBed -i filtered_regions.bed";
    #my $merged = `bash -O extglob -c \'$queryString\'`;
    my $merged = Scoring::backticks_bash($queryString);
    return $merged;
}

sub processVar {
    my $variants   = $_[0]; # List of all overlapping variants
    my $parameters = $_[1]; # All the submitted parameters to filter variants.

    my %hash = (); # SNP info container.
    my %genotypeContainer = (); # genotype information container.

    # Which build are we using?
    my $build = "GRCh".$parameters->{"build"};

    print "[Info] Filtering variants:" if $verbose;

    my @total_vars = split("\n", $variants);
    printf "[Info] Total number of overlapping variants: %s\n", scalar(@total_vars) if $verbose;

    # Looping through all variants that satisfy the submitted criteria:
    while (@total_vars){
        # Removing one variant at a time:
        my $variant = shift @total_vars;

        # line: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EGAN00001033155
        my ($chr, $pos, $id, $a1, $a2, $qual, $filter, $info, $format, @genotypes) = split(/\t/, $variant);

        # Generating variant name (Sometimes the long allele names cause problems):
        my $short_a1 = length $a1 > 5 ? substr($a1,0,4) : $a1;
        my $short_a2 = length $a2 > 5 ? substr($a2,0,4) : $a2;

        my $SNPID = sprintf("%s_%s_%s_%s", $chr, $pos, $short_a1, $short_a2);

        # Parsing info field for relevant information:
        (my $ac )= $info =~ /AC=(.+?)[;\b]/;
        (my $an )= $info =~ /AN=(.+?)[;\b]/;
        (my $consequence ) = $info =~ /consequence=(.+?)[;\b]/;
        # If AN and AC values are not found we skip the variant:
	if (! $an){
	    print "[Warning] $SNPID has no AN; skipping";
	    next;
	}

	if (! $ac){
	    print "[Warning] $SNPID has no AC; skipping";
	    next;
	}
	
        # We add NA if consequence was not found:
        $consequence = "NA" unless $consequence;

        # We don't consider multialleleic sites this time.
        if ( $a2 =~ /,/){
            print  "[Warning] $SNPID will be omitted because it's multiallelic! ($a2).";
            next;
        }

        print "WARNING: $variant has no AC" if $ac eq "NA";
        print "WARNING: $variant has no AN" if $an eq "NA";

        # Calculating values for filtering:
        my $missingness = (scalar(@genotypes)*2 - $an)/(scalar(@genotypes)*2);
        my $MAF = $ac/$an;

        # This flag shows if the non-ref allele is the major one:
        my $genotypeFlip = 0; #
        if ( $MAF > 0.5 ){
            print "[Info] MAF of $SNPID is $MAF is greater then 0.5, genotype is flipped.";
            $MAF = 1 - $MAF;
            $genotypeFlip = 1;
        }

        my $MAC = $ac;
        $MAC = $an - $ac if $MAC > $an / 2;

        # Looking at the memory usage:
        # print STDERR total_size(\%hash), "\n";

        # We don't consider indels if weights are used:
        if (( length($a2) > 1 or length($a1) > 1 ) && $parameters->{"score"} ne "NA"){
            print  "[Warning] $SNPID will be omitted because it' an indel! ($a1/$a2).";
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
            printf "[Warning] $SNPID will be omitted because of consequence is not lof (%s)\n", $consequence;
            next;
        }
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
        $hash{$SNPID}{"alleles"} = [$a1, $a2];
        $hash{$SNPID}{$build} = [$chr, $pos - 1, $pos];
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

    printf  ( "[Info] Number of filtered variants: %s\n", scalar(keys %hash));
    return (\%hash, \%genotypeContainer);
}

sub print_SNPlist {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];
    my $outputhandle = $_[2];

    # Get the list of variants:
    my @snpIDs = keys %hash;

    # Check if we have scores as well:
    my $flag = $hash{$snpIDs[0]}{"score"} ? 1 : 0;

    # The line with the snps will be written any ways:
    print $outputhandle "$gene_name\t$flag\t", join("\t", @snpIDs),"";

    # We return if the flag is 0, so there is no scores given:
    return 1 if $flag == 0;

    # If scores are given, we add an extra line:
    print $outputhandle "$gene_name\t0\t";
    my @scores = ();
    foreach my $snpid (keys %hash){
        # print  $outputhandle "$hash{$snpid}{score}\t"
        # push (@scores, $hash{$snpid}{"score"}{"processed_eigen"}) if $score eq "Eigen";
        push (@scores, $hash{$snpid}{"score"});
    }
    print $outputhandle join("\t", @scores),"";

    return 1;
}

sub print_genotypes {

    my %genotype      = %{$_[0]};
    my $outputhandler = $_[1];
    my $parameters    = $_[2];
    my $gene_counter  = $_[3];

    my $vcfChrFile = $parameters->{"vcfFile"};
    $vcfChrFile =~ s/\%/11/g; # TODO ???

    # Checking if the file exists or not:
    if (! -e $vcfChrFile ){
        print STDERR "[Error] The vcf file does not exists. $vcfChrFile";
    }

    # Assembling header for the first gene:
    if ( $gene_counter == 0){
        # Get list of sample IDs:
        my $samples = `zgrep -m1 "#CHROM"  $vcfChrFile | cut -f10-`;
        print $outputhandler "0\t$samples";
    }

    # Saving data:
    for my $var (keys %genotype){
        print $outputhandler "$var\t", join("\t", @{$genotype{$var}});
    }
}

sub print_SNP_info {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];
    my $outfilehandle = $_[2];
    my $gene_counter = $_[3];
    my $build = $_[4];

    # Extract variant names:
    my @variants = keys %hash;

    # Assembling header for the first gene:
    if ( $gene_counter == 0){
        my $header = "Gene_name\tSNPID(GRCh${build})\trsID\tAlleles\tconsequences\tMAF\tAC\tmissigness";
        $header .= "\tScore" if exists $hash{$variants[0]}{"score"};
        print $outfilehandle $header;
    }

    # printing out info for each gene and variant:
    foreach my $variant (values %hash){

        # generating all the required fields:
        my $alleleString = join("/", @{$variant->{'alleles'}});
        my $rsID = $variant->{'rsID'};
        my $MAF = sprintf "%.4f", $variant->{'frequencies'}->[2];
        my $missingness = sprintf "%.4f", $variant->{'missingness'};
        my $weight = $variant->{'score'} if exists $variant->{"score"};
        my $consequence = $variant->{'consequence'};
        my $AC = $variant->{'frequencies'}->[0];
        my $SNPID = $variant->{'GRCh'.$build}->[0].":".$variant->{'GRCh'.$build}->[2];
        my $line = join ("\t", $gene_name, $SNPID, $rsID, $alleleString, $consequence, $MAF, $AC, $missingness);
        $line .= "\t$weight" if exists $variant->{"score"};
        print $outfilehandle $line;
    }

}
###

sub usage {
    print "Usage::";
    print("      Required:");
    print("          --input <input file>");
    print("          --output <output prefix>");
    print("          --vcfFile <input VCF>");
    print("          --gencode-file <GENCODE file>");
    print("          --tempdir <TEMP DIR>");
    print("          --configFile <config file>");
    print("      Optional:");    
#    print("          --build <genome build; default: 38>");
    print("          --GENCODE <comma separated list of GENCODE features (gene, exon, transcript, CDS or UTR)>");
    print("          --GTEx <comma separated list of GTEx features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)>");
    print("          --overlap <comma separated list of overlap features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)>");
    print("          --extend <by how many basepairs the GENCODE features should be extended.>");
    print("          --MAF <MAF upper threshold>");
    print("          --MAC <MAC lower threshold>");
    print("          --maxVars <max number of variants in a gene (default: 1000)>");
    print("          --SkipMinor <skip minor transcripts by APPRIS>");
    print("          --verbose <increase verbosity>");
    print("          --score <which score to use to weight variants; one of: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred, Linsight, Mixed >");
    print("          --lof <only select high impact variants: transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification, inframe_insertion, inframe_deletion, splice_region_variant>");
    print("          --loftee <only select high and low confident loss of function variants>");
    print("          --lofteeHC <only select high confident loss of function variants>");
    print("          --missingness <missingness upper threshold>");
    print("          --shift <shift scores by this value>");
    print("          --cutoff <score threshold below which the variants will be removed>");
    print("          --floor <scores below this threshold will be set to this value>");
    print("          --chromosome-prefix <chromosome prefix in VCF files; default: \"chr\">");
    print("          --help <this help>");
}

sub checkVCFs{
    my $inputvcfpattern=$_;

    for (1 .. 22){
	my $fname=$inputvcfpattern;
	$fname=~s/\%/$_/;
	return 1 if -e $fname;
    }

    return undef;
}
