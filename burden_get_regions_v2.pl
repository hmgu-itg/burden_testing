#!/usr/bin/env perl


# Version information:
our $version = "v3.3 Last modified: 2017.01.05";

# Accepted feature names:
    # GENCODE: gene, exon, transcript, CDS, UTR
    # GTEx: promoter, enhancer, promoterFlank, openChrom, TF_bind, allreg
    # overlap: promoter, enhancer, promoterFlank, openChrom, TF_bind, allreg

# Only a handful packages were used:
use strict;
use warnings;
use Data::Dumper; # used for diagnostic purposes.
use Pod::Usage qw(pod2usage); # Used for providing useful error/warning and user messages.
use JSON;
use DateTime;
use File::Basename;
use HTTP::Tiny; # used to query Ensembl to get consequences.
use Getopt::Long qw(GetOptions);

# Status report:
print "[Info] The script was called with the following parameters:\n", join(" ", $0, @ARGV), "\n";
print "[Info] Script version: $version\n";
printf "[Info] Run date: %s\n", DateTime->now->strftime("%Y. %b %d %H:%M");

# Storing script directory:
our $scriptdir = dirname(__FILE__);

# Initializing command line parameters with default values if required:
our $score = "NA";
my $shift = 0;
my $k = 50; # Weighting constant.
my $extend = 0;
my $MAF = 0.05; # 5%
my $MAC = 0;
my $missingthreshold = 0.001; # 0.1%
my $configFileName = "config.txt"; #
my $floor = "NA"; # Variants will not be removed below the threshold, but assigned with this value
my $cutoff = 0; # That is the default value, below which every variant will be removed.


# Command line options with no default values:
my ($inputFile, $outputFile, $GENCODE, $GTEx, $verbose, $overlap, $help,
    $minor, $tissue_list, $MAF_weight);

# A command line switch. If turned on, only loss of function variants will be included in the test.
our ($lof, $loftee);

GetOptions(
    # Input/Output:
    'input|i=s'  => \$inputFile,
    'output|o=s' => \$outputFile,

    # select source (arguments are the class of the features.):
    "GENCODE|G=s" => \$GENCODE,
    "GTEx|E=s"    => \$GTEx,
    "overlap|L=s" => \$overlap,

    # Extend regions with a defined length:
    "extend=s" => \$extend,

    # Appris:
    "SkipMinor" => \$minor,

    # Only considering regulatory elements associated in the
    # following tissues (this only filters based on ensembl tissues):
    "tissues|T=s" => \$tissue_list,

    # Variant features:
    'MAF=s' => \$MAF,
    'MAC=s' => \$MAC,

    # Verbose output:
    'verbose|v' => \$verbose,

    # specifying config file:
    'conf|c=s' => \$configFileName,

    # Which score we need:
    'score|s=s' => \$score,

    # A flag to weight by MAF:
    'MAFweight|w'  => \$MAF_weight,
    'ExpConst|k=s' => \$k,

    # Do we need only loss of function:
    'lof' => \$lof,
    'loftee' => \$loftee, # Filters in only high confident loss of function variants

    # Accepting missingness filter:
    'missingness|m=s' => \$missingthreshold,

    # Changing scores that will be used for weights:
    'shift=s'  => \$shift,  # The value used to shift eigen scores:
    'cutoff=s' => \$cutoff, # hard threshold that will be applied on scores:
    'floor=s'  => \$floor,  # How do we want to floor Eigen values:

    # Asking for help:
    'help|h' => \$help
);

# Open config file:
my %ConfFiles = %{&read_param($configFileName)};

our $geneBedFile      = $ConfFiles{geneBedFile};
our $vcfFile          = $ConfFiles{vcfFile};
our $liftoverPath     = $ConfFiles{liftoverPath};
our $EigenPath        = $ConfFiles{EigenPath};
our $caddPath         = $ConfFiles{caddPath};
our $temporaryBedFile = $ConfFiles{temporaryBedFile};
our $GENCODEFile      = $ConfFiles{GENCODEFile};
our $GWAVA_DIR        = $ConfFiles{GWAVA_DIR};
our $bigWigToolDir    = $ConfFiles{bigWigTools};
our $linsigthFile     = $ConfFiles{Linsight};

# All files must exist otherwise the scrip quits:
foreach my $sourcefile ($geneBedFile, $vcfFile, $liftoverPath, $EigenPath, $caddPath){
    my $file = sprintf($sourcefile, "chr12");
    pod2usage({-verbose => 99, -message => "[Error] One of the requested file does not exits ($file). Exiting.\n", -sections => "FILES|SYNOPSIS" }) unless ( -e $file);
}

# Checking essential parameters:
pod2usage({-verbose => 2}) if $help;
pod2usage({-message => "[Error] Input file has to be specified with the -i switch!\n", -verbose => 0}) unless $inputFile;
pod2usage({-message => "[Error] Output file prefix has to be specified with the -o switch!\n", -verbose => 0}) unless $outputFile;

# Parsing command line arguments:
our %GENCODE = %{&parseGENCODE($GENCODE)} if $GENCODE;
our %GTEx    = %{&parseRegulation($GTEx)} if $GTEx;
our %overlap = %{&parseRegulation($overlap)} if $overlap;
our %Variant = ('MAF' => $MAF, 'MAC' => $MAC);

# Optional parameters regarding the genes:
$GENCODE{'minor'}  = 1 if $minor;
$GENCODE{'extend'} = $extend;

# Reporting score shifts:
unless ($score eq "NA"){
    my $tag = $floor eq "NA" ? "removed." : "assigned to $floor.";
    printf "[Info] %s scores are shifted by %s.\n", $score, $shift;
    printf "[Info] Variants with score below %s, will be %s.\n", $cutoff, $tag;
}

# Missingness:
printf  "[Info] Missingness threshold: %s (variants above %s%% missing genotypes will be excluded).\n", $missingthreshold, $missingthreshold * 100;

# Reading genomic coordinates:
our ($GENCODE_name, $GENCODE_ID) = &read_GENCODE($scriptdir."/gencode_genes_V25.tsv.gz");

# The submitted scoring method might not be supported! Check for it now!
my %acceptedScores = ("GWAVA" => 1,
                      "CADD" => 1,
                      "Eigen" => 1,
                      "EigenPC" => 1,
                      "EigenPhred" => 1,
                      "EigenPCPhred" => 1,
                      "Linsight" => 1,
                      "Mixed" => 1,
                      "NA" => 1);
unless ( exists $acceptedScores{$score}){
    pod2usage({-message => sprintf('[Warning] %s is not a supported scoring method.
[Warning] Implemented scoring methods: %s.
[Warning] No weights will be used.\n', $score, join(",", keys %acceptedScores)), -verbose => 0, -exitval => "NOEXIT"})
}

# This is the list of those consequences that will be retained upon switchin --lof
our %lof_cons = (
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
    "missense"                 => 1,
    "splice_region_variant"    => 1,
    "coding_sequence_variant"  => 1
);

###
### Reading file, processing the list of genes, line by line.
###
open(my $variant_output, ">", $outputFile."_variants") or die "[Error] Output file ($outputFile\_variants) could not opened.\n";

my $genotypes = {};
my $ID = $inputFile;
chomp $ID;

#########
######### The following block will be done for the gene.
#########
my ($chr, $start, $end, $stable_ID, $name, $CollapsedBed);
# If the input is not a region a few extra steps will be taken:
unless ($ID =~ /chr(\d+)_(\d+)-(\d+)/i){
    ($chr, $start, $end, $stable_ID, $name) = &GetCoordinates($ID);

    # if a false gene name or ID was given, we record it and kill the process:
    unless ($name){
        printf STDERR "[Error] Gene %s was not found in the GENCODE data! Exiting.\n    ", $ID;
        exit;
    }

    # Removing non-alphanumeric characters from gene names:
    $name =~ s/[^0-9a-z]//gi;

    # Skipping genes that were not found in the GENCODE dataset.
    if ($start eq "NA") {
        print STDERR "[Warning] Gene $ID was not found in the GENCODE data. Is it a walid gene name? This gene will be skipped!\n";
        exit;
    }

    print "\n\n[Info] Queried gene: $name (Ensembl ID: $stable_ID), Genomic location: chr$chr:$start-$end (Input: $ID)\n";
    my $bedlines = &BedToolsQuery($chr, $start, $end, $stable_ID);

    $CollapsedBed = &FilterLines($bedlines, $stable_ID);
}
# If the
else {
    ($chr, $start, $end) = $ID =~ /(chr\d+)_(\d+)-(\d+)/i;
    $CollapsedBed = join("\t", $chr, $start, $end);
    $name = $ID;
    printf "\n\n[Info] Queried region: %s:%s-%s\n", $chr, $start, $end;
}
#########
#########
#########

# Retrieve overlapping variations:
die "[Error] No genomic region was selected. Exiting.\n" unless $CollapsedBed;
my $variants = &GetVariants($CollapsedBed);

# Process variant list:
my $hash;
($hash, $genotypes) = &processVar($variants, $genotypes);

# We exit if there is not variant for the gene:
die "[Warning] Less than two variants remained in the gene, will be skipped.\n" unless scalar(keys %{$hash}) > 1;

# If we apply scores, coordinates have to lifted over:
if ( $score ne "NA" ) {

    # Saving temporary bedfile for liftover:
    my $flag;
    ($hash, $flag) = &liftover($hash, $name);

    # Exit if there is no lifted over variant:
    die "[Warning] None of the variants of this gene was lifted over! So no scoring is doable. Check if liftOver works, and the filtering parameters.\n" if ($flag == 0);

    # Returning Eigen scores if that's the case:
    $hash = &get_Eigen_Score($hash) if $score =~ /Eigen/i;

    # Returning CADD scores if that's the case:
    $hash = &get_CADD_GERP($hash, "CADD") if $score eq "CADD";

    # Returning Insight scores if that's required:
    $hash = &get_linsight($hash, $name) if $score eq "Linsight";

    # Returning GWAVA scores if that's the case:
    $hash = &get_GWAVA($hash, $name) if $score eq "GWAVA";

    # Returning mixed Eigen/CADD scores if that's the case:
    $hash = &get_mixed($hash, $name) if $score eq "Mixed";

    # Processing scores based on the flooring method and the cutoff:
    $hash = &process_score($hash, $floor, $shift, $cutoff);
}

# Saving SNP data into file:
&print_SNPlist($hash, $variant_output, $name);
&print_SNP_info($hash, $name);

# Diagnostic dumper of the hash:
# print Dumper $hash;

# print "Genotypes:", Dumper ($genotypes);
open(my $genotype_output, ">", $outputFile."_genotype") or die "[Error] Output file ($outputFile\_genotype) could not opened.\n";
&print_genotypes($genotypes, $genotype_output);

##
## Saving SNP file with consequences and scores
##
sub print_SNP_info {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];

    # print Dumper %hash;## Diagnostic line: what do we have sotred in the hash at the last stage?
    # Looping through all the variants and save the data:
    open(my $outfilehandle, ">", $outputFile."_SNP_info") or die "[Error] Output file ($outputFile\_SNP_info) could not opened.\n";

    # Printing out header:
    print $outfilehandle "gene_name\tSNPID\trsID\tallele_string\tconsequence\tMAF\tmissingness\tscore\n";

    foreach my $variant (values %hash){

        # generating all the required fields:
        my $alleleString = join("/", @{$variant->{'alleles'}});
        my $rsID = $variant->{'rsID'};
        my $MAF = sprintf "%.4f", $variant->{'frequencies'}->[2];
        my $missingness = sprintf "%.4f", $variant->{'missingness'};
        my $weight = exists $variant->{'score'} ? $variant->{'score'} : "NA";
        my $consequence = $variant->{'consequence'};
        my $SNPID = $variant->{'GRCh38'}->[0].":".$variant->{'GRCh38'}->[1];

        $SNPID .= "(".$variant->{'GRCh37'}->[0].":".$variant->{'GRCh37'}->[2]."[b37])" if exists $variant->{'GRCh37'};


        print $outfilehandle join ("\t", $gene_name, $SNPID, $rsID, $alleleString, $consequence, $MAF, $missingness, $weight), "\n";
    }

}

##
## Function to parsing input parameters
##
sub parseGENCODE {
    my %AcceptedFeatures = ( "gene" => 1, "exon" => 1, "transcript" => 1, "CDS" => 1, "UTR" => 1);
    my $gencodeString = $_[0];
    my %hash = ();
    foreach my $feature (split(",", $gencodeString)){

        # Checking if the provided feature is exists or not:
        if (exists $AcceptedFeatures{$feature}) {
            $hash{$feature} = 1;
        }
        else {
            printf STDERR "[Error] The provided GENCODE feature name is not supported: %s. Use these: %s\n", $feature, join(", ", keys %AcceptedFeatures);
        }


    }
    return \%hash
}
sub parseRegulation {
    # Accepted features: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg
    my %AcceptedFeatures = ( "promoter" => 1, "CTCF" => 1, "enhancer" => 1, "promoterFlank" => 1, "openChrom" => 1, "TF_bind" => 1, "allreg" => 1 );
    my $regString = $_[0];
    my %hash = ();
    foreach my $feature (split(",", $regString)){

        unless (exists $AcceptedFeatures{$feature}) {
            printf STDERR "[Error] The provided regulatory feature name is not supported: %s. Use these: %s\n", $feature, join(", ", keys %AcceptedFeatures);
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
    return \%hash
}

##
## Saving genotype information:
##
sub print_genotypes {
    my %genotype = %{$_[0]};
    my $outputhandler = $_[1];
    my $vcfChrFile = sprintf($vcfFile, "chr11");

    # Get list of sample IDs:
    my $samples = `zgrep -m1 "#CHROM"  $vcfChrFile | cut -f10-`;

    # Saving data:
    print $outputhandler "0\t$samples";
    for my $var (keys %genotype){
        print $outputhandler "$var\t", join("\t", @{$genotype{$var}}), "\n";
    }
}

##
## Weight scores with MAF. (lower MAF higher weight.)
##
sub weight_score {
    my %hash = %{$_[0]};
    my $k = $_[1];

    foreach my $snpID (keys %hash){
        my $AC = $hash{$snpID}{'frequencies'}[0];
        my $AN = $hash{$snpID}{'frequencies'}[1];
        my $MAF = $hash{$snpID}{'frequencies'}[2];
        $hash{$snpID}{score} = $hash{$snpID}{score} * exp( - ($AC - 1)/$AN * $k);
    }
    print "[Info]";
    return \%hash;
}

##
## Saving variants in SNP list:
##
sub print_SNPlist {
    my %hash = %{$_[0]};
    my $outputhandle = $_[1];
    my $gene_name = $_[2];

    # Get the list of variants:
    my @snpIDs = keys %hash;

    # Check if we have scores as well:
    my $flag = $hash{$snpIDs[0]}{score} ? 1 : 0;

    # The line with the snps will be written any ways:
    print $outputhandle "$gene_name\t$flag\t", join("\t", @snpIDs),"\n";

    # We return if the flag is 0, so there is no scores given:
    return 1 if $flag == 0;

    # If scores are given, we add an extra line:
    print $outputhandle "$gene_name\t0\t";
    my @scores = ();
    foreach my $snpid (keys %hash){
        # print  $outputhandle "$hash{$snpid}{score}\t"
        # push (@scores, $hash{$snpid}{"score"}{"processed_eigen"}) if $score eq "Eigen";
        push (@scores, $hash{$snpid}{"score"}) if $score;
    }
    print $outputhandle join("\t", @scores),"\n";

    return 1;
}


##
## Flooring Eigen scores:
##
sub process_score {

    my %hash   = %{$_[0]}; # hash with data
    my $floor  = $_[1]; # What to do with values above the cutoff. By default, leave as it is.
    my $shift  = $_[2]; # With how much we should push eigen scores.
    my $cutoff = $_[3]; # Below which variants will be excluded.

    printf "[Info] Processing %s scores... ", $score if $verbose;

    # Looping through all variants and modify scores:
    foreach my $var ( keys %hash){

        # Shifting scores:
        $hash{$var}{"score"} = $hash{$var}{"score"} + $shift;

        # How to deal with variants that are below the cutoff:
        if ($floor ne "NA" && $hash{$var}{"score"} <= $cutoff ) {
            $hash{$var}{"score"} = $floor;
            printf "[Warning] score of %s is set to %s, because %s score is below threshold: %s!\n", $var, $floor, $score, $hash{$var}{"score"};
        }
        elsif ($floor eq "NA" && $hash{$var}{"score"} <= $cutoff){
            printf "[Warning] %s is deleted, because %s score is below threshold: %s!\n", $var, $score, $hash{$var}{"score"};
            delete $hash{$var};
        }
        # If the score is above cutoff, we don't do anything.
    }
    print "Done.\n" if $verbose;
    return \%hash;
}


##
## Get Eigen scores
##
sub get_Eigen_Score {
    my %hash = %{$_[0]};

    # just sayin'
    print  "[Info] Adding Eigen scores for variants.\n" if $verbose;

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){
        my $EigenFile = $EigenPath;
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;

        # Two tabix queries will be submitted regardless of the output...
        my $tabix_query = sprintf("tabix %s %s:%s-%s | grep %s", $EigenFile, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2], $hash{$var}{alleles}[1]);
        print "$tabix_query\n" if $verbose;
        my $lines = `bash -O extglob -c \'$tabix_query\'`;
        $hash{$var}{"score"} = "NA"; # Initialize Eigen score.

        foreach my $line (split("\n", $lines)){
            chomp $line;
            # chr     pos     a1      a2      Eigen   EigenPC EigenPhred      EigenPCPhred
            my @array = split(" ", $line);

            # Testing alleles:
            if (($array[2] eq $hash{$var}{alleles}[0] and $array[3] eq $hash{$var}{alleles}[1]) or
                ($array[2] eq $hash{$var}{alleles}[1] and $array[3] eq $hash{$var}{alleles}[0])) {

                $hash{$var}{"score"} = $array[4] || "NA" if $score eq "Eigen"; # Un-scaled raw Eigen score.
                $hash{$var}{"score"} = $array[5] || "NA" if $score eq "EigenPC"; # Un-scaled Eigen PC score.
                $hash{$var}{"score"} = $array[6] || "NA" if $score eq "EigenPhred"; # Phred scaled Eigen score.
                $hash{$var}{"score"} = $array[7] || "NA" if $score eq "EigenPCPhred"; # Phred scaled Eigen PC score.
            }
        }

        # Reporting if no score has been found for the variant:
        if ( $hash{$var}{"score"} eq "NA") {
            printf ( "[Warning] %s score was not found for variant %s in the non-coding set! Removing variant.\n", $score, $var);
            delete $hash{$var};
        }
    }

    printf "[Info] Eigen scores have been added to variants (Number of variants: %s).\n\n", scalar keys %hash if $verbose;
    return \%hash
}


###
### Retrieve and assign linsight scores to variants:
###
sub get_linsight {
    my %hash = %{$_[0]};
    my $gene_name = $_[1];

    # just sayin'
    print  "[Info] Adding Linsight scores for variants.\n" if $verbose;

    # Calling BigWig tools to get the scores:
    #./bigWigAverageOverBed /lustre/scratch115/projects/t144_helic_15x/analysis/HA/weights/LINSIGHT/LINSIGHT.bw APOC3_GRCh37.bed  out.tab
    my $linsightOut = sprintf("%s_GRCh37_linsight.tab", $gene_name);
    my $bigwigQuery = sprintf("%s/bigWigAverageOverBed %s %s_GRCh37.bed %s", $bigWigToolDir, $linsigthFile, $gene_name, $linsightOut);
    print $bigwigQuery,"\n";
    `$bigwigQuery`;

    # Reading linsight scores from the output file:
    open( my $LIN,"<",$linsightOut) or warn "[Warning] Linsight scores were not generated. Check bigWitTools!\n";
    my %linsightScore = ();
    while ( my $line = <$LIN>) {
        chomp $line;
        my @a = split("\t", $line);
        $linsightScore{$a[0]} = $a[3];
    }

    # Adding linsight scores to the hash:
    foreach my $var (keys %hash){

        if (exists $linsightScore{$var}) {
            $hash{$var}{"score"} = $linsightScore{$var};
        }
        else {
            printf ( "[Warning] linsight score was not found for variant %s! Removing variant.\n", $var);
            delete $hash{$var};
        }
    }

    printf "[Info] linsight scores have been added to variants (Number of variants: %s).\n\n", scalar keys %hash if $verbose;
    return \%hash
}

##
## Get CADD scores
##
sub get_mixed {
    my %hash = %{$_[0]};
    my $GeneName = $_[1];

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){

        # Extracting consequence:
        my $consequence = $hash{$var}{'consequence'};
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;

        # Initialize score:
        $hash{$var}{"score"} = "NA";

        # If the consequence is severe, we use cadd:
        if ( exists $lof_cons{$consequence} ){

            printf "[Info] %s is a %s so CADD-Phred scores are used.\n", $var, $consequence;
            # Tabixing for CADD score:
            my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f1-5,25-28,115-", $caddPath, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2]);
            print "[Info] $tabix_query\n";
            my $lines = `bash -O extglob -c \'$tabix_query\'`;

            # Looping through all ouput lines:
            foreach my $line (split("\n", $lines)){

                # Line: 12	56482614	C	C	G	5.52	4.63	510	1.76195e-61	4.931423	25.0
                chomp $line;
                my ($Chrom, $Pos, $ref, $anc, $alt, $GerpN, $GerpS, $GerpRS, $GerpRSpval, $RawScore, $PHRED) = split("\t", $line);

                # Testing alleles:
                if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
                    ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {
                        $hash{$var}{"score"} = $PHRED;
                }
            }
        }
        else {
            my $EigenFile = $EigenPath;
            my $tabix_query = sprintf("tabix %s %s:%s-%s | grep %s", $EigenFile, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2], $hash{$var}{alleles}[1]);
            printf "[Info] %s is a %s so EigenPhred scores are used.\n", $var, $consequence;
            print "$tabix_query\n" if $verbose;
            my $lines = `bash -O extglob -c \'$tabix_query\'`;

            foreach my $line (split("\n", $lines)){
                chomp $line;
                # chr     pos     a1      a2      Eigen   EigenPC EigenPhred      EigenPCPhred
                my @array = split(" ", $line);

                # Testing alleles:
                if (($array[2] eq $hash{$var}{alleles}[0] and $array[3] eq $hash{$var}{alleles}[1]) or
                    ($array[2] eq $hash{$var}{alleles}[1] and $array[3] eq $hash{$var}{alleles}[0])) {

                    $hash{$var}{"score"} = $array[6]; # Phred scaled Eigen score.
                }
            }

        }

        # Reporting if no score has been found for the variant:
        if ( $hash{$var}{"score"} eq "NA") {
            printf ( "[Warning] %s score was not found for variant %s! Removing variant.\n", $score, $var);
            delete $hash{$var};
        }
    }
    print "[Info] Mixed EigenPhred/CADDPhred scores have been added to variants.\n" if $verbose;
    return \%hash
}


##
## Get CADD scores
##
sub get_CADD_GERP {
    my %hash = %{$_[0]};
    my $score = $_[1];

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;
        my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f1-5,25-28,115-", $caddPath, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2]);
        print "[Info] $tabix_query\n";
        my $lines = `bash -O extglob -c \'$tabix_query\'`;
        #$hash{$var}{CADD} = ["NA", "NA"]; # Initialize CADD scores.
        #$hash{$var}{GERP} = ["NA", "NA", "NA", "NA"]; # Initalize GERP scores
        $hash{$var}{"score"}{${score}."_score"} = "NA";

        foreach my $line (split("\n", $lines)){
            # Line: 12	56482614	C	C	G	5.52	4.63	510	1.76195e-61	4.931423	25.0
            chomp $line;
            my ($Chrom, $Pos, $ref, $anc, $alt, $GerpN, $GerpS, $GerpRS, $GerpRSpval, $RawScore, $PHRED) = split("\t", $line);

            # Testing alleles:
            if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
                ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {

                    #$hash{$var}{CADD} = [$RawScore, $PHRED]; # Initialize CADD scores.
                    #$hash{$var}{GERP} = [$GerpN, $GerpS, $GerpRS, $GerpRSpval]; # Initalize GERP scores
                    $hash{$var}{"score"} = $PHRED if $score eq "CADD";
                    $hash{$var}{"score"} = $GerpS if $score eq "GERP";
            }
        }
    }
    print "[Info] $score scores have been added to variants.\n" if $verbose;
    return \%hash
}
##
## Get GWAVA scores as weights:
##
sub get_GWAVA {

    my %hash = %{$_[0]};
    my $gene_name = $_[1];

    # Report progress:
    print "[Info] Run GWAVA annotation... ";

    my %gwava_scores = ();

    # before running gwava we have to check if the directory and the scripts are there:
    if ( -e "$GWAVA_DIR/src/gwava_annotate.py") {

        # Running GWAVA annotation:
        `export GWAVA_DIR=${GWAVA_DIR}; python ${GWAVA_DIR}/src/gwava_annotate.py ${gene_name}_GRCh37.bed ${gene_name}_GRCh37.annot`;

        # Running GWAVA prediction:
        `export GWAVA_DIR=${GWAVA_DIR}; python ${GWAVA_DIR}/src/gwava.py tss ${gene_name}_GRCh37.annot ${gene_name}_GRCh37.gwava`;

        # Now let's test if the file is exists, if not, all weights are 1:
        if ( -e $gene_name."_GRCh37.gwava" ) {
            open(my $GWF, "<", $gene_name."_GRCh37.gwava" ) or warn "[Warning] ${gene_name}_GRCh37.gwava could not be opened!\n";
            while ( my $line = <$GWF>) {
                chomp $line;
                my ($chr, $start, $end, $ID, $gwava) = split("\t", $line);
                $gwava_scores{$ID} = $gwava;
            }
        } else { print "[Warning] GWAVA file could not be opened!\n" }
    }
    else {
        print "\n[Warning] GWAVA executables were not found in the specified path: $GWAVA_DIR/src/gwava_annotate.py\n";
        print "[Warning] Please update config file!\n[Info] Adding NA-s as scores. "
    }
    # Looping through all variants in the hash:
    foreach my $var (keys %hash){
        $hash{$var}{"score"} = exists $gwava_scores{$var} ? $gwava_scores{$var} : "NA"; # Initialize GWAVA score.
    }

    # Report progress:
    print "Done.\n";
    return \%hash;
}


##
## Lifting over the coordinates of the variants to the older build:
##
sub liftover {
    my $gene_name = $_[1];
    my %hash = %{$_[0]};

    my $tempFileName = sprintf($temporaryBedFile, $gene_name);
    open( my $tempbed, ">", $tempFileName) or die "[Error] Temporary bedfile could not be opened for writing: $tempFileName\n";
    foreach my $variant (keys %hash){
        printf $tempbed "%s\t%s\t%s\t%s\n", $hash{$variant}{GRCh38}[0], $hash{$variant}{GRCh38}[1] - 1, $hash{$variant}{GRCh38}[1], $variant;
    }

    # Liftover query:
    my $liftover_query = sprintf("%s %s %s/hg38ToHg19.over.chain.gz %s_GRCh37.bed %s_unmapped.bed  2> /dev/null", $liftoverPath, $tempFileName, $scriptdir, $gene_name, $gene_name);

    # Calling liftover:
    `bash -O extglob -c \'$liftover_query\'`;

    # Reading mapped file:
    my $lifted_file = sprintf("%s_GRCh37.bed", $gene_name);
    open(my $lifted, "<", $lifted_file) or die "[Error] After liftover run, the mapped file could not be opened.\n";
    my $liftedVarNo  = 0;
    while (my $line = <$lifted>) {
        chomp $line;
        my ($chr, $start, $end, $SNPID) = split("\t", $line);
        $hash{$SNPID}{"GRCh37"} = [$chr, $start, $end];
        $liftedVarNo ++;
    }

    print  "[Info] Number of variants successfully lifted over: $liftedVarNo\n\n" if $verbose;
    my $flag = $liftedVarNo == 0 ? 0 : 1;
    return (\%hash, $flag);
}


##
## Filtering variation list based on MAF and MAC
##
sub processVar {
    my $variants = $_[0]; # List of all overlapping variants
    my %genotypeContainer = %{$_[1]}; # hash to contain genotype information
    my $total = 0;
    my %hash = ();

    print  "[Info] Filtering variants:\n" if $verbose;

    my @total_vars = split("\n", $variants);
    printf  "[Info] Total number of overlapping variants: %s\n", scalar(@total_vars) if $verbose;

    foreach my $variant (@total_vars){
        $total ++;

        #line: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EGAN00001033155
        my ($chr, $pos, $id, $a1, $a2, $qual, $filter, $info, $format, @genotypes) = split(/\t/, $variant);

        # Generating variant name:
        my $SNPID = sprintf("%s_%s_%s_%s", $chr, $pos, $a1, $a2);

        # Parsing info field for relevant information:
        $info =~ /AC=(.+?);.+;AN=(.+?);.+consequence=(\S+?);/;
        my ($ac, $an) = $2 ? ($1, $2) : ("NA", "NA");

        my $consequence = $3 ? $3 : "NA";

        # We don't consider multialleleic sites this time.
        if ( $a2 =~ /,/){
            print  "[Warning] $SNPID will be omitted because of multiallelic! ($a2).\n";
            next;
        }

        print "\n\nWARNING: $variant\n\n" if $ac eq "NA";

        # Calculating values for filtering:
        my $missingness = (scalar(@genotypes)*2 - $an)/(scalar(@genotypes)*2);
        my $MAF = $ac/$an;

        # This flag shows if the non-ref allele is the major:
        my $genotypeFlip = 0; #
        if ( $MAF > 0.5 ){
            print "[Info] $MAF is greater then 0.5, genotype is flipped.\n";
            $MAF = 1 - $MAF;
            $genotypeFlip = 1;
        }

        my $MAC = $ac;
        $MAC = $an - $ac if $MAC > $an / 2;

        # We don't consider indels this time.
        if (( length($a2) > 1 or length($a1) > 1 )  && $score ne "NA"){
            print  "[Warning] $SNPID will be omitted because indel! ($a1/$a2).\n";
            next;
        }
        # Filter out variant because of high missingness:
        if ( $missingness > $missingthreshold && $score ne "NA"){
            print  "[Warning] $SNPID will be omitted because of high missingness ($missingness).\n";
            next;
        }
        # Filter out variant because of high MAF (regardless of the applied weight):
        if ( $MAF > $Variant{'MAF'} ){
            printf "[Warning] $SNPID will be omitted because of high MAF (%.3f, cutoff: %s).\n", $MAF, $Variant{'MAF'};
            next;
        }
        # Filter out variant because of high MAF:
        if ( $ac < $Variant{'MAC'} ){
            printf "[Warning] $SNPID will be omitted because of low minor allele count ($ac, cutoff: %s).\n", $Variant{'MAC'};
            next;
        }
        # If loss of function variants are required, we skipp althose variants that are not LOf:
        if ($lof && not exists $lof_cons{$consequence}) {
            printf "[Warning] $SNPID will be omitted because of consequence is not lof (%s).\n", $consequence;
            next;
        }
        # If loftee is enabled, only high-confidence loss-of-function variants will be selected:
        if ($loftee && $info =~ /LoF_conf\=-/ ) {
            printf "[Warning] $SNPID will be omitted because it is not a high-confidence loss of function variant.\n";
            next;
        }

        # If everything went fine, initializing variant by adding missingness to the hash:
        $hash{$SNPID}{"missingness"} = $missingness;


        #Storing variant data for all variant:
        $hash{$SNPID}{"alleles"} = [$a1, $a2];
        $hash{$SNPID}{"GRCh38"} = [$chr, $pos];
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


##
## bcftools query:
##
sub GetVariants {
    my $merged = $_[0];
    my $distance = 0;

    # Finding out which chromosome are we on:
    my ($chr) = $merged =~ /(chr.+?)\t/;

    # Print info:
    print  "\n[Info] Extracting variants from vcf files:\n" if $verbose;
    my $chrSpecVcfFile = sprintf("$vcfFile", $chr);
    my $bcftoos_query = sprintf("tabix %s ", $chrSpecVcfFile);

    # looping through all lines:
    foreach my $line (split("\n", $merged)){
        my ($chr, $start, $end) = split("\t", $line);
        $distance += $end - $start;
        $bcftoos_query .= sprintf(" %s:%s-%s", $chr, $start, $end);
    }

    # Extract overlapping variants:
    #foreach my $line (split("\n", $merged)){
    #    my ($chr, $start, $end) = split("\t", $line);
    #    $distance += $end - $start;
    #    my $bcftoos_query = sprintf("tabix %s %s:%s-%s", $vcfFile, $chr, $start, $end);
    print  "$bcftoos_query\n" if $verbose;
    my $variations = `bash -O extglob -c \'$bcftoos_query\'`;

    print sprintf("[Info] Total covered genomic regions: %s bp\n", $distance) if $verbose;
    print  "\n" if $verbose;

    return $variations;

}

##
## Filtering lines based on submitted criteria
##
sub FilterLines {
    my $lines = $_[0];
    my $stable_ID = $_[1];
    my @output_lines = ();
    my %hash = ();

    # The following hashes have to be checked:
    # %GENCODE
    # %GTEx
    # %overlap

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

            push (@output_lines, formatLines(\%annot_hash, $GENCODE{'extend'}));
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
            print  "\n[Info] From the GENCODE dataset the following features were extracted: $report\n" if $verbose;
        }
        if (exists $hash{GTEx}) {
            my $report = '';
            foreach my $feature (keys %{$hash{GTEx}}){
                $report .= "$feature: $hash{GTEx}{$feature}, ";
            }
            print  "[Info] The following regulatory features were linked to the gene: $report\n" if $verbose;
        }
        if (exists $hash{overlap}) {
            my $report = '';
            foreach my $feature (keys %{$hash{overlap}}){
                $report .= "$feature: $hash{overlap}{$feature}, ";
            }
            print  "[Info] The following overlapping regulatory features were extracted: $report\n"  if $verbose;
        }
    }

    my $filteredLines = join("\n", @output_lines);
    print "[Info] Selected lines:\n", $filteredLines,"\n" if $verbose;

    # Collapsing overlapping features:
    my $queryString = sprintf("mergeBed -i <(echo -e \"%s\" | sort -k1,1 -k2,2n )", $filteredLines);
    my $merged = `bash -O extglob -c \'$queryString\'`;
    return $merged;
}

##
## This function prepares the second bedfile that will be used for extract variants from
## the vcf files.
##
sub formatLines {
    my %hash = %{$_[0]};
    my $ext = $_[1] // 0;
    return sprintf("%s\t%s\t%s\t%s", $hash{"chr"}, $hash{"start"} - $ext, $hash{"end"} + $ext, encode_json(\%hash))
}


##
## This function returns all lines corresponding to the given gene using bedtools.
##
sub BedToolsQuery {
    my ($chr, $start, $end, $stable_ID) = @_;
    my $queryString = sprintf("intersectBed -wb -a <(echo -e \"chr%s\\t%s\\t%s\\t%s\") -b %s -sorted | cut -f9-",
                                $chr, $start, $end, $stable_ID, $geneBedFile, $stable_ID);
    print "[Info] IntersectBed query string: $queryString\n" if $verbose;
    my $query = `bash -O extglob -c \'$queryString\'`;
    return $query;

}

##
## This function returns genomic coordinates of a gene given its stable ID or external name.
##
sub GetCoordinates {
    my $ID = $_[0];

    # If Stable ID is given:
    my ($chr, $start, $end, $stable_ID, $name) = ('NA' x 5);
    if ( exists $GENCODE_ID->{$ID} ) {
        $chr        = $GENCODE_ID->{$ID}->{chr};
        $start      = $GENCODE_ID->{$ID}->{start};
        $end        = $GENCODE_ID->{$ID}->{end};
        $stable_ID  = $GENCODE_ID->{$ID}->{ID};
        $name       = $GENCODE_ID->{$ID}->{name};
    }
    elsif (exists $GENCODE_name->{$ID}) { # Assuming gene name is given:
        $chr        = $GENCODE_name->{$ID}->{chr};
        $start      = $GENCODE_name->{$ID}->{start};
        $end        = $GENCODE_name->{$ID}->{end};
        $stable_ID  = $GENCODE_name->{$ID}->{ID};
        $name       = $GENCODE_name->{$ID}->{name};
    }
    else {
        print  "[Warning] $ID was not found in the GENCODE database. Only GENCODE gene names and stable Ensembl IDs are accepted. Gene will be skipped!\n";
    }


    return ($chr, $start, $end, $stable_ID, $name);
}

##
## Function for reading config file.
##
sub read_param {
    my $filename = $_[0]; # This value is given from the command line.

    # Create absolute path:
    my $absConfigFile = $filename ? `readlink -f $filename` : "NA";

    chomp $absConfigFile;

    # if any problem occur, we use the default value:
    unless ( -e $absConfigFile ) {
        $absConfigFile = $scriptdir."/config.txt";
        $absConfigFile = `readlink -f $absConfigFile`;
        chomp $absConfigFile;
    }

    print  "[Info] Config file: $absConfigFile\n";

    # Reading file:
    open(my $CONF, "cat $absConfigFile | grep -v \"#\" | ") or die "[Error] Config file could not be oppended. Exiting.\n";

    # Reading file:
    my %hash;
    while ( my $line = <$CONF>) {
        chomp $line;
        next unless $line;

        my @pair = split(/=/, $line);
        $hash{$pair[0]} = $pair[1];
    }

    return \%hash;
}

##
## Adding consequences:
##
sub get_Consequences {
    print "[Info] Returning variant consequences." if $verbose;
    my %source_hash = %{$_[0]};
    # hash with all
    my @array = ();
    foreach my $var (keys(%source_hash)){
        my ($chr, $pos, $ref, $alt) = split("_", $var);
        $chr =~ s/^chr//i;

        # Building hash:
        my $string = sprintf("%s %s %s %s %s . . .", $chr, $pos, $var, $ref, $alt);
        push(@array, $string);

        # once we reach the limit we have to process:
        if (scalar(@array) == 200) {
            print "." if $verbose;
            %source_hash = %{&submitEnsembl(\%source_hash, \@array)};
            @array = ();
        }

    }

    # Submitting the remaining variants:
    %source_hash = %{&submitEnsembl(\%source_hash, \@array)};

    print " Done.\n" if $verbose;
    return \%source_hash;
}

sub submitEnsembl {

    my %source_hash = %{$_[0]};

    # Generate json formatted query string:
    my $jsonstring = encode_json({"variants" => $_[1]});

    # Initialize http request:
    my $http = HTTP::Tiny->new();
    my $server = 'http://rest.ensembl.org';
    my $ext = '/vep/homo_sapiens/region';
    my $response = $http->request('POST', $server.$ext, {
      headers => {
        'Content-type' => 'application/json',
        'Accept' => 'application/json'
      },
      content => $jsonstring
    });

    # Did everything went smoothly?
    warn sprintf("[Error] Returning consequences from Ensembl has failed!\nResponse status: %s\nReason: %s\nContent: %s",
                $response->{status}, $response->{reason}, $response->{status})unless $response->{success};

    # Decoding response content:
    my @VEP = @{decode_json($response->{content})};

    # Looping through the returned array:
    foreach my $var (@VEP){
        my ($input) = $var->{"input"} =~ / (chr\S+) /;
        my $most_severe_consequece = $var->{"most_severe_consequence"};

        # Adding consequence to hash, if exsists:
        $source_hash{$input}{"Consequence"} = $most_severe_consequece ? $most_severe_consequece : "NA";

        # Adding colocated known variant if exists:
        $source_hash{$input}{"rsID"} = exists $var->{colocated_variants}->[0]->{id} ? $var->{colocated_variants}->[0]->{id} : "NA";
    }

    return \%source_hash;
}

##
## Reading GENCODE gene coordinates:
##
sub read_GENCODE {
    my $gencodeFile = $_[0];

    # Reading the file and create and return two hashes. Each contain the chromosome and
    # coordinate of the genes. Keys:
    my %gene_names = ();
    my %gene_IDs = ();
    open(my $FILE, "zcat $gencodeFile | ") or die "[Error] GENCODE file with gene coordinates could not be opened. Exiting.\n";
    while (my $line = <$FILE>) {
        next if $line =~ /^#/; # skipping header
        chomp $line;
        my ($chr, $start, $end, $name, $ID) = split("\t", $line);
        $chr =~ s/^chr//i;

        # generate hash ref:
        my $ref = {"chr" => $chr,
                   "start" => $start,
                   "end" => $end,
                   "name" => $name,
                   "ID" => $ID};

        # adding to gene names:
        $gene_names{$name} = $ref;
        $gene_IDs{$ID} = $ref;
    }

    return \%gene_IDs, \%gene_names;
}


##
## Below is the documentation for the script:
##
=pod

=head1 DESCRIPTION

This script was written to select variants for collapsing tests (MONSTER). The script allows to specify
a custom set of criteria to define genomic regions of interest, scoring and weighting methods.

In more details: the selection is perfomed in a gene basis: the user can specify which functional
part of the gene is of interest (eg. exon, CDS stc). Besides the GENCODE elements, Ensembl regulatory
features are also available to be selected: if a feature overlaps with the queried gene, or
overlaps with with an eQTL signal which is linked to the gene.

Overlapping variants from the 15X helic data are then selected using intersectbed. After applying
a various set of filters, scores are assigned to the variants to use them as weights by the
collapsing test. Optionally weigh can be adjusted for allele frequencies to add
more weight to rare variants.

=head1 VERSION

v.2.2. Last modified: 2016.08.11 by Daniel Suveges, mail: ds26@sanger.ac.uk

=head1 SYNOPSIS

perl $0
    --input|i <gene name>
    --output|o <output prefix>
    --GENCODE|G <gencode features>
    --GTEx|E <regulatory features>
    --overlap|L <regulatory features>
    --extend <region extension in bp>
    --SkipMinor
    --tissues|T <list of tissues>
    --MAF <upper MAF threshold>
    --MAC <lower minor allele count>
    --verbose|v
    --conf|c <config file name>
    --score|s <scoring method>
    --MAFweight|w <weighting scores based on MAF>
    --ExpConst|k <weighting const, float>
    --lof
    --loftee
    --missingness|m <upper missingness threshold, float>
    --shift <by which the scores are shifted, float>
    --cutoff <score cutoff, float>
    --floor <lower threshold of scores>
    --help|h


B<Input gene name:>

Single Ensembl gene name or stable gene IDs. Required parameter!

B<Output files:>

Output prefix ${output}_genotype and ${output}_variants files are saved. Genotype
and SNP files are suitable input for MONSTER (more information of the format:
I<www.stat.uchicago.edu/~mcpeek/software/MONSTER/MONSTER_v1.2_doc.pdf>)
Output prefix is a required parameter!

B<Feature selection:>

=over 4

B<GENCODE>: accepted gencode parameters: comma separated list of GENCODE feature
names: gene, exon, transcript, CDS, UTR (without spaces!)

B<GTEx>: specifying the types of Ensembl regulatory features that were linked to
the gene by overlapping GTEx eQTL signal.

B<Overlap>: specifying the types of Ensembl regulatory features that were
overlapping the queried gene.

For GTEx, Overlap the following regulatory feature names are accepted: promoter,
enhancer, promoterFlank, openChrom, TF_bind, allreg (where allreg equal to selecting
all regulatory features.) (without spaces!)

B<tissues>: expect a comma separated list of tissues. Will use only those features
that were found to be active in the listed features.

B<extend>: the selected GENCODE regions will be extended with the submitted
value (default: 0bp)

=back

B<Varinat filters:>

=over 4

B<MAF>/B<MAC> for setting the minimal allele count and maximal allele frequencies.
Both MAF/MAC values are calculated from the corresponding vcf line. (default MAF: 0.05, MAC: 0)

B<SkipMinor>: any GENCODE feature will be excluded that belongs to a transcript
annotated as minor in the APPRIS database.

B<MAF_weight>: MAF weighting will be applied. Exponential scoring the lower MAF
regions. C<s_w = s * e^(-k * ($AC - 1)/$AN)>
    where s_w: weighted score, s: score, AC: allele count, AN: allele number,
    k: exponential constant

=back

B<Scoring options:>

=over 4

B<score>: available scoring meghods are: CADD, Eigen, GWAVA are supported currently
(default: no weight)

B<ExpConst>: describes the steepness of the weighting function. (default: 50)

B<conf>: the specified that will used to read parameters instead of the default config.txt

B<lof>: only variants with high predicted consequences were included in the variant set.
Consequences are expected to be stored in the INFO field of the vcf file after the "consequence" tag.

B<loftee>: Only high and low confidence loss of function variants are included based
on loftee prediction. This information is expected to be sotred in the vcf file after
the loftee tag.

B<missingness>: upper threshold of missingness. Variants with missingness higher
than this threshold will be excluded. (calculated based on the vcf file AN and AF field)

B<shift>: a float describing how much the scores are shifted. By default this value is 0, so
the scores are not sifted.

B<cutoff>: variants with lower score than the cutoff will be excluded. by default it is 0, as
MONSTER does not accept negaive weights.

B<floor>: this threshold score will be added to every variants below this threshold. Default
value is 0.

=back

B<Other options:>

=over 4

B<Verbose>: prints out more reports during the progression of the script.

B<help>: prints out help message.

=back


=head1 FILES:

To properly run the script all of these files have to be found. Otherwise the
script exits. The path of these files have to be specified in the configuration
file (config.txt) in the script directory, or has to be specified using --config.


=over 4

B<geneBedFile>: GENCODE based genome-wide annotation file. Mapped to GRCh38 coordinates.

B<vcfFile>: input vcf files

B<liftoverPath>: As scores are not available to the newer b38 build, the selected variants
have to be mapped to the b37 build.

B<chainFile>: chain file to map GRCh38 coordinates to GRCh37 build.

B<EigenPathNonCoding/EigenPathCoding>: the genome wide set of eigen scores. The newer
release (v1.1) of the Eigen scores are divided into coding and non-coding files.

B<caddPath>: The genome wide set of CADD scores.

B<GENCODEFile>: A list of gencode gene names/gene IDs and coordinates to quickly gent the
coordinates.

=back

=cut
