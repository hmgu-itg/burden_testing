#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
use lib dirname(__FILE__);
use Scoring;

my $gene;
my $variant;

my %C = (
    "transcript_ablation"      => 36,
    "splice_acceptor_variant"  => 35,
    "splice_donor_variant"     => 34,
    "stop_gained"              => 33,
    "frameshift_variant"       => 32,
    "stop_lost"                => 31,
    "start_lost"               => 30,
    "transcript_amplification" => 29,
    "inframe_insertion"        => 28,
    "inframe_deletion"         => 27,
    "missense_variant" => 26,
    "protein_altering_variant" => 25,
    "splice_region_variant" => 24,
    "incomplete_terminal_codon_variant" => 23,
    "start_retained_variant" => 22,
    "stop_retained_variant" => 21,
    "synonymous_variant" => 20,
    "coding_sequence_variant" => 19,
    "mature_miRNA_variant" => 18,
    "5_prime_UTR_variant" => 17,
    "3_prime_UTR_variant" => 16,
    "non_coding_transcript_exon_variant" => 15,
    "intron_variant" => 14,
    "NMD_transcript_variant" => 13,
    "non_coding_transcript_variant" => 12,
    "upstream_gene_variant" => 11,
    "downstream_gene_variant" => 10,
    "TFBS_ablation" => 9,
    "TFBS_amplification" => 8,
    "TF_binding_site_variant" => 7,
    "regulatory_region_ablation" => 6,
    "regulatory_region_amplification" => 5,
    "feature_elongation" => 4,
    "regulatory_region_variant" => 3,
    "feature_truncation" => 2,
    "intergenic_variant" => 1
);

    
sub getVariantType{
    my ($ref,$alt)=@_;

    return "SNP" if length($ref)==1 && length($alt)==1;
    return "DEL" if length($ref)>1 && length($alt)==1;
    return "INS" if length($ref)==1 && length($alt)>1;

    return "NA";    
}

sub getConsequences{
    my $variant=$_[0];
    my $lof=$_[1];
    my $stable_ID=$_[2];

    my $vepin;
    my $vepout;
    
    my $cons;

    local $\="\n";
    local $,="\t";

    my $fname1="temp_vep_input.txt";

    open ($vepin, ">", $fname1) or die "[Error] Input file for VEP could not be opened.";
    
    my ($chr, $pos, $ref, $alt) = split(/:/, $variant);
    (my $c = $chr ) =~ s/chr//i;
    my $varID=$c."_".$pos."_".$ref."_".$alt;

    # skip multiallelics
    if ($alt=~/,/){
	print "Warning: $variant is multiallelic";
	return "NA";
    }

    my $vtype=getVariantType($ref,$alt);
    if ($vtype eq "SNP"){
	print $vepin $c,$pos,$pos,$ref."/".$alt,"+",$varID;
    }
    elsif($vtype eq "DEL"){
	my $r=substr($ref,1,length($ref)-1);
	print $vepin $c,$pos+1,$pos+length($ref)-1,$r."/-","+",$varID;
    }
    elsif($vtype eq "INS"){
	my $a=substr($alt,1,length($alt)-1);
	print $vepin $c,$pos+1,$pos,"-/".$a,"+",$varID;
    }
    else{
	print "[Error] could not determine variant type of $variant";
	return "NA";
    }
    
    close($vepin);

    my $queryString="vep -i ".$fname1." --dir /usr/local/bin/.vep --dir_cache /usr/local/bin/.vep -o STDOUT --offline --no_stats | grep -v \"^#\" | awk -v g=".$stable_ID." 'BEGIN{FS=\"\\t\";}\$4==g{print \$0;}' | cut -f 1,7";
    my $query =Scoring::backticks_bash($queryString);
    my $max_severity=0;
    foreach my $line (split (/\n/, $query)){
	my @a=split(/\t/,$line);
	my $ID=$a[0];
	my $effs=$a[1];

	foreach my $e (split(/,/,$effs)){
	    if (! exists($lof->{$e})){ # low severity
		    $cons=$e if ($max_severity==0);
	    }
	    else{
		if ($lof->{$e} > $max_severity){
		    $max_severity=$lof->{$e};
		    $cons=$e;
		}
	    }
	}

    }

    return $cons;
}

GetOptions('gene=s' => \$gene,'variant=s' => \$variant);

my $c=getConsequences($variant,\%C,$gene);
$,="\t";
$\="\n";

print $variant,$gene,$c;



    

    
