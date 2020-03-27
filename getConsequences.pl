#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use lib dirname(__FILE__);
use Scoring;

my $gene;
my $variant;

%C = {
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
    my $fname2="temp_vep_output.txt";

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
	$count++;
    }
    elsif($vtype eq "DEL"){
	my $r=substr($ref,1,length($ref)-1);
	print $vepin $c,$pos+1,$pos+length($ref)-1,$r."/-","+",$varID;
	$count++;
    }
    elsif($vtype eq "INS"){
	my $a=substr($alt,1,length($alt)-1);
	print $vepin $c,$pos+1,$pos,"-/".$a,"+",$varID;
	$count++;
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



    

    
