#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my ($inputFile,$outputFile,$excludeFile);

GetOptions('input=s' => \$inputFile,'output=s' => \$outputFile,'exclude=s' => \$excludeFile);

my %vars_to_exclude;
my $fh;

$\="\n";

open (my $fh, "<", $excludeFile) or die "[Error] File with variants to exclude ($excludeFile) could not be opened. Exiting.";
while(my $v=<$fh>){
    chomp $v;
    $vars_to_exclude{$v}=1;
}
close($fh);

open (my $fh, "<", $inputFile) or die "[Error] Input file ($inputFile) could not be opened. Exiting.";
my %H;
while(my $line=<$fh>){
    chomp $line;
    my @a=split(/\t/,$line);

    my @ind=();
    my $flag=0;
    
    if ($a[1]==1){
	# this line contains variant IDs
	# next line contains weights
	$flag=1;
	$H{$a[0]}->{"vars"}=join("\t",@a[2..$#a]);
    }
    elsif($a[1]==0){
	if ($flag==1){
	    $H{$a[0]}->{"wights"}=join("\t",@a[2..$#a]);
	}
	else{
	    $H{$a[0]}->{"vars"}=join("\t",@a[2..$#a]);
	    $H{$a[0]}->{"weights"}="NA";
	}

	$flag=0;
    }
    else{
	die "[Error] The second field in $line is neither 0 nor 1. Exiting.";
    }

}
close($fh);

my %mono;
foreach my $gene (keys(%H)){
    my $vars=$H{$gene}->{"vars"};
    my @a=split(/\t/,$vars);
    my $c=0;
    for (my $i=0;$i<scalar(@a);$i++){
	if exists($vars_to_exclude{$a[$i]}){
	    $a[$i]="NA";
	    $c++;
	}
    }

    my $rem=scalar(@a)-$c;
    $mono{$gene}=1 if($rem<2);

    $H{$gene}->{"vars"}=join("\t",@a);
}

$,="\t";
open (my $fh, ">", $outputFile) or die "[Error] Output ($outputFile) could not be opened. Exiting.";
foreach my $gene (keys(%H)){
    if (exists($mono{$gene})){
	print STDERR $gene;
    }
    else{
	my $vars=$H{$gene}->{"vars"};
	my $weights=$H{$gene}->{"weights"};
	my @a=split(/\t/,$vars);
	my @c=();
	my @d=();

	for (my $i=0;$i<scalar(@a);$i++){
	    push(@c,$a[$i]) unless($a[$i] eq "NA");
	}
	
	if ($weights eq "NA"){
	    print $fh $gene,"0",join("\t",@c);
	}
	else{
	    my @b=split(/\t/,$weights);
	    die "[Error] Numbers of variants and weights for $gene are different. Exiting." if scalar(@b)!=scalar(@b);
	    for (my $i=0;$i<scalar(@b);$i++){
		push(@d,$b[$i]) unless($a[$i] eq "NA");
	    }
	    print $fh $gene,"1",join("\t",@c);
	    print $fh $gene,"0",join("\t",@d);
	}
    }
}
close($fh);

