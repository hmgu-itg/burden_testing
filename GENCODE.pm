package GENCODE;

## This class was created to return genomic coordinates for a gene
## given its name or stable ID on Ensembl.

## Usage:

# Initialize:
#
# use GENCODE;
# my $G = GENCODE->new($GENCODE_fileName);

# Extract coordinates:
# ($chr, $start, $end, $stable_ID, $name) = $G->GetCoordinates($ID);

use strict;
use warnings;

sub new {
    my ( $class, $parameters) = @_;
    my $self = {};
    
    bless( $self, $class );
    
    my $GENCODE_filename = $parameters->{"gencode_file"};
    $self->_initialize($GENCODE_filename);
    return $self;
}

# reading gzipped gencode file
sub _initialize {
    my $self = shift;
    my $gencodeFile = shift;
    my $FILE;

    $self->{"failed"}=0;

    if (! -e $gencodeFile){
        print "[Error] GENCODE file ($gencodeFile) could not be opened\n";
        $self->{"failed"}=1;
        return;
    }
    
    my %counts;
    open($FILE, "zcat $gencodeFile | ");
    while (my $line = <$FILE>) {
        next if $line =~ /^#/;
	
        chomp $line;
	
	#1-based coordinates
        my ($chr, $start, $end, $name, $ID) = split("\t", $line);
        $chr =~ s/^chr//i;
        if ( $chr eq "Y" && (($start > 10001 && $end < 2781479) || ($start > 56887903 && $end < 57217415))){
		# we are in PAR Y. The genes are also on PAR X so we skip
		next;
        } 
        my $ref = {"chr" => $chr,
                   "start" => $start,
                   "end" => $end,
                   "name" => $name,
                   "ID" => $ID};

	# to detect (name) duplicates
	$counts{$name}++;
	
	# by name
	# array can contain more than one ref (in case o duplicate names)
	push @{$self->{"gene_names"}->{$name}}, $ref;

	# by ID
        # shouldn't happen as IDs are supposed to be unique
	if (exists($self->{"gene_names"}->{$ID})){
	    print "[Error] GENCODE::_initialize : duplicate ID $ID";
	    $self->{"failed"}=1;
	    return;
	}
	else{
	    push @{$self->{"gene_names"}->{$ID}}, $ref;
	}

	# also use ID prefix as key
	if ($ID=~/(ENSG\d+)\./){ # should always be the case
	    my $x=$1;
	    push @{$self->{"gene_names"}->{$x}}, $ref;
	}
    }

    foreach my $n (keys %counts){
	if ($counts{$n}>1){
	    $self->{"duplicate_names"}->{$n}=1;
	}
    }

    if (exists($self->{"duplicate_names"})){
	print "[Warning]: duplicate gene names detected in $gencodeFile:";
	foreach my $n (keys %{$self->{"duplicate_names"}}){
	    print $n;
	}
    }
    
    close($FILE);
}

# checks if provided gene name occurs more than once in GENCODE data
sub isDuplicateName {
    my $self = shift;
    my $name = shift;

    if (exists($self->{"duplicate_names"})){
	if (exists($self->{"duplicate_names"}->{$name})){
	    return 1;
	}
	else{
	    return 0;
	}
    }
    else{
	return 0;
    }
}

# get gene coordinates based on the gene name, stable ID or stable ID prefix.
sub GetCoordinates {
    my $self = shift;
    my $ID = shift; # name , full ID or ID prefix
    my $ret=undef;

    if ( exists $self->{"gene_names"}->{$ID} ) {
	foreach my $rec (@{$self->{"gene_names"}->{$ID}}){
	    my $stID=$rec->{ID};
	    $ret->{$stID}->{chr}        = $rec->{chr};
	    $ret->{$stID}->{start}      = $rec->{start};
	    $ret->{$stID}->{end}        = $rec->{end};
	    $ret->{$stID}->{name}       = $rec->{name};
	}
    }
    elsif ($ID=~/(ENSG\d+)\./){
	my $x=$1;
	if ( exists $self->{"gene_names"}->{$x} ) {
	    foreach my $rec (@{$self->{"gene_names"}->{$x}}){
		my $stID=$rec->{ID};
		$ret->{$stID}->{chr}        = $rec->{chr};
		$ret->{$stID}->{start}      = $rec->{start};
		$ret->{$stID}->{end}        = $rec->{end};
		$ret->{$stID}->{name}       = $rec->{name};
	    }
	}
    }
    else {
        print  "[Warning] GENCODE::GetCoordinates: $ID was not found\n";
    }

    return $ret;
}
1;
