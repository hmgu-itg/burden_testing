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

# reading the data stored in the gzipped gencode file:
sub _initialize {
    my $self = shift;
    my $gencodeFile = shift;
    my $FILE;
    
    die "[Error] GENCODE file with gene coordinates could not be opened\n" unless -e $gencodeFile;
    $self->{"failed"}=0;

    open($FILE, "zcat $gencodeFile | ");
    while (my $line = <$FILE>) {
        next if $line =~ /^#/;
	
        chomp $line;
	
	#1-based coordinates
        my ($chr, $start, $end, $name, $ID) = split("\t", $line);
        $chr =~ s/^chr//i;

        my $ref = {"chr" => $chr,
                   "start" => $start,
                   "end" => $end,
                   "name" => $name,
                   "ID" => $ID};

	# by name
	push @{$self->{"gene_names"}->{$name}}, $ref;

	# by ID
        # shouldn't happen as IDs are supposed to be unique
	if (exists($self->{"gene_names"}->{$ID})){
	    print "[Error] GENCODE::_initialize : $ID is already in the hash";
	    $self->{"failed"}=1;
	    return;
	}
	else{
	    $self->{"gene_names"}->{$ID} = $ref;
	}

	# also use ID prefix as key
	if ($ID=~/(ENSG\d+)\./){
	    my $x=$1;
	    push @{$self->{"gene_names"}->{$x}}, $ref;
	}
    }
    
    close($FILE);
}

# get gene coordinates based on the gene name, stable ID or stable ID prefix.
sub GetCoordinates {
    my $self = shift;
    my $ID = shift;
    my $ret=undef;

    if ( exists $self->{"gene_names"}->{$ID} ) {
	my $stID=$self->{"gene_names"}->{$ID}->{ID};
        $ret->{$stID}->{chr}        = $self->{"gene_names"}->{$ID}->{chr};
        $ret->{$stID}->{start}      = $self->{"gene_names"}->{$ID}->{start};
        $ret->{$stID}->{end}        = $self->{"gene_names"}->{$ID}->{end};
        $ret->{$stID}->{name}       = $self->{"gene_names"}->{$ID}->{name};
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
