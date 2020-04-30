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

# Simply reading the data stored in the gzipped gencode file:
sub _initialize {
    my $self = shift;
    my $gencodeFile = shift;
    my $FILE;
    
    die "[Error] GENCODE file with gene coordinates could not be opened. Exiting.\n" unless -e $gencodeFile;
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

        # adding to hash, unless already exists:
	if (exists($self->{"gene_names"}->{$name})){
	    print "[Warning] GENCODE::_initialize : $name is already in the hash; skipping $name";
	}
	else{
	    $self->{"gene_names"}->{$name} = $ref;
	}

        # adding to hash, unless already exists:
	if (exists($self->{"gene_names"}->{$ID})){
	    print "[Warning] GENCODE::_initialize : $ID is already in the hash; skipping $ID";
	}
	else{
	    $self->{"gene_names"}->{$ID} = $ref;
	}

    }
    close($FILE);
}

# A method to extract gene coordinates based on the gene name or stable ID.
sub GetCoordinates {
    my $self = shift;
    my $ID = shift; # stable ID or gene name
    my ($chr, $start, $end, $stable_ID, $name) = ('NA') x 5;

    if ( exists $self->{"gene_names"}->{$ID} ) {
        $chr        = $self->{"gene_names"}->{$ID}->{chr};
        $start      = $self->{"gene_names"}->{$ID}->{start};
        $end        = $self->{"gene_names"}->{$ID}->{end};
        $stable_ID  = $self->{"gene_names"}->{$ID}->{ID};
        $name       = $self->{"gene_names"}->{$ID}->{name};
    }
    else {
        print  "[Warning] GENCODE::GetCoordinates: $ID was not found\n";
    }

    return ($chr, $start, $end, $stable_ID, $name);
}
1;
