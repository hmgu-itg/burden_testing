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
    
    
    my $GENCODE_filename = sprintf("%s/gencode_genes_V25_GRCh%s.tsv.gz", $parameters->{"scriptDir"}, $parameters->{"build"});
    $self->_initialize($GENCODE_filename);
    return $self;
}

# Simply reading the data stored in the gzipped gencode file:
sub _initialize {
    my $self = shift;
    my $gencodeFile = shift;

    # Reading the file and create and return two hashes. Each contain the chromosome and
    # coordinate of the genes. Keys:
    
    die "[Error] GENCODE file with gene coordinates could not be opened. Exiting.\n" unless -e $gencodeFile;
    open(my $FILE, "zcat $gencodeFile | ");
    while (my $line = <$FILE>) {
        next if $line =~ /^#/;
        chomp $line;
        my ($chr, $start, $end, $name, $ID) = split("\t", $line);
        $chr =~ s/^chr//i;

        # generate hash ref for each gene:
        my $ref = {"chr" => $chr,
                   "start" => $start,
                   "end" => $end,
                   "name" => $name,
                   "ID" => $ID};
        # adding to hash:
        $self->{"gene_names"}->{$name} = $ref;
        $self->{"gene_names"}->{$ID} = $ref;
    }
}

# A method to extract the coordinates of any gene based on the gene name or stable ID.
sub GetCoordinates {
    my $self = shift;
    my $ID = shift;
    my ($chr, $start, $end, $stable_ID, $name) = ('NA') x 5;

    # If Stable ID is given:
    if ( exists $self->{"gene_names"}->{$ID} ) {
        $chr        = $self->{"gene_names"}->{$ID}->{chr};
        $start      = $self->{"gene_names"}->{$ID}->{start};
        $end        = $self->{"gene_names"}->{$ID}->{end};
        $stable_ID  = $self->{"gene_names"}->{$ID}->{ID};
        $name       = $self->{"gene_names"}->{$ID}->{name};
    }
    else {
        print  "[Warning] $ID was not found in the GENCODE database. Only HGNC names and stable Ensembl IDs are accepted. 'NA'-s will be returned!\n";
    }


    return ($chr, $start, $end, $stable_ID, $name);
}
1;