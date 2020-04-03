#!/usr/bin/perl

use  strict;
use Getopt::Long;

my @input;

GetOptions('input=s' => \@input,'output=s' => \$outputFile);

