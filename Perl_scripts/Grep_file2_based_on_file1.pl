#!/usr/bin/perl -w

# Script used to collect all lines in file2 which match one query in file1

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use List::Util 'first';

# Define the different input/output files
my $file1;    # Input file listing all queries
my $file2;  # Input file containing the information to retrieve
my $output; # Output file containing for each query (of file1) its associated information (of file2)

# Define the parameter in order to submit input files to this script
&GetOptions (
    'file1=s' => \$file1,
    'file2=s' => \$file2,
    'output=s' => \$output,
);

# Open the input file1
unless ($file1) {
    die "Please specify the file containing the list of query via -file1 parameter!\n";
}
open (FILE1, "<$file1") || die "Cannot open $file1: $!\n"; $_="1";
my @file1 = <FILE1>;
close (FILE1);

# Open the input file2
unless ($file2) {
    die "Please specify the file containing the information to retrieve via -file2 parameter!\n";
}
open (FILE2, "<$file2") || die "Cannot open $file2: $!\n"; $_="1";

# Determine the format of input file2
my $format;
if ($file2 =~ /\.gtf/) {
    $format = "gtf";
}
elsif ($file2 =~ /\.fa/) {
    $format = "fasta";
}
else {
    die "The input file2 is not in one of the accepted format, please seek advice with this script!\n";
}

# Open the output file
unless ($output) {
    die "Please specify the output file via -output parameter!\n";
}
if (-e $output) {
    die "This file: $output already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Define variables required for reading the file1
my $query;

# Define variables required for reading the input file
my $line;        # Scalar containing line from the input file2
my $line_count = 0; # Number of line already read in file2
my $match_line = 0;  # Number of line in file2 matching a query in file1
my $matching = "no";    # Scalar which contains yes or no depending if file2 current value matches in file1

# Process input files
while(1) {
    chomp($line = <FILE2>);    # Read lines one by one from the information file2
    $line_count ++;
    if ($format eq "gtf") {
        if ($line =~ /transcript_id \"(.*?)\";/){   # Identify the transcript id from line of file2
            my $transcript = $1;
            if (first { /$transcript/ } @file1) {
                print OUTPUT "$line\n";
                $match_line ++;
            }
        }
        else {
            die "The file $file2 does not contain a transcript id in gtf file at line $line_count, please check your file!\n";
        }
    }
    elsif ($format eq "fasta") {
        if ($line =~ /^>(ENSP\d{11})/){   # Identify the protein id from line of file2
            my $prot = $1;
            if (first { /$prot/ } @file1) {
                print OUTPUT "$line\n";
                $match_line ++;
                $matching = "yes";
            }
            else {
                $matching = "no";
            }
        }
        elsif ($matching eq "yes") {
            print OUTPUT "$line\n";
        }
        elsif ($line =~ /^>/) {
            die "The file $file2 does not contain a protein id in fasta file at line $line_count, please check your file!\n";
        }
    }
    last if eof (FILE2);  # If the file2 was fully read, then exit reading the file
}
close (FILE2);    # Close second file
close (OUTPUT); # Close the output file

print STDERR "There were $match_line line from file $file2 which matched a query of file $file1!\n";

__END__