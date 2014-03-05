#!/usr/bin/perl -w

# Script used on sequence fasta file for parsing into smaller fasta files

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the different input/output files such as fastq sequence file, indices list file, each index output fastq file
my $fasta;    # Input file listing all indices used in the pool RNA-seq library
my $fh;         # Filehandle output fastq file for each index
my $parsing;    # Number of sequence per parsed file

# Define the parameter in order to submit input files to this script
&GetOptions (
    'fasta=s' => \$fasta,
    'parsing=i' => \$parsing,
);

# Open the input fasta file
unless ($fasta) {
    die "Please specify the fasta file via -fasta parameter!\n";
}
open (FASTA, "<$fasta") || die "Cannot open $fasta: $!\n"; $_="1";

# Check that parsing number has been provided
unless ($parsing) {
    die "Please specify the number of sequence per parsed file via -parsing parameter!\n";
}

# Define variables required for reading the input file
my $name = $fasta;
$name =~ s/.fa//;
my $index = 0;
my %indices;

# Define variables required for reading the input fastq sequence file
my $line;        # Scalar containing line from the fasta file
my $count = 0;  # Scalar containing number of transcript per output fasta file

# Read in the fastq sequence file
while(1) {
    chomp($line = <FASTA>);    # Read lines one by one from the fasta file
    if ($line =~ /^>/){   # Identify the id from line of fasta file
        $count ++;
        if ($count == 1) {
            $index ++;
            $indices{$index} = $count;  # Put all indices sequence into hash
            my $outfile = "${name}_${index}.fa";
            open_fh($index, $outfile);  # Create and open an outfile for each individual fasta
            print_to_fh ($index, "$line\n");
        }
        elsif ($count < $parsing) {
            $indices{$index} = $count;  # Put all indices sequence into hash
            print_to_fh ($index, "$line\n");
        }
        elsif ($count == $parsing) {
            $indices{$index} = $count;  # Put all indices sequence into hash
            print_to_fh ($index, "$line\n");
            $count = 0;
        }
    }
    else {
        print_to_fh ($index, "$line\n");
    }
    last if eof (FASTA);  # If the sequence fastq file was fully read, then exit reading the fastq file
}
close (FASTA);    # Close sequence fasta file

print STDERR "Fasta_number\tNumber_of_sequences\n";
foreach my $keys (keys %indices){
    close_fh ($keys);
    print STDERR "$keys\t$indices{$keys}\n";
}

################
# Sub-routines #
################

# Subroutine to open a file handle outfile for each sample
sub open_fh {
    local *FH = shift;
    my $file = shift;
    open (*FH, ">>$file") || die "opening $file: $!"; $_="1" ;
}

# Subroutine to print output into file handle outfile for each sample
sub print_to_fh {
    local *FH = shift;
    my $string = shift;
    print FH $string;
}

# Subroutine to close a file handle outfile for each sample
sub close_fh {
    local *FH = shift;
    my $file = shift;
    close (*FH);
}

__END__