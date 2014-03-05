#!/usr/bin/perl -w

# Script used to reformat (and compile) report results from STAR run

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $report;    # Input file containing the final report from STAR run

# Define the parameter in order to submit input files to this script
&GetOptions (
    'report=s' => \$report,
);

# Open the input fasta file
unless ($report) {
    die "Please specify the report file via -report parameter!\n";
}
open (REPORT, "<$report") || die "Cannot open $report: $!\n"; $_="1";

# Find the name of the sample on which the report is based
my $name = ${report};
$name =~ s/_?Log\.final\.out//;
$name =~ s/.*\///;

# Define variables required for reading the input file
my $header;         # Scalar containing first column of file
my $value;          # Scalar containing second column of file
my $input = "Number of input reads";        # Scalar containing header of interest from file
my $input_number;                           # Scalar containing value of interest from file
my $map = "Uniquely mapped reads number";   # Scalar containing header of interest from file
my $map_number;                             # Scalar containing value of interest from file
my $map_percent = "Uniquely mapped reads %";# Scalar containing header of interest from file
my $map_percent_number;                     # Scalar containing value of interest from file
my $average_length = "Average mapped length";# Scalar containing header of interest from file
my $map_length;                     # Scalar containing value of interest from file
my $multihit = "Number of reads mapped to multiple loci";# Scalar containing header of interest from file
my $multihit_number;                     # Scalar containing value of interest from file
my $multihit_percent = "% of reads mapped to multiple loci";# Scalar containing header of interest from file
my $multihit_percent_number;                     # Scalar containing value of interest from file
my $multihit_excluded = "Number of reads mapped to too many loci";# Scalar containing header of interest from file
my $multihit_excluded_number;                     # Scalar containing value of interest from file
my $multihit_excluded_percent = "% of reads mapped to too many loci";# Scalar containing header of interest from file
my $multihit_excluded_percent_number;                     # Scalar containing value of interest from file
my $unmapped = "Unmapped reads";    # Scalar containing header of interest from file
my $unmapped_number;                # Scalar containing value of interest from file
my $unmapped_percent = "% Unmapped reads";    # Scalar containing header of interest from file
my $unmapped_percent_number;                # Scalar containing value of interest from file

# Read in the fastq sequence file
while (<REPORT>){
    ($header, $value) = (split(/\t/));  # Read index sequence and sample name in a tab separated file
    chomp($header);
    if ($value) {
        chomp($value);
        if ($header =~ /${input}/) {
            $input_number = $value;
        }
        elsif ($header =~ /${map}/) {
            $map_number = $value;
        }
        elsif ($header =~ /${map_percent}/) {
            $map_percent_number = $value;
        }
        elsif ($header =~ /${average_length}/) {
            $map_length = $value;
        }
        elsif ($header =~ /${multihit}/) {
            $multihit_number = $value;
        }
        elsif ($header =~ /${multihit_percent}/) {
            $multihit_percent_number = $value;
        }
        elsif ($header =~ /${multihit_excluded}/) {
            $multihit_excluded_number = $value;
        }
        elsif ($header =~ /${multihit_excluded_percent}/) {
            $multihit_excluded_percent_number = $value;
        }
    }
    last if (defined $multihit_excluded_percent_number);  # If the report file was fully read, then exit reading it
}
close (REPORT);    # Close report file

$unmapped_number = ($input_number-($map_number+$multihit_number+$multihit_excluded_number));
$unmapped_percent_number = ($unmapped_number*100/$input_number);

# Open the output file and print out results
my $outfile = "All_star_log_final_out.txt";

if (-e $outfile) {
    open (OUTFILE, ">>$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">>$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
    print OUTFILE "Sample_id\t$input\t$map\t$map_percent\t$average_length\t$multihit\t$multihit_percent\t$multihit_excluded\t$multihit_excluded_percent\t$unmapped\t$unmapped_percent\n";
}

print OUTFILE "$name\t$input_number\t$map_number\t$map_percent_number\t$map_length\t$multihit_number\t$multihit_percent_number\t$multihit_excluded_number\t$multihit_excluded_percent_number\t$unmapped_number\t$unmapped_percent_number\n";

close (OUTFILE);   # Close output file

__END__