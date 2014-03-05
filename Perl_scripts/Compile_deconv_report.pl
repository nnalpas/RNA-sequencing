#!/usr/bin/perl -w

# Script used to compile

use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $list_reports;    # Input file containing the list of all reports from decconvolution perl script
my $output;     # Output file which will contain all combined reports

# Define the parameter in order to submit input files to this script
&GetOptions (
    'list_reports=s' => \$list_reports,
    'output=s' => \$output,
);

# Open the input fasta file
unless ($list_reports) {
    die "Please specify the file containing the list of report to compile in a single report via -list_reports parameter!\n";
}
open (LIST_REPORTS, "<$list_reports") || die "Cannot open $list_reports: $!\n"; $_="1";

unless ($output) {
    die "Please specify the output file which will contain compiled reports via -output parameter!\n";
}
if (-e $output) {
    die "The output file which you specified already exists, this perl script cannot overwrite an existing file, please specify another output file name or delete the existing file!\n";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Define all necessary variables
my %deconvoluted_reads_counts;  # Hash which will contain number of reads attributed to individual samples also showing their index sequence
my $index;
my $sample;
my $number_reads;
my $number_reads_kept;
my $percent1;
my $number_reads_adapter;
my $percent2;
my $number_reads_quality;
my $percent3;

while (1) {
    chomp (my $report = <LIST_REPORTS>);
    open (REPORT, "<$report") || die "Cannot open $report: $!\n"; $_="1";       # Open all report listed in the list_report file one after the other to collect data from them
    PROCESS_REPORT: while (1) {
        chomp (my $line = <REPORT>);
        if  ($line =~ /Index\ssequence/) {
            while (1) {
                chomp (my $sample_line = <REPORT>);
                if ($sample_line eq "") {       # Close the report file and exit the current loop if all indices from the report have been processed
                    close (REPORT);
                    last PROCESS_REPORT;
                }
                else {
                    ($index, $sample, $number_reads, $number_reads_kept, $percent1, $number_reads_adapter, $percent2, $number_reads_quality, $percent3) = (split(/\t/, $sample_line));    # Split the sample line into the different variables
                    if (exists ($deconvoluted_reads_counts{$sample})){       # If the index is already into the hash add up the number to previous value
                        $deconvoluted_reads_counts{$sample}{number_reads} += $number_reads;
                        $deconvoluted_reads_counts{$sample}{number_reads_kept} += $number_reads_kept;
                        $deconvoluted_reads_counts{$sample}{number_reads_adapter} += $number_reads_adapter;
                        $deconvoluted_reads_counts{$sample}{number_reads_quality} += $number_reads_quality;
                    }
                    else {      # Otherwise the index is entered into the hash for the first time and the sample id and current number values are defined as well
                        $deconvoluted_reads_counts{$sample}{index} = $index;
                        $deconvoluted_reads_counts{$sample}{number_reads} = $number_reads;
                        $deconvoluted_reads_counts{$sample}{number_reads_kept} = $number_reads_kept;
                        $deconvoluted_reads_counts{$sample}{number_reads_adapter} = $number_reads_adapter;
                        $deconvoluted_reads_counts{$sample}{number_reads_quality} = $number_reads_quality;
                    }
                }
            }
        }
    }
    last if eof (LIST_REPORTS);
}
close (LIST_REPORTS);

print OUTPUT "\nSample id\tIndex sequence\tNo. reads pre-filtering\tNo. reads post-filtering\tPercentage reads post-filtering\tNo. reads discarded due to adapter sequence\tPercentage reads with adapter sequence\tNo. reads discarded due to poor overall quality\tPercentage reads with poor overall quality\n";
foreach my $key (sort keys %deconvoluted_reads_counts) {
    if ($key eq "Ignore") {
        $deconvoluted_reads_counts{$key}{number_reads_kept} = 0;
    }
    unless ($deconvoluted_reads_counts{$key}{number_reads_adapter}) {
        $deconvoluted_reads_counts{$key}{number_reads_adapter} = 0;
    }
    unless ($deconvoluted_reads_counts{$key}{number_reads_quality}) {
        $deconvoluted_reads_counts{$key}{number_reads_quality} = 0;
    }
    my $number_reads_kept_percent = ($deconvoluted_reads_counts{$key}{number_reads_kept}*100/$deconvoluted_reads_counts{$key}{number_reads});
    my $number_reads_adapter_percent = ($deconvoluted_reads_counts{$key}{number_reads_adapter}*100/$deconvoluted_reads_counts{$key}{number_reads});
    my $number_reads_quality_percent = ($deconvoluted_reads_counts{$key}{number_reads_quality}*100/$deconvoluted_reads_counts{$key}{number_reads});
    print OUTPUT "$key\t$deconvoluted_reads_counts{$key}{index}\t$deconvoluted_reads_counts{$key}{number_reads}\t$deconvoluted_reads_counts{$key}{number_reads_kept}\t$number_reads_kept_percent\t$deconvoluted_reads_counts{$key}{number_reads_adapter}\t$number_reads_adapter_percent\t$deconvoluted_reads_counts{$key}{number_reads_quality}\t$number_reads_quality_percent\n";
}
print OUTPUT "\n";

__END__