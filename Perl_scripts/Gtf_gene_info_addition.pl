#!/usr/bin/perl -w

# Copyright (C) November/2013 by Nicolas C Nalpas

# Script use to create gtf annotation files for sense, antisense and novel genes

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use List::Util 'first';

# Define the different input/output files
my $gtf;    # Input gtf file
my $chrom_length;    # Input chromosome length file
my $sense_gtf; # Output gtf file which will contain all the different sense genes info
my $antisense_gtf; # Output gtf file which will contain all the different antisense genes info
my $novel_gtf; # Output gtf file which will contain all the different novel genes info
my $all_gtf; # Output gtf file which will contain all the different genes info
my $gene_info;  # Output file which will contain all the genes info

# Define the parameter in order to submit input files to this script
&GetOptions (
    'gtf=s' => \$gtf,
    'chrom_length=s' => \$chrom_length,
    'sense_gtf=s' => \$sense_gtf,
    'antisense_gtf=s' => \$antisense_gtf,
    'novel_gtf=s' => \$novel_gtf,
    'all_gtf=s' => \$all_gtf,
    'gene_info=s' => \$gene_info,
);

my $start_date = localtime;
print STDERR ("Start = $start_date\n\n");

# Check that the input gtf file has been provided
unless ($gtf) {
    die "Please specify the gtf file containing the gene annotation via -gtf parameter!\n";
}

# Check that the input chromosome length file has been provided
unless ($chrom_length) {
    die "Please specify the chromosome length file containing the gene annotation via -chrom_length parameter!\n";
}

# Check that the sense gene output file has been provided
unless ($sense_gtf) {
    die "Please specify the output file via -sense_gtf parameter!\n";
}

# Check that the antisense gene output file has been provided
unless ($antisense_gtf) {
    die "Please specify the output file via -antisense_gtf parameter!\n";
}

# Check that the novel gene output file has been provided
unless ($novel_gtf) {
    die "Please specify the output file via -novel_gtf parameter!\n";
}

# Check that the all gene output file has been provided
unless ($all_gtf) {
    die "Please specify the output file via -all_gtf parameter!\n";
}

# Check that the gene info output file has been provided
unless ($gene_info) {
    die "Please specify the output file via -gene_info parameter!\n";
}

# Open the input gtf file and store the data in an array
open (GTF, "<$gtf") || die "Cannot open $gtf: $!\n"; $_="1";
my @gtf = <GTF>;
close (GTF);

# Open the input chromosome length file
open (CHROM_LENGTH, "<$chrom_length") || die "Cannot open $chrom_length: $!\n"; $_="1";

# Open the sense gene output file
if (-e $sense_gtf) {
    die "This file: $sense_gtf already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (SENSE_GTF, ">$sense_gtf") || die "Cannot open $sense_gtf: $!\n"; $_="1";
}

# Open the antisense gene output file
if (-e $antisense_gtf) {
    die "This file: $antisense_gtf already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (ANTISENSE_GTF, ">$antisense_gtf") || die "Cannot open $antisense_gtf: $!\n"; $_="1";
}

# Open the novel gene output file
if (-e $novel_gtf) {
    die "This file: $novel_gtf already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (NOVEL_GTF, ">$novel_gtf") || die "Cannot open $novel_gtf: $!\n"; $_="1";
}

# Open the all gene output file
if (-e $all_gtf) {
    die "This file: $all_gtf already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (ALL_GTF, ">$all_gtf") || die "Cannot open $all_gtf: $!\n"; $_="1";
}

# Open the gene info output file
if (-e $gene_info) {
    die "This file: $gene_info already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (GENE_INFO, ">$gene_info") || die "Cannot open $gene_info: $!\n"; $_="1";
}

# Print the header line into the gene info output file
print GENE_INFO "bostaurus_gene_id\tgene_biotype\tstart_position\tend_position\tstrand\tchromosome_name\tbostaurus_transcript_id\n";

# Process input chromosome length file
my %chrom_length;   # Hash containing the chromosome length
while(1) {
    chomp(my $line = <CHROM_LENGTH>);    # Read lines one by one from the chromosome length input file
    my ($chrom, $length) = (split(/\t/, $line));
    $chrom_length{$chrom} = $length;
    last if eof CHROM_LENGTH;
}

# Define variables required for reading the input file
my $line_count = 0; # Number of line already read in input file
my $gene_count = 0;  # Number of gene in gtf file
my $current_gene = "start"; # Scalar containing the current gene id being processed
my %gene;   # Hash containing all the genes information
my @current_gene;   # Array containing all information for the current gene
my %overlap;    # Hash containing all genes on a chromosome basis for overlap assessment
my $count_total = 0; # Scalar containing the total number of lines in the all gene gtf output file
my $gene_info_count = 0;    # Scalar containing the number of lines in the all genes info file
my $prom_term_count = 0;    # Scalar conataining the number of lines corresponding to promotor or terminator
my $cds_utr_count = 0;    # Scalar containing the number of lines corresponding to new CDS and UTR
my %cds;    # Hash containing the transcript ID and a value for either 5' UTR or 3' UTR

# Process input gtf file
while(1) {
    my ($line, $line_utr, $gtf_line, $chromosome, $software, $feature, $start, $end, $score, $strand, $frame, $attribute, $gene_id, $transcript_id, $gene_name, $gene_biotype);
    my @lines;  # Array containing the current lines being processed
    LINE_TYPE: {
        chomp($line = $gtf[0]);    # Read lines one by one from the gtf input file in array
        shift @gtf; # Remove the current line from the gtf array
        $line_count ++;
        if ($line =~ /\texon.*?BTAG\d{11}/) {
            ($chromosome, $software, $feature, $start, $end, $score, $strand, $frame, $attribute) = (split(/\t/, $line));
            ($gene_id = $attribute) =~ s/.*?gene_id \"(.{3}BTAG\d{11}).*/$1/;
            ($transcript_id = $attribute) =~ s/.*?transcript_id \"(.{3}BTAT\d{11})\".*/$1/;
            ($gene_biotype = $attribute) =~ s/.*?gene_biotype \"(.*?)\".*/$1/;
            if ($attribute =~ /.*?gene_name \"(.*?)\".*/) {
                $gene_name = $1;
            }
            else {
                $gene_name = "empty";
            }
            push (@lines, $line);
            if ($software eq "protein_coding") {
                unless (exists $cds{$transcript_id}{utr}) {
                    $cds{$transcript_id}{utr} = "five_prime_UTR";
                }
                UTR: {
                    undef $line_utr;
                    if ($cds{$transcript_id}{utr} eq "five_prime_UTR") {
                        if ($gtf[0] =~ /\texon\t.*?$gene_id.*?$transcript_id/) {
                            ($line_utr = $line) =~ s/(.*?\t)exon(\t.*?\;) exon_number \"\d*?\"\;(.*?\;) exon_id \".*?\"\;/$1$cds{$transcript_id}{utr}$2$3/;
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && ((($gtf[1] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && $gtf[2] =~ /\texon\t/) || ($gtf[1] =~ /\texon\t/ && $gtf[1] !~ /$transcript_id/)) || $gtf[1] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/ || ($gtf[1] =~ /\tstart_codon\t.*?$gene_id.*?$transcript_id/ && $gtf[2] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/) || ($gtf[1] =~ /\tstart_codon\t.*?$gene_id.*?$transcript_id/ && (($gtf[2] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && $gtf[3] =~ /\texon\t/) || ($gtf[2] =~ /\texon\t/ && $gtf[2] !~ /$transcript_id/))))) {
                            if ($strand eq "+") {
                                (my $cds_start = $gtf[0]) =~ s/^.*?\tCDS\t(\d*)\t.*\s/$1/;
                                unless ($start == $cds_start) {
                                    my $utr_end = ($cds_start-1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1$cds{$transcript_id}{utr}$2$utr_end$3/;
                                }
                            }
                            elsif ($strand eq "-") {
                                (my $cds_end = $gtf[0]) =~ s/^.*?\tCDS\t\d*\t(\d*)\t.*\s/$1/;
                                unless ($end == $cds_end) {
                                    my $utr_start = ($cds_end+1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1$cds{$transcript_id}{utr}\t$utr_start$2/;
                                }
                            }
                            $cds{$transcript_id}{utr} = "CDS";
                            $cds{$transcript_id}{cds} = $gtf[0];
                            if ($line_utr) {
                                push (@lines, $line_utr);
                                $cds_utr_count ++;
                            }
                            redo UTR;
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/) {
                            if ($strand eq "+") {
                                (my $cds_start = $gtf[0]) =~ s/^.*?\tCDS\t(\d*)\t.*\s/$1/;
                                unless ($start == $cds_start) {
                                    my $utr_end = ($cds_start-1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1$cds{$transcript_id}{utr}$2$utr_end$3/;
                                }
                            }
                            elsif ($strand eq "-") {
                                (my $cds_end = $gtf[0]) =~ s/^.*?\tCDS\t\d*\t(\d*)\t.*\s/$1/;
                                unless ($end == $cds_end) {
                                    my $utr_start = ($cds_end+1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1$cds{$transcript_id}{utr}\t$utr_start$2/;
                                }
                            }
                            $cds{$transcript_id}{utr} = "CDS";
                            $cds{$transcript_id}{cds} = $gtf[0];
                        }
                        else {
                            die "At line $line_count the exon for transcript $transcript_id is not followed by either an exon or a CDS, please check your file!\n";
                        }
                        if ($line_utr) {
                            push (@lines, $line_utr);
                            $cds_utr_count ++;
                        }
                    }
                    elsif ($cds{$transcript_id}{utr} eq "CDS") {
                        if ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && $gtf[1] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/) {
                            (my $stop_codon_end = $gtf[1]) =~ s/^.*?\tprotein_coding\tstop_codon\t\d*\t(\d*)\t.*\s/$1/;
                            (my $stop_codon_start = $gtf[1]) =~ s/^.*?\tprotein_coding\tstop_codon\t(\d*)\t\d*\t.*\s/$1/;
                            if (($stop_codon_end-$stop_codon_start) == 2) {
                                if ($strand eq "+") {
                                    $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t\d*\t)\d*(\t.*)\s/$1$stop_codon_end$2/;
                                    my $cds_end = $stop_codon_end;
                                    unless ($end == $cds_end) {
                                        my $utr_start = ($cds_end+1);
                                        ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1three_prime_UTR\t$utr_start$2/;
                                    }
                                }
                                elsif ($strand eq "-") {
                                    $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t)\d*(\t\d*\t.*)\s/$1$stop_codon_start$2/;
                                    my $cds_start = $stop_codon_start;
                                    unless ($start == $cds_start) {
                                        my $utr_end = ($cds_start-1);
                                        ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1three_prime_UTR$2$utr_end$3/;
                                    }
                                }
                                $cds{$transcript_id}{utr} = "three_prime_UTR";
                                chomp($gtf_line = $gtf[0]);
                                push (@lines, $gtf_line);
                                chomp($gtf_line = $gtf[1]);
                                push (@lines, $gtf_line);
                                shift @gtf;
                                shift @gtf;
                                $line_count += 2;
                            }
                            else {
                                die "The stop_codon after the $line_count line in $gtf for transcript $transcript_id does not equal 3 bases length, please check your file!\n";
                            }
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && $gtf[1] =~ /\tstart_codon\t.*?$gene_id.*?$transcript_id/ && $gtf[2] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/) {
                            (my $stop_codon_end = $gtf[2]) =~ s/^.*?\tprotein_coding\tstop_codon\t\d*\t(\d*)\t.*\s/$1/;
                            (my $stop_codon_start = $gtf[2]) =~ s/^.*?\tprotein_coding\tstop_codon\t(\d*)\t\d*\t.*\s/$1/;
                            if (($stop_codon_end-$stop_codon_start) == 2) {
                                if ($strand eq "+") {
                                    $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t\d*\t)\d*(\t.*)\s/$1$stop_codon_end$2/;
                                    my $cds_end = $stop_codon_end;
                                    unless ($end == $cds_end) {
                                        my $utr_start = ($cds_end+1);
                                        ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1three_prime_UTR\t$utr_start$2/;
                                    }
                                }
                                elsif ($strand eq "-") {
                                    $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t)\d*(\t\d*\t.*)\s/$1$stop_codon_start$2/;
                                    my $cds_start = $stop_codon_start;
                                    unless ($start == $cds_start) {
                                        my $utr_end = ($cds_start-1);
                                        ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1three_prime_UTR$2$utr_end$3/;
                                    }
                                }
                                $cds{$transcript_id}{utr} = "three_prime_UTR";
                                chomp($gtf_line = $gtf[0]);
                                push (@lines, $gtf_line);
                                chomp($gtf_line = $gtf[1]);
                                push (@lines, $gtf_line);
                                chomp($gtf_line = $gtf[2]);
                                push (@lines, $gtf_line);
                                shift @gtf;
                                shift @gtf;
                                shift @gtf;
                                $line_count += 3;
                            }
                            else {
                                die "The stop_codon after the $line_count line in $gtf for transcript $transcript_id does not equal 3 bases length, please check your file!\n";
                            }
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && ((($gtf[1] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && $gtf[2] =~ /\texon\t/) || ($gtf[1] =~ /\texon\t/ && $gtf[1] !~ /$transcript_id/)) || ($gtf[1] =~ /\tstart_codon\t.*?$gene_id.*?$transcript_id/ && (($gtf[2] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && $gtf[3] =~ /\texon\t/) || ($gtf[2] =~ /\texon\t/ && $gtf[2] !~ /$transcript_id/))))) {
                            if ($strand eq "+") {
                                (my $cds_end = $gtf[0]) =~ s/^.*?\tprotein_coding\tCDS\t\d*\t(\d*)\t.*\s/$1/;
                                unless ($end == $cds_end) {
                                    my $utr_start = ($cds_end+1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1three_prime_UTR\t$utr_start$2/;
                                }
                            }
                            elsif ($strand eq "-") {
                                (my $cds_start = $gtf[0]) =~ s/^.*?\tprotein_coding\tCDS\t(\d*)\t\d*\t.*\s/$1/;
                                unless ($start == $cds_start) {
                                    my $utr_end = ($cds_start-1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1three_prime_UTR$2$utr_end$3/;
                                }
                            }
                            $cds{$transcript_id}{utr} = "three_prime_UTR";
                            chomp($gtf_line = $gtf[0]);
                            push (@lines, $gtf_line);
                            if ($gtf[1] =~ /\tstart_codon\t.*?$gene_id.*?$transcript_id/) {
                                chomp($gtf_line = $gtf[1]);
                                push (@lines, $gtf_line);
                                shift @gtf;
                                $line_count ++;
                            }
                            shift @gtf;
                            $line_count ++;
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && $gtf[1] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && $gtf[2] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/ && $gtf[3] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/) {
                            (my $stop_codon_end = $gtf[2]) =~ s/^.*?\tprotein_coding\tstop_codon\t\d*\t(\d*)\t.*\s/$1/;
                            (my $stop_codon_start = $gtf[2]) =~ s/^.*?\tprotein_coding\tstop_codon\t(\d*)\t\d*\t.*\s/$1/;
                            if ($strand eq "+" && $start <= $stop_codon_start && $end == $stop_codon_end) {
                                $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t\d*\t)\d*(\t.*)\s/$1$stop_codon_end$2/;
                            }
                            elsif ($strand eq "-" && $start == $stop_codon_start && $end >= $stop_codon_end) {
                                $gtf[0] =~ s/^(.*?\tprotein_coding\tCDS\t)\d*(\t\d*\t.*)\s/$1$stop_codon_start$2/;
                            }
                            chomp($gtf_line = $gtf[0]);
                            push (@lines, $gtf_line);
                            chomp($gtf_line = $gtf[2]);
                            push (@lines, $gtf_line);
                            $gtf[2] = $gtf[1];
                            shift @gtf;
                            shift @gtf;
                            $line_count += 2;
                        }
                        elsif ($gtf[0] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/) {
                            my $new_cds = $gtf[0];
                            chomp($new_cds);
                            $cds{$transcript_id}{cds} =~ /.*( protein_id \"ENSBTAP\d{11}\"\;)/;
                            $new_cds = $new_cds . $1;
                            $new_cds =~ s/^(.*?\t)stop_codon\t(\d*)\t(\d*)(\t.*)/$1CDS\t$2\t$3$4/;
                            my $cds_start = $2;
                            my $cds_end = $3;
                            if ($strand eq "+") {
                                unless ($end == $cds_end) {
                                    my $utr_start = ($cds_end+1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon\t\d*(\t\d*\t.*)/$1three_prime_UTR\t$utr_start$2/;
                                }
                            }
                            elsif ($strand eq "-") {
                                unless ($start == $cds_start) {
                                    my $utr_end = ($cds_start-1);
                                    ($line_utr = $line) =~ s/^(.*?\t)exon(\t\d*\t)\d*(\t.*)/$1three_prime_UTR$2$utr_end$3/;
                                }
                            }
                            $cds{$transcript_id}{utr} = "three_prime_UTR";
                            push (@lines, $new_cds);
                            $cds_utr_count ++;
                            chomp($gtf_line = $gtf[0]);
                            push (@lines, $gtf_line);
                            shift @gtf;
                            $line_count ++;
                        }
                        elsif ($gtf[0] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ && $gtf[1] =~ /\texon\t.*?$gene_id.*?$transcript_id/ && ($gtf[2] =~ /\tCDS\t.*?$gene_id.*?$transcript_id/ || $gtf[2] =~ /\tstop_codon\t.*?$gene_id.*?$transcript_id/)) {
                            chomp($gtf_line = $gtf[0]);
                            push (@lines, $gtf_line);
                            shift @gtf;
                            $line_count ++;
                        }
                        else {
                            die "At line $line_count the exon for transcript $transcript_id is not followed by a CDS, please check your file!\n$line\n\n";
                        }
                        if ($line_utr) {
                            push (@lines, $line_utr);
                            $cds_utr_count ++;
                        }
                    }
                    elsif ($cds{$transcript_id}{utr} eq "three_prime_UTR") {
                        if ($gtf[0] =~ /\texon\t.*?$gene_id.*?$transcript_id/) {
                            ($line_utr = $line) =~ s/(.*?\t)exon(\t.*?\;) exon_number \"\d*?\"\;(.*?\;) exon_id \".*?\"\;/$1$cds{$transcript_id}{utr}$2$3/;
                        }
                        elsif ($gtf[0] =~ /\texon\t/ && $gtf[0] !~ /$transcript_id/) {
                            ($line_utr = $line) =~ s/(.*?\t)exon(\t.*?\;) exon_number \"\d*?\"\;(.*?\;) exon_id \".*?\"\;/$1$cds{$transcript_id}{utr}$2$3/;
                            $cds{$transcript_id}{utr} = "next_transcript";
                        }
                        else {
                            die "At line $line_count the exon for transcript $transcript_id is not followed by an exon, please check your file!\n";
                        }
                        if ($line_utr) {
                            push (@lines, $line_utr);
                            $cds_utr_count ++;
                        }
                    }
                    else {
                        die "At line $line_count the transcript $transcript_id has unacceptable value $cds{$transcript_id}{utr}, please check your file!\n";
                    }
                }
            }
        }
        elsif ($line =~ /$current_gene/) {
            push (@current_gene, $line);
            redo LINE_TYPE;
        }
        else {
            die "Cannot find a gene ID in $gtf file at line $line_count, please seek advice!\n";
        }
    }
    if ($gene_id eq $current_gene) {
        @current_gene = (@current_gene, @lines);
        if ($start < $gene{$current_gene}{start}) {
            $gene{$current_gene}{start} = $start;
        }
        if ($end > $gene{$current_gene}{end}) {
            $gene{$current_gene}{end} = $end;
        }
        if ($gene{$current_gene}{transcript} !~ /$transcript_id/) {
            $gene{$current_gene}{transcript} = $gene{$current_gene}{transcript} . ', ' . $transcript_id;
        }
    }
    elsif ($gene_id ne $current_gene) {
        unless ($current_gene eq "start") {
            GENE: {
                my $gene_line;
                print GENE_INFO "$current_gene\t$gene{$current_gene}{gene_biotype}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{chromosome}\t$gene{$current_gene}{transcript}\n"; # Print in the gene info output file the info for sense and novel genes
                $gene_info_count ++;
                if ($gene{$current_gene}{gene_name} ne "empty") {
                    $gene_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\t$gene{$current_gene}{feature}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
                else {
                    $gene_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\t$gene{$current_gene}{feature}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
                unshift (@current_gene, $gene_line);
            }
            PROMOTER: {
                unless ($current_gene =~ /NOVBTAG\d{11}/) {
                    my $current_chromo = $gene{$current_gene}{chromosome};
                    if ($gene{$current_gene}{strand} eq "+" || $gene{$current_gene}{strand} eq "-") {
                        if ($gene{$current_gene}{start} == 1) {
                            $gene{$current_gene}{prom_start} = "no_prom";
                            $gene{$current_gene}{prom_end} = "no_prom";
                        }
                        elsif (($gene{$current_gene}{start}-2000) < 1) {
                            $gene{$current_gene}{prom_start} = 1;
                            $gene{$current_gene}{prom_end} = ($gene{$current_gene}{start}-1);
                        }
                        else {
                            $gene{$current_gene}{prom_start} = ($gene{$current_gene}{start}-2000);
                            $gene{$current_gene}{prom_end} = ($gene{$current_gene}{start}-1);
                        }
                        if ($gene{$current_gene}{end} == $chrom_length{$current_chromo}) {
                            $gene{$current_gene}{term_start} = "no_term";
                            $gene{$current_gene}{term_end} = "no_term";
                        }
                        elsif (($gene{$current_gene}{end}+2000) > $chrom_length{$current_chromo}) {
                            $gene{$current_gene}{term_start} = ($gene{$current_gene}{end}+1);
                            $gene{$current_gene}{term_end} = $chrom_length{$current_chromo};
                        }
                        else {
                            $gene{$current_gene}{term_start} = ($gene{$current_gene}{end}+1);
                            $gene{$current_gene}{term_end} = ($gene{$current_gene}{end}+2000);
                        }
                    }
                    else {
                        die "The gene $current_gene at line has a non accepted value for strand: $gene{$current_gene}{strand}, please check your file!\n";
                    }
                    my ($prom_line, $term_line);
                    if ($gene{$current_gene}{strand} eq "+") {
                        if ($gene{$current_gene}{gene_name} ne "empty") {
                            $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                            $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                        }
                        else {
                            $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                            $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                        }
                    }
                    elsif ($gene{$current_gene}{strand} eq "-") {
                        if ($gene{$current_gene}{gene_name} ne "empty") {
                            $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                            $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                        }
                        else {
                            $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                            $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                        }
                    }
                    unless ($gene{$current_gene}{prom_start} eq "no_prom") {
                        unshift (@current_gene, $prom_line);
                        $prom_term_count ++;
                    }
                    unless ($gene{$current_gene}{term_start} eq "no_term") {
                        push (@current_gene, $term_line);
                        $prom_term_count ++;
                    }
                }
            }
            if ($current_gene =~ /ENSBTAG/) {
                foreach (@current_gene) {
                    $count_total ++;
                    print SENSE_GTF "$_\n";
                    print ALL_GTF "$_\n";
                }
            }
            elsif ($current_gene =~ /NOVBTAG/) {
                foreach (@current_gene) {
                    $count_total ++;
                    print NOVEL_GTF "$_\n";
                    print ALL_GTF "$_\n";
                }
            }
            else {
                die "The gene id $current_gene at line $line_count of the input file $gtf is not a correct gene ID, please check your file!\n";
            }
        }
        if (exists $gene{$gene_id}) {
            die "The gene $gene_id is already present in the hash generated by this perl script, each gene present in the gtf file should have unique ID, please seek advice!\n";
        }
        else {
            $current_gene = $gene_id;
            undef @current_gene;
            @current_gene = @lines;
            $gene_count ++;
            $gene{$current_gene}{chromosome} = $chromosome;
            $gene{$current_gene}{software} = $software;
            $gene{$current_gene}{feature} = "gene";
            $gene{$current_gene}{start} = $start;
            $gene{$current_gene}{end} = $end;
            $gene{$current_gene}{score} = $score;
            $gene{$current_gene}{strand} = $strand;
            $gene{$current_gene}{frame} = $frame;
            $gene{$current_gene}{transcript} = $transcript_id;
            $gene{$current_gene}{gene_name} = $gene_name;
            $gene{$current_gene}{gene_biotype} = $gene_biotype;
            unless (exists $overlap{$chromosome}{$current_gene} || $current_gene =~ /NOVBTAG\d{11}/) {
                $overlap{$chromosome}{$current_gene} = 0;
            }
        }
    }
    unless (@gtf) {
        my $gene_line;
        print GENE_INFO "$current_gene\t$gene{$current_gene}{gene_biotype}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{chromosome}\t$gene{$current_gene}{transcript}\n"; # Print in the gene info output file the info for sense and novel genes
        $gene_info_count ++;
        if (exists $gene{$current_gene}{gene_name}) {
            $gene_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\t$gene{$current_gene}{feature}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";    
        }
        else {
            $gene_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\t$gene{$current_gene}{feature}\t$gene{$current_gene}{start}\t$gene{$current_gene}{end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
        }
        unshift (@current_gene, $gene_line);
        unless ($current_gene =~ /NOVBTAG\d{11}/) {
            my $current_chromo = $gene{$current_gene}{chromosome};
            if ($gene{$current_gene}{strand} eq "+" || $gene{$current_gene}{strand} eq "-") {
                if ($gene{$current_gene}{start} == 1) {
                    $gene{$current_gene}{prom_start} = "no_prom";
                    $gene{$current_gene}{prom_end} = "no_prom";
                }
                elsif (($gene{$current_gene}{start}-2000) < 1) {
                    $gene{$current_gene}{prom_start} = 1;
                    $gene{$current_gene}{prom_end} = ($gene{$current_gene}{start}-1);
                }
                else {
                    $gene{$current_gene}{prom_start} = ($gene{$current_gene}{start}-2000);
                    $gene{$current_gene}{prom_end} = ($gene{$current_gene}{start}-1);
                }
                if ($gene{$current_gene}{end} == $chrom_length{$current_chromo}) {
                    $gene{$current_gene}{term_start} = "no_term";
                    $gene{$current_gene}{term_end} = "no_term";
                }
                elsif (($gene{$current_gene}{end}+2000) > $chrom_length{$current_chromo}) {
                    $gene{$current_gene}{term_start} = ($gene{$current_gene}{end}+1);
                    $gene{$current_gene}{term_end} = $chrom_length{$current_chromo};
                }
                else {
                    $gene{$current_gene}{term_start} = ($gene{$current_gene}{end}+1);
                    $gene{$current_gene}{term_end} = ($gene{$current_gene}{end}+2000);
                }
            }
            else {
                die "The gene $current_gene at line has a non accepted value for strand: $gene{$current_gene}{strand}, please check your file!\n";
            }
            my ($prom_line, $term_line);
            if ($gene{$current_gene}{strand} eq "+") {
                if ($gene{$current_gene}{gene_name} ne "empty") {
                    $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                    $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
                else {
                    $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                    $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
            }
            elsif ($gene{$current_gene}{strand} eq "-") {
                if ($gene{$current_gene}{gene_name} ne "empty") {
                    $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                    $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_name \"$gene{$current_gene}{gene_name}\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
                else {
                    $prom_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tpromotor\t$gene{$current_gene}{term_start}\t$gene{$current_gene}{term_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                    $term_line = "$gene{$current_gene}{chromosome}\t$gene{$current_gene}{software}\tterminator\t$gene{$current_gene}{prom_start}\t$gene{$current_gene}{prom_end}\t$gene{$current_gene}{score}\t$gene{$current_gene}{strand}\t$gene{$current_gene}{frame}\tgene_id \"$current_gene\"\; gene_biotype \"$gene{$current_gene}{gene_biotype}\"\;";
                }
            }
            unless ($gene{$current_gene}{prom_start} eq "no_prom") {
                unshift (@current_gene, $prom_line);
                $prom_term_count ++;
            }
            unless ($gene{$current_gene}{term_start} eq "no_term") {
                push (@current_gene, $term_line);
                $prom_term_count ++;
            }
        }
        if ($current_gene =~ /ENSBTAG/) {
            foreach (@current_gene) {
                $count_total ++;
                print SENSE_GTF "$_\n";
                print ALL_GTF "$_\n";
            } 
        }
        elsif ($current_gene =~ /NOVBTAG/) {
            foreach (@current_gene) {
                $count_total ++;
                print NOVEL_GTF "$_\n";
                print ALL_GTF "$_\n";
            }
        }
        last;  # If the gtf input file was fully read, then exit reading the file
    }
}
close (SENSE_GTF);  # Close the sense gene gtf file
close (NOVEL_GTF);  # Close the novel gene gtf file
close (ALL_GTF);  # Close the all gene gtf file (containing sense and novel genes info)

# Check genes for overlap in order to define antisense genes
my %antisense;  # Hash containing the antisense gene info
foreach my $gene1 (keys %gene) {
    unless ($gene1 =~ /NOVBTAG\d{11}/) {
        my $chromo = $gene{$gene1}{chromosome};
        foreach my $gene2 (keys %{$overlap{$chromo}}) {
            unless ($gene1 eq $gene2) {
                my ($start1, $start2, $end1, $end2);
                if (($gene{$gene1}{start}-5000) < 1) {
                    $start1 = 1;
                }
                else {
                    $start1 = ($gene{$gene1}{start}-5000);
                }
                if (($gene{$gene2}{start}-5000) < 1) {
                    $start2 = 1;
                }
                else {
                    $start2 = ($gene{$gene2}{start}-5000);
                }
                if (($gene{$gene1}{end}+5000) > $chrom_length{$chromo}) {
                    $end1 = $chrom_length{$chromo};
                }
                else{
                    $end1 = ($gene{$gene1}{end}+5000);
                }
                if (($gene{$gene2}{end}+5000) > $chrom_length{$chromo}) {
                    $end2 = $chrom_length{$chromo};
                }
                else {
                    $end2 = ($gene{$gene2}{end}+5000);
                }
                if ($gene{$gene2}{start} >= $start1 && $gene{$gene2}{start} <= $end1) {
                    $overlap{$chromo}{$gene1} ++;
                }
                elsif ($start2 <= $gene{$gene1}{start} && $end2 >= $gene{$gene1}{start}) {
                    $overlap{$chromo}{$gene1} ++;
                }
            }
        }
    }
}

# Define the info for antisense gene (passing overlap filtering)
my $numb_antisense = 0; 
foreach my $chromo (keys %overlap) {
    foreach my $key (keys %{$overlap{$chromo}}) {
        if ($overlap{$chromo}{$key} == 0) {
            $numb_antisense ++;
            $antisense{$key}{chromosome} = $gene{$key}{chromosome};
            $antisense{$key}{software} = "antisense";
            $antisense{$key}{feature} = "gene";
            if (($gene{$key}{start}-2000) < 1) {
                $antisense{$key}{start} = 1;
            }
            else {
                $antisense{$key}{start} = ($gene{$key}{start}-2000);
            }
            if (($gene{$key}{end}+2000) > $chrom_length{$chromo}) {
                $antisense{$key}{end} = $chrom_length{$chromo};
            }
            else {
                $antisense{$key}{end} = ($gene{$key}{end}+2000);
            }
            $antisense{$key}{score} = $gene{$key}{score};
            if ($gene{$key}{strand} eq "+") {
                $antisense{$key}{strand} = "-";
            }
            elsif ($gene{$key}{strand} eq "-") {
                $antisense{$key}{strand} = "+";
            }
            $antisense{$key}{frame} = $gene{$key}{frame};
            $antisense{$key}{gene_biotype} = "antisense";
            print ANTISENSE_GTF "$antisense{$key}{chromosome}\t$antisense{$key}{software}\t$antisense{$key}{feature}\t$antisense{$key}{start}\t$antisense{$key}{end}\t$antisense{$key}{score}\t$antisense{$key}{strand}\t$antisense{$key}{frame}\tgene_id \"NAT$key\"\; gene_biotype \"$antisense{$key}{gene_biotype}\"\;\n";
            print GENE_INFO "NAT$key\t$antisense{$key}{gene_biotype}\t$antisense{$key}{start}\t$antisense{$key}{end}\t$antisense{$key}{strand}\t$antisense{$key}{chromosome}\n";
            $gene_info_count ++;
        }
    }
}
close (GENE_INFO); # Close the output gene info file

print STDERR "There were $line_count lines from file $gtf, corresponding to $gene_count number of sense and novel genes and $numb_antisense number of antisense genes!\n\n";
print STDERR "Antisense were determined based on full length sense gene information plus an extra 2000 bases prior 5'-end and also after 3'-end (to be consider as promoter and terminator region), antisense gene strand was also define as opposite to his sense strand; any neighbours genes overlapping or within 5000 bases of each other were ignored for the antisense genes identification!\n\n";
if ($gene_info_count == ($gene_count+$numb_antisense) && $count_total == ($line_count+$gene_count+$prom_term_count+$cds_utr_count)) {
    print STDERR "Process completed successfully!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR ("Finish = $finish_date\n\n");

__END__