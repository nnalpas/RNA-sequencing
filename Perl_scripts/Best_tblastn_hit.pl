#!/usr/bin/perl -w

# Script used to find best hit from a tblastn run within BLAST, and then identify RBH

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $blast;    # Input file containing tblastn results
my $homo_gtf;    # Input file containing gtf annotation from Homo sapiens
my $bos_gtf;    # Input file containing gtf annotation from Bos taurus
my $rbh;    # Output file containing RBH results

# Define the parameter in order to submit input files to this script
&GetOptions (
    'blast=s' => \$blast,
    'homo_gtf=s' => \$homo_gtf,
    'bos_gtf=s' => \$bos_gtf,
    'rbh=s' => \$rbh,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input tblastn file
unless ($blast) {
    die "Please specify the input tblastn file in csv format via -tblastn parameter!\n";
}
open (BLAST, "<$blast") || die "Cannot open $blast: $!\n"; $_="1";

# Open the input gtf annotation file for Homo sapiens
unless ($homo_gtf) {
    die "Please specify the input Homo sapiens annotation file in gtf format via -homo_gtf parameter!\n";
}
open (HOMOGTF, "<$homo_gtf") || die "Cannot open $homo_gtf: $!\n"; $_="1";

# Open the input gtf annotation file for Bos taurus
unless ($bos_gtf) {
    die "Please specify the input Bos taurus annotation file in gtf format via -bos_gtf parameter!\n";
}
open (BOSGTF, "<$bos_gtf") || die "Cannot open $bos_gtf: $!\n"; $_="1";

# Open the ouput file which will contain the reciprocal best hit results
if (-e $rbh) {
    die "This file: $rbh already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (RBH, ">$rbh") || die "Cannot open $rbh: $!\n"; $_="1";
}

# Define output file of the best hits for tblastn
my $name = `basename $blast`;
chomp ($name);
$name =~ s/(.*)\.tblastn/$1/;

# Open the ouput file which will contain the best hit blast data for each novel transcript
my $outfile = "${name}_besthit.csv";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

# Open the Homo sapiens gtf file and for each protein ID get its associated gene ID
my %protein;
while (1) {
    chomp (my $line = <HOMOGTF>);
    if ($line =~ /^(.*?)\t.*?gene_id \"(ENSG\d{11})\";.*?gene_name \"(.*?)\";.*?protein_id \"(ENSP\d{11})\";/) {
        unless (exists $protein{$3}) {
            $protein{$4}{gene} = $2;
            $protein{$4}{chromosome} = $1;
            $protein{$4}{name} = $3;
        }
    }
    last if eof HOMOGTF;
}
close (HOMOGTF);

# Open the Bos taurus gtf file and for each transcript ID get its associated gene ID
my %transcript;
while (1) {
    chomp (my $line = <BOSGTF>);
    if ($line =~ /^(.*?)\t.*?gene_id \"(NOVBTAG\d{11}|ENSBTAG\d{11})\";.*?transcript_id \"(NOVBTAT\d{11}|ENSBTAT\d{11})\";/) {
        unless (exists $transcript{$3}) {
            $transcript{$3}{gene} = $2;
            $transcript{$3}{chromosome} = $1;
            $transcript{$3}{name} = "none";
        }
    }
    if ($line =~ /transcript_id \"(ENSBTAT\d{11})\";.*?gene_name \"(.*?)\";/) {
        $transcript{$1}{name} = $2;
    }
    last if eof BOSGTF;
}
close (BOSGTF);

# Define variables required for reading input file
my $line;
my %best_blast;
my $ref_line = 0;
my $novel_number = 0;
my $ambiguous = 0;
my $better = 0;
my $badder = 0;
my $multihit = 0;

# Read the input file one line at a time and find the hit with the lowest e-value, it will be the best hit
while (1) {
    chomp (my $line = <BLAST>);
    my ($query_id, $subject_id, $identity_perc, $length, $substitution_num, $insertion_num, $query_start, $query_end, $subject_start, $subject_end, $e_value, $score) = (split(/,/, $line));
    $ref_line ++;
    if (exists $best_blast{$query_id}) {
        if (($e_value < $best_blast{$query_id}{e_value})  || ($e_value == $best_blast{$query_id}{e_value} && $score > $best_blast{$query_id}{score}) || ($e_value == $best_blast{$query_id}{e_value} && $score == $best_blast{$query_id}{score} && $identity_perc > $best_blast{$query_id}{identity})) {
            $best_blast{$query_id}{subject} = $subject_id;
            $best_blast{$query_id}{identity} = $identity_perc;
            $best_blast{$query_id}{length} = $length;
            $best_blast{$query_id}{substitution} = $substitution_num;
            $best_blast{$query_id}{insertion} = $insertion_num;
            $best_blast{$query_id}{query_start} = $query_start;
            $best_blast{$query_id}{query_end} = $query_end;
            $best_blast{$query_id}{subject_start} = $subject_start;
            $best_blast{$query_id}{subject_end} = $subject_end;
            $best_blast{$query_id}{e_value} = $e_value;
            $best_blast{$query_id}{score} = $score;
            $best_blast{$query_id}{bos_gene} = $transcript{$subject_id}{gene};
            $best_blast{$query_id}{bos_chromosome} = $transcript{$subject_id}{chromosome};
            $better ++;
        }
        elsif ($e_value == $best_blast{$query_id}{e_value} && $score == $best_blast{$query_id}{score} && $identity_perc == $best_blast{$query_id}{identity}) {
            if ($transcript{$subject_id}{gene} eq $best_blast{$query_id}{bos_gene}) {
                if ($best_blast{$query_id}{same_gene} eq "none") {
                    $best_blast{$query_id}{same_gene} = $subject_id;
                }
                else {
                    $best_blast{$query_id}{same_gene} = $best_blast{$query_id}{same_gene} . "," . $subject_id;   
                }
                $multihit ++;
            }
            elsif ($transcript{$subject_id}{chromosome} =~ /^GJ/ && $best_blast{$query_id}{bos_chromosome} !~ /^GJ/) {
                $badder ++;
            }
            elsif ($transcript{$subject_id}{chromosome} !~ /^GJ/ && $best_blast{$query_id}{bos_chromosome} =~ /^GJ/) {
                $best_blast{$query_id}{subject} = $subject_id;
                $best_blast{$query_id}{identity} = $identity_perc;
                $best_blast{$query_id}{length} = $length;
                $best_blast{$query_id}{substitution} = $substitution_num;
                $best_blast{$query_id}{insertion} = $insertion_num;
                $best_blast{$query_id}{query_start} = $query_start;
                $best_blast{$query_id}{query_end} = $query_end;
                $best_blast{$query_id}{subject_start} = $subject_start;
                $best_blast{$query_id}{subject_end} = $subject_end;
                $best_blast{$query_id}{e_value} = $e_value;
                $best_blast{$query_id}{score} = $score;
                $best_blast{$query_id}{bos_gene} = $transcript{$subject_id}{gene};
                $best_blast{$query_id}{bos_chromosome} = $transcript{$subject_id}{chromosome};
                $best_blast{$query_id}{bos_name} = $transcript{$subject_id}{name};
                $best_blast{$query_id}{same_gene} = "none";
                $better ++;
            }
            elsif ($transcript{$subject_id}{name} ne "none" || $best_blast{$query_id}{bos_name} ne "none") {
                if ($transcript{$subject_id}{name} eq $best_blast{$query_id}{bos_name} && (($transcript{$subject_id}{chromosome} =~ /^GJ/ && $best_blast{$query_id}{bos_chromosome} =~ /^GJ/) || ($transcript{$subject_id}{chromosome} eq $best_blast{$query_id}{bos_chromosome}))) {
                    if ($best_blast{$query_id}{same_gene} eq "none") {
                        $best_blast{$query_id}{same_gene} = $subject_id;
                    }
                    else {
                        $best_blast{$query_id}{same_gene} = $best_blast{$query_id}{same_gene} . "," . $subject_id;
                    }
                    $multihit ++;
                }
                else {
                    $ambiguous ++;
                }
            }
            else {
                $best_blast{$query_id}{subject} = "ambiguous";
                $ambiguous ++;
            }
        }
        else {
            $badder ++;
        }
    }
    else {
        $best_blast{$query_id}{subject} = $subject_id;
        $best_blast{$query_id}{identity} = $identity_perc;
        $best_blast{$query_id}{length} = $length;
        $best_blast{$query_id}{substitution} = $substitution_num;
        $best_blast{$query_id}{insertion} = $insertion_num;
        $best_blast{$query_id}{query_start} = $query_start;
        $best_blast{$query_id}{query_end} = $query_end;
        $best_blast{$query_id}{subject_start} = $subject_start;
        $best_blast{$query_id}{subject_end} = $subject_end;
        $best_blast{$query_id}{e_value} = $e_value;
        $best_blast{$query_id}{score} = $score;
        $best_blast{$query_id}{bos_gene} = $transcript{$subject_id}{gene};
        $best_blast{$query_id}{bos_chromosome} = $transcript{$subject_id}{chromosome};
        $best_blast{$query_id}{bos_name} = $transcript{$subject_id}{name};
        $best_blast{$query_id}{same_gene} = "none";
        $novel_number ++;
    }
    last if eof BLAST;
}
close (BLAST);      # Close the input file

# Print out the results in the tblastn output files
my $outfile_count = 0;
my $ambiguous_novel = 0;
foreach my $keys (keys %best_blast) {
    if ($best_blast{$keys}{subject} ne "ambiguous") {
        print OUTFILE "$keys,$best_blast{$keys}{subject},$best_blast{$keys}{identity},$best_blast{$keys}{e_value},$best_blast{$keys}{score},$best_blast{$keys}{same_gene}\n";
        $outfile_count ++;
    }
    else {
        $ambiguous_novel ++;
    }
}
close (OUTFILE);    # Close tblastn output file containing the best match blast data

# Find the Reciprocal Best Hits
my %rbh;
my $non_reciprocal = 0;
foreach my $keys (keys %best_blast) {
    unless ($best_blast{$keys}{subject} eq "ambiguous") {
        my $homo_gene = $protein{$keys}{gene};
        unless (exists $rbh{$homo_gene}) {
            $rbh{$homo_gene}{homo_protein} = $keys;
            $rbh{$homo_gene}{identity} = $best_blast{$keys}{identity};
            $rbh{$homo_gene}{e_value} = $best_blast{$keys}{e_value};
            $rbh{$homo_gene}{score} = $best_blast{$keys}{score};
            my $value = $best_blast{$keys}{subject};
            my @value;
            unless ($best_blast{$keys}{same_gene} eq "none") {
                @value = (split(/,/, $best_blast{$keys}{same_gene}));
            }
            push @value, $value;
            foreach my $value (@value) {
                unless ($rbh{$homo_gene}{bos_gene}) {
                    $rbh{$homo_gene}{bos_gene} = $transcript{$value}{gene};
                    $rbh{$homo_gene}{bos_transcript} = $value;
                }
                elsif ($rbh{$homo_gene}{bos_gene} eq $transcript{$value}{gene}) {
                    $rbh{$homo_gene}{bos_transcript} = $rbh{$homo_gene}{bos_transcript} . "," . $value;
                }
                else {
                    $non_reciprocal ++;
                    $rbh{$homo_gene}{bos_gene} = "ambiguous";
                }
            }
        }
        else {
            if ($best_blast{$keys}{e_value} < $rbh{$homo_gene}{e_value} || ($best_blast{$keys}{e_value} == $rbh{$homo_gene}{e_value} && $best_blast{$keys}{score} > $rbh{$homo_gene}{score}) || ($best_blast{$keys}{e_value} == $rbh{$homo_gene}{e_value} && $best_blast{$keys}{score} == $rbh{$homo_gene}{score} && $best_blast{$keys}{identity} > $rbh{$homo_gene}{identity})) {
                $rbh{$homo_gene}{homo_protein} = $keys;
                $rbh{$homo_gene}{identity} = $best_blast{$keys}{identity};
                $rbh{$homo_gene}{e_value} = $best_blast{$keys}{e_value};
                $rbh{$homo_gene}{score} = $best_blast{$keys}{score};
                my $value = $best_blast{$keys}{subject};
                my @value;
                unless ($best_blast{$keys}{same_gene} eq "none") {
                    @value = (split(/,/, $best_blast{$keys}{same_gene}));
                }
                push @value, $value;
                foreach my $value (@value) {
                    unless ($rbh{$homo_gene}{bos_gene}) {
                        $rbh{$homo_gene}{bos_gene} = $transcript{$value}{gene};
                        $rbh{$homo_gene}{bos_transcript} = $value;
                    }
                    elsif ($rbh{$homo_gene}{bos_gene} eq $transcript{$value}{gene}) {
                        $rbh{$homo_gene}{bos_transcript} = $rbh{$homo_gene}{bos_transcript} . "," . $value;
                    }
                    else {
                        $non_reciprocal ++;
                        $rbh{$homo_gene}{bos_gene} = "ambiguous";
                    }
                }
            }
            elsif ($best_blast{$keys}{e_value} == $rbh{$homo_gene}{e_value} && $best_blast{$keys}{score} == $rbh{$homo_gene}{score} && $best_blast{$keys}{identity} == $rbh{$homo_gene}{identity}) {
                my $value = $best_blast{$keys}{subject};
                my @value;
                unless ($best_blast{$keys}{same_gene} eq "none") {
                    @value = (split(/,/, $best_blast{$keys}{same_gene}));
                }
                push @value, $value;
                foreach my $value (@value) {
                    unless ($rbh{$homo_gene}{bos_gene}) {
                        $rbh{$homo_gene}{bos_gene} = $transcript{$value}{gene};
                        $rbh{$homo_gene}{bos_transcript} = $value;
                    }
                    elsif ($rbh{$homo_gene}{bos_gene} eq $transcript{$value}{gene}) {
                        $rbh{$homo_gene}{homo_protein} = $rbh{$homo_gene}{homo_protein} . "," . $keys;
                        $rbh{$homo_gene}{bos_transcript} = $rbh{$homo_gene}{bos_transcript} . "," . $value;
                    }
                    else {
                        $non_reciprocal ++;
                        $rbh{$homo_gene}{bos_gene} = "ambiguous";
                    }
                }
            }
        }
    }
}

# Print out the results in the tblastn output files
my $rbh_count = 0;
my $rbh_ambiguous = 0;
my $known_gene = 0;
foreach my $keys (keys %rbh) {
    if ($rbh{$keys}{bos_gene} =~ /NOVBTAG/) {
        print RBH "$keys\t$rbh{$keys}{homo_protein}\t$rbh{$keys}{bos_gene}\t$rbh{$keys}{bos_transcript}\t$rbh{$keys}{identity}\t$rbh{$keys}{e_value}\t$rbh{$keys}{score}\n";
        $rbh_count ++;
    }
    elsif ($rbh{$keys}{bos_gene} =~ /ENSBTAG/) {
        $known_gene ++;
    }
    elsif ($rbh{$keys}{bos_gene} eq "ambiguous") {
        $rbh_ambiguous ++;
    }
    else {
        die "$keys\t$rbh{$keys}{homo_protein}\t$rbh{$keys}{bos_gene}\t$rbh{$keys}{bos_transcript}\t$rbh{$keys}{identity}\t$rbh{$keys}{e_value}\t$rbh{$keys}{score}\n";
    }
}
close (RBH);    # Close RBH output file containing the reciprocal best hits


print STDERR "This perl script processed $ref_line number of blast matches for a total of $novel_number number of proteins, including $ambiguous_novel number of best ids ties!\n\n";
print STDERR "This perl script processed $ref_line number of blast matches for a total of $rbh_count number of identified RBH, as well as $rbh_ambiguous number of ambiguous RBH and $known_gene number of known Bos taurus genes!\n\n";
if ($ref_line == ($novel_number+$better+$ambiguous+$badder+$multihit) && $outfile_count == ($novel_number-$ambiguous_novel)) {
    print STDERR "Process completed successfully!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__