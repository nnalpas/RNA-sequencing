#!/usr/bin/perl -w

# Script used to select best hit from a blastx run within BLAST, in order to then perform the reciprocal blast

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $blast;    # Input file containing blast results
my $homo_gtf;    # Input file containing gtf annotation from Homo sapiens

# Define the parameter in order to submit input files to this script
&GetOptions (
    'blast=s' => \$blast,
    'homo_gtf=s' => \$homo_gtf,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input blastx file
unless ($blast) {
    die "Please specify the input blast file in csv format via -blast parameter!\n";
}
open (BLAST, "<$blast") || die "Cannot open $blast: $!\n"; $_="1";

# Open the input gtf annotation file for Homo sapiens
unless ($homo_gtf) {
    die "Please specify the input Homo sapiens annotation file in gtf format via -homo_gtf parameter!\n";
}
open (HOMOGTF, "<$homo_gtf") || die "Cannot open $homo_gtf: $!\n"; $_="1";

# Define output file
my $name = `basename $blast`;
chomp ($name);
$name =~ s/(.*)\.blastx/$1/;

# Open the ouput file which will contain the best hit blast data for each novel transcript
my $outfile = "${name}_besthit.csv";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

# Open the ouput file which will contain the best hit id for each novel transcript
my $out_protein = "${name}_hitid.txt";
if (-e $out_protein) {
    die "This file: $out_protein already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPROTEIN, ">$out_protein") || die "Cannot open $out_protein: $!\n"; $_="1";
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
            $best_blast{$query_id}{homo_gene} = $protein{$subject_id}{gene};
            $best_blast{$query_id}{homo_chromosome} = $protein{$subject_id}{chromosome};
            $best_blast{$query_id}{homo_name} = $protein{$subject_id}{name};
            $better ++;
        }
        elsif ($e_value == $best_blast{$query_id}{e_value} && $score == $best_blast{$query_id}{score} && $identity_perc == $best_blast{$query_id}{identity}) {
            if ($protein{$subject_id}{gene} eq $best_blast{$query_id}{homo_gene}) {
                if ($best_blast{$query_id}{same_gene} eq "none") {
                    $best_blast{$query_id}{same_gene} = $subject_id;
                }
                else {
                    $best_blast{$query_id}{same_gene} = $best_blast{$query_id}{same_gene} . "," . $subject_id;   
                }
                $multihit ++;
            }
            elsif ($protein{$subject_id}{chromosome} =~ /^(HG|HS|GL)/ && $best_blast{$query_id}{homo_chromosome} !~ /^(HG|HS|GL)/) {
                $badder ++;
            }
            elsif ($protein{$subject_id}{chromosome} !~ /^(HG|HS|GL)/ && $best_blast{$query_id}{homo_chromosome} =~ /^(HG|HS|GL)/) {
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
                $best_blast{$query_id}{homo_gene} = $protein{$subject_id}{gene};
                $best_blast{$query_id}{homo_chromosome} = $protein{$subject_id}{chromosome};
                $best_blast{$query_id}{homo_name} = $protein{$subject_id}{name};
                $best_blast{$query_id}{same_gene} = "none";
                $better ++;
            }
            elsif ($protein{$subject_id}{name} eq $best_blast{$query_id}{homo_name} && (($protein{$subject_id}{chromosome} =~ /^(HG|HS|GL)/ && $best_blast{$query_id}{homo_chromosome} =~ /^(HG|HS|GL)/) || ($protein{$subject_id}{chromosome} eq $best_blast{$query_id}{homo_chromosome}))) {
                if ($best_blast{$query_id}{same_gene} eq "none") {
                    $best_blast{$query_id}{same_gene} = $subject_id;
                }
                else {
                    $best_blast{$query_id}{same_gene} = $best_blast{$query_id}{same_gene} . "," . $subject_id;  
                }
                $multihit ++;
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
        $best_blast{$query_id}{homo_gene} = $protein{$subject_id}{gene};
        $best_blast{$query_id}{homo_chromosome} = $protein{$subject_id}{chromosome};
        $best_blast{$query_id}{homo_name} = $protein{$subject_id}{name};
        $best_blast{$query_id}{same_gene} = "none";
        $novel_number ++;
    }
    last if eof BLAST;
}
close (BLAST);      # Close the input file

# Remove duplicate values from the list of best protein hits
my %uniq_best_ids;
my $ambiguous_novel = 0;
foreach my $keys (keys %best_blast) {
    my $value = $best_blast{$keys}{subject};
    if ($value ne "ambiguous") {
        my @value;
        unless ($best_blast{$keys}{same_gene} eq "none") {
            @value = (split(/,/, $best_blast{$keys}{same_gene}));
        }
        push @value, $value;
        foreach my $value (@value) {
            if (exists $uniq_best_ids{$value}) {
                $uniq_best_ids{$value} += 1;
            }
            else {
                $uniq_best_ids{$value} = 1;
            }
        }
    }
}

# Print out the results in the output files
my $outfile_count = 0;
my $protein_count = keys %uniq_best_ids;
foreach my $keys (keys %best_blast) {
    if ($best_blast{$keys}{subject} ne "ambiguous") {
        print OUTFILE "$keys,$best_blast{$keys}{subject},$best_blast{$keys}{identity},$best_blast{$keys}{e_value},$best_blast{$keys}{score},$best_blast{$keys}{same_gene}\n";
        $outfile_count ++;
    }
    else {
        $ambiguous_novel ++;
    }
}
foreach my $keys (keys %uniq_best_ids) {
    print OUTPROTEIN "$keys\n";
}
close (OUTFILE);    # Close output file containing  the best match blast data
close (OUTPROTEIN); # Close the output file containing the match protein id

print STDERR "This perl script processed $ref_line number of blast matches for a total of $novel_number number of novel transcripts, including $ambiguous_novel number of best ids ties; while there were $protein_count number of unique best match IDs!\n\n";
if ($ref_line == ($novel_number+$better+$ambiguous+$badder+$multihit) && $outfile_count == ($novel_number-$ambiguous_novel)) {
    print STDERR "Process completed successfully!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__