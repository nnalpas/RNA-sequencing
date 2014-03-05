#!/usr/bin/perl -w

# Script used to modify the gtf file generated after Cufflinks

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $newgtf;    # Input file containing annotation from Cufflinks
my $type;  # Input parameter to define the type of feature within the gtf file
my $refgtf;  # Input file containing reference annotation in gtf format

# Define the parameter in order to submit input files to this script
&GetOptions (
    'newgtf=s' => \$newgtf,
    'type=s' => \$type,
    'refgtf=s' => \$refgtf,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input gtf file
unless ($newgtf) {
    die "Please specify the gtf file to edit via -newgtf parameter!\n";
}
open (NEWGTF, "<$newgtf") || die "Cannot open $newgtf: $!\n"; $_="1";

# Read the type of gtf file provided by user
unless ($type) {
    die "Please specify the type of gtf file that you have provided for editing via -type parameter (authorised value are novel or nat)!\n";
}

# Define output file
my $name = `basename $newgtf`;
chomp ($name);
$name =~ s/(.*)\.gtf/$1/;

# Open the ouput file which will be in gtf format
my $outfile = "${name}_edited.gtf";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

# Define variables required for reading files
my %reference_id;
my $ref_line = 0;
my $gtf_line = 0;
my %gene;
my %transcript;
my %exon;
my $num_gene = 0;
my $num_transcript = 0;
my $num_exon = 0;

# According to the type of feature in the new gtf file open or not the reference gtf file
if ($type eq "nat") {
    unless ($refgtf) {
        die "Please specify the reference gtf file via -refgtf parameter!\n";
    }
    open (REFGTF, "<$refgtf") || die "Cannot open $refgtf: $!\n"; $_="1";
    print STDERR "Working with natural antisense transcript annotation, reference annotation is required and has been opened!\n";
    while (1) {
        chomp (my $line = <REFGTF>);
        $ref_line ++;
        if ($line =~ /gene_id \"(.*?)\"; transcript_id \"(.*?)\";/) {
            $reference_id{$2} = $1;
        }
        else {
            die "The reference gtf file does not have a gene ID and/or transcript ID at line $ref_line!\n";
        }
        last if eof REFGTF;
    }
    close (REFGTF);
}
elsif ($type eq "novel") {
    print STDERR "Working with novel transcript annotation, reference annotation not required!\n";
}
else {
    die "The value of type parameter is non authorised, note that authorised value are novel or nat!\n";
}

# Read in the gtf file to edit
if ($type eq "novel") {
    NOVEL_EDIT: while (1) {
        my $new_gene_id;
        my $new_transcript_id;
        my $new_exon_id;
        chomp (my $line = <NEWGTF>);
        $gtf_line ++;
        my ($chromosome, $rna_type, $feature, $start, $end, $score, $strand, $frame, $attribute) = (split(/\t/, $line));
        my ($gene_id, $transcript_id, $exon_num, $oid, $class_code, $tss_id) = (split(/;/, $attribute));
        if (exists $gene{$gene_id}) {
            $new_gene_id = $gene{$gene_id};
        }
        else {
            $num_gene ++;
            $new_gene_id = "gene_id \"NOVBTAG" . sprintf("%011d", $num_gene) . "\";";
            $gene{$gene_id} = $new_gene_id;
        }
        if (exists $transcript{$transcript_id}) {
            $new_transcript_id = $transcript{$transcript_id};
        }
        else {
            $num_transcript ++;
            $new_transcript_id = "transcript_id \"NOVBTAT" . sprintf("%011d", $num_transcript) . "\";";
            $transcript{$transcript_id} = $new_transcript_id;
        }
        my $exon_id = $transcript_id;
        $exon_id =~ s/ transcript_id \"(.*)\"/$1 . \. . $exon_num/;
        if (exists $exon{$exon_id}) {
            $new_exon_id = $exon{$exon_id};
        }
        else {
            $num_exon ++;
            $new_exon_id = "exon_id \"NOVBTAE" . sprintf("%011d", $num_exon) . "\";";
            $exon{$exon_id} = $new_exon_id;
        }
        my $gene_biotype = "gene_biotype \"novel_transcript\";";
        $exon_num =~ s/^ (exon_number \"\d*\")/$1;/;
        $class_code =~ s/^ (class_code \".\")/$1;/;
        print OUTFILE "$chromosome\t$rna_type\t$feature\t$start\t$end\t$score\t$strand\t$frame\t$new_gene_id $new_transcript_id $exon_num $gene_biotype $new_exon_id $class_code\n";
        last if eof NEWGTF;
    }
}
elsif ($type eq "nat") {
    NAT_EDIT: while (1) {
        my $new_gene_id;
        my $new_transcript_id;
        my $new_exon_id;
        chomp (my $line = <NEWGTF>);
        $gtf_line ++;
        my ($chromosome, $rna_type, $feature, $start, $end, $score, $strand, $frame, $attribute) = (split(/\t/, $line));
        my ($gene_id, $transcript_id, $exon_num, $gene_name, $oid, $nearest_ref, $class_code, $tss_id) = (split(/;/, $attribute));
        $nearest_ref =~ s/ nearest_ref \"(.*)\"/$1/;
        if (exists $reference_id{$nearest_ref}) {
            $new_gene_id = "gene_id \"NAT-" . $reference_id{$nearest_ref} . "\";";
        }
        else {
            die "Check your gtf file to edit at line $gtf_line, it appears not to have a correct nearest_ref!\n";
        }
        if (exists $transcript{$transcript_id}) {
            $new_transcript_id = $transcript{$transcript_id};
        }
        else {
            $num_transcript ++;
            $new_transcript_id = "transcript_id \"NATBTAT" . sprintf("%011d", $num_transcript) . "\";";
            $transcript{$transcript_id} = $new_transcript_id;
        }
        my $exon_id = $transcript_id;
        $exon_id =~ s/ transcript_id \"(.*)\"/$1 . \. . $exon_num/;
        if (exists $exon{$exon_id}) {
            $new_exon_id = $exon{$exon_id};
        }
        else {
            $num_exon ++;
            $new_exon_id = "exon_id \"NATBTAE" . sprintf("%011d", $num_exon) . "\";";
            $exon{$exon_id} = $new_exon_id;
        }
        my $gene_biotype = "gene_biotype \"NAT\";";
        $exon_num =~ s/^ (exon_number \"\d*\")/$1;/;
        $class_code =~ s/^ (class_code \".\")/$1;/;
        print OUTFILE "$chromosome\t$rna_type\t$feature\t$start\t$end\t$score\t$strand\t$frame\t$new_gene_id $new_transcript_id $exon_num $gene_biotype $new_exon_id $class_code\n";
        last if eof NEWGTF;
    }
}
close (NEWGTF);
close (OUTFILE);   # Close output file

print STDERR "This perl script is designed to edit gtf annotation file generated by Cufflinks for novel or NAT transcripts!\n\n";

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__