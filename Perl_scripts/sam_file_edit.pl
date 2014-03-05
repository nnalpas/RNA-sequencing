#!/usr/bin/perl -w

# Script used to add optional field to aligned reads in sam file from STAR run

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $sam;    # Input file containing aligned reads from STAR run
my $sense;  # Input file containing reads with assigned sense gene from featureCounts
my $antisense;  # Input file containing reads with assigned antisense gene from featureCounts
my $novel;  # Input file containing reads with assigned novel gene from featureCounts

# Define the parameter in order to submit input files to this script
&GetOptions (
    'sam=s' => \$sam,
    'sense=s' => \$sense,
    'antisense=s' => \$antisense,
    'novel=s' => \$novel,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input sam file
unless ($sam) {
    die "Please specify the sam file via -sam parameter!\n";
}
open (SAM, "<$sam") || die "Cannot open $sam: $!\n"; $_="1";

# Open the input assigned sense read file
unless ($sense) {
    die "Please specify the file containing reads with assigned sense gene via -sense parameter!\n";
}
open (SENSE, "<$sense") || die "Cannot open $sense: $!\n"; $_="1";

# Open the input assigned antisense read file
unless ($antisense) {
    die "Please specify the file containing reads with assigned antisense gene via -antisense parameter!\n";
}
open (ANTISENSE, "<$antisense") || die "Cannot open $antisense: $!\n"; $_="1";

# Open the input assigned novel gene read file
unless ($novel) {
    die "Please specify the file containing reads with assigned novel gene via -novel parameter!\n";
}
open (NOVEL, "<$novel") || die "Cannot open $novel: $!\n"; $_="1";

# Get the name of the samples
my $name = `basename $sam`;
chomp ($name);
$name =~ s/(.*)(_Aligned.out)?.sam/$1/;

# Open the ouput file which will be in SAM format
my $outfile = "${name}_assigned.sam";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

# Define variables required for reading the featureCounts reads file
my %feature;    # Hash containing the read assignment for sense, antisense or novel featureCounts results
my $count_line_feature = 0;
my $assigned_read = 0;
my $unassigned_read = 0;
my $ambiguous_read = 0;
my $multihit_read = 0;

# Read the featureCounts for sense gene reads file and create a hash of assigned or unassigned reads
while(1) {
    chomp (my $feature_line = <SENSE>);
    $count_line_feature ++;
    (my $read = $feature_line) =~ s/(.*?)\t.*/$1/;
    unless (exists $feature{$read}) {
        if ($feature_line =~ /\tAssigned\t(ENSBTAG\d{11})(.*)?/) {
            my $gene = $1;
            $feature{$read} = $gene;
            $assigned_read ++;
        }
        elsif ($feature_line =~ /\tUnassigned_No_Features/) {
            $feature{$read} = "Unassigned";
            $unassigned_read ++;
        }
        elsif ($feature_line =~ /\tUnassigned_Ambiguous/) {
            $feature{$read} = "Ambiguous";
            $ambiguous_read ++;
        }
        elsif ($feature_line =~ /\tUnassigned_Multimapping/) {
            $feature{$read} = "Multihit";
            $multihit_read ++;
        }
        else {
            die "The featureCounts read file have unexpected value at $count_line_feature, please check your file $sense!\n";
        }
    }
    else {
        if ($feature_line =~ /\tUnassigned_Multimapping/ && $feature{$read} eq "Multihit") {
            $multihit_read ++;
        }
        else {
            die "The featureCounts read file have reads at $count_line_feature previously encountered which are not multihits, this is not allowed, please check your file $sense!\n";
        }
    }
    
    last if eof (SENSE);
}
close (SENSE);  # Close featureCounts sense gene read file

# Read the featureCounts for antisense gene reads file and add info to hash of assigned or unassigned reads
$count_line_feature = 0;
while(1) {
    chomp(my $feature_line = <ANTISENSE>);
    $count_line_feature ++;
    (my $read = $feature_line) =~ s/(.*?)\t.*/$1/;
    if (exists $feature{$read}) {
        if ($feature_line =~ /\tAssigned\t(NATENSBTAG\d{11})(.*)?/ && $feature{$read} eq "Unassigned") {
            my $gene = $1;
            $feature{$read} = $gene;
            $assigned_read ++;
            $unassigned_read --;
        }
        elsif ($feature_line =~ /\tUnassigned_No_Features/ && $feature{$read} eq "Unassigned") {
            unless (eof ANTISENSE) {
                redo;
            }
        }
        elsif ($feature_line =~ /\tUnassigned_Ambiguous/ && $feature{$read} eq "Unassigned") {
            $feature{$read} = "Ambiguous";
            $ambiguous_read ++;
            $unassigned_read --;
        }
        else {
            die "The featureCounts read files $sense and $antisense both have assigned gene for the read $read, please check your files!\n";
        }
    }
    else {
        die "The read $read in file $antisense does not exists in the $sense featureCounts file at line $count_line_feature, please check your file!\n";
    }
    last if eof (ANTISENSE);
}
close (ANTISENSE);  # Close featureCounts antisense gene read file

# Read the featureCounts for novel gene reads file and add info to hash of assigned or unassigned reads
$count_line_feature = 0;
while(1) {
    chomp(my $feature_line = <NOVEL>);
    $count_line_feature ++;
    (my $read = $feature_line) =~ s/(.*?)\t.*/$1/;
    if (exists $feature{$read}) {
        if ($feature_line =~ /\tAssigned\t(NOVBTAG\d{11})(.*)?/ && $feature{$read} eq "Unassigned") {
            my $gene = $1;
            $feature{$read} = $gene;
            $assigned_read ++;
            $unassigned_read --;
        }
        elsif ($feature_line =~ /\tUnassigned_No_Features/ && $feature{$read} eq "Unassigned") {
            unless (eof NOVEL) {
                redo;
            }
        }
        elsif ($feature_line =~ /\tUnassigned_Ambiguous/ && $feature{$read} eq "Unassigned") {
            $feature{$read} = "Ambiguous";
            $ambiguous_read ++;
            $unassigned_read --;
        }
        else {
            die "The featureCounts read files $sense and/or $antisense and/or $novel have assigned gene for the read $read, please check your files!\n";
        }
    }
    else {
        die "The read $read in file $novel does not exists in the $sense featureCounts file at line $count_line_feature, please check your file!\n";
    }
    last if eof (NOVEL);
}
close (NOVEL);  # Close featureCounts novel gene read file

# Define variables required for reading the input file
my $header_line = '@HD';
my $ref_seq = '@SQ';
my $read_group = '@RG';
my $program = '@PG';
my $comment = '@CO';
my $read_count = 0;
my $count_line_input = 0;
my $count_line_output = 0;
my $strand;
my $header_count = 0;

# Read in the sam file and the sense and antisense read files
SAM_PROCESS: while (1){
    chomp (my $line = <SAM>);
    $count_line_input ++;
    if ($line =~ /^$header_line/) {
        print OUTFILE "$line\n";
        $count_line_output ++;
        $header_count ++;
    }
    elsif ($line =~ /^$ref_seq/) {
        if ($line =~ /\tAS:UMD3.1.71\tSP:Bos taurus/) {
            print OUTFILE "$line\n";
            $count_line_output ++;
            $header_count ++;
        }
        else {
            print OUTFILE "$line\tAS:UMD3.1.71\tSP:Bos taurus\n";
            $count_line_output ++;
            $header_count ++;
        }
    }
    elsif ($line =~ /^$read_group/) {
        if ($line =~ /\@RG\tID:$name\tCN:BGI\tDS:Bovine alveolar macrophages challenge project\tPL:Illumina HiSeq 2000\tSM:$name/) {
            print OUTFILE "$line\n";
            $count_line_output ++;
            $header_count ++;
        }
        else {
            print OUTFILE ('@RG', "\tID:$name\tCN:BGI\tDS:Bovine alveolar macrophages challenge project\tPL:Illumina HiSeq 2000\tSM:$name\n");
            $count_line_output ++;
            $header_count ++;
        }
    }
    elsif ($line =~ /^$program/) {
        print OUTFILE "$line\n";
        $count_line_output++;
        $header_count ++;
    }
    elsif ($line =~ /^$comment/) {
        print OUTFILE "$line\n";
        $count_line_output++;
        $header_count ++;
    }
    elsif ($line =~ /^(FC.{5}ACXX\:\d\:\d+\:\d+\:\d+\#\w{6})/) {
        my $read_id = $1;
        (my $flag_mate = $line) =~ s/^$read_id\t(\d+)\t.*/$1/;
        if ($flag_mate == 163 || $flag_mate == 83 || $flag_mate == 419 || $flag_mate == 339) {
            $strand = "-";
        }
        elsif ($flag_mate == 99 || $flag_mate == 147 || $flag_mate == 355 || $flag_mate == 403) {
            $strand = "+";
        }
        else {
            die "The flag of read: $flag_mate at line $count_line_input is not an authorised flag value, please check input file!\n";
        }
        if (exists $feature{$read_id}) {
            print OUTFILE "$line\tXS:A:$strand\tXG:Z:$feature{$read_id}\n";
            $count_line_output ++;
            $read_count ++;
        }
        else{
            die "Read ID $read_id at line $count_line_input from $sam file is missing from the featureCounts file $sense, please check your files!\n";
        }
    }
    else {
        die "The line does not correspond to any fields of SAM file format at line number $count_line_input for file $sam, please check your SAM file!\n"
    }
    last if eof SAM;
}
close (SAM);    # Close SAM file
close (SENSE);  # Close assigned sense read file
close (ANTISENSE);  # Close assigned antisense read file
close (NOVEL);  # Close assigned novel read file
close (OUTFILE);   # Close output file

print STDERR "This perl script is designed to work specifically with STAR aligned reads SAM files, with the Bos taurus UMD3.1.71 genome and with a featureCounts reads assignment file; script will need to be modified for any variations!\n\n";
print STDERR "There were $count_line_input lines in the input SAM file $sam, and $count_line_output lines in the ouput SAM file $outfile of which $read_count read fragments!\n\n";
if ($count_line_output == ($read_count + $header_count) && ($read_count/2) == ($assigned_read+$unassigned_read+$ambiguous_read+$multihit_read)) {
    print STDERR "Process completed correctly!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__