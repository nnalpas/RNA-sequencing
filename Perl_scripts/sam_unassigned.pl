#!/usr/bin/perl -w

# Script used to select all unassigned reads from SAM files

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $sam;    # Input file containing aligned reads from STAR run
my $featurecount;  # Input file containing reads with assigned gene from featureCounts
my $custom;   # Custom name for unassigned reads

# Define the parameter in order to submit input files to this script
&GetOptions (
    'sam=s' => \$sam,
    'featurecount=s' => \$featurecount,
    'custom=s' => \$custom,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input sam file
unless ($sam) {
    die "Please specify the sam file via -sam parameter!\n";
}
open (SAM, "<$sam") || die "Cannot open $sam: $!\n"; $_="1";

# Open the input assigned sense read file
unless ($featurecount) {
    die "Please specify the file containing reads with assigned sense, antisense or novel gene via -featurecount parameter!\n";
}
open (FEATURECOUNT, "<$featurecount") || die "Cannot open $featurecount: $!\n"; $_="1";

# Check that a custom name for the reads output file was provided by user
unless ($custom) {
    die "Please specify the custom name to give the unassigned reads output file via -custom parameter!\n";
}

# Get the name of the samples
my $name = `basename $sam`;
chomp ($name);
$name =~ s/(.*).sam/$1/;

# Open the ouput file which will be in SAM format
my $outfile = "${name}_${custom}_unassigned.sam";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

# Define variables required for reading the featureCounts reads file
my %feature;
my $count_line_feature = 0;
my $unassigned_count = 0;

# Read the featureCounts reads file and create a hash of unassigned reads
while(1) {
    chomp (my $feature_line = <FEATURECOUNT>);
    $count_line_feature ++;
    if ($feature_line =~ /(.*)\tUnassigned_No_Features/) {
        my $unassigned = $1;
        $feature{$unassigned} = "Unassigned";
        $unassigned_count ++;
    }
    elsif ($feature_line !~ /[Assigned]|[Unassigned_Ambiguous]|[Unassigned_Multimapping]/) {
        die "The featureCounts read file have unexpected value at $count_line_feature, please check your file!\n";
    }
    last if eof FEATURECOUNT;
}
close (FEATURECOUNT);  # Close assigned read file

# Define variables required for reading the input file
my $header_line = '@HD';
my $ref_seq = '@SQ';
my $read_group = '@RG';
my $program = '@PG';
my $comment = '@CO';
my $read_count = 0;
my $count_line_input = 0;
my $count_line_output = 0;
my $header_count = 0;

# Read in the sam file and select only unassigned reads
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
        chomp (my $line2 = <SAM>);
        $count_line_input ++;
        if ($line2 =~ /^$read_id/) {
            if (exists $feature{$read_id}) {
                print OUTFILE "$line\n";
                print OUTFILE "$line2\n";
                $count_line_output += 2;
                $read_count ++;
            }
        }
        else {
            die "The reads are not paired at line $count_line_input, please check your SAM file!\n";
        }
    }
    else {
        die "The line does not correspond to any fields of SAM file format at line number $count_line_input for file $sam, please check your SAM file!\n"
    }
    last if eof SAM;
}
close (SAM);    # Close SAM file
close (OUTFILE);   # Close output file

print STDERR "This perl script is designed to work specifically with STAR aligned reads SAM files, with the Bos taurus UMD3.1.71 genome and with a featureCounts reads assignment file; script will need to be modified for any variations!\n\n";
print STDERR "There were $count_line_input lines in the input SAM file $sam, and $count_line_output lines in the ouput SAM file $outfile of which $read_count read fragments, corresponding to $unassigned_count unassigned reads!\n\n";
if ($count_line_output == ($read_count * 2 + $header_count) && $read_count == $unassigned_count) {
    print STDERR "Process completed correctly!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__