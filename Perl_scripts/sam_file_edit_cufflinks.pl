#!/usr/bin/perl -w

# Script used to add optional field to aligned reads in sam file from STAR run

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $sam;    # Input file containing aligned reads from STAR run

# Define the parameter in order to submit input files to this script
&GetOptions (
    'sam=s' => \$sam,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input sam file
unless ($sam) {
    die "Please specify the sam file via -sam parameter!\n";
}
open (SAM, "<$sam") || die "Cannot open $sam: $!\n"; $_="1";
print STDERR "Processing the SAM file: $sam\n\n";

# Get the name of the samples
my $name = `basename $sam`;
chomp ($name);
$name =~ s/(.*)\.sam/$1/;
$name =~ s/(.*)_Aligned.out/$1/;

# Open the ouput files which will be in SAM format
my $outfile = "${name}_novel.sam";
if (-e $outfile) {
    die "This file: $outfile already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTFILE, ">$outfile") || die "Cannot open $outfile: $!\n"; $_="1";
}

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
    elsif ($line =~ /^(FC.1\w{3}ACXX\:\d\:\d+\:\d+\:\d+\#\w{6})/) {
        my $read_id = $1;
        chomp (my $line2 = <SAM>);
        $count_line_input ++;
        my $flag_mate1;
        my $flag_mate2;
        if ($line2 =~ /^$read_id/) {
            if ($line =~ /^$read_id\t(\d+)\t/) {
                $flag_mate1 = $1;
            }
            if ($line2 =~ /^$read_id\t(\d+)\t/) {
                $flag_mate2 = $1;
            }
            if ($flag_mate1 == 163 && $flag_mate2 == 83) {
                print OUTFILE "$line\tRG:Z:$name\tXS:A:-\n";
                print OUTFILE "$line2\tRG:Z:$name\tXS:A:-\n";
            }
            elsif ($flag_mate1 == 99 && $flag_mate2 == 147) {
                print OUTFILE "$line\tRG:Z:$name\tXS:A:+\n";
                print OUTFILE "$line2\tRG:Z:$name\tXS:A:+\n";
            }
            elsif ($flag_mate1 == 419 && $flag_mate2 == 339) {
                print OUTFILE "$line\tRG:Z:$name\tXS:A:-\n";
                print OUTFILE "$line2\tRG:Z:$name\tXS:A:-\n";
            }
            elsif ($flag_mate1 == 355 && $flag_mate2 == 403) {
                print OUTFILE "$line\tRG:Z:$name\tXS:A:+\n";
                print OUTFILE "$line2\tRG:Z:$name\tXS:A:+\n";
            }
            else {
                die "The flag of read mate1: $flag_mate1 and/or the flag of read mate2: $flag_mate2, do not match at line $count_line_input, please check input file!\n";
            }
            $count_line_output += 2;
            $read_count ++;
        }
        else {
            die "The reads are not paired at line $count_line_input, please check your SAM file!\n";
        }
    }
    else {
        die "The line does not correspond to any fields of SAM file format at line number $count_line_input, please check your SAM file!\n"
    }
    last if eof SAM
}
close (SAM);    # Close SAM file
close (OUTFILE);   # Close output file

print STDERR "This perl script is designed to work specifically with STAR aligned reads SAM files, with the Bos taurus UMD3.1.71 genome and with a featureCounts reads assignment file; script will need to be modified for any variations!\n\n";
print STDERR "There were $count_line_input lines in the input SAM file $sam, and $count_line_output lines in the ouput SAM file $outfile of which $read_count read fragments!\n\n";
if ($count_line_output == ($read_count * 2 + $header_count)) {
    print STDERR "Process completed correctly!\n\n";
}
else {
    print STDERR "Process completed uncorrectly!\n\n";
}

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__