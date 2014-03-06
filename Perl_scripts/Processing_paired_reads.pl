#!/usr/bin/perl -w

# Copyright (C) November/2013 by Nicolas C Nalpas

# Script used on sequence fastq file for deconvolution of reads from each individual samples based on unique index barcode, and for adapter filtering, and for overall read quality and for 3'-end read trimming

# Define all modules to be used in script
use lib '/home/people/09155562/myperl_lib/x86_64-linux-thread-multi/'; # location of String::Approx module on Phaeton server
use lib 'C:/Perl/lib/String-Approx-3.26';   # location of String::Approx on windows laptop
use strict;
use warnings;
use Getopt::Long;
use String::Approx qw(amatch);
use String::Approx qw(aindex);
use IO::File;
use File::Basename;

# Define the different input/output files such as fastq sequence file, indices list file, each index output fastq file
my $seqfile1;    # Input fastq sequence paired1 file
my $seqfile2;    # Input fastq sequence paired2 file
my $indices;    # Input file listing all indices used in the pool RNA-seq library
my $trim_mate1;      # Input number of bases to be trimmed at 3'-end of reads mate1
my $trim_mate2;      # Input number of bases to be trimmed at 3'-end of reads mate2
my $adapter_filt;    # Variable to determine if adapter filtering need to be performed
my $amatch_param;       # Variable to determine the parameters to provide amatch for adpater filtering
my $quality_filt;     # Variable to determine if overall quality read filtering need to be performed
my $quality_bases;      # Variable to determine the percentage of bases required to be above the quality threshold
my $illumina_version;   # Variable to determine the illumina version used for phred score encoding
my $fh;         # Filehandle output fastq file for each index

# Define the parameter in order to submit input files to this script
&GetOptions (
    'seqfile1=s' => \$seqfile1,
    'seqfile2=s' => \$seqfile2,
    'indices=s' => \$indices,
    'trim_mate1=s' => \$trim_mate1,
    'trim_mate2=s' => \$trim_mate2,
    'adapter_filt=s' => \$adapter_filt,
    'amatch_param=s' => \$amatch_param,
    'quality_filt=i' => \$quality_filt,
    'quality_bases=i' => \$quality_bases,
    'illumina_version=s' => \$illumina_version,
);

my $start_date = localtime;
print STDERR ("Start = $start_date\n\n");

# Open the input file containing fastq sequence paired1 if it has been submitted by user, otherwise ask for submission of this file
unless ($seqfile1) {
    die "Please specify the sequence paired1 file (either in gz compressed or uncompressed format) via -seqfile1 parameter!\n";
}
if ($seqfile1 =~ /\.gz$/) {
    open(SEQFILE1, "gunzip -c $seqfile1 |") || die "Cannot umcompress $seqfile1 in pipe: $!\n"; $_="1";
}
else {
    open (SEQFILE1, "<$seqfile1") || die "Cannot open $seqfile1: $!\n"; $_="1";
}
print STDERR "Processing $seqfile1\n\n";

# Open the input file containing fastq sequence paired1 if it has been submitted by user, otherwise ask for submission of this file
unless ($seqfile2) {
    die "Please specify the sequence paired2 file (either in gz compressed or uncompressed format) via -seqfile2 parameter!\n";
}
if ($seqfile2 =~ /\.gz$/) {
    open(SEQFILE2, "gunzip -c $seqfile2 |") || die "Cannot umcompress $seqfile2 in pipe: $!\n"; $_="1";
}
else {
    open (SEQFILE2, "<$seqfile2") || die "Cannot open $seqfile2: $!\n"; $_="1";
}
print STDERR "Processing $seqfile2\n\n";

# Open the input file containing index list used for library preparation if it has been submitted by user, otherwise ask for submission of this file
unless ($indices) {
    die "Please specify the indices file via -indices parameter!\n";
}
if ($indices eq "no") {
    print STDERR "As you requested, deconvolution was not perfomed!\n\n";
}
else {
    open (INDICES, "<$indices") || die "Cannot open $indices: $!\n"; $_="1";
}

# Process the input number of reads to trim at 3'-end of reads if it has been submitted by user, otherwise ask for submission of this number
my $trim5_mate1;    # Scalar containing number of bases to trim off at 5'-end of read mate1
my $keep_mate1;     # Scalar containing number of bases to keep in read mate1
my $trim3_mate1;    # Scalar containing number of bases to trim off at 3'-end of read mate1
my $trim5_mate2;    # Scalar containing number of bases to trim off at 5'-end of read mate2
my $keep_mate2;     # Scalar containing number of bases to keep in read mate2
my $trim3_mate2;    # Scalar containing number of bases to trim off at 3'-end of read mate2
my $length_mate1;
my $length_mate2;

unless ($trim_mate1) {
    die "Please specify the number of bases to trim at 5'-end and 3'-end of reads mate1 if trimming requested via the -trim_mate1 parameter, if trimming not required state zero where appropriate; note that parameters need to look like this e.g.: 5end_10_keep_70_3end_10!\n";
}
if ($trim_mate1 =~ /5end_(\d*)_keep_(\d*)_3end_(\d*)/) {     # Split the trimming parameters for read mate1 into several scalar
    $trim5_mate1 = $1;
    $keep_mate1 = $2;
    $trim3_mate1 = $3;
    $length_mate1 = ($trim5_mate1+$keep_mate1+$trim3_mate1);
    print STDERR "Parameters provided for trimming process of mate1 are correctly formatted!\n\n";
    if ($trim5_mate1 == 0 && $trim3_mate1 == 0) {
        $trim_mate1 = "no";
        print STDERR "As you requested, trimming of read mate1 was not performed!\n\n";
    }
}
else {
    die "Parameters provided for trimming process of mate1 are uncorrectly formatted, note that parameters need to look like this e.g.: 5end_10_keep_70_3end_10!\n";
}

unless ($trim_mate2) {
    die "Please specify the number of bases to trim at 5'-end and 3'-end of reads mate2 if trimming requested via the -trim_mate2 parameter, if trimming not required state zero where appropriate; note that parameters need to look like this e.g.: 5end_10_keep_70_3end_10!\n";
}
if ($trim_mate2 =~ /5end_(\d*)_keep_(\d*)_3end_(\d*)/) {     # Split the trimming parameters for read mate2 into several scalar
    $trim5_mate2 = $1;
    $keep_mate2 = $2;
    $trim3_mate2 = $3;
    $length_mate2 = ($trim5_mate2+$keep_mate2+$trim3_mate2);
    print STDERR "Parameters provided for trimming process of mate2 are correctly formatted!\n\n";
    if ($trim5_mate2 == 0 && $trim3_mate2 == 0) {
        $trim_mate2 = "no";
        print STDERR "As you requested, trimming of read mate2 was not performed!\n\n";
    }
}
else {
    die "Parameters provided for trimming process of mate2 are uncorrectly formatted, note that parameters need to look like this e.g.: 5end_10_keep_70_3end_10!\n";
}

if ($keep_mate1 == 0 || $keep_mate2 == 0) {
    die "The script was terminated because your trimming parameters are not allowed, you cannot keep zero bases for mate1 and/or mate2, please change your parameters!\n";
}

# Process the input file containing the adapter to filter out of reads if it has been submitted by user, otherwise ask for submission of this file
unless ($adapter_filt) {
    die "Please specify if you request adapter filtering to be performed by providing a file containing adapter sequence or stating no via -adapter_filt parameter!\n";
}
if ($adapter_filt eq "no") {
    print STDERR "As you requested, adapter filtering was not perfomed!\n\n";
}
else {
    open (ADAPTER_FILT, "<$adapter_filt") || die "Cannot open $adapter_filt: $!\n"; $_="1";
}

unless ($amatch_param) {
    die "Please specify the parameters for matching adapter sequence during adapter filtering process otherwise state no via -amatch_param parameter; note that parameters need to look like this e.g.: i_3_S3_D2_I2 means case-insensitive match, up to 3 base changes, of which a max 3 substitutions and a max 2 insertions or deletions!\n";
}
if ($amatch_param eq "no") {
    print STDERR "As you requested, adapter filtering was not perfomed!\n\n";
    if ($adapter_filt ne "no") {
        die "Please specify if adapter filtering is required, if yes both parameters -adapter_filt and -amatch_param must have values different from no, else state no for both parameters!\n";
    }
}
elsif ($amatch_param =~ /i_\d_S\d_D\d_I\d/) {
    $amatch_param =~ s/_/ /g;
    $amatch_param =~ s/(.*)/[\"$1\"]/;
    print STDERR "Parameters provided for adapter sequence matching during filtering process are correctly formatted!\n\n";
}
else {
    die "Parameters provided for adapter sequence matching during filtering process are uncorrectly formatted; note that parameters need to look like this e.g.: i_3_S3_D2_I2 means case-insensitive match, up to 3 base changes, of which a max 3 substitutions and a max 2 insertions or deletions!\n"
}

# Process the input number of overall quality threshold for reads if it has been submitted by user, otherwise ask for submission of this number
unless ($quality_filt) {
    print STDERR "Please specify the overall quality of reads via the -quality_filt parameter, default is currently zero!\n";
    $quality_filt = 0;
    print STDERR "As you requested, the overall read quality filtering was not performed!\n\n";
}
chomp($quality_filt);

unless ($quality_bases) {
    print STDERR "Please specify the maximum percentage of bases per read mate (note: read mate are processed as different reads) allowed to be below quality threshold via the -quality_bases parameter, default is currently zero!\n";
    $quality_bases = 0;
    print STDERR "As you requested, the overall read quality filtering was not performed!\n\n";
    if ($quality_filt != 0) {
        die "Please specify if overall read quality filtering is required, if yes both parameters -quality_filt and -quality_bases must have values different from 0, else state 0 for both parameters or do not even write these two parameters!\n";
    }
}
chomp($quality_bases);
my $quality_bases_mate1 = ($quality_bases*$keep_mate1/100);
my $quality_bases_mate2 = ($quality_bases*$keep_mate2/100);

unless ($illumina_version) {
    die "Please specify the illumina version of reads phred score encoding via the -illumina_version parameter!\n";
}
chomp($illumina_version);

# Define the phred score encoding based on illumina versions
my $quality_threshold;
if ($illumina_version >= 1.0 && $illumina_version < 1.8) {
    $quality_threshold = ($quality_filt+64);
}
elsif ($illumina_version >= 1.8) {
    $quality_threshold = ($quality_filt+33);
}
else {
    die "This illumina version: $illumina_version is not allowed; please check your illumina version for phred score encoding!\n";
}

# Check that at least one process will be performed otherwise terminate script
if ($indices eq "no" && $trim_mate1 eq "no" && $trim_mate2 eq "no" && $adapter_filt eq "no" && $quality_filt == 0) {
    die "You have not requested any processes to be performed, please check values of parameters provided!\n";
}

# Define the seqfile basename according to data generated by BGI or MSU or Conway and open the report file, which will contain all stats on individual samples and tags occurence
my ($file_basename,$file_path) = fileparse($seqfile1);
my $data_centre;    # Scalar containing the name of the data centre which generated the data
if ($file_basename =~ /(Raw_.*)_1\..*$/) {
    $file_basename = $1;
    print STDERR "Data seems to have been generated by BGI, if this is not the case please seek advice about this perl script!\n\n";
    $data_centre = "BGI";
    # Open the report file, which will contain all stats on individual samples and tags occurence
    my $report = "${file_basename}_1_2_deconv.report";
    open (REPORT, ">$report") || die "Cannot open $report: $!\n"; $_="1";
}
elsif ($file_basename =~ /(.*?)_R1_(\d{3})\..*$/) {
    $file_basename = $1;
    my $lane = $2;
    print STDERR "Data seems to have been generated by MSU, if this is not the case please seek advice about this perl script!\n\n";
    $data_centre = "MSU";
    # Open the report file, which will contain all stats on individual samples and tags occurence
    my $report = "${file_basename}_1_2_${lane}_deconv.report";
    open (REPORT, ">$report") || die "Cannot open $report: $!\n"; $_="1";
    if ($indices ne "no") {
        die "Your data seems to have been generated by MSU, deconvolution is not allowed on such data type yet; please seek advice about this perl script!";
    }
}
elsif ($file_basename =~ /(.*_sequence)\.txt\.gz/) {
    $file_basename = $1;
    print STDERR "Data seems to have been generated by Conway, if this is not the case please seek advice about this perl script!\n\n";
    $data_centre = "Conway";
    # Open the report file, which will contain all stats on individual samples and tags occurence
    my $report = "${file_basename}_deconv.report";
    open (REPORT, ">$report") || die "Cannot open $report: $!\n"; $_="1";
}
else {
    die "Your data seems to have been generated by an unknown sequencing centre; please seek advice about this perl script!";
}

# Open two output files, which will contain all excluded paired1 reads and paired2 reads respectively based on their unassigned barcode
if ($indices ne "no") {
    open (EXCLUDED1, ">${file_basename}_1.excluded") || die "Cannot open ${file_basename}1.excluded: $!\n"; $_="1";
    open (EXCLUDED2, ">${file_basename}_2.excluded") || die "Cannot open ${file_basename}2.excluded: $!\n"; $_="1";
}

# Define variables required for reading the input indices file
my $index;  # Scalar containing first index sequence
my $sample; # Scalar containing first sample name
my @indices;    # Array containing all indices sequence
my %indices;    # Hash containing all indices sequence associated with their sample name

# Read in the input file containing the indices used for library preparation and open filehandle output files for each index
unless ($indices eq "no") {
    while (<INDICES>){
        ($index, $sample) = (split(/\t/));  # Read index sequence and sample name in a tab separated file
        chomp($index);
        chomp($sample);
        push @indices, "$index\n";  # Put all indices sequence into array
        next if $index =~ /^\s*$/;
        $indices{$index} = $sample; # Create hash of indices sequence with corresponding samples name
        my $outfile1 = "${sample}_pe1.fastq";
        my $outfile2 = "${sample}_pe2.fastq";
        open_fh("${index}1", $outfile1);  # Create and open an outfile for each paired1 reads individual samples
        open_fh("${index}2", $outfile2);  # Create and open an outfile for each paired2 reads individual samples
        if ($adapter_filt ne "no" || $quality_filt != 0) {
            my $excluded1 = "${sample}_pe1.excluded";
            my $excluded2 = "${sample}_pe2.excluded";
            open_fh("${index}1_excluded", $excluded1);  # Create and open an outfile for each paired1 reads individual samples
            open_fh("${index}2_excluded", $excluded2);  # Create and open an outfile for each paired2 reads individual samples
        }
        last if eof (INDICES);
    }
    close (INDICES);
}
else {
    $index = "$file_basename";
    $indices{$index} = $index; # Create hash of indices sequence with corresponding samples name
    my $outfile1 = "${file_basename}_pe1_filt.fastq";
    my $outfile2 = "${file_basename}_pe2_filt.fastq";
    open_fh("${index}1", $outfile1);  # Create and open an outfile for each paired1 reads individual samples
    open_fh("${index}2", $outfile2);  # Create and open an outfile for each paired2 reads individual samples
    if ($adapter_filt ne "no" || $quality_filt != 0) {
        my $excluded1 = "${file_basename}_pe1.excluded";
        my $excluded2 = "${file_basename}_pe2.excluded";
        open_fh("${index}1_excluded", $excluded1);  # Create and open an outfile for each paired1 reads individual samples
        open_fh("${index}2_excluded", $excluded2);  # Create and open an outfile for each paired2 reads individual samples
    }
}

# Define variables required for reading the input adapter sequence file
my %adap_filter; # Hash containing all adapter sequence, adapter name and its associated read mate
my $adap_seq;    # Scalar containing adapter sequence
my $adap_id;     # Scalar containing adapter name
my $mate;   # Scalar containing read mate number
my $mate1;
my $mate2;
my $adap_seq1;
my $adap_id1;
my $adap_seq2;
my $adap_id2;
my $i = 0;  # Scalar to check correct formatting of adapter file

unless ($adapter_filt eq "no") {
    while (<ADAPTER_FILT>){
        ($adap_seq, $adap_id, $mate) = (split(/\t/));    # Read adapter sequence, adapter name and mate number separated by tab
        $i++;
        chomp($adap_seq);
        chomp($adap_id);
        chomp($mate);
        unless (exists $adap_filter{$adap_id}) {  # Create hash of adapter name with corresponding adapter sequence
            $adap_filter{$adap_id}{sequence} = $adap_seq;
            $adap_filter{$adap_id}{mate} = $mate;
        }
        if ($mate =~ /Read1/i) {
            $adap_id1 = $adap_id;
            $adap_seq1 = $adap_seq;
            $mate1 = $mate;
        }
        elsif ($mate =~ /Read2/i) {
            $adap_id2 = $adap_id;
            $adap_seq2 = $adap_seq;
            $mate2 = $mate;
        }
        elsif ($i > 2) {
            die "The adapter filtering file is not formatted correctly, this file should contain between one and two adapter sequences, adapter ids and read mates to filter (tab separated)!"
        }
        else {
            die "The adapter filtering file is not formatted correctly, this file should contain between one and two adapter sequences, adapter ids and read mates to filter (tab separated)!"
        }
        last if eof (ADAPTER_FILT);
    }
    close (ADAPTER_FILT);
}

# Define variables required for reading the input fastq sequence file
my $id_line1_1;        # Scalar containing first id line from the pe1 fastq sequence file
my $sequence_line_1;  # Scalar containing sequence line from the pe1 fastq sequence file
my $id_line2_1;       # Scalar containing second id line from the pe1 fastq sequence file
my $quality_line_1;   # Scalar containing quality line from the pe1 fastq sequence file
my $id_line1_2;        # Scalar containing first id line from the pe2 fastq sequence file
my $sequence_line_2;  # Scalar containing sequence line from the pe2 fastq sequence file
my $id_line2_2;       # Scalar containing second id line from the pe2 fastq sequence file
my $quality_line_2;   # Scalar containing quality line from the pe2 fastq sequence file
my $untrim_sequence_line_1; # Scalar containing untrimmed read sequence from the pe1 fastq sequence file
my $untrim_quality_line_1; # Scalar containing untrimmed read quality from the pe1 fastq sequence file
my $untrim_sequence_line_2; # Scalar containing untrimmed read sequence from the pe2 fastq sequence file
my $untrim_quality_line_2; # Scalar containing untrimmed read quality from the pe2 fastq sequence file
my $tag;            # Scalar containing the first 6 bases from the sequence which corresponds to the tag (for de-convoluting of reads per sample)
my %matches;        # Hash containing all tags matching exactly to an index and their number of occurence
my %mismatch_allowed;   # Hash containing all mismatch tags which could be attributed to a unique index with their number of occurence in the fastq file and the index sequence they have been attributed to
my %ignore;         # Hash containing all tags which are to be ignored because they are one mismatch away from several indices with their number of occurence in the fastq file and the indices sequence they are ambiguous to
my %unassigned;     # Hash containing all tags not matching to any indices (based on current criterion) with their number of occurence in the fastq file
my %deconvoluted_reads_counts;  # Hash containing number of reads attributed to individual samples also showing their index sequence
foreach my $keys (keys %indices) {
    $deconvoluted_reads_counts{$keys}{sample} = $indices{$keys};
}
my %adapter_per_sample;  # Hash containing number of adapter per individual samples
foreach my $keys (keys %indices) {
    $adapter_per_sample{$keys}{sample} = $indices{$keys};
    $adapter_per_sample{$keys}{adap_id1} = $adap_id1;
    $adapter_per_sample{$keys}{adap_id2} = $adap_id2;
}
my %quality_bases;       # Hash containing the number of mate reads below the quality threshold
foreach my $keys (keys %indices) {
    $quality_bases{$keys}{sample} = $indices{$keys};
}
my $input_count = 0;    # Scalar to count total number of read in input fastq file

# Define the matches hash once and for all
foreach my $keys (keys %indices) {
    $matches{$keys}{occurence}=0;
    $matches{$keys}{index}="Perfect match";
}

# Read in the fastq sequence file
PROCESSING: while(1) {
    chomp($id_line1_1 = <SEQFILE1>);    # Read 4 lines one by one at a time from the pe1 fastq file, with first line being the pe1 id line
    chomp($sequence_line_1 = <SEQFILE1>);  # Second line being the pe1 sequence line
    chomp($id_line2_1 = <SEQFILE1>);   # Third line being the pe1 id line again
    chomp($quality_line_1 = <SEQFILE1>);   # Fourth line being the pe1 phred score quality line
    chomp($id_line1_2 = <SEQFILE2>);    # Read 4 lines one by one at a time from the pe2 fastq file, with first line being the pe2 id line
    chomp($sequence_line_2 = <SEQFILE2>);  # Second line being the pe2 sequence line
    chomp($id_line2_2 = <SEQFILE2>);   # Third line being the pe2 id line again
    chomp($quality_line_2 = <SEQFILE2>);    # Fourth line being the pe2 phred score quality line
    $input_count ++;    # Count number of reads in input fastq file
    my $keep = "yes";   # Scalar to determine if a deconvoluted read shall be kept post-filtering
    PAIRED_CHECK: {
        my $id1 = $id_line1_1;
        my $id2 = $id_line1_2;
        if ($data_centre eq "BGI") {
            $id1 =~ s/\/1$//;
            $id2 =~ s/\/2$//;
        }
        elsif ($data_centre eq "MSU") {
            $id1 =~ s/(\d )1(\:N\:0\:\w{6})$/"$1" . "$2"/;
            $id2 =~ s/(\d )2(\:N\:0\:\w{6})$/"$1" . "$2"/;
        }
        if ($id1 ne $id2) { # Test if read id of first and second fastq file are identical otherwise stop the script
            die "The read id of first fastq file do not match the read id of second fastq file for read number $input_count, note that data must be paired and in the same order between first and second fastq file\n";
        }
    }
    DECONV: if ($indices ne "no" && $id_line1_1 =~ /\#(\w{6})/){   # Remove the tag (barcode) present in the id line at the index location
        $tag = $1;  # Put the removed tag sequence into a scalar
        $index = "";    # Start fresh with scalar index as null
        if (exists ($matches{$tag})){   # If the tag sequence matches exactly one of the indices used for library preparation then count it as one more occurence and define scalar index as the tag sequence
            $index = $tag;
            $matches{$tag}{occurence}++;
        }
        elsif (exists ($mismatch_allowed{$tag})) {
            $index = $mismatch_allowed{$tag}{index};
            $mismatch_allowed{$tag}{occurence}++;
        }
        elsif (exists ($ignore{$tag}) ){    # If the tag sequence matches exactly one of the tag to ignore then count it as one more occurence and define scalar index as null
            $index="";
            $ignore{$tag}{occurence}++;
        }
        elsif (exists ($unassigned{$tag})) {    # If the tag sequence matches exactly one of the tag to leave unassigned then count it as one more occurence and define scalar index as null
            $index="";
            $unassigned{$tag}++;
        }
        else {  # When tag does not match any of the indices or any of the tag to ignore it will undergo test to determine if tag is one mismatch away from zero, one or several indices
            my $barcode;
            my %mismatches;   # Hash containing current tag which is tested for mismatch against all indices, with the number of indices it is one mismatch away from and these indices sequence
            my @test_indices = @indices;    # Define into an array the list of indices against which the tag will be tested
            foreach my $test_indices (@test_indices) {    # Perform the test against one index at a time by removing first index on top of array and defining this index as a scalar
                chomp ($test_indices);
                if (amatch($test_indices, ["S1 i D0"], $tag)) { # If the tag is one mismatch away from the currently tested index then count one more potentially matching index and store matching index sequence
                    $barcode = $test_indices;
                    if (!exists $mismatches{$tag}) {
                        $mismatches{$tag}{occurence} = 1;
                        $mismatches{$tag}{index} = $barcode;
                    }        
                    else {
                        $mismatches{$tag}{occurence}++;
                        $mismatches{$tag}{index} = join (", ", $mismatches{$tag}{index}, $barcode);
                    }
                }
            }
            if (!exists $mismatches{$tag}) {    # If the tag was not one mismatch away from any indices count zero potentially matching index and attribute this tag to unassigned hash and count one occurence
                $mismatches{$tag}{occurence} = 0;
                $index = "";
                $unassigned{$tag} = 1;
            }
            elsif ($mismatches{$tag}{occurence} == 1) {    # If the tag was found to be one mismatch away from a unique index then count one more occurence and store its matching index sequence and define scalar index as the tag sequence
                $index = $barcode;
                $mismatch_allowed{$tag}{occurence} = 1;
                $mismatch_allowed{$tag}{index} = $index;
            }
            elsif ($mismatches{$tag}{occurence} > 1) {  # If the tag was found to be one mismatch away from several indices then count one more occurence and store all potentially matching indices sequence and update the list of tag to ignore with this one and define scalar index as null
                $index = "";
                $ignore{$tag}{occurence} = 1;
                $ignore{$tag}{index} = $mismatches{$tag}{index};
            }
        }
    }
    TRIMMING: {
        $untrim_sequence_line_1 = $sequence_line_1;
        $untrim_quality_line_1 = $quality_line_1;
        $untrim_sequence_line_2 = $sequence_line_2;
        $untrim_quality_line_2 = $quality_line_2;
        if ($length_mate1 == length $sequence_line_1 && $length_mate2 == length $sequence_line_2) {
            if ($trim_mate1 ne "no") {
                $sequence_line_1 =~ s/^.{$trim5_mate1}(.{$keep_mate1}).{$trim3_mate1}$/$1/;
                $quality_line_1 =~ s/^.{$trim5_mate1}(.{$keep_mate1}).{$trim3_mate1}$/$1/;
            }
            if ($trim_mate2 ne "no") {
                $sequence_line_2 =~ s/^.{$trim5_mate2}(.{$keep_mate2}).{$trim3_mate2}$/$1/;
                $quality_line_2 =~ s/^.{$trim5_mate2}(.{$keep_mate2}).{$trim3_mate2}$/$1/;
            }
        }
        else {
            die "The script was terminated because the read mate1 length: ", length $sequence_line_1, " and/or read mate2 length: ", length $sequence_line_2, " do not match the trimming parameters which you provided!\n";
        }
    }
    ADAP_FILT: if ($adapter_filt ne "no" && $index ne "" && $keep ne "no_bc_quality") {
        if ($mate1 ne "" && amatch($adap_seq1, ${amatch_param}, $sequence_line_1)){  # If the adapter sequence has fuzzy matching as provided with parameters, then read will be output into excluded pe1
            $keep = "no_bc_adapter";
            $adapter_per_sample{$index}{occurence1}++;
        }
        if ($mate2 ne "" && amatch($adap_seq2, ${amatch_param}, $sequence_line_2)){  # If the adapter sequence has fuzzy matching as provided with parameters, then read will be output into excluded pe2
            $keep = "no_bc_adapter";
            $adapter_per_sample{$index}{occurence2}++;
        }
    }
    QUAL_FILT: if ($quality_filt != 0 && $index ne "" && $keep ne "no_bc_adapter") {
        my @quality_read1 = reverse(unpack ("C*", $quality_line_1));
        my @quality_read2 = reverse(unpack ("C*", $quality_line_2));
        my $quality_bases1 = 0;
        my $quality_bases2 = 0;
        QUAL_READ1: foreach my $base (@quality_read1) {
            if ($base < $quality_threshold) {
                $quality_bases1 ++;
                if ($quality_bases1 > $quality_bases_mate1) {
                    $quality_bases{$index}{occurence1} ++;
                    $keep = "no_bc_quality";
                    last QUAL_READ1;
                }
            }
        }
        QUAL_READ2: foreach my $base (@quality_read2) {
            if ($base < $quality_threshold) {
                $quality_bases2 ++;
                if ($quality_bases2 > $quality_bases_mate2) {
                    $quality_bases{$index}{occurence2} ++;
                    $keep = "no_bc_quality";
                    last QUAL_READ2;
                }
            }
        }
    }
    PRINT: {
        if ($index ne "") {  # If scalar index is not null then add the tag sequence to the header pe1 and header pe2, and print id lines, sequence line and quality line to the pe1 and pe2 outfiles of the corresponding index sequence sample
            if ($indices ne "no") {
                $id_line1_1 =~ s/\#\w{8}/"\#" . $tag/e;
                $id_line1_2 =~ s/\#\w{8}/"\#" . $tag/e;
            }
            $deconvoluted_reads_counts{$index}{number_reads}++; # Count number of reads added to each individual samples prior-filtering
            if ($keep eq "yes") {
                $deconvoluted_reads_counts{$index}{number_reads_kept}++; # Count number of reads added to each individual samples post-filtering
                print_to_fh ("${index}1", join ("\n", $id_line1_1, $sequence_line_1,  $id_line2_1 , $quality_line_1) . "\n");
                print_to_fh ("${index}2", join ("\n", $id_line1_2, $sequence_line_2,  $id_line2_2 , $quality_line_2) . "\n");
            }
            elsif ($keep eq "no_bc_adapter") {
                $deconvoluted_reads_counts{$index}{number_reads_adapter}++; # Count number of filtered out reads                    
                print_to_fh ("${index}1_excluded", join ("\n", $id_line1_1, $untrim_sequence_line_1,  $id_line2_1 , $untrim_quality_line_1) . "\n");
                print_to_fh ("${index}2_excluded", join ("\n", $id_line1_2, $untrim_sequence_line_2,  $id_line2_2 , $untrim_quality_line_2) . "\n");
            }
            elsif ($keep eq "no_bc_quality") {
                $deconvoluted_reads_counts{$index}{number_reads_quality}++; # Count number of filtered out reads
                print_to_fh ("${index}1_excluded", join ("\n", $id_line1_1, $untrim_sequence_line_1,  $id_line2_1 , $untrim_quality_line_1) . "\n");
                print_to_fh ("${index}2_excluded", join ("\n", $id_line1_2, $untrim_sequence_line_2,  $id_line2_2 , $untrim_quality_line_2) . "\n");
            }
        }
        elsif ($indices ne "no" && $index eq "") {  # If scalar index is null then add id lines, sequence line and quality line into the excluded pe1 and pe2 read file
            print EXCLUDED1 join ("\n", $id_line1_1, $untrim_sequence_line_1,  $id_line2_1 , $untrim_quality_line_1) . "\n";
            print EXCLUDED2 join ("\n", $id_line1_2, $untrim_sequence_line_2,  $id_line2_2 , $untrim_quality_line_2) . "\n";
            $index = "Excluded";
            $deconvoluted_reads_counts{$index}{sample} = "Ignore";
            $deconvoluted_reads_counts{$index}{number_reads}++; # Count number of reads added to excluded file
        }
    }
    last if eof (SEQFILE1);  # If the sequence fastq file was fully read, then exit reading the fastq file
}
close (SEQFILE1);    # Close sequence pe1 fastq file
close (SEQFILE2);    # Close sequence pe2 fastq file

if ($indices ne "no") {
    close (EXCLUDED1);   # Close excluded pe1 reads fastq file
    close (EXCLUDED2);   # Close excluded pe2 reads fastq file
}
    
foreach $index (keys %indices){
    close_fh ("${index}1");
    close_fh ("${index}2");
    if ($adapter_filt ne "no" || $quality_filt != 0) {
        close_fh ("${index}1_excluded");
        close_fh ("${index}2_excluded");
    }
}

# Print into report the details of the input fastq file
print REPORT "\n\nNumber of reads in input pool fastq file: $input_count\n\n";

# Print into report the details of adapter filtering per individual samples
print REPORT "Number of reads per deconvoluted samples prior- and post-filtering:\n";
print REPORT "Index sequence\tSample id\tNo. reads pre-filtering\tNo. reads post-filtering\tPercentage reads post-filtering\tNo. reads discarded due to adapter sequence\tPercentage reads with adapter sequence\tNo. reads discarded due to poor overall quality\tPercentage reads with poor overall quality\n";
foreach my $key (sort keys %deconvoluted_reads_counts) {
    if ($key eq "Excluded") {
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
    print REPORT "$key\t$deconvoluted_reads_counts{$key}{sample}\t$deconvoluted_reads_counts{$key}{number_reads}\t$deconvoluted_reads_counts{$key}{number_reads_kept}\t$number_reads_kept_percent\t$deconvoluted_reads_counts{$key}{number_reads_adapter}\t$number_reads_adapter_percent\t$deconvoluted_reads_counts{$key}{number_reads_quality}\t$number_reads_quality_percent\n";
}
print REPORT "\n";

# Print into report the details of reads trimming
print REPORT "There were $trim5_mate1 and $trim5_mate2 bases trimmed at 5'-end of read mate1 and mate2 respectively; with $trim3_mate1 and $trim3_mate2 bases trimmed at 3'-end of read mate1 and mate2 respectively; with $keep_mate1 and $keep_mate2 bases kept for read mate1 and mate2 respectively.\n";
print REPORT "\n";

# Print into report the details of adapter filtered out
if ($adapter_filt ne "no" && $amatch_param ne "no") {
    print REPORT "Adapter filtering was performed with fuzzy matching determine by these parameters: $amatch_param; with i meaning case-insensitive match, up to ... base changes, of which a max S... substitutions and a max D... deletions or I... insertions.\n";
    print REPORT "Adapter id\tAdapter sequence\tAdapter present in read\n";
    foreach my $key (sort keys %adap_filter) {
        print REPORT "$key\t$adap_filter{$key}{sequence}\t$adap_filter{$key}{mate}\n";
    }
    print REPORT "\n";
    print REPORT "Number of adapter filtered out per individual samples:\n";
    print REPORT "Index sequence\tSample id\tAdapter1 id\tNo. reads containing Adapter1\tAdapter2 id\tNo. reads containing Adapter2\n";
    foreach my $key (sort keys %adapter_per_sample) {
        print REPORT "$key\t$adapter_per_sample{$key}{sample}\t$adapter_per_sample{$key}{adap_id1}\t$adapter_per_sample{$key}{occurence1}\t$adapter_per_sample{$key}{adap_id2}\t$adapter_per_sample{$key}{occurence2}\n";
    }
    print REPORT "\n";
}
else {
    print REPORT "As you requested, adapter filtering was not perfomed!\n\n";
}

# Print into report the details of overall quality filtering
if ($quality_filt != 0) {
    print REPORT "The overall quality filtering using Illumina version $illumina_version was set up such that for any read mates with more than $quality_bases% of bases having a phred score below $quality_filt, the whole paired-end read was filtered out.\n";
    print REPORT "\n";
    print REPORT "Number of filtered out read mates below quality threshold per individual samples:\n";
    print REPORT "Index sequence\tSample id\tNo. reads mate1 below quality threshold\tNo. reads mate2 below quality threshold\n";
    foreach my $key (sort keys %quality_bases) {
        print REPORT "$key\t$quality_bases{$key}{sample}\t$quality_bases{$key}{occurence1}\t$quality_bases{$key}{occurence2}\n";
    }
    print REPORT "\n";
}
else {
    print REPORT "As you requested, the overall read quality filtering was not performed!\n\n";
}

# Print into report the number of tags encountered in the fastq file
if ($indices ne "no") {
    my $processed_tag = ((keys %matches)+(keys %mismatch_allowed)+(keys %ignore)+(keys %unassigned)); # Scalar containing the total number of matching tags to accept or to ignore
    print REPORT "Number of different tags encountered in the fastq file is: $processed_tag\n\n";
    print REPORT "Tags which match exactly to a unique index:\n";   # Print into report the occurence of tags matching unique index
    print REPORT "Matching tag\tNo. occurence\tMatched index\n";
    foreach my $key (keys %matches) {
        print REPORT "$key\t$matches{$key}{occurence}\t$matches{$key}{index}\n";
    }
    print REPORT "\n";
    print REPORT "Tags which are one mismatch away from a unique index:\n"; # Print into report the occurence of tags one mismatch away from unique index
    print REPORT "Matching tag\tNo. occurence\tMatched index\n";
    foreach my $key (keys %mismatch_allowed) {
        print REPORT "$key\t$mismatch_allowed{$key}{occurence}\t$mismatch_allowed{$key}{index}\n";
    }
    print REPORT "\n";
    print REPORT "Tags which were ambiguous (and so ignored) because they are one mismatch away from several indices:\n";   # Print into report the occurence of tags ambiguous to several indices
    print REPORT "Ambiguous tag\tNo. occurence\tAmbiguous index\n";
    foreach my $key (keys %ignore) {
        print REPORT "$key\t$ignore{$key}{occurence}\t$ignore{$key}{index}\n";
    }
    print REPORT "\n";
    print REPORT "Tags which were unassigned because they are more than one mismatch away from any indices:\n"; # Print into report the occurence of tags unassigned to any indices
    print REPORT "Unassigned tag\tNo. occurence\n";
    foreach my $key (keys %unassigned) {
        print REPORT "$key\t$unassigned{$key}\n";
    }
    print REPORT "\n";
}
else {
    print REPORT "As you requested, deconvolution was not perfomed!\n\n";
}

close (REPORT); # Close report file

# Print note to user on the use of the perl script
print STDERR "NOTE TO USER:\nThis perl script is optimised to work specifically with paired-end reads with index barcode of 6 bases length positioned at the tag position of the read header, also allowing adapter sequence filtering on paired-end reads, overall phred score quality filtering for paired-end reads and trimming of paired-end reads!\n\n";
print STDERR "Please note that overall read quality filtering is performed after adapter filtering, in consequence reads excluded based on adapter filtering will not be quality filtered!\n\n";
print STDERR "Deconvolution, adapter filtering, overall quality filtering and read trimming processes completed\n\n";
my $finish_date = localtime;
print STDERR ("Finish = $finish_date\n\n");

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