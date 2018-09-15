#!/usr/bin/perl -w
# the first line of perl code has to be as above
#
# Author: Alejandro Schaffer
#
# Code to compare annotations  of Norovirus entries from dnaorg_annotate and blastx
# Usage: compare_predictions.pl \
#        --input_nucleotide_table <input nucleotide prediction table> \
#        --input_blastx_directory <directory in which to find blastx summaries> \
#        --input_accessions <list of accessions for which there are nucleotide predictions> \
#        --input_protein_names <names of proteins> \
#        --alignment_tolerance <number of nucleotides>
#        --indel_tolerance <number of nucleotides> 
#        --output_table <combined output table> \
#        --output_directory <output directory for single accession summaries>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


# input/output file names
my $input_predict_file;        #input file that summarizes predictions
my $input_names_file;          #input file with protein names
my $input_accessions_file;     #input file that lists accessions to check
my $output_table;              #output file that summarizes predictions
my $one_blastx_summary_file;   #name of one file for a blastx summary
my $one_validation_file;       #name of one file for a validation report on one accession

# input/output directories
my $input_blastx_directory;    #directory in which blastx results are stored
my $output_directory;          #directory in which individual reports are stored  

# variables used in processing input
my $nextline;                  #one line of the input summary file
my @fields;                    #entries in one row
my @protein_names;             #array of protein names in the proteome of interest
my @protein_accessions;        #array of protein accessions in the proteome of interest, indices of protein_names and Protein_accessions should match up
my @accessions;                #array of accessions for which nucleotide predictions are being validated
my @accessions_no_version;     #array of accessions without the version numbers
my $num_proteins;              #number of proteins in this proteome
my $num_headers_input;         #number of column headers in input table
my $num_headers_output;        #number of column headers in output table
my @column_headers_input;      #array of column headers in the input predictions table
my @column_headers_output;     #array of column headers in the output table
my $num_accessions;            #number of accesions to consider
my $num_values_on_row;         #number of values on a single input row 
my $i;                         #loop index
my $j;                         #loop index
my $p;                         #loop index
my $input_header_offset;       #index for columns inside loop
my $output_header_offset;      #index for columns inside loop
my @all_input_data;            #all the values in the input file in a two-dimensional array, where the first dimension is by accession
my @all_output_data;           #all the values in the output table in a two-dimensional array, where the first dimension is by accession
my $nucleotide_protein_tolerable_discrepancy;  #upper bound on the discrepancy between the nucleotide alignment endpoint and protein alignment endpoint tolerable for foosh
my $indel_tolerance;            #how large can an in-frame indel be before the sequence becomes ineligible for foosh


#some constants
my $num_generic_input_headers = 3;    #number of headers in the input table that are not specific to a protein
my $num_specific_input_headers = 2;   #number of headers in the input table that are specific to each protein
my $num_generic_output_headers = 4;   #number of headers in the output table that are not specific to a protein
my $num_specific_output_headers = 12; #number of headers in the output table that are specific to each protein
my $dividing_line = "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\";    #used to make one output file easier to read
my $FOOSH_COLUMN_INDEX = 2;                 # the index of the output column that has the foosh decision
my $INDEL_COLUMN_INDEX = 3;                 # the index of the output column that has the longest indel
my $DEBUG = 1;

#for each possible protein the output table has seven columns and these are the header index offsets
my $HEADER_OFFSET_START_NUC  = 0;
my $HEADER_OFFSET_START_PROT = 1;
my $HEADER_OFFSET_STOP_NUC   = 2;
my $HEADER_OFFSET_STOP_PROT  = 3;
my $HEADER_OFFSET_INS        = 4;
my $HEADER_OFFSET_DEL        = 5;
my $HEADER_OFFSET_MAX_INDEL  = 6;
my $HEADER_OFFSET_TRC        = 7;
my $HEADER_OFFSET_FRAME      = 8;   #what is the frame of the alignment?
my $HEADER_OFFSET_NTERM      = 9;   #does the alignment extend to the protein N-terminus? 
my $HEADER_OFFSET_CTERM      = 10;   #does the alignment extend to the protein C-terminus?
my $HEADER_OFFSET_STOP_CODON = 11;  #if the difference between nucleotide end and protein end is exactly 3, then is that because the three nucleotides encode a stop codon?



# variables for dealing with options (epn-options.pm)
my %opt_HH           = ();
my @opt_order_A      = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
$opt_group_desc_H{"1"} = "basic options";
#       option                         type      default group   requires incompat preamble-outfile                                                                help-outfile
opt_Add("-h",                          "boolean",0,          0,    undef, undef,   undef,                                          "display this help",                         \%opt_HH, \@opt_order_A);
opt_Add("--input_nucleotide_table",             "string", undef,      1,    undef, undef,   "input file of nucleotide annotations",   "File name <s> of nucleotide predictions",        \%opt_HH, \@opt_order_A);
opt_Add("--input_blastx_directory",             "string", undef,      1,    undef, undef,   "input directory with blastx summaries",    "Directory name <s> of blastx summaries",        \%opt_HH, \@opt_order_A);
opt_Add("--input_accessions",             "string", undef,      1,    undef, undef,   "input list of accessions to validate",          "File name <s> of accesions to validate",        \%opt_HH, \@opt_order_A);
opt_Add("--input_protein_names",             "string", undef,      1,    undef, undef,   "input file of protein names",                "File name <s> of protein list",        \%opt_HH, \@opt_order_A);
opt_Add("--alignment_tolerance",            "integer", undef,      1,    undef, undef,   "tolerance for alignment discrepancy",                "Tolerance <n> of protein nucleotide discrepancy",        \%opt_HH, \@opt_order_A);
opt_Add("--indel_tolerance",            "integer", undef,      1,    undef, undef,   "tolerance for in-frame indels",                "Tolerance <n> for in-frame indels",        \%opt_HH, \@opt_order_A);
opt_Add("--output_table",            "string", undef,      1,    undef, undef,   "combined output table",                               "Name <s> of combined output table to create",         \%opt_HH, \@opt_order_A);
opt_Add("--output_directory",            "string", undef,      1,    undef, undef,   "output directory with one output file per accession",  "Name <s> of directory with one output file per accession",         \%opt_HH, \@opt_order_A);

my $synopsis = "compare_predictions.pl: compare nucleotide predictions and protein alignments \n";
my $usage    = "Usage:\n\n";
$usage      .= "\tcompare_predictions.pl \\\n";
$usage      .= "\t--input_nucleotide_table <input nucleotide prediction table>  \\\n";
$usage      .= "\t--input_blastx_directory <directory in which to find blastx summaries> \\\n";
$usage      .= "\t--input_accessions <list of accessions for which there are nucleotide predictions> \\\n";
$usage      .= "\t--input_protein_names <names of proteins> \\\n";
$usage      .= "\t--alignment_tolerance <protein nucleotide discrepancy tolerance> \\\n";
$usage      .= "\t--indel_tolerance <indel length tolerance> \\\n";
$usage      .= "\t--output_table <combined output table>  \\\n";
$usage      .= "\t--output_directory <output directory for single accession summaries> \n\n";
$usage      .= "\nFor example:\n";
$usage      .= "compare_predictions.pl --input_nucleotide_table nucleotide_annotations.txt ";
$usage      .= "--input_blastx_directory ../blastx ";
$usage      .= "--input_accessions accessions.txt ";
$usage      .= "--input_protein_names protein_names.txt ";
$usage      .= "--alignment_tolerance 5 ";
$usage      .= "--indel_tolerance 27 ";
$usage      .= "--output_table nucleotide_annotations_validation_table.txt ";
$usage      .= "--output_directory single_accession_validation_files ";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay =
    &GetOptions('h'                            => \$GetOptions_H{"-h"},
		'input_nucleotide_table=s'         => \$GetOptions_H{"--input_nucleotide_table"},
		'input_blastx_directory=s'     => \$GetOptions_H{"--input_blastx_directory"},				
		'input_accessions=s'      => \$GetOptions_H{"--input_accessions"},		
		'input_protein_namese=s'           => \$GetOptions_H{"--input_protein_names"},
		'alignment_tolerance=n'           => \$GetOptions_H{"--alignment_tolerance"},
		'indel_tolerance=n'           => \$GetOptions_H{"--indel_tolerance"},						
		'output_table=s'                => \$GetOptions_H{"--output_table"},
		'output_directory=s'           => \$GetOptions_H{"--output_directory"});

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status
}

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);


# retrieve options
$input_predict_file      = opt_Get("--input_nucleotide_table", \%opt_HH);
$input_blastx_directory  = opt_Get("--input_blastx_directory", \%opt_HH);
$input_accessions_file   = opt_Get("--input_accessions", \%opt_HH);
$input_names_file        = opt_Get("--input_protein_names", \%opt_HH);
$nucleotide_protein_tolerable_discrepancy = opt_Get("--alignment_tolerance", \%opt_HH);
$indel_tolerance = opt_Get("--indel_tolerance", \%opt_HH);
$output_table            = opt_Get("--output_table", \%opt_HH);
$output_directory        = opt_Get("--output_directory", \%opt_HH);

# die if any of the required options were not used
my $errmsg = undef;
if(! defined $input_predict_file)      { $errmsg .= "ERROR, --input_nucleotide_table option not used.\n"; }
if(! defined $input_blastx_directory)  { $errmsg .= "ERROR, --input_blastx_directory option not used.\n"; }
if(! defined $input_accessions_file)   { $errmsg .= "ERROR, --input_accessions_file option not used.\n"; }
if(! defined $input_names_file)        { $errmsg .= "ERROR, --input_names_file option not used.\n"; }
if(! defined $nucleotide_protein_tolerable_discrepancy)        { $errmsg .= "ERROR, --alignment_tolerance option not used.\n"; }
if(! defined $indel_tolerance)        { $errmsg .= "ERROR, --indel_tolerance option not used.\n"; }
if(! defined $output_table)            { $errmsg .= "ERROR, --output_file option not used.\n"; }
if(! defined $output_directory)        { $errmsg .= "ERROR, --output_directory option not used.\n"; }

if(defined $errmsg) {
    die $errmsg . "\n$usage\n";
}

# open output files
open(SUMMARY, "<", $input_predict_file)         or die "Cannot open 1 $input_predict_file for input\n";
open(PROTEINS, "<", $input_names_file)          or die "Cannot open 2 $input_names_file for input\n";
open(ACCESSIONS,  "<", $input_accessions_file)        or die "Cannot open 3 $input_accessions_file for output\n";
open(OUTPUT,  ">", $output_table)                or die "Cannot open 4 $output_table for output\n";



$num_proteins = read_protein_names();
close(PROTEINS);
($num_headers_input, $num_headers_output) = make_column_headers($num_proteins);
if ($DEBUG) {
    print "Numbers of headers are $num_headers_input $num_headers_output\n";
}

$num_accessions = 0;
while(defined($nextline = <ACCESSIONS>)) {
    chomp($nextline);
    $accessions[$num_accessions] = $nextline;
    if ($nextline =~ m/\./) {
	($accessions_no_version[$num_accessions]) = ($nextline =~ m/^(\S+)\./);
    }
    else {
	$accessions_no_version[$num_accessions] = $accessions[$num_accessions];
    }
    if ($DEBUG) {
	print "Accesion $num_accessions is $accessions[$num_accessions] and without the version is $accessions_no_version[$num_accessions]\n";
    }
    $num_accessions++;
}
close(ACCESSIONS);

#initialize output data
for ($i = 0; $i < $num_accessions; $i++) {
    for ($j = 0; $j < $num_headers_output; $j++) {    
	$all_output_data[$i][$j] = "NP";
    }
    $all_output_data[$i][$FOOSH_COLUMN_INDEX] = "Yes";    
    $all_output_data[$i][$HEADER_OFFSET_STOP_CODON] = "MOOT";
    $all_output_data[$i][$INDEL_COLUMN_INDEX] = 0;
}

$i = 0;
$nextline = <SUMMARY>;
while(defined($nextline = <SUMMARY>)) {
    chomp($nextline);
    @fields = split /\t/, $nextline;
    for ($j = 0; $j < $num_headers_input; $j++) {
	$all_input_data[$i][$j] = $fields[$j];
    }
    $all_output_data[$i][0] = $all_input_data[$i][1]; #copy accession
    $all_output_data[$i][1] = $all_input_data[$i][2]; #copy length
    for ($p = 0; $p < $num_proteins; $p++) {
	$input_header_offset = $num_generic_input_headers + ($p * $num_specific_input_headers);
	$output_header_offset = $num_generic_output_headers + ($p * $num_specific_output_headers);
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_START_NUC] = $all_input_data[$i][$input_header_offset];
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_STOP_NUC] = $all_input_data[$i][$input_header_offset+1];
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_START_PROT] = "NP";
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_STOP_PROT] = "NP";
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_INS] = "";
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_DEL] = "";
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_MAX_INDEL] = 0;	
	$all_output_data[$i][$output_header_offset + $HEADER_OFFSET_TRC] = "";				     
    }
    $i++;
}
if ($DEBUG) {
    print "Initial values of output table\n";
    for ($j = 0; $j < $num_headers_output; $j++) {
	print "$all_output_data[0][$j]\t"
    }
    print "\n";    
}

for ($i = 0; $i < $num_accessions; $i++) {
   process_blastx_summary_file($i, $input_blastx_directory);    
}    

output_decisions($num_headers_output, $num_accessions);
close(OUTPUT);

################################################
# List of subroutines:
#
# read_protein_names():          Read the file of protein names and store them in teh array protein_names
# make_column_headers();         Generate the values for the array column_headers
# output_decisions();           Output the large table of all predictions and the single files for each accession
#
################################################

# Subroutine: read_protein_names()
# Synopsis: Read the file of protein names and store them in teh array protein_names
#
# Args: none
#
# Returns: number of proteins

sub read_protein_names {

    my $sub_name = "read_protein_names()";
    my $nargs_exp = 0;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my $local_nextline;
    my $local_num_proteins;
    my $local_i;
    
    $local_nextline = <PROTEINS>;
    chomp($local_nextline);
    $local_num_proteins = $local_nextline;
    for($local_i = 0; $local_i < $local_num_proteins; $local_i++) {
	$local_nextline = <PROTEINS>;
	chomp($local_nextline);
	$protein_names[$local_i] = $local_nextline;
    }

    return $local_num_proteins;
}

# Subroutine: make_column_headers()
# Synopsis: Generate the values for the array column_headers
#
# Args: local_num_proteins  #number of distinct proteins in this genome
#
# Returns: number of headers

sub make_column_headers {

    my $sub_name = "make_column_headers()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_num_proteins) = @_;
    my $local_num_headers_input = 0;
    my $local_num_headers_output = 0;    
    my $local_i;

    $column_headers_input[$local_num_headers_input] = "index";
    $local_num_headers_input++;
    $column_headers_input[$local_num_headers_input] = "accession";    
    $column_headers_output[$local_num_headers_output] = "accession";
    $local_num_headers_input++;
    $local_num_headers_output++;    
    $column_headers_input[$local_num_headers_input] = "totlen";
    $column_headers_output[$local_num_headers_output] = "totlen";
    $local_num_headers_input++;
    $local_num_headers_output++;
    $column_headers_output[$local_num_headers_output] = "foosh?";
    $local_num_headers_output++;
    $column_headers_output[$local_num_headers_output] = "max_indel";
    $local_num_headers_output++;    
    for($local_i = 0; $local_i < $local_num_proteins; $local_i++) {
	$column_headers_input[$local_num_headers_input] = $protein_names[$local_i] . "_start_nuc";
	$local_num_headers_input++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_start_nuc";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_start_prot";
	$local_num_headers_output++;
	$column_headers_input[$local_num_headers_input] = $protein_names[$local_i] . "_stop_nuc";
	$local_num_headers_input++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_stop_nuc";
	$local_num_headers_output++;	
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_stop_prot";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_ins_prot";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_del_prot";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_max_indel";
	$local_num_headers_output++;				
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_trc_prot";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_frame";
	$local_num_headers_output++;
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_N_terminus";
	$local_num_headers_output++;				
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_C_terminus";
	$local_num_headers_output++;				
	$column_headers_output[$local_num_headers_output] = $protein_names[$local_i] . "_last_codon_stop";
	$local_num_headers_output++;					
    }
    return $local_num_headers_input, $local_num_headers_output;
}


# Subroutine: output_decisions()
# Synopsis: Output the large table of all predictions and the single files for each accession 
#
# Args: $num_column_headers, $num_accessions
#
# Returns: nothing

sub output_decisions {
    my $sub_name = "output_decisions()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_num_column_headers, $local_num_accessions) = @_;    
    my $local_i; #loop index
    my $local_j; #loop index
    my $local_one_accession; #a single accession to summarize
    my $local_summary_filename; #name of file with the summary for one accession
    my $local_input_header_offset;
    my $local_output_header_offset;
    my $local_first_stop; #nucleotide position for first stop codon
    my $local_last_stop; #nucleotide position for first stop codon
    my $local_nucleotide_start;
    my $local_nucleotide_stop;
    my $local_valid;  #is the nucleotide prediction validated
    

    for ($local_i = 0; $local_i < $local_num_accessions; $local_i++) {    
	$local_one_accession = $accessions_no_version[$local_i];
	$local_summary_filename = $output_directory . $local_one_accession . ".validation.txt";
	open(VALIDATION,  ">", $local_summary_filename)                or die "Cannot open 5 $local_summary_filename for output\n";
	for ($local_j = 0; $local_j < $num_proteins; $local_j++) {
	    if (0 != $local_j) {
		print VALIDATION "$dividing_line\n";
	    }
	    $local_input_header_offset = $num_generic_input_headers + ($local_j * $num_specific_input_headers);
	    $local_output_header_offset = $num_generic_output_headers + ($local_j * $num_specific_output_headers);
	    if (($all_input_data[$local_i][$local_input_header_offset] =~ /\d+/) || ($all_input_data[$local_i][$local_input_header_offset+1] =~ /\d+/)) { #is there a nucleotide prediction?
		$local_nucleotide_start = $all_input_data[$local_i][$local_input_header_offset];
		$local_nucleotide_stop = $all_input_data[$local_i][$local_input_header_offset+1];
		$local_valid = 1; #start with the presumption that prediction is valid, look to invalidate
		print  VALIDATION "$protein_names[$local_j]: testing nucleotide annotation with protein alignment\n";
		print VALIDATION "Is there a significant blastx alignment between the nucleotide accession $accessions_no_version[$local_i] and this protein $protein_accessions[$local_j]?\n";
		if (($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT] =~ /\d+/) && ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT] =~ /\d+/)) {
		    print VALIDATION "Answer: Yes\n";
		    if ($all_input_data[$local_i][$local_input_header_offset] =~ /\d+/) {
			print VALIDATION "Is the nucleotide 5\' endpoint $all_input_data[$local_i][$local_input_header_offset] at least as far as the 5' extent $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT] of the blastx alignment?\n";
			if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT] <  $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT]) {
			    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_NUC] <= $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT]) {
				print VALIDATION "Answer: Yes\n";
			    }
			    else {
				print VALIDATION "Answer: No\n";
				$all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";				
				$local_valid = 0;
			    }
			}
			else {
			    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_NUC] >= $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT]) {
				print VALIDATION "Answer: Yes\n";
			    }
			    else {
				print VALIDATION "Answer: No\n";
				$all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";				
				$local_valid = 0;
			    }
			}
			if (abs($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_NUC] - $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT]) >
			     $nucleotide_protein_tolerable_discrepancy) {  
			    $all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";
			}
		    }
		    if ($all_input_data[$local_i][$local_input_header_offset+1] =~ /\d+/) {
			print VALIDATION "Is the nucleotide 3\' endpoint $all_input_data[$local_i][$local_input_header_offset + 1] at least as far as the 3' extent $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT] of the blastx alignment?\n";
			if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_START_PROT] <  $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT]) {
			    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_NUC] >= $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT]) {
				print VALIDATION "Answer: Yes\n";
			    }
			    else {
				print VALIDATION "Answer: No\n";
				$all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";								
				$local_valid = 0;
			    }
			}
			else {
			    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_NUC] <= $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT]) {
				print VALIDATION "Answer: Yes\n";
			    }
			    else {
				print VALIDATION "Answer: No\n";
				$all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";								
				$local_valid = 0;
			    }
			}
			if (abs($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_NUC] - $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_PROT]) >
			     $nucleotide_protein_tolerable_discrepancy) {  
			    $all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";
			    $local_valid = 0;
			}			
		    }
                    print VALIDATION "What is the frame of the blastx alignment?\n";
		    print VALIDATION "Answer: $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_FRAME]\n"; 
		    print VALIDATION "Does the blastx alignment extend to the N-terminus of $protein_names[$local_j]?\n";
		    print VALIDATION "Answer: $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_NTERM]\n"; 
		    print VALIDATION "Does the blastx alignment extend to the C-terminus of $protein_names[$local_j]?\n";
		    print VALIDATION "Answer: $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_CTERM]\n";
		    if ((!("MOOT" eq $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_CODON])) && 
                        (!("NP" eq $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_CODON]))) { 
			print VALIDATION "The predicted nucleotide stop is exactly 3 positions away from the predicted protein stop in nucleotide positions. Is this because blastx omitted the stop codon?\n";
			print VALIDATION "Answer: $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_STOP_CODON]\n"; 
		    }
                    print VALIDATION "Does the blastx alignment have a predicted truncation before the predicted nucleotide end?\n";
		    if (($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_TRC] =~ /\d/) && ($local_nucleotide_start =~ /\d/) && ($local_nucleotide_stop =~ /\d/)) {
			($local_first_stop) = ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_TRC] =~ /^(\d+)/);
			($local_last_stop) = ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_TRC] =~ /(\d+)$/);
			if ((($local_nucleotide_start <= $local_first_stop) && ($local_first_stop <= $local_nucleotide_stop)) || 
			    (($local_nucleotide_start <= $local_last_stop) && ($local_last_stop <= $local_nucleotide_stop)) ||
			    (($local_nucleotide_stop <= $local_first_stop) && ($local_first_stop <= $local_nucleotide_start)) || 
			    (($local_nucleotide_stop <= $local_last_stop) && ($local_last_stop <= $local_nucleotide_start))) {
			    print VALIDATION "Answer: Yes, protein prediction truncations $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_TRC] include at least one between $local_nucleotide_start and $local_nucleotide_stop\n";
			    $local_valid = 0;
   			    $all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";								
			}
			else {
			    print VALIDATION "Answer: No\n";						
			}
		    }
		    else {
			print VALIDATION "Answer: No\n";						
		    }
		    print VALIDATION "Does the query have in-frame insertions w.r.t. the RefSeq protein?\n";
		    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_INS] =~ /\d/) {
			print VALIDATION "Answer: Yes, $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_INS]\n";
		    }
		    else {
			print VALIDATION "Answer: No\n";			
		    }
		    print VALIDATION "Does the query have in-frame deletions w.r.t. the RefSeq protein?\n";
		    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_DEL] =~ /\d/) {
			    print VALIDATION "Answer: Yes, $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_DEL]\n";
		    }
		    else {
			print VALIDATION "Answer: No\n";			
		    }
		    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_MAX_INDEL] > $all_output_data[$local_i][$INDEL_COLUMN_INDEX]) {
			$all_output_data[$local_i][$INDEL_COLUMN_INDEX] = $all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_MAX_INDEL];
		    }
		    if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_MAX_INDEL] > 0) {
			print VALIDATION "Is the length of the longest in-frame indel for this RefSeq protein > $indel_tolerance ?\n";
			if ($all_output_data[$local_i][$local_output_header_offset + $HEADER_OFFSET_MAX_INDEL] > $indel_tolerance) {
			    print VALIDATION "Answer: Yes\n";
			    $all_output_data[$local_i][$FOOSH_COLUMN_INDEX] = "No";
			}
			else {
			    print VALIDATION "Answer: No\n";
			}
		    }
		}
		else {
		    print VALIDATION "Answer: No\n";					    
		}
		print VALIDATION "Does the protein alignment validate the nucleotide annotation of $protein_names[$local_j]?\n";
		if ($local_valid) {
		    print VALIDATION "Answer: Yes\n";
		}
		else {
		    print VALIDATION "Answer: No\n";
		}
	    }
	    else {
		print  VALIDATION "$protein_names[$local_j]: no prediction\n";
	    }
	}
	print VALIDATION "What is the length of the longest relevant indel? $all_output_data[$local_i][$INDEL_COLUMN_INDEX]\n";
	print VALIDATION "Do the protein alignments allow the nucleotide sequence $accessions_no_version[$local_i] to retain eligibility for foosh? $all_output_data[$local_i][$FOOSH_COLUMN_INDEX]\n";
	close(VALIDATION);
    }
    for($local_j = 0; $local_j < $local_num_column_headers; $local_j++) {
	print OUTPUT "$column_headers_output[$local_j]\t";
    }
    print OUTPUT "\n";
    for ($local_i = 0; $local_i < $local_num_accessions; $local_i++) {
	for($local_j = 0; $local_j < $local_num_column_headers; $local_j++) {
	    print OUTPUT "$all_output_data[$local_i][$local_j]\t";
	}
	print OUTPUT "\n";
    }

}

# Subroutine: process_blastx_summary_file()
# Synopsis: Find the alignments summary for one accession and process into the ith row of all_output_data
#
# Args: $accession_index, $blastx_directory
#
# Returns: nothing

sub process_blastx_summary_file {
    my $sub_name = "process_blastx_summary_file()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_accession_index, $local_blastx_directory) = @_;    
    my $local_accession_name;           #accession without version
    my $local_blastx_filename;          #name of file that has the blastx summary for this accession
    my $local_nextline;                 #one line of blastx summary file
    my @local_fields;                   #array to split local_nextline
    my $local_matching_protein;         #the name of the protein for hich we have a match
    my $local_protein_index;            #the index in @protein_names of the protein for hich we have a match
    my $local_protein_accession;        #the accession of the protein for which we have a match
    my $local_i = -1;                   #the index of the output data to update for this accession
    my $local_j;                        #loop index
    my $local_hsp = 0;                  #count of hsps in blastx summary output for one matching protein
    my $local_qstart;                   #start of an alignment in the query
    my $local_qend;                     #end of an alignment in the query
    my $local_sstart;                   #start of an alignment in the subject
    my $local_send;                     #end of an alignment in the subject
    my $local_slen;                     #length of subject in amino acids
    
    
    $local_accession_name = $accessions_no_version[$local_accession_index];
    $local_blastx_filename = $local_blastx_directory . $local_accession_name . ".blastx.summary.txt";
    if ($DEBUG) {
	print "Called process_blastx_summary_file for $local_blastx_filename\n";
    }
    open(BLASTX,  "<", $local_blastx_filename)                or die "Cannot open 6 $local_blastx_filename for input\n";
    while(defined($local_nextline = <BLASTX>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
	if ($local_fields[0] eq "HDEF") {
	    ($local_protein_accession,$local_matching_protein) = ($local_fields[1] =~ m/^(\S+)\s+(\S+)/);  #match first and second tokens
	    if ($DEBUG) {
		print "Found HDEF line for $local_matching_protein with accession $local_protein_accession\n";
	    } 
	    $local_protein_index = -1;
	    for ($local_j = 0; $local_j < $num_proteins; $local_j++)   {
		if ($protein_names[$local_j] eq $local_matching_protein) {
		    $local_protein_index = $local_j;
		    $protein_accessions[$local_j] = $local_protein_accession;
		}
	    }
	    if (-1 == $local_protein_index) {
		print "ERROR: Unrecognized matching protein $local_matching_protein";
	    }
	    $local_hsp = 0;
	}
	if ($local_fields[0] eq "HSP") {
	    $local_hsp = $local_fields[1];
	}
	if ($local_fields[0] eq "SLEN") {
	    ($local_slen) = ($local_fields[1] =~ m/^(\d+)/);
	}
	if (1 == $local_hsp) {
	    if ($local_fields[0] eq "INS") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_INS;
		$all_output_data[$local_accession_index][$local_i] = $local_fields[1];
	    }
	    if ($local_fields[0] eq "DEL") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_DEL;
		$all_output_data[$local_accession_index][$local_i] = $local_fields[1];
	    }
	    if ($local_fields[0] eq "MAXIN") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_MAX_INDEL;
		if ($local_fields[1] > $all_output_data[$local_accession_index][$local_i]) {
		    $all_output_data[$local_accession_index][$local_i] = $local_fields[1];
		}
	    }
	    if ($local_fields[0] eq "MAXDE") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_MAX_INDEL;
		if ($local_fields[1] > $all_output_data[$local_accession_index][$local_i]) {
		    $all_output_data[$local_accession_index][$local_i] = $local_fields[1];
		}
	    }	    
	    if ($local_fields[0] eq "STOP") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_TRC;
		$all_output_data[$local_accession_index][$local_i] = $local_fields[1];
	    }
	    if ($local_fields[0] eq "FRAME") {
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_FRAME;
		$all_output_data[$local_accession_index][$local_i] = $local_fields[1];
	    }
	    if ($local_fields[0] eq "QRANGE") {
		($local_qstart) = ($local_fields[1] =~ m/^(\d+)/);
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_START_PROT;
		$all_output_data[$local_accession_index][$local_i] = $local_qstart;
		($local_qend) = ($local_fields[1] =~ m/(\d+)$/);
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_PROT;
		$all_output_data[$local_accession_index][$local_i] = $local_qend;
	    }
	    #SRANGE is the last field for the hsp, so we can take cross-filed decisions after seeing SRANGE
	    if ($local_fields[0] eq "SRANGE") {
		($local_sstart) = ($local_fields[1] =~ m/^(\d+)/);
		($local_send) = ($local_fields[1] =~ m/(\d+)$/);
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_NTERM;
		if ((1 == $local_sstart) || (1 == $local_send)) {
		    $all_output_data[$local_accession_index][$local_i] = "Yes"
		}
		else {
		    $all_output_data[$local_accession_index][$local_i] = "No"
		}
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_CTERM;
		if (($local_slen == $local_sstart) || ($local_slen == $local_send)) {
		    $all_output_data[$local_accession_index][$local_i] = "Yes"
		}
		else {
		    $all_output_data[$local_accession_index][$local_i] = "No"		    
		}
		$local_i = $num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_CODON;
		if ((!("NP" eq $all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_NUC])) &&
		    (!("?" eq $all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_NUC])) &&
		    (!("NP" eq $all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_PROT])) &&
		    (!("?" eq $all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_PROT])) &&
			     3 == abs($all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_NUC] -
			  $all_output_data[$local_accession_index][$num_generic_output_headers + ($local_protein_index * $num_specific_output_headers) + $HEADER_OFFSET_STOP_PROT])) {
		    if (($local_slen == $local_sstart) || ($local_slen == $local_send)) {
			$all_output_data[$local_accession_index][$local_i] = "Yes";
		    }
		    else {
			$all_output_data[$local_accession_index][$local_i] = "No";			
		    }
		}
		else {
		    $all_output_data[$local_accession_index][$local_i] = "MOOT";					    
		}
	    }	    	    
	}
    }
    close(BLASTX);
}

