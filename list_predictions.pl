#!/usr/bin/perl -w
# the first line of perl code has to be as above
#
# Author: Alejandro Schaffer
#
# Code to simplify the table of annotations to extract some rows and make them tab delimited
# Usage: list_predictions.pl \
#        --input_predict_table <input prediction table> \
#        --input_names_table <input protein names table> \
#        --output_table <simplified output table> \
#        --output_accessions <output file of accessions to check>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


# input/output file names
my $input_predict_file;        #input file that summarizes predictions
my $input_names_file;          #input file with protein names
my $output_file;               #output file that summarizes predictions
my $accessions_file;           #output file that lists accessions to check

# variables used in processing input
my $nextline;                  #one line of the input summary file
my @fields;                    #entries in one row
my @protein_names;             #array of protein names
my $num_proteins;              #number of proteins in this genome
my $num_headers;               #number of column headers
my @column_headers;            #array of column headers
my $first_index;               #first index used in reading input table
my $second_index;              #second index used in reading input table
my $accession_index;           #index on accessions
my $num_accessions;            #number of accesions to consider
my $num_values_on_row;         #number of values on a single input row 
my $i;                         #loop index
my $j;                         #loop index
my @all_values;                #all the values in the input file in a two-dimensional array, where the first dimension is by accession
my @exists_prediction;         #is there at least one start or stop prediction for this accession
my $short_accession;           #accession without version
my $num_generic_headers = 3;   #number of columns not specific to any protein
my $DEBUG = 0;                 #for debugging

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
opt_Add("-h",                          "boolean",0,          0,    undef, undef,   undef,                                                                          "display this help",                         \%opt_HH, \@opt_order_A);
opt_Add("--input_predict_file",             "string", undef,      1,    undef, undef,   "input file of annotations",                                              "File name <s> of predictions",        \%opt_HH, \@opt_order_A);
opt_Add("--input_names_file",             "string", undef,      1,    undef, undef,   "input file of names",                                              "File name <s> of protein list",        \%opt_HH, \@opt_order_A);
opt_Add("--output_file",            "string", undef,      1,    undef, undef,   "output file to create",                                                  "Name <s> of output file to create",         \%opt_HH, \@opt_order_A);
opt_Add("--output_accessions",            "string", undef,      1,    undef, undef,   "output file with accessions tp check",           "Name <s> of output file with accessions to check",         \%opt_HH, \@opt_order_A);

my $synopsis = "list_predictions.pl: simplify annotation towards protein validation\n";
my $usage    = "Usage:\n\n";
$usage      .= "\tlist_predictions.pl \\\n";
$usage      .= "\t--input_predict_file <input table file> \\\n";
$usage      .= "\t--input_names_file <input names file> \\\n";
$usage      .= "\t--output_file <output file>\\\n";
$usage      .= "\t--output_accessions <output file with accessions to check>\n\n";
$usage      .= "\nFor example:\n";
$usage      .= "list_predictions.pl --input_predict_file all_annotations.txt ";
$usage      .= "--input_names_file protein_names.txt ";
$usage      .= "--output_file simplified_annotations.txt ";
$usage      .= "--output_accessions accessions_to_check.txt ";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay =
    &GetOptions('h'                            => \$GetOptions_H{"-h"},
		'input_predict_file=s'           => \$GetOptions_H{"--input_predict_file"},
		'input_names_file=s'              => \$GetOptions_H{"--input_names_file"},		
		'output_file=s'                  => \$GetOptions_H{"--output_file"},
		'output_accessions=s'                  => \$GetOptions_H{"--output_accessions"});

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


# define file names
$input_predict_file        = opt_Get("--input_predict_file", \%opt_HH);
$input_names_file        = opt_Get("--input_names_file", \%opt_HH);
$output_file               = opt_Get("--output_file", \%opt_HH);
$accessions_file         = opt_Get("--output_accessions", \%opt_HH);

# die if any of the required options were not used
my $errmsg = undef;
if(! defined $input_predict_file)        { $errmsg .= "ERROR, --input_predict_file option not used.\n"; }
if(! defined $input_names_file)        { $errmsg .= "ERROR, --input_names_file option not used.\n"; }
if(! defined $output_file)               { $errmsg .= "ERROR, --output_file option not used.\n"; }
if(! defined $accessions_file)               { $errmsg .= "ERROR, --output_accessions option not used.\n"; }
if(defined $errmsg) {
    die $errmsg . "\n$usage\n";
}

# open output files
open(SUMMARY, "<", $input_predict_file)         or die "Cannot open 1 $input_predict_file for input\n";
open(PROTEINS, "<", $input_names_file)          or die "Cannot open 2 $input_names_file for input\n";
open(OUTPUT,  ">", $output_file)                or die "Cannot open 3 $output_file for output\n";
open(ACCESSIONS,  ">", $accessions_file)        or die "Cannot open 4 $accessions_file for output\n";


$num_proteins = read_protein_names();
close(PROTEINS);
$num_headers = make_column_headers($num_proteins);
$first_index = 0;
$second_index = 0;
$accession_index  =0;
while(defined($nextline = <SUMMARY>)) {
    chomp($nextline);
    $second_index = $first_index % ($num_headers + 1);
    if (($DEBUG) && ($first_index > 49)) {
	print "first $first_index  second $second_index $nextline\n";
    } 
    if ($second_index != ($num_headers)) {
	@fields = split /\t/, $nextline;
	$num_values_on_row = @fields;
	if (0 == $second_index) {
	    for($i = 0; $i < ($num_values_on_row-1); $i++) {
		$exists_prediction[$accession_index + $i]  = 0;
		if (($DEBUG) && ($first_index > 49)) {
		    my $temp = $accession_index + $i; 
		    print "Setting exists_prediction index $temp to zero\n";
		} 
	    }
	}
	for($i = 1; $i < $num_values_on_row; $i++) {
	    $all_values[$accession_index + $i -1][$second_index] = $fields[$i];
	    
	}
    }
    else {
	if (($DEBUG) && ($first_index > 49)) {
	    print "testing for predictions on row $first_index\n";
	}
	for($i = 0; $i < ($num_values_on_row-1); $i++) {
	    for ($j = 0; $j < $num_proteins; $j++) {
	       if ((!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)] eq "NP")) && 
		   (!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)] eq "?")) &&
		   (!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)] =~m/\[/)) &&
		   (!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)+1] =~m/\[/)) &&
		   (!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)+1] eq "NP")) && 
		   (!($all_values[$accession_index + $i][$num_generic_headers + (2 * $j)+1] eq "?"))) {
		   $exists_prediction[$accession_index + $i] = 1;
		   if (($DEBUG) && ($first_index > 49)) {
		    my $temp2 = $accession_index + $i; 
		    print "Setting exists_prediction index $temp2 to one\n";
		   }
	       }
	    }
	}
	$accession_index = $accession_index + $num_values_on_row -1;
    }
    $first_index++;
}
$num_accessions = $accession_index;
for($i = 0; $i < $num_headers; $i++) {
    print OUTPUT "$column_headers[$i]\t";
}
print OUTPUT "\n";
for ($j = 0; $j < $num_accessions; $j++) {
    if (1 == $exists_prediction[$j]) {
	for($i = 0; $i < $num_headers; $i++) {
	    print OUTPUT "$all_values[$j][$i]\t";
	}
	print OUTPUT "\n";
	if ($all_values[$j][1] =~m/\.\d+/) {
	    ($short_accession) = ($all_values[$j][1] =~m/(\S+)\.\d+/);
	}
	else {
	    $short_accession = $all_values[$j][1];
	}
	print ACCESSIONS "$short_accession\n";
    }
}
close(OUTPUT);
close(ACCESSIONS);

################################################
# List of subroutines:
#
# read_protein_names();    Read the file of protein names
# make_column_headers();   Generate the values for the array column_headers
#
################################################

# Subroutine: read_protein_names()
# Synopsis: Read the file of protein names
#
# Args: none
#
# Returns: number of proteins

sub read_protein_names {

    my $sub_name = "read_protein_names()";
    my $nargs_exp = 0;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my $local_nextline;       #one line of the file
    my $local_num_proteins;   # number of proteins seen so far
    my $local_i;              # loop index
    
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
    my $local_num_headers = 0;     # number of columns
    my $local_i;                   # loop index
    
    
    $column_headers[$local_num_headers] = "index";
    $local_num_headers++;
    $column_headers[$local_num_headers] = "accession";
    $local_num_headers++;
    $column_headers[$local_num_headers] = "totlen";
    $local_num_headers++;
    for($local_i = 0; $local_i < $local_num_proteins; $local_i++) {
	$column_headers[$local_num_headers] = $protein_names[$local_i] . "_start";
	$local_num_headers++;
	$column_headers[$local_num_headers] = $protein_names[$local_i] . "_stop";
	$local_num_headers++;	
    }
    return $local_num_headers;
}
