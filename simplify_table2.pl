#!/usr/bin/perl -w
# the first line of perl code has to be as above
#
# Author: Alejandro Schaffer
#
# Code to simplify the table of annotations to extract some rows and make them tab delimited
# Usage: simplify_table2.pl \
#        --input_table <input table> \
#        --output_table <simplified output table>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


# input/output file names
my $input_file;        #input file that summarizes predictions
my $output_file;       #output file that summarizes predictions

# variables used in processing input
my $nextline;  #one line of the input summary file
my @fields;    #entries in one row


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
opt_Add("--input_file",             "string", undef,      1,    undef, undef,   "input file of annotations",                                              "File name <s> of vecscreen summary",        \%opt_HH, \@opt_order_A);
opt_Add("--output_file",            "string", undef,      1,    undef, undef,   "output file to create",                                                  "Name <s> of output file to create",         \%opt_HH, \@opt_order_A);

my $synopsis = "simplify_table2.pl: simplify annotation towards protein validation\n";
my $usage    = "Usage:\n\n";
$usage      .= "\tsimplify_table2.pl.pl \\\n";
$usage      .= "\t--input_file <input table file> \\\n";
$usage      .= "\t--output_file <output file>\n\n";
$usage      .= "\nFor example:\n";
$usage      .= "simplify_table2.pl --input_summary all_annotations.txt ";
$usage      .= "--output_file simplified_annotations.txt\n";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay =
    &GetOptions('h'                            => \$GetOptions_H{"-h"},
		'input_file=s'              => \$GetOptions_H{"--input_file"},
		'output_file=s'                    => \$GetOptions_H{"--output_file"});

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
$input_file        = opt_Get("--input_file", \%opt_HH);
$output_file               = opt_Get("--output_file", \%opt_HH);

# die if any of the required options were not used
my $errmsg = undef;
if(! defined $input_file)        { $errmsg .= "ERROR, --input_file option not used.\n"; }
if(! defined $output_file)               { $errmsg .= "ERROR, --output_file option not used.\n"; }
if(defined $errmsg) {
    die $errmsg . "\n$usage\n";
}

# open output files
open(SUMMARY, "<", $input_file) or die "Cannot open 1 $input_file for input\n";
open(OUTPUT,  ">", $output_file)        or die "Cannot open 2 $output_file for output\n";

while(defined($nextline = <SUMMARY>)) {
    chomp($nextline);
    if ($nextline =~m/idx/) {
	digest_line($nextline);
    }     
    if ($nextline =~m/accession/) {
	digest_line($nextline);
    }
    if ($nextline =~m/totlen/) {
	digest_line($nextline);
    }    
    if (($nextline =~m/start1/) && ($nextline =~m/VP1/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VP1\:start1/VP1_start/;
	digest_line($nextline);
    }
    if (($nextline =~m/stop1/) && ($nextline =~m/VP1/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VP1\:stop1/VP1_stop/;
	digest_line($nextline);
    }
    if (($nextline =~m/start1/) && ($nextline =~m/VP2/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VP2\:start1/VP2_start/;
	digest_line($nextline);
    }
    if (($nextline =~m/stop1/) && ($nextline =~m/VP2/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VP2\:stop1/VP2_stop/;
	digest_line($nextline);
    }
    if (($nextline =~m/start1/) && ($nextline =~m/VF1/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VF1\:start1/VF1_start/;
	digest_line($nextline);
    }
    if (($nextline =~m/stop1/) && ($nextline =~m/VF1/)) {
	$nextline =~s/CDS \#\d \[single exon\; \+\]\:VF1\:stop1/VF1_stop/;
	digest_line($nextline);
    }                 
    if (($nextline =~m/polyprotein\:start/) && (! ($nextline =~m/polyprotein\:startc/))) {
	$nextline =~s/CDS\(MP\) \#1 \[6 mature_peptide\(s\)\]\:nonstructural polyprotein:start/polyprotein_start/;
	digest_line($nextline);
    }
    if (($nextline =~m/polyprotein\:stop/) && (! ($nextline =~m/polyprotein\:stopc/))) {
	$nextline =~s/CDS\(MP\) \#1 \[6 mature_peptide\(s\)\]\:nonstructural polyprotein:stop/polyprotein_stop/;	
	digest_line($nextline);
    }
}

################################################
# List of subroutines:
# 
# digest_line();   prints a line of input to the output with some replacements
#
################################################

# Subroutine: digest_line()
# Synopsis: prints a line of input to the output with some replacements
#
# Args: $input_line
#
# Returns: nothing

sub digest_line {

    my $sub_name = "digest_line()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_line) = @_;
    my @local_fields;      #entries on the line
    my $local_i;           #loop index
    my $local_num_fields;  #number of fields in the tab-delimited @local_line
    
    @local_fields = split /\s+/, $local_line;
    $local_num_fields = @local_fields;

    for($local_i = 0; $local_i < $local_num_fields; $local_i++) {
	print OUTPUT "$local_fields[$local_i]\t";
    }
    print OUTPUT "\n";

    return;
}
