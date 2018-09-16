#!/usr/bin/perl -w
# the first line of perl code has to be above
# Author: Alejandro Schaffer with help from Eric Nawrocki
#
# Code to do Norovirus protein validation in several high-level steps
#
#
# Usage: NPV.pl --input_fasta <input fasta file> \                  [REQUIRED]
#               --input_dir <directory for dnaorg_classify>  \      [REQUIRED]
#               --refseq    <six digits of refseq>           \      [REQUIRED]
#               --protein_names <names in this genome>       \      [REQUIRED]
#               --align_tolerance <difference between alignments \  [REQUIRED]
#               --indel_tolerance <allowed size of indel>         \ [REQUIRED]
#               --output_dir <output directory>                   \ [REQUIRED]
#               --verbose                                         \ [OPTIONAL]
#               --keep                                              [OPTIONAL]

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday); # for timings

require "epn-options.pm";
require "dnaorg.pm";

#Variables for handling command-line parameters
my $output_dir       = undef; # directory for output files
my $input_fasta_file = undef; # input fasta file
my $input_dir        = undef; # directory for input files
my $refseq           = undef; # six digits of refseq
my $protein_names    = undef; # names of proteins in this genome
my $align_tolerance  = undef; # tolerance for difference between endpoints of nucleotide and protein aligments
my $indel_tolerance  = undef; # tolerance for largest indel
my $verbose_mode;             # did user specify verbose mode, in which case extra columns are printed
my $verbose_string;           # string to add as argument to called programs depending on $verbose_mode 
my $keep_mode;                # did user specify keep mode
my $keep_string;              # string to add as argument to called programs depending on $verbose_mode
my $NPV_exec_dir  = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/2018.02/Norovirus_proteins/repository/";  #directory with code files

my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();


# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option              type       default group requires incompat preamble-outfile                                       help-outfile
opt_Add("-h",               "boolean", 0,         0,    undef, undef,  undef,                                                 "display this help",                                         \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input_fasta",     "string",  undef,     1,    undef, undef,  "input fasta file",                    "REQUIRED: file name <s> with sequences in fasta format",    \%opt_HH, \@opt_order_A);
opt_Add("--input_dir",       "string",  undef,     1,    undef, undef,  "input directory",                     "REQUIRED: input directory <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--refseq",          "string",  undef,     1,    undef, undef,  "refseq",                              "REQUIRED: refseq in six digits <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--protein_names",   "string",  undef,     1,    undef, undef,  "protein names",                       "REQUIRED: file of protein names in genome <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--align_tolerance", "integer", undef,     1,    undef, undef,  "alignment tolerance",                 "REQUIRED: difference in alignment lengths tolerated <n>",    \%opt_HH, \@opt_order_A);
opt_Add("--indel_tolerance", "integer", undef,     1,    undef, undef,  "indel tolerance",                     "REQUIRED: largest indel tolerated <n>",    \%opt_HH, \@opt_order_A);
opt_Add("--output_dir",      "string",  undef,     1,    undef, undef,  "output directory",                    "REQUIRED: output directory <s>",    \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options (not required)";
opt_Add("--verbose",        "boolean", 0,         2,    undef, undef,  "output more intermediate information",                      "output more intermediate information",                            \%opt_HH, \@opt_order_A);
opt_Add("--keep",           "boolean", 0,         2,    undef, undef,  "keep all intermediate files (e.g. protein predictions)", "keep all intermediate files (e.g. protein predictions)",       \%opt_HH, \@opt_order_A);


my $synopsis = "NPV.pl :: run protein validation for Norovirus";
my $usage    = "Usage:\n";
$usage      .= "\tNPV.pl --input_fasta <input fasta file> \\\n";
$usage      .= "\t--input_dir <directory for dnaorg_classify>  \\\n";
$usage      .= "\t--refseq    <six digits of refseq>           \\\n";
$usage      .= "\t--protein_names <names in this genome>       \\\n";
$usage      .= "\t--align_tolerance <difference between alignments \\\n";
$usage      .= "\t--indel_tolerance <allowed size of indel>        \\\n";
$usage      .= "\t--output_dir <output directory>                  \\\n";
$usage      .= "\t--verbose <boolean, default 0> \\\n";
$usage      .= "\t--keep    <boolen, default 0>  \\\n";
$usage      .= "\nFor example:\n";
$usage      .= "NPV.pl --input_fasta seqs.fa ";
$usage      .= "--input_dir classify-dir1 ";
$usage      .= "--refseq 029646 ";
$usage      .= "--protein_names protein_names3.txt ";
$usage      .= "--align_tolerance 5 ";
$usage      .= "--indel_tolerance 27 ";
$usage      .= "--output_dir . ";
$usage      .= "--keep\n ";                  

my $total_seconds = -1 * seconds_since_epoch(); # by multiplying by -1, we can just add another seconds_since_epoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.01";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions('h'                  => \$GetOptions_H{"-h"},
		'input_fasta=s'      => \$GetOptions_H{"--input_fasta"},
		'input_dir=s'        => \$GetOptions_H{"--input_dir"},
		'refseq=s'           => \$GetOptions_H{"--refseq"},
		'protein_names=s'    => \$GetOptions_H{"--protein_names"},						
		'align_tolerance=n'  => \$GetOptions_H{"--align_tolerance"},
		'indel_tolerance=n'  => \$GetOptions_H{"--indel_tolerance"},		
		'output_dir=s'       => \$GetOptions_H{"--output_dir"},
		'verbose'            => \$GetOptions_H{"--verbose"},
		'keep'               => \$GetOptions_H{"--keep"});

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# retrieve options
$output_dir       = opt_Get("--output_dir", \%opt_HH);
$input_fasta_file = opt_Get("--input_fasta", \%opt_HH);
$input_dir        = opt_Get("--input_dir", \%opt_HH);
$refseq           = opt_Get("--refseq", \%opt_HH);
$protein_names    = opt_Get("--protein_names", \%opt_HH);
$align_tolerance  = opt_Get("--align_tolerance", \%opt_HH);
$indel_tolerance  = opt_Get("--indel_tolerance", \%opt_HH);
$verbose_mode     = opt_Get("--verbose", \%opt_HH);         
$keep_mode        = opt_Get("--keep", \%opt_HH);

# exit if necessary (if options were messed up)
# first, determine if all required options were actually used
my $reqopts_errmsg = "";
if(! defined $output_dir)        { $reqopts_errmsg .= "ERROR, --output_dir option not used. It is required.\n"; }
if(! defined $input_fasta_file)  { $reqopts_errmsg .= "ERROR, --input_fasta option not used. It is required.\n"; }
if(! defined $input_dir)         { $reqopts_errmsg .= "ERROR, --input_dir option not used. It is required.\n"; }
if(! defined $refseq)            { $reqopts_errmsg .= "ERROR, --refseq option not used. It is required.\n"; }
if(! defined $protein_names)     { $reqopts_errmsg .= "ERROR, --protein_names option not used. It is required.\n"; }
if(! defined $align_tolerance)   { $reqopts_errmsg .= "ERROR, --align_tolerance option not used. It is required.\n"; }
if(! defined $indel_tolerance)   { $reqopts_errmsg .= "ERROR, --indel_tolerance option not used. It is required.\n"; }

if(($reqopts_errmsg ne "") || (! $all_options_recognized) || ($GetOptions_H{"-h"})) {
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if   ($reqopts_errmsg ne "")     { die $reqopts_errmsg; }
  elsif(! $all_options_recognized) { die "ERROR, unrecognized option;"; }
  else                             { exit 0; } # -h, exit with 0 status
}

$verbose_string = ($verbose_mode) ? "--verbose" : "";
$keep_string    = ($keep_mode)    ? "--keep"    : "";

# make sure the executable files we need exist 
my %execs_H = (); # hash with paths to all required executables
$execs_H{"break_fasta"}                  = $NPV_exec_dir . "break_fasta.pl";
$execs_H{"simplify_table2"}              = $NPV_exec_dir . "simplify_table2.pl";
$execs_H{"list_predictions"}             = $NPV_exec_dir . "list_predictions.pl";
$execs_H{"run_blast_validation_scripts"} = $NPV_exec_dir . "run_blast_validation_scripts.pl";
$execs_H{"compare_predictions"}         = $NPV_exec_dir . "compare_predictions.pl";
validate_executable_hash(\%execs_H);

# set output file names

my $nextline;  #one line of input
my $num_proteins; #number of proteins in this genome

#find number of proteins
open(PROTEINS, "<", $protein_names)          or die "Cannot open 1 $protein_names for input\n";
$nextline = <PROTEINS>;
chomp($nextline);
$num_proteins = $nextline;
close(PROTEINS);

###############################
# Step 1. Call break_fasta.pl #
###############################

my $progress_w = 51; # the width of the left hand column in our progress output, hard-coded
my $start_secs = output_progress_prior("Running break_fasta.pl", $progress_w, undef, *STDOUT);
my $cwd =getcwd();
chdir ".." ;
run_command($execs_H{"break_fasta"} . " " . $cwd . "/" . "$input_fasta_file", 0); # 0: don't echo command to STDOUT
chdir $cwd;
my $desc_str = sprintf("split input FASTA into individual files");
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 2. Simplify predictions #
##################################

$start_secs = output_progress_prior("Simplifying predictions", $progress_w, undef, *STDOUT);
my $file_substring1 = $input_dir ."-NC_" . $refseq .".pass";
my $input_file_call1 = $input_dir . "/" . $file_substring1 . "/" . $file_substring1 . ".dnaorg_annotate.tbl";
my $temp_output_file1 = "NC_" . $refseq . "_table_simplified.txt";
run_command($execs_H{"simplify_table2"} . " --input_file $input_file_call1  --output_file $temp_output_file1 ", 0); # 0: don't echo command to STDOUT
    $desc_str = sprintf("output saved as $temp_output_file1 %s", $keep_mode ? "]" : " (temporarily)");
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 3. List predictions #
##################################

$start_secs = output_progress_prior("Listing predictions", $progress_w, undef, *STDOUT);
my $temp_output_file2 = "NC_" . $refseq ."_listed_predictions.txt";
my $temp_output_file3 = "NC_" . $refseq . "_accessions_to_test.txt";
run_command($execs_H{"list_predictions"} . " --input_predict_file $temp_output_file1  --input_names_file $protein_names --output_file $temp_output_file2 --output_accessions $temp_output_file3 ", 0); # 0: don't echo command to STDOUT
    $desc_str = sprintf("output saved as $temp_output_file2 and $temp_output_file3 %s", $keep_mode ? "]" : " (temporarily)");
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 4. Run blastx #
##################################

$start_secs = output_progress_prior("Running blastx", $progress_w, undef, *STDOUT);
my $blast_database = "../../" . $refseq . "db";
run_command($execs_H{"run_blast_validation_scripts"} . " --input $temp_output_file3  --seq_directory ../ --num_proteins 3 --db $blast_database ", 0); # 0: don't echo command  to STDOUT
unlink "validation.qsub.script";
chdir ".." ;
unlink <*blastx.*sh>;
if (!$keep_mode) {
    unlink <*blastx.*out>;    
}
chdir $cwd;    
unlink <blastrun*>;
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 5. Compare predictions #
##################################

$start_secs = output_progress_prior("Listing predictions", $progress_w, undef, *STDOUT);
my $output_table_name = $refseq . "_nucleotide_annotations_validation_table.txt";
run_command($execs_H{"compare_predictions"} . " --input_nucleotide_table $temp_output_file2  --input_blastx_directory ../ --input_accessions $temp_output_file3 --input_protein_names $protein_names --alignment_tolerance $align_tolerance --indel_tolerance $indel_tolerance --output_table $output_table_name --output_directory ./ ", 0); # 0: don't echo command to STDOUT
    $desc_str = sprintf("output saved as $temp_output_file2 and $temp_output_file3 %s", $keep_mode ? "]" : " (temporarily)");
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 6. Clean up #
##################################
if (!$keep_mode) {
    unlink $temp_output_file1;
    unlink $temp_output_file2;
    unlink $temp_output_file3;
    chdir ".." ;
    unlink <*.na>;
    chdir $cwd;
}


#####################################################################
# SUBROUTINES
#####################################################################
# List of subroutines:
#
#
# Miscellaneous functions:
# output_progress_prior:      output routine for a step, prior to running the step
# output_progress_complete:   output routine for a step, after the running the step
# run_command:              run a command using system()
# seconds_since_epoch:      number of seconds since the epoch, for timings

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016
#
# Purpose:     Runs a command using system() and exits in error
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
    my $sub_name = "run_command()";
    my $nargs_expected = 2;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($cmd, $be_verbose) = @_;

    if($be_verbose) {
	print ("Running cmd: $cmd\n");
    }

    my ($seconds, $microseconds) = gettimeofday();
    my $start_time = ($seconds + ($microseconds / 1000000.));

    system($cmd);

    ($seconds, $microseconds) = gettimeofday();
    my $stop_time = ($seconds + ($microseconds / 1000000.));

    if($? != 0) {
	die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine : seconds_since_epoch()
# Incept:      EPN, Sat Feb 13 06:17:03 2016
#
# Purpose:     Return the seconds and microseconds since the
#              Unix epoch (Jan 1, 1970) using
#              Time::HiRes::gettimeofday().
#
# Arguments:   NONE
#
# Returns:     Number of seconds and microseconds
#              since the epoch.
#
#################################################################
sub seconds_since_epoch {
  my $nargs_expected = 0;
  my $sub_name = "seconds_since_epoch()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($seconds, $microseconds) = gettimeofday();
  return ($seconds + ($microseconds / 1000000.));
}

#################################################################
# Subroutine : validate_executable_hash()
#
# Purpose:     Given a reference to a hash in which the
#              values are paths to executables, validate
#              that those files are executable.
#
# Arguments:
#   $execs_HR: REF to hash, keys are short names to executable
#              e.g. "vecscreen", values are full paths to that
#              executable, e.g. "/usr/local/infernal/1.1.1/bin/vecscreen"
#
# Returns:     void
#
# Dies:        if one or more executables does not exist
#
#################################################################
sub validate_executable_hash {
  my $nargs_expected = 1;
  my $sub_name = "validate_executable_hash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  my ($execs_HR) = (@_);

  my $fail_str = undef;
  foreach my $key (sort keys %{$execs_HR}) {
    if(! -e $execs_HR->{$key}) {
      $fail_str .= "\t$execs_HR->{$key} does not exist.\n";
    }
    elsif(! -x $execs_HR->{$key}) {
      $fail_str .= "\t$execs_HR->{$key} exists but is not an executable file.\n";
    }
  }

  if(defined $fail_str) {
    die "ERROR in $sub_name(),\n$fail_str";
  }

  return;
}	


#################################################################
# Subroutine : output_progress_prior()
# Incept:      EPN, Fri Feb 12 17:22:24 2016
#
# Purpose:      Output to $FH1 (and possibly $FH2) a message indicating
#               that we're about to do 'something' as explained in
#               $outstr.
#
#               Caller should call *this* function, then do
#               the 'something', then call output_progress_complete().
#
#               We return the number of seconds since the epoch, which
#               should be passed into the downstream
#               output_progress_complete() call if caller wants to
#               output running time.
#
# Arguments:
#   $outstr:     string to print to $FH
#   $progress_w: width of progress messages
#   $FH1:        file handle to print to, can be undef
#   $FH2:        another file handle to print to, can be undef
#
# Returns:     Number of seconds and microseconds since the epoch.
#
#################################################################
sub output_progress_prior {
    my $nargs_expected = 4;
    my $sub_name = "output_progress_prior()";
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
    my ($outstr, $progress_w, $FH1, $FH2) = @_;

    if(defined $FH1) { printf $FH1 ("# %-*s ... ", $progress_w, $outstr); }
    if(defined $FH2) { printf $FH2 ("# %-*s ... ", $progress_w, $outstr); }

    return seconds_since_epoch();
}

#################################################################
# Subroutine : output_progress_complete()
# Incept:      EPN, Fri Feb 12 17:28:19 2016 
#
# Purpose:     Output to $FH1 (and possibly $FH2) a
#              message indicating that we've completed
#              'something'.
#
#              Caller should call *this* function,
#              after both a call to output_progress_prior()
#              and doing the 'something'.
#
#              If $start_secs is defined, we determine the number
#              of seconds the step took, output it, and
#              return it.
#
# Arguments:
#   $start_secs:    number of seconds either the step took
#                   (if $secs_is_total) or since the epoch
#                   (if !$secs_is_total)
#   $extra_desc:    extra description text to put after timing
#   $FH1:           file handle to print to, can be undef
#   $FH2:           another file handle to print to, can be undef
#
# Returns:     Number of seconds the step took (if $secs is defined,
#              else 0)
#
#################################################################
sub output_progress_complete {
    my $nargs_expected = 4;
    my $sub_name = "output_progress_complete()";
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
    my ($start_secs, $extra_desc, $FH1, $FH2) = @_;

    my $total_secs = undef;
    if(defined $start_secs) {
	$total_secs = seconds_since_epoch() - $start_secs;
    }

    if(defined $FH1) { printf $FH1 ("done."); }
    if(defined $FH2) { printf $FH2 ("done."); }

    if(defined $total_secs || defined $extra_desc) {
	if(defined $FH1) { printf $FH1 (" ["); }
	if(defined $FH2) { printf $FH2 (" ["); }
    }
    if(defined $total_secs) {
	if(defined $FH1) { printf $FH1 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
	if(defined $FH2) { printf $FH2 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
    }
    if(defined $extra_desc) {
	if(defined $FH1) { printf $FH1 $extra_desc };
	if(defined $FH2) { printf $FH2 $extra_desc };
    }
    if(defined $total_secs || defined $extra_desc) {
	if(defined $FH1) { printf $FH1 ("]"); }
	if(defined $FH2) { printf $FH2 ("]"); }
    }

    if(defined $FH1) { printf $FH1 ("\n"); }
    if(defined $FH2) { printf $FH2 ("\n"); }

    return (defined $total_secs) ? $total_secs : 0.;
}
    
