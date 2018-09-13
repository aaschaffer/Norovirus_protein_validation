#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to convert kraken classification of reads to a summary table
# Usage: break_fasta.pl <fasta_file> 

use strict;

my $infile; #input FASTA file with many sequences
my $outfile; #output FASTA file with one sequence
my $nofile = 1;
my $nextline;
my $same_sequence = 2;
my $state = $nofile;
my $id;  #accession of sequence
my $id_no_version;  #accession of sequence without the version
my $new_file_name;

$infile = $ARGV[0];

open(SEQS, "<$infile") or die "Cannot open $infile\n";
while(defined($nextline = <SEQS>)) {
    chomp($nextline);
    if ($nofile == $state) {
	($id) =($nextline=~m/>(\S+)/);
	if ($id =~ m/\.\d+/) {
	    ($id_no_version) = ($id =~ m/^(\S+)\.\d+/);
	}
	else {
	    $id_no_version = $id;
	}
	$new_file_name = $id_no_version . "." . "na";
	open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
	$state = $same_sequence;
	print ONEFILE "$nextline";
	print ONEFILE "\n";
    }
    else {
	if ($same_sequence == $state) {
	    if ($nextline=~m/>/ ) {
		close(ONEFILE);
		($id) =($nextline=~m/>(\S+)/);
		if ($id =~ m/\.\d+/) {
		    ($id_no_version) = ($id =~ m/^(\S+)\.\d+/);
		}
		else {
		    $id_no_version = $id;
		}		
		$new_file_name = $id_no_version . "." . "na";
		open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
		$state = $same_sequence;
	    }
	    print ONEFILE "$nextline";
	    print ONEFILE "\n";
	}
    }
}
close(ONEFILE);
close(SEQS);
