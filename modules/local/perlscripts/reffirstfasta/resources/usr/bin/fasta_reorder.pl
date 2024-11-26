#!/usr/bin/env perl
use strict;
use warnings;

sub croak { die " $0: @_: $!\n" }

sub usage {
    my $usage =<<'END_USAGE';
Usage: perl fasta_reorder.pl <file>.aln <seq_name> <order>

Takes in a multi-fasta file and moves <seq_name> to the to first or last sequence
 in the <out>.aln.

Arguments:
  <file>.aln    original fasta
  <seq_name>    name of sequence to move
  <order>       1: first place 2: last place
Output:
  reordered_<file>.aln  reordered fasta

Example:
  perl fasta_reoder.pl test.aln Reference 1 #moves '>Reference' to first sequence
  perl fasta_reoder.pl test.aln Sample2 2 #moves '>Sample2' to last sequence
END_USAGE

    return $usage;
}

# Check command line arguments
if (@ARGV != 3) {
    print usage();
    croak "Incorrect number of arguments.\n";
    exit 1;
}

our ($filename, $seq_name, $order) = @ARGV;

# Check order param
if ($order != 1 and $order != 2) {
    print usage();
    croak "Invalid <order>\n.";
    exit 1;
}

# Open the FASTA file
open(my $fh, '<', $filename) or croak "Could not open file '$filename'.";

my $query_seq = '';
my $other_seqs = '';
my $current_seq = '';
my $in_query = 0;
my $seq_count = 0;

while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>/) {
        $seq_count++;
        # If starting a new sequence header
        if ($line =~ /^>$seq_name/) {
            # Found the query sequence
            $in_query = 1;
            # Save any previously accumulated non-query sequence data
            $other_seqs .= $current_seq if $current_seq;
            $current_seq = $line . "\n";  # Start new sequence
        } else {
            # Save any previously accumulated sequence data
            if ($in_query){                
                $query_seq .= $current_seq if $current_seq;
            } else {
                $other_seqs .= $current_seq if $current_seq;
            }

            # Not in query sequence
            $in_query = 0;
            $current_seq = $line . "\n";  # Start new sequence
        }
    } else {
        # Continue accumulating sequence data
        $current_seq .= $line . "\n";
    }
}

# Append the last processed sequence to the appropriate variable
if ($in_query) {
    $query_seq .= $current_seq;
} else {
    $other_seqs .= $current_seq;
}


close($fh);

# Output the reordered sequences

# Generate the output filename
my $prefix = 'reordered_';
my $outfile = $prefix . $filename;

# Open for writing
open(my $out_fh, '>', $outfile) or croak "Could not write to file '$outfile'.";
if ($order == 1){
    print $out_fh $query_seq;
    print $out_fh $other_seqs;
} elsif ($order == 2){
    print $out_fh $other_seqs;
    print $out_fh $query_seq;
}

close($out_fh);

if ($order == 1){
    print "Reordering completed. Moved $seq_name to first sequence in $outfile.\n";
} elsif ($order == 2){
    print "Reordering completed. Moved $seq_name to last sequence in $outfile.\n";
}