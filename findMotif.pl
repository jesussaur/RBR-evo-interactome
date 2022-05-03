#!/usr/bin/perl

# findMotif.pl
# Juan Caballero <juan.caballero.perez@gmail.com>

use strict;
use warnings;

$ARGV[2] or die "usage: perl findMotif.pl MOTIF FASTA_IN FASTA_OUT\n";
my $motif  = shift @ARGV;
my $fa_in  = shift @ARGV;
my $fa_out = shift @ARGV;

my $pattern = $motif;
$pattern =~ s/x/./gi;

open (my $in,  "<", $fa_in)  or die "cannot open $fa_in\n";
open (my $out, ">", $fa_out) or die "cannot open $fa_out\n";
$/="\n>";
my $tot = 0;
my $mat = 0;
while (<$in>) {
    $tot++;
    s/>//;
    my ($id, @seq) = split (/\n/, $_);
    my $seq = join ("", @seq);

    if ($seq =~ /$pattern/) {
        print $out ">$_";
        $mat++;
    }
}
close $in;
close $out;

warn "Searched motif $motif in $tot sequences, $mat sequences matched\n";
