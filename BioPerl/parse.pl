# Benchmark the parsing of a PDB file given as an argument

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $pdb_filepath = $ARGV[0];

sub parse {
    my $structio = Bio::Structure::IO->new(-file => $pdb_filepath);
    return $structio->next_structure;
}

my $start = time();
parse();
my $end = time();

print $end - $start, "\n";
