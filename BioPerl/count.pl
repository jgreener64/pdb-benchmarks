# Benchmark the counting of alanine residues in a PDB file

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $pdb_filepath = "data/1AKE.pdb";
my $structio = Bio::Structure::IO->new(-file => $pdb_filepath);
my $struc = $structio->next_structure;

sub count {
    my $c = 0;
    for my $chain ($struc->get_chains) {
        for my $res ($struc->get_residues($chain)) {
            if (substr($res->id, 0, 3) eq "ALA") {
                $c++;
            }
        }
    }
    return $c;
}

my $start = time();
count();
my $end = time();

print $end - $start, "\n";
