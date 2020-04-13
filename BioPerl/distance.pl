# Benchmark the calculation of a distance in a PDB file
# The distance is the closest distance between any atoms of residues 50 and 60
#   of chain A in 1AKE

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $pdb_filepath = "data/1AKE.pdb";
my $structio = Bio::Structure::IO->new(-file => $pdb_filepath);
my $struc = $structio->next_structure;

sub distance {
    my @coords_50 = ();
    my @coords_60 = ();
    for my $chain ($struc->get_chains) {
        if ($chain->id eq "A") {
            for my $res ($struc->get_residues($chain)) {
                if (substr($res->id, -3, 3) eq "-50") {
                    for my $atom ($struc->get_atoms($res)) {
                        push @coords_50, [$atom->xyz];
                    }
                } elsif (substr($res->id, -3, 3) eq "-60") {
                    for my $atom ($struc->get_atoms($res)) {
                        push @coords_60, [$atom->xyz];
                    }
                }
            }
        }
    }
    my $min_sq_dist = "Infinity";
    for (my $i = 0; $i < scalar(@coords_50); $i++) {
        for (my $j = 0; $j < scalar(@coords_60); $j++) {
            my $sq_dist = ($coords_50[$i][0]-$coords_60[$j][0]) ** 2 + ($coords_50[$i][1]-$coords_60[$j][1]) ** 2 + ($coords_50[$i][2]-$coords_60[$j][2]) ** 2;
            if ($sq_dist < $min_sq_dist) {
                $min_sq_dist = $sq_dist;
            }
        }
    }
    return sqrt($min_sq_dist);
}

my $start = time();
distance();
my $end = time();

print $end - $start, "\n";
