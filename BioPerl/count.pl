# Benchmark the counting of alanine residues in a PDB file

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $pdb_filepath = "pdbs/1AKE.pdb";
my $runs = 1;

my $structio = Bio::Structure::IO->new(-file => $pdb_filepath);
my $struc = $structio->next_structure;
my $times = 0.0;

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

for (my $i = 0; $i < $runs; $i++) {
    my $start_run = time();
    my $c = count();
    my $end_run = time();
    my $run_time = $end_run - $start_run;
    $times = $times + $run_time;
}

print "Average time per run: ", $times / $runs, "\n";
