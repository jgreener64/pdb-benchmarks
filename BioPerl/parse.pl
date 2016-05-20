# Benchmark the parsing of a PDB file given as an argument

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $pdb_filepath = $ARGV[0];
my $runs = 1;

my $times = 0.0;

sub parse {
    my $structio = Bio::Structure::IO->new(-file => $pdb_filepath);
    return $structio->next_structure;
}

for (my $i = 0; $i < $runs; $i++) {
    my $start_run = time();
    my $struc = parse();
    my $end_run = time();
    my $run_time = $end_run - $start_run;
    $times = $times + $run_time;
}

print "Average time per run: ", $times / $runs, "\n";
