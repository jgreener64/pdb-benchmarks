# Benchmark the parsing of a PDB file given as an argument

use Bio::Structure::IO;
use Time::HiRes qw(time);
use strict;

my $file = $ARGV[0];
my $runs = 5;

my $times = 0.0;

for( $a = 0; $a < $runs; $a = $a + 1 ){
    my $start_run = time();
    my $structio = Bio::Structure::IO->new(-file => $file);
    my $struc = $structio->next_structure;
    my $end_run = time();
    my $run_time = $end_run - $start_run;
    $times = $times + $run_time;
}

print "Average time per run: ", $times / $runs, "\n";
