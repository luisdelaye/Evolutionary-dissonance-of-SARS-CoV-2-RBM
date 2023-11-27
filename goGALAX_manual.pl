#!/usr/bin/perl

use strict;

my @mbf = glob("../MrBayes/*.run2.t");

for (my $i = 0; $i <= $#mbf; $i++){
    $mbf[$i] =~ s/\.\.\/MrBayes\///;
    print ("$mbf[$i]\n");
    for (my $j = $i; $j <= $#mbf; $j++){
        $mbf[$j] =~ s/\.\.\/MrBayes\///;
        print ("\t$mbf[$j]\n");
        open (ROB, ">infile") or die ("no puedo abrir infile\n");
        print ROB ("/Users/jose/Documents/Proyectos/Entropy/spike14_manual/MrBayes/$mbf[$i]\n");
        print ROB ("/Users/jose/Documents/Proyectos/Entropy/spike14_manual/MrBayes/$mbf[$j]\n");
        close (ROB);
        system ("cat infile");
        my $outfile1 = $mbf[$i];
        my $outfile2 = $mbf[$j];
        $outfile1 =~ s/.nex.run2.t//;
        $outfile2 =~ s/.nex.run2.t//;
        my $outfile = $outfile1.'_AND_'.$outfile2.'_run2';
        print ("\t\toutfile: $outfile\n");
        system ("galax1 --listfile infile --skip 200 --outgroup 5 --details");
        my $outf = $outfile.'_output-galax.txt';
        system ("mv output-galax.txt $outf");
        $outf = $outfile.'_output-galax-details.txt';
        system ("mv output-galax-details.txt $outf");
        $outf = $outfile.'_output-galax-merged.tre';
        system ("mv output-galax-merged.tre $outf");
        system ("ls -l");
        my $pausa = <STDIN>;
    }
}
