#!/usr/bin/perl

use strict;

my $ML;
my %hashML;

open (MIA, "MLmanuales.txt") or die ("No puedo abrir MLmanuales.txt\n");
print ("frgmt_1s\tlnL1\tfrgmt_2s\tlnL2\tfrgmt_1y2\tlnM0\n");
while (my $linea = <MIA>){
    if ($linea =~ /^\s+Mean:   ([\d\.-]+)/){
      $ML = $1;
    }
    if ($linea =~ /More statistics on stepping-stone sampling is dumped to ([\w\d\.]+) file./){
      my $file = $1;
      if ($file =~ /spike14_([\d]+)_[\d]+.nex.ss/){
        my $and = $1.'_AND_'.$1;
        $hashML{$1} = $ML;
        print ("$1\t$ML\t$1\t$ML\t$and\tNA\n");
      } elsif ($file =~ /spike14_([\d]+)_to_[\d]+_AND_([\d]+)_to_[\d]+.nex.ss/) {
        if ($1 < $2){
          my $and = $1.'_AND_'.$2;
          print ("$1\t$hashML{$1}\t$2\t$hashML{$2}\t$and\t$ML\n");
        } else {
          my $and = $2.'_AND_'.$1;
          print ("$2\t$hashML{$2}\t$1\t$hashML{$1}\t$and\t$ML\n");
        }
      } else {
        die ("el patron fallo: $linea\n");
      }
    }
}
close (MIA);
