#!/usr/bin/perl

use strict;

my $file = $ARGV[0];
my $s    = $ARGV[1];
my $e    = $ARGV[2];

my $outfile = 'outfile_'.$s.'_'.$e.'.fasta';

my $r = 0;
my $sec;
my $name;

my @IDs;
my %genes;

open (MIA, "$file") or die ("No puedo abrir $file\n");
while (my $linea = <MIA>){
	chomp ($linea);
	if ($linea =~ />/ && $r == 1){
    $genes{$name} = $sec;
    push (@IDs, $name);
		$r = 0;
		$name = $sec = ();
	}
	if ($linea !~ />/ && $r == 1){
		$sec = $sec.$linea;
	}
	if ($linea =~ />/ && $r == 0){
		$name = $linea;
		$r = 1;
	}
}
if ($r == 1){
  $genes{$name} = $sec;
  push (@IDs, $name);
	$r = 0;
	$name = $sec = ();
}
close (MIA);

open (ROB, ">$outfile") or die ("No puedo abrir $outfile\n");
for (my $i = 0; $i <= $#IDs; $i++){
    my $sec = $genes{$IDs[$i]};
    $sec = '0'.$sec;
    my @asec = split (//, $sec);
    my @subs = @asec[$s..$e];
    print     ("$IDs[$i]\n");
    print ROB ("$IDs[$i]\n");
    my $c = 0;
    for (my $j = 0; $j <= $#subs; $j++){
        $c++;
        if ($c < 70){
            print     ("$subs[$j]");
            print ROB ("$subs[$j]");
        } else {
            print     ("$subs[$j]\n");
            print ROB ("$subs[$j]\n");
            $c = 0;
        }
    }
    print     ("\n");
    print ROB ("\n");
}
close (ROB);
