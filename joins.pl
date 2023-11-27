#!/usr/bin/perl

use strict;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my $outfile = 'outfile_'.$file1.'_'.$file2.'.fasta';

my $r = 0;
my $sec;
my $name;

my @IDs;
my %genes1;

open (MIA, "$file1") or die ("No puedo abrir $file1\n");
while (my $linea = <MIA>){
	chomp ($linea);
	if ($linea =~ />/ && $r == 1){
    $genes1{$name} = $sec;
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
  $genes1{$name} = $sec;
  push (@IDs, $name);
	$r = 0;
	$name = $sec = ();
}
close (MIA);

my %genes2;

open (MIA, "$file2") or die ("No puedo abrir $file2\n");
while (my $linea = <MIA>){
	chomp ($linea);
	if ($linea =~ />/ && $r == 1){
    $genes2{$name} = $sec;
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
  $genes2{$name} = $sec;
	$r = 0;
	$name = $sec = ();
}
close (MIA);

open (ROB, ">$outfile") or die ("No puedo abrir $outfile\n");
for (my $i = 0; $i <= $#IDs; $i++){
    print     ("$IDs[$i]\n");
    print ROB ("$IDs[$i]\n");
    my $concat = $genes1{$IDs[$i]}.$genes2{$IDs[$i]};
    my @concat = split (//, $concat);
    my $c = 0;
    for (my $j = 0; $j <= $#concat; $j++){
        $c++;
        if ($c < 70){
            print     ("$concat[$j]");
            print ROB ("$concat[$j]");
        } else {
            print     ("$concat[$j]\n");
            print ROB ("$concat[$j]\n");
            $c = 0;
        }
    }
    print     ("\n");
    print ROB ("\n");
}
close (ROB);
