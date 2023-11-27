#!/usr/bin/perl

# This script fragments a multiple sequence alignment into smaller alignments.

# use
# perl fragment.pl sequence_alignment.fasta window_size steps

# out: fragments/fragment_1_to_n.fasta ...

# See https://github.com/luisdelaye/X/ for more details.

# Author
# Luis Jose Delaye Arredondo
# Laboratorio de Genomica Evolutiva
# Departamento de Ingenieria Genetica
# CINVESTAV - Irapuato
# Mexico
# luis.delaye@cinvestav.mx
# Copy-Left  = : - )

# beta.1.0 version

#-------------------------------------------------------------------------------

use strict;

my $file  = $ARGV[0];
my $wsize = $ARGV[1];
my $steps = $ARGV[2];

my %secs;
my @head;

my $r = 0;
my $sec;
my $name;

#-------------------------------------------------------------------------------
# Check if the directory fragments/ exists

if (-e "fragments"){
	if (!-d "fragments"){
		die ("There is a file called 'fragments' please rename it before proceeding\n");
	}
}
if (!-d "fragments"){
	system ("mkdir fragments");
} else {
	#print ("A directory named 'fragments' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
		system ("rm fragments/*");
	#}
}

#-------------------------------------------------------------------------------
# Read the alignment

open (MIA, "$file") or die ("No puedo abrir $file\n");
while (my $linea = <MIA>){
	chomp ($linea);
	if ($linea =~ />/ && $r == 1){
    $secs{$name} = $sec;
		$r = 0;
		$name = $sec = ();
	}
	if ($linea !~ />/ && $r == 1){
		$sec = $sec.$linea;
	}
	if ($linea =~ />/ && $r == 0){
    $linea =~ s/>//;
		$name = $linea;
    push (@head, $name);
		$r = 1;
	}
}
if ($r == 1){
  $secs{$name} = $sec;
	$r = 0;
	$name = $sec = ();
}
close (MIA);

#-------------------------------------------------------------------------------
# Fragment the multiple sequence alignment

my $asize = length($secs{$head[0]});
my $Nfrag = 0;

for (my $i = 0; $i <= $asize; $i = $i + $steps){
	my $col1 = $i + 1;
	##my $col2 = $i + $steps;
	my $col2 = $i + $wsize; # nuevo
	my $outfile = 'fragment_'.$col1.'_to_'.$col2.'.fasta';
	if ($i + $wsize <= $asize) {
		$Nfrag++;
		open (ROB, ">fragments/$outfile") or die ("Can't open fragments/$outfile\n");
		print ("\n-----\n");
		for (my $j = 0; $j <= $#head; $j++){
			#print (">$head[$j] (first sequence in alignment)\n") if ($j == 0);
			my $subs = substr($secs{$head[$j]}, $i, $wsize);
			#print ("$subs\n") if ($j == 0);
			# Print the fragments to the file
			print ROB (">$head[$j]\n");
			print ROB ("$subs\n");
		}
		# Borrar este for y el print
		#for (my $j = 0; $j < $wsize; $j++){
		#	print ("-");
		#}
		my $tmp = $i + $wsize;
		print ("fragment number: $Nfrag; ($i + $wsize == $tmp) <= $asize\n");
		print ("outfile: fragments/$outfile\n");
	} else {
		# Borrar esta parte
		#for (my $j = 0; $j <= $#head; $j++){
		#	print ("\n>$head[$j]\n") if ($j == 0);
		#	my $subs = substr($secs{$head[$j]}, $i, $wsize);
		#	print ("$subs\n") if ($j == 0);
		#}
		#for (my $j = 0; $j < $wsize; $j++){
		#	print ("-");
		#}
		#print ("\n$asize:$i\n");
		#print ("STOP\n");
	}
	close (ROB);
}

#-------------------------------------------------------------------------------
# Make statistics

$Nfrag = 0;
my $extra = 0;

for (my $i = 0; $i <= $asize; $i = $i + $steps){
	if ($i + $wsize <= $asize) {
		$Nfrag++;
	} else {
		$extra = length(substr($secs{$head[0]}, $i, $wsize));
	}
}

print ("\nStatistics:\n");
print ("Alignment length....: $asize\n");
print ("Windows size........: $wsize\n");
print ("Sliding positions...: $steps\n");
print ("Number of fragments.: $Nfrag\n");
print ("Columns not analyzed: $extra\n");
print ("\n");
