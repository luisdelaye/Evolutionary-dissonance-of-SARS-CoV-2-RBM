#!/usr/bin/perl

# This script concatenate sequence fasta files for ss analysis in MrBayes.

# use
# perl concatena.pl

# out: concatena2/fragments_X_to_Y_AND_Z_to_W.fasta ...

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

my $r = 0;
my $sec;
my $name;

#-------------------------------------------------------------------------------
# Check if the directory fragments/ exists

if (-e "concatena2"){
	if (!-d "concatena2"){
		die ("There is a file called 'concatena2' please rename it before proceeding\n");
	}
}
if (!-d "concatena2"){
	system ("mkdir concatena2");
} else {
	#print ("A directory named 'concatena2' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
		system ("rm concatena2/*");
	#}
}

#-------------------------------------------------------------------------------
# Read the files to concatenate

my @files = sort {$a<=>$b} glob("fragments/*.fasta");

my %ffiles;

foreach my $file (@files){
  if ($file =~ /(fragments[_\d]*\/fragment_\d+_to_\d+\.fasta)/){
    my $f = $1;
    $file =~ /fragment_(\d+)_to_\d+\.fasta/;
    $ffiles{$f} = $1;
  }
}

my @ffiles = reverse ordena (\%ffiles);

#for (my $i = 0; $i <= $#ffiles; $i++){
#	print ("$ffiles[$i]\n");
#}

#-------------------------------------------------------------------------------
# Concatenate files

for (my $i = 0; $i <= $#ffiles; $i++){
	print ("$ffiles[$i]\n");
	my $X = ();
	my $Y = ();
	if ($ffiles[$i] =~ /fragment_([\d]+)_to_([\d]+).fasta/){
		$X = $1;
		$Y = $2;
	} else {
		die ("pattern 1 failed: $ffiles[$i]\n");
	}
	for (my $j = $i+1; $j <= $#ffiles; $j++){
		print ("\t$ffiles[$j]\n");
		my $Z = ();
		my $W = ();
		if ($ffiles[$j] =~ /fragment_([\d]+)_to_([\d]+).fasta/){
			$Z = $1;
			$W = $2;
			my $outfile = 'fragments_'.$X.'_to_'.$Y.'_AND_'.$Z.'_to_'.$W;
			print ("\t->$outfile\n");
			concatena ($ffiles[$i], $ffiles[$j], $outfile);
		} else {
			die ("pattern 2 failed: $ffiles[$j]\n");
		}
	}
}


#-------------------------------------------------------------------------------
# Subrutines

sub concatena {
	my $file1   = $_[0];
	my $file2   = $_[1];
	my $outfile = 'concatena2/'.$_[2].'.fasta';
	my $X = ();
	my $Z = ();
	if ($_[2] =~ /fragments_(\d+)_to_\d+_AND_(\d+)_to_\d+/){
		$X = $1;
		$Z = $2;
	} else {
		die ("pattern inside sub failed: $outfile\n");
	}
	my %secs1;
	my %secs2;

	my $r = 0;
	my $sec;
	my $name;
	my @head;

	open (MIA, "$file1") or die ("No puedo abrir $file1\n");
	while (my $linea = <MIA>){
		chomp ($linea);
		if ($linea =~ />/ && $r == 1){
	    $secs1{$name} = $sec;
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
	  $secs1{$name} = $sec;
		$r = 0;
		$name = $sec = ();
	}
	close (MIA);

	open (MIA, "$file2") or die ("No puedo abrir $file2\n");
	while (my $linea = <MIA>){
		chomp ($linea);
		if ($linea =~ />/ && $r == 1){
	    $secs2{$name} = $sec;
			$r = 0;
			$name = $sec = ();
		}
		if ($linea !~ />/ && $r == 1){
			$sec = $sec.$linea;
		}
		if ($linea =~ />/ && $r == 0){
	    $linea =~ s/>//;
			$name = $linea;
			$r = 1;
		}
	}
	if ($r == 1){
	  $secs2{$name} = $sec;
		$r = 0;
		$name = $sec = ();
	}
	close (MIA);

	open (ROB, ">$outfile") or die ("Can't open $outfile\n");
	for (my $i = 0; $i <= $#head; $i++){
		if (exists $secs2{$head[$i]}){
			if ($X + 500 == $Z){
				my $con = substr($secs1{$head[$i]},0,500).$secs2{$head[$i]};
				print ROB (">$head[$i]\n$con\n");
			} elsif ($X + 1000 == $Z){
				my $con = substr($secs1{$head[$i]},0,1000).$secs2{$head[$i]};
				print ROB (">$head[$i]\n$con\n");
			} elsif ($X + 1500 == $Z){
				my $con = substr($secs1{$head[$i]},0,1500).$secs2{$head[$i]};
				print ROB (">$head[$i]\n$con\n");
			} else {
				my $con = $secs1{$head[$i]}.$secs2{$head[$i]};
				print ROB (">$head[$i]\n$con\n");
			}
		} else {
			die ("The sequence $head[$i] is lacking in $file2\n");
		}
	}
	close (ROB);
}

sub ordena {

	# Recibe una referencia a un hash
	my $hashref = $_[0];
	my @ordenados = ();

	my @a = keys (%{$hashref});
	my $N_de_elementos = $#a +1;

	my $suma = 0;
	until ($suma == $N_de_elementos){
		for (my $i = 0; $i <= $#a; $i++){
			#print ("\nOrdenando: $a[$i] ($hash{$a[$i]})\n");
			my $elmayor = 'si';
			for (my $j = 0; $j <= $#a; $j++){
				if ($i != $j){
					#print ("\t$a[$j] ($hash{$a[$j]})\n");
					if ($hashref->{$a[$i]} < $hashref->{$a[$j]}){
						$elmayor = 'no';
					}
				}
			}
			#print ("es el mayor: $elmayor\n");
			if ($elmayor eq 'si'){
				push (@ordenados, $a[$i]);
				my $salta = $i;
				my $ordenado = splice (@a, $salta, 1);
				#print ("\nordenado: $ordenado\n");
				$suma++;
			}
			#$pausa = <STDIN>;
		}

	}

	return (@ordenados);
}
