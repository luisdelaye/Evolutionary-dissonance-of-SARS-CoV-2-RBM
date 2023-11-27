#!/usr/bin/perl

# This script estimates the information content by using GALAX.

# use
# perl goGALAX_local.pl N
# out:
#  information/
#   fragment_X_to_Y.output-galax-details.txt
#   fragment_X_to_Y.output-galax-merged.txt
#   fragment_X_to_Y.output-galax.txt
#   fragments_X_to_Y_AND_Z_to_W.run1.output-galax-details.txt
#   fragments_X_to_Y_AND_Z_to_W.run1.output-galax-merged.txt
#   fragments_X_to_Y_AND_Z_to_W.run1.output-galax.txt
#   fragments_X_to_Y_AND_Z_to_W.run2.output-galax-details.txt
#   fragments_X_to_Y_AND_Z_to_W.run2.output-galax-merged.txt
#   fragments_X_to_Y_AND_Z_to_W.run2.output-galax.txt

# Where N refers to the number of trees to skip (burnin)

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

my $skip = $ARGV[0];

my @files = sort {$a<=>$b} glob("phylogenies/fragment_*.t");

my %tfiles;

foreach my $file (@files){
  if ($file =~ /(fragment_\d+_to_\d+\.mcmc)\.run[12].t/){
    my $f = $1;
    $file =~ /fragment_(\d+)_to_\d+\.mcmc\.run[12].t/;
    $tfiles{$f} = $1;
  }
}

my @tfiles = reverse ordena (\%tfiles);

#foreach my $f (@tfiles){
#  print ("$f\n");
#}
#die ("bien!\n");

#-------------------------------------------------------------------------------
# Check if the directory information/ exists

if (-e "information"){
	if (!-d "information"){
		die ("There is a file called 'information' please rename it before proceeding\n");
	}
}
if (!-d "information"){
	system ("mkdir information");
} else {
	#print ("A directory named 'information' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
		system ("for f in information/* ; do rm \"\$f\"; done");
	#}
}

#-------------------------------------------------------------------------------
# Calculate information content and phylogenetic disonance

for(my $i = 0; $i <= $#tfiles; $i++){
  my $tfile = $tfiles[$i];
  print ("\n-----\n");
  print ("$tfile\n");
  my $f1 = 'phylogenies/'.$tfile.'.run1.t';
  my $f2 = 'phylogenies/'.$tfile.'.run2.t';
  open (ROB, ">infile") or die ("Can't open infile\n");
  print ROB ("$f1\n");
  print ROB ("$f2\n");
  close (ROB);
  system ("galax1 --listfile infile --skip $skip --outgroup 5 --details");
  my $outfile = $tfile;
  $outfile =~ s/\.mcmc//;
  $outfile = 'information/'.$outfile.'.output-galax.txt';
  system ("mv output-galax.txt $outfile");
  $outfile = $tfile;
  $outfile =~ s/\.mcmc//;
  $outfile = 'information/'.$outfile.'.output-galax-details.txt';
  system ("mv output-galax-details.txt $outfile");
  $outfile = $tfile;
  $outfile =~ s/\.mcmc//;
  $outfile = 'information/'.$outfile.'.output-galax-merged.txt';
  system ("mv output-galax-merged.tre $outfile");
  for(my $j = $i+1; $j <= $#tfiles; $j++){
    my $tfile2 = $tfiles[$j];
    # run1.t
    my $fs1 = 'phylogenies/'.$tfile2.'.run1.t';
    open (ROB, ">infile") or die ("Can't open infile\n");
    print ROB ("$f1\n");
    print ROB ("$fs1\n");
    close (ROB);
    system ("galax1 --listfile infile --skip $skip --outgroup 5 --details");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run1.output-galax.txt';
    system ("mv output-galax.txt $outfile");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run1.output-galax-details.txt';
    system ("mv output-galax-details.txt $outfile");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run1.output-galax-merged.txt';
    system ("mv output-galax-merged.tre $outfile");
    # run2.t
    my $fs2 = 'phylogenies/'.$tfile2.'.run2.t';
    open (ROB, ">infile") or die ("Can't open infile\n");
    print ROB ("$f2\n");
    print ROB ("$fs2\n");
    close (ROB);
    system ("galax1 --listfile infile --skip $skip --outgroup 5 --details");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run2.output-galax.txt';
    system ("mv output-galax.txt $outfile");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run2.output-galax-details.txt';
    system ("mv output-galax-details.txt $outfile");
    $outfile = $tfile;
    #$outfile =~ s/\.mb$/.galax/;
    $outfile =~ s/fragment_//;
    $tfile =~ /_(\d+_to_\d+)/;
    my $iindex = $1;
    $tfile2 =~ /_(\d+_to_\d+)/;
    my $jindex = $1;
    $outfile = 'information/fragments_'.$iindex.'_AND_'.$jindex.'.run2.output-galax-merged.txt';
    system ("mv output-galax-merged.tre $outfile");
  }
  #print ("Press enter to run the next phylogenetic analysis\n");
  #my $pausa = <STDIN>;
}
system ("rm infile");

#-------------------------------------------------------------------------------
# Subrutines

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
