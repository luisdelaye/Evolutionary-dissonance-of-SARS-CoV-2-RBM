#!/usr/bin/perl

# This script format fasta files and infer marginal likelihoods with MrBayes.

# use
# perl goMrBayesss_local_threads.pl start end id
# out:
#  marginalL/
#   fragment_X_to_Y.ss.ckp
#   fragment_X_to_Y.ss.ckp~
#   fragment_X_to_Y.ss.mcmc
#   fragment_X_to_Y.ss.out.txt
#   fragment_X_to_Y.ss.run1.p
#   fragment_X_to_Y.ss.run1.t
#   fragment_X_to_Y.ss.run2.p
#   fragment_X_to_Y.ss.run2.t
#   fragment_X_to_Y.ss.ss

# Where:
# start indicates the first file to analyze in the list list_for_MrBayes.txt
# end indicates the last file to analyze in the list list_for_MrBayes.txt
# id sets the id number of the analysis

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

my $start = $ARGV[0];
my $end   = $ARGV[1];
my $id    = $ARGV[2];

# Con el siguiente codigo ordeno los archivos con el codigo de summary_H.pl
my @ffiles;
open (MIA, "list_for_MrBayes_ss.txt") or die ("Can't open list_for_MrBayes_ss.txt\n");
while (my $linea = <MIA>){
  chomp ($linea);
  push (@ffiles, $linea);
}
close (MIA);

if ($end -1 > $#ffiles){
  die ("Choose a number smaller than $#ffiles for the last file to analyze\n");
}

#my @files = sort {$a<=>$b} glob("phylogenies/fragment_*.nex");

#my %ffiles;

#foreach my $file (@files){
#  if ($file =~ /(phylogenies\/fragment_\d+_to_\d+\.nex)/){
#    my $f = $1;
#    $file =~ /fragment_(\d+)_to_\d+.nex/;
#    $ffiles{$f} = $1;
#  }
#}

#my @ffiles = reverse ordena (\%ffiles);
#foreach my $file (@ffiles){
#  print ("$file\n");
#}

#-------------------------------------------------------------------------------
# Check if the directory marginalL/ exists

if (-e "marginalL"){
	if (!-d "marginalL"){
		die ("There is a file called 'marginalL' please rename it before proceeding\n");
	}
}
if (!-d "marginalL"){
	system ("mkdir marginalL");
} else {
	#print ("A directory named 'marginalL' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
  #  system ("rm marginalL/*");
	#}
}

my $processed = 'processed_goMrBayesss_local_threads_'.$id.'.txt';
if (-e "$processed"){
  print ("A file named '$processed' already exists, if you continue");
	print ("its contents will be deleted: [y|n]\n");
	my $answer = <STDIN>;
	if ($answer =~ /n/i){
		die;
	} else {
    system ("rm $processed");
  }
}

#-------------------------------------------------------------------------------
# Make the phylogenetic analyses

print ("Processing files:\n");
for (my $i = ($start -1); $i <= ($end -1); $i++){
  margen($ffiles[$i], $id, $processed);
}

#-------------------------------------------------------------------------------
# Subrutines

sub margen {

  my $file  = $_[0];
  my $id    = $_[1];
  my $pfile = $_[2];

  print ("$file\n");
  # $file: phylogenies/fragment_X_to_Y.nex
  #--------------------------------
  # Prepare infile and run MrBayes
  my $infile = $file;
  $infile =~ s/phylogenies\///;

  #print ("\$infile: $infile\n");
  # $infile: fragment_X_to_Y.nex
  my $outss = $file;
  $outss =~ s/phylogenies\///;
  $outss =~ s/.nex//;
  $outss = $outss.'.ss';
  #print ("\$outss: $outss\n");
  # $outss: fragment_X_to_Y.ss
  my $infile_ss = 'infile_ss'.$id;
  open (ROB, ">$infile_ss") or die ("Can't open $infile_ss\n");
  print ROB ("exe phylogenies/$infile;\n");
  print ROB ("lset  nst=6 rates=invgamma ngammacat=4;\n");
  print ROB ("prset brlenspr=unconstrained:GammaDir(1.0, 0.01, 1.0, 1.0) shapepr=exp(1.0) statefreqpr=dirichlet(1.0,1.0,1.0,1.0) revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0);\n");
  print ROB ("ss ngen=1000000 samplefreq=1000 printfreq=10000 filename=marginalL/$outss;\n");
  close (ROB);
  $outss =~ s/.ss/.ss.out.txt/;
  system ("mb $infile_ss > marginalL/$outss");
  open (ROB, ">>$pfile") or die ("Can't open $pfile\n");
  print ROB ("marginalL/$outss\n");
  close (ROB);
  #print ("presiona enter para continuar\n");
  #my $pausa = <STDIN>;
  system ("rm $infile_ss");
}
