#!/usr/bin/perl

# This script format fasta files and infer marginal likelihoods with MrBayes.

# use
# perl goMrBayesss_con_local_threads.pl start end id
# out:
#  marginalLcon/
#   fragments_X_to_Y_AND_Z_to_W.ss.ckp
#   fragments_X_to_Y_AND_Z_to_W.ss.ckp~
#   fragments_X_to_Y_AND_Z_to_W.ss.mcmc
#   fragments_X_to_Y_AND_Z_to_W.ss.out.txt
#   fragments_X_to_Y_AND_Z_to_W.ss.run1.p
#   fragments_X_to_Y_AND_Z_to_W.ss.run1.t
#   fragments_X_to_Y_AND_Z_to_W.ss.run2.p
#   fragments_X_to_Y_AND_Z_to_W.ss.run2.t
#   fragments_X_to_Y_AND_Z_to_W.ss.ss

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
open (MIA, "list_for_MrBayes_ss_con.txt") or die ("Can't open list_for_MrBayes_ss_con.txt\n");
while (my $linea = <MIA>){
  chomp ($linea);
  push (@ffiles, $linea);
}
close (MIA);

if ($end -1 > $#ffiles){
  die ("Choose a number smaller than $#ffiles for the last file to analyze\n");
}

#-------------------------------------------------------------------------------
# Check if the directory marginalLcon/ exists

if (-e "marginalLcon"){
	if (!-d "marginalLcon"){
		die ("There is a file called 'marginalLcon' please rename it before proceeding\n");
	}
}
if (!-d "marginalLcon"){
	system ("mkdir marginalLcon");
} else {
	#print ("A directory named 'marginalLcon' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
	#	system ("rm marginalLcon/*");
	#}
}

my $processed = 'processed_goMrBayesss_con_local_threads_'.$id.'.txt';
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
  # $file: concatena2/fragments_X_to_Y_AND_Z_to_W.fasta
  #--------------------------------
  # Format the fasta file to nexus
  my $outfile = $file;

  $outfile =~ s/fasta/nex/;
  #print ("\$outfile: $outfile\n");
  # $outfile: concatena2/fragments_X_to_Y_AND_Z_to_W.nex
  # Falta ver como hacer para que emita el formato correcto: DNA o proteina
  system ("readal -in $file -out $outfile -nexus");
  my $tmpfile = 'tmp_sscon'.$id;
  open (MIA, "$outfile") or die ("Can't open $outfile\n");
  open (ROB, ">$tmpfile") or die ("Can't open $tmpfile\n");
  while (my $linea = <MIA>){
    if ($linea =~ /FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-;/){
      print ROB ("FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-;");
    } else {
      print ROB ("$linea");
    }
  }
  close (ROB);
  close (MIA);
  system ("mv $tmpfile $outfile");
  #print ("\$outfile: $outfile\n");
  # $outfile: concatena2/fragments_X_to_Y_AND_Z_to_W.nex
  #--------------------------------
  # Prepare infile and run MrBayes
  my $infile = $file;
  $infile =~ s/concatena2\///;
  $infile =~ s/.fasta/.nex/;
  #print ("\$infile: $infile\n");
  # $infile: fragments_X_to_Y_AND_Z_to_W.nex
  my $outss = $file;
  $outss =~ s/concatena2\///;
  $outss =~ s/.fasta//;
  $outss = $outss.'.ss';
  #print ("\$outss: $outss\n");
  # $outss: fragments_X_to_Y_AND_Z_to_W.ss
  my $infile_sscon = 'infile_sscon'.$id;
  open (ROB, ">$infile_sscon") or die ("Can't open $infile_sscon\n");
  print ROB ("exe concatena2/$infile;\n");
  print ROB ("lset  nst=6 rates=invgamma ngammacat=4;\n");
  print ROB ("prset brlenspr=unconstrained:GammaDir(1.0, 0.01, 1.0, 1.0) shapepr=exp(1.0) statefreqpr=dirichlet(1.0,1.0,1.0,1.0) revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0);\n");
  print ROB ("ss ngen=1000000 samplefreq=1000 printfreq=10000 filename=marginalLcon/$outss;\n");
  close (ROB);
  $outss =~ s/.ss/.ss.out.txt/;
  system ("mb $infile_sscon > marginalLcon/$outss");
  open (ROB, ">>$pfile") or die ("Can't open $pfile\n");
  print ROB ("marginalLcon/$outss\n");
  close (ROB);
  #print ("presiona enter para continuar\n");
  #my $pausa = <STDIN>;
  system ("rm $infile_sscon");
}
