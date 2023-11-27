#!/usr/bin/perl

# This script format fasta files and infer phylogenies with MrBayes.

# use
# perl goMrBayes_local_threads.pl start end id
# out:
#  phylogenies/
#   fragment_X_to_Y.mcmc.ckp
#   fragment_X_to_Y.mcmc.ckp~
#   fragment_X_to_Y.mcmc.con.tre
#   fragment_X_to_Y.mcmc.lstat
#   fragment_X_to_Y.mcmc.mcmc
#   fragment_X_to_Y.mcmc.parts
#   fragment_X_to_Y.mcmc.pstat
#   fragment_X_to_Y.mcmc.run1.p
#   fragment_X_to_Y.mcmc.run1.t
#   fragment_X_to_Y.mcmc.run2.p
#   fragment_X_to_Y.mcmc.run2.t
#   fragment_X_to_Y.mcmc.trprobs
#   fragment_X_to_Y.mcmc.tstat
#   fragment_X_to_Y.mcmc.vstat
#   fragment_X_to_Y.mcmc.out.txt
#   fragment_X_to_Y.nex

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
open (MIA, "list_for_MrBayes_mcmc.txt") or die ("Can't open list_for_MrBayes_mcmc.txt\n");
while (my $linea = <MIA>){
  chomp ($linea);
  push (@ffiles, $linea);
}
close (MIA);

if ($end -1 > $#ffiles){
  die ("Choose a number smaller than $#ffiles for the last file to analyze\n");
}

#-------------------------------------------------------------------------------
# Check if the directory phylogenies/ exists

if (-e "phylogenies"){
	if (!-d "phylogenies"){
		die ("There is a file called 'phylogenies' please rename it before proceeding\n");
	}
}
if (!-d "phylogenies"){
	system ("mkdir phylogenies");
} else {
	#print ("A directory named 'phylogenies' already exists, if you continue");
	#print ("its contents will be deleted: [y|n]\n");
	#my $answer = <STDIN>;
	#if ($answer =~ /n/i){
	#	die;
	#} else {
	#	system ("rm phylogenies/*");
	#}
}

my $processed = 'processed_goMrBayes_local_threads_'.$id.'.txt';
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
  # $file: fragments/fragment_X_to_Y.fasta
  #--------------------------------
  # Format the fasta file to nexus
  my $outfile = $file;
  $outfile =~ s/fragments/phylogenies/;
  $outfile =~ s/fasta/nex/;
  #print ("\$outfile: $outfile\n");
  # $outfile: phylogenies/fragment_X_to_Y.nex
  # Falta ver como hacer para que emita el formato correcto: DNA o proteina
  system ("readal -in $file -out $outfile -nexus");
  my $tmpfile = 'tmp_mcmc'.$id;
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
  # $outfile: phylogenies/fragment_X_to_Y.nex
  #-------------------------------
  # Prepare infile and run MrBayes
  my $outmb = $outfile;
  $outfile =~ s/phylogenies\///;
  #print ("\$outfile: $outfile\n");
  # $outfile: fragment_X_to_Y.nex
  $outmb =~ s/.nex/.mcmc/;
  #print ("\$outmb: $outmb\n");
  # $outmb: phylogenies/fragment_X_to_Y.mcmc
  my $infile_mcmc = 'infile_mcmc'.$id;
  open (ROB, ">$infile_mcmc") or die ("Can't open $infile_mcmc\n");
  print ROB ("exe phylogenies/$outfile;\n");
  print ROB ("lset  nst=6 rates=invgamma ngammacat=4;\n");
  print ROB ("prset brlenspr=unconstrained:GammaDir(1.0, 0.01, 1.0, 1.0) shapepr=exp(1.0) statefreqpr=dirichlet(1.0,1.0,1.0,1.0) revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0);\n");
  print ROB ("mcmc ngen=1000000 samplefreq=500 printfreq=10000 burninfrac=0.25 nchains=4 nruns=2 filename=$outmb;\n");
  print ROB ("sump filename=$outmb;\n");
  print ROB ("sumt filename=$outmb;\n");
  close (ROB);
  $outmb =~ s/.mcmc/.mcmc.out.txt/;
  #print ("\$outmb: $outmb\n");
  # $outmb: phylogenies/fragment_X_to_Y.mcmc.out.txt
  system ("mb $infile_mcmc > $outmb");
  open (ROB, ">>$pfile") or die ("Can't open $pfile\n");
  print ROB ("$outmb\n");
  close (ROB);
  #print ("presiona enter para continuar\n");
  #my $pausa = <STDIN>;
  system ("rm $infile_mcmc");
}
