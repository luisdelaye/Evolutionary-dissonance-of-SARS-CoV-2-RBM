#!/usr/bin/perl

# This script prepare a list of files for goMrBayes_local_threads.pl.

# use
# perl prepare_list_mcmc.pl
# out: list_for_MrBayes_mcmc.txt

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

# Con el siguiente codigo ordeno los archivos con el codigo de summary_H.pl
my @filesFs1;
my %filesFs1;
my %temporal;
my @temporal;
my @temp = sort {$a<=>$b} glob("fragments/fragment_*.fasta");
my %hash0;
for (my $i = 0; $i <= $#temp; $i++){
  if ($temp[$i] =~ /fragment_(\d+)_to_(\d+).fasta/){
    $hash0{$1} = $temp[$i];
  }
}
@temp = sort {$a<=>$b} keys (%hash0);
my $menor  = 1000000;
my $mayor  = 0;
my $step   = 1000000;
my $window = 0;
my $v = 0;
foreach my $file (@temp){
  $v++;
  #print ("->$file\t$hash0{$file}\n");
  if ($hash0{$file} =~ /fragment_(\d+)_to_(\d+).fasta/){
    push (@temporal, $hash0{$file});
    #print ("->$v\t$file\n");
    if ($1 < $menor){
      $menor = $1;
    }
    if ($v == 2){
        #print ("-->\$v = $v\n");
        #print ("-->\$1 = $1\n");
        #print ("-->\$2 = $2\n");
        #print ("-->\$menor = $menor\n");
        #print ("-->\$mayor = $mayor\n");
        $step = $1 - $menor;
        #print ("-->\$step = $step\n");
    }
    if ($2 > $mayor){
      $mayor = $2;
    }
    $window = $2 - $1 + 1;
  }
}
print ("numero menor: $menor\n");
print ("numero mayor: $mayor\n");
print ("intervalo...: $step\n");
print ("ventana.....: $window\n");
# Create the names/index of the files
for (my $i = $menor; $i <= $mayor - $window +1; $i = $i + $step){
  my $f = 'fragment_'.$i.'_to_'.($i+$window-1);
  #print ("$f\n");
  push (@filesFs1, $f);
  my $l = 'fragment_'.$i.'_to_'.($i+$window-1);;
  $filesFs1{$f} = $l;
  $temporal{$f} = 0;
}
# Check whether all names correspond to a real file
foreach my $f (@filesFs1){
  my $fi = 'fragments/'.$f.'.fasta';
  if (!-e $fi){
    die ("file does not exists in folder fragments/: $fi\n");
  }
}
# Check whether all files have a name in @filesFs1
foreach my $f (@temporal){
  $f =~ s/fragments\///;
  $f =~ s/.fasta//;
  if (!exists $temporal{$f}){
    die ("file does not exists in \@temporal: $f\n");
  }
}
my $n = 0;
my @ffiles;
open (ROB, ">list_for_MrBayes_mcmc.txt") or die ("Can't open list_for_MrBayes_mcmc.txt\n");
for (my $i = 0; $i <= $#filesFs1; $i++){
  $ffiles[$i] = 'fragments/'.$filesFs1[$i].'.fasta';
  print ROB ("$ffiles[$i]\n");
  $n++;
}
close (ROB);
print ("Number of files: $n\n");
