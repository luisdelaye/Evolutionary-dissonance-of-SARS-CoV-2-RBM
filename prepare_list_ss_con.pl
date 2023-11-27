#!/usr/bin/perl

# This script prepare a list of files for goMrBayes_con_local_threads.pl.

# use
# perl prepare_list.pl
# out: list_for_MrBayes_ss_con.txt

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
my @temporal = sort {$a<=>$b} glob("concatena2/fragments_*.fasta");
my $menor  = 1000000;
my $mayor  = 0;
my $step   = 1000000;
my $window = 0;
foreach my $file (@temporal){
  if ($file =~ /fragments_(\d+)_to_(\d+)_AND_(\d+)_to_(\d+).fasta/){
    if ($1 < $menor){
      $menor = $1;
    }
    if ($4 > $mayor){
      $mayor = $4;
    }
    if ($step > ($3 - $1)){
      $step = $3 - $1;
      #$step = $1 - $menor;
    }
    $window = $2 - $1 + 1;
  }
}
#print ("numero menor: $menor\n");
#print ("numero mayor: $mayor\n");
#print ("intervalo...: $step\n");
#print ("ventana.....: $window\n");
# Create the names/index of the files
for (my $i = $menor; $i <= $mayor - $window; $i = $i + $step){
  for (my $j = $i+$step; $j <= $mayor - $window +1; $j = $j + $step){
    my $f = 'fragments_'.$i.'_to_'.($i+$window-1).'_AND_'.$j.'_to_'.($j+$window-1);
    #print ("$f\n");
    push (@filesFs1, $f);
    my $l = $i.'_to_'.($i+$window-1).'_AND_'.$j.'_to_'.($j+$window-1);
    $filesFs1{$f} = $l;
    $temporal{$f} = 0;
  }
}
# Check whether all names correspond to a real file
foreach my $f (@filesFs1){
  my $fi = 'concatena2/'.$f.'.fasta';
  if (!-e $fi){
    die ("file does not exists in folder concatena2/: $fi\n");
  }
}
# Check whether all files have a name in @filesFs1
foreach my $f (@temporal){
  $f =~ s/concatena2\///;
  $f =~ s/.fasta//;
  if (!exists $temporal{$f}){
    die ("file does not exists in \@temporal: $f\n");
  }
}
my $n = 0;
my @ffiles;
open (ROB, ">list_for_MrBayes_ss_con.txt") or die ("Can't open list_for_MrBayes_ss_con.txt\n");
for (my $i = 0; $i <= $#filesFs1; $i++){
  $ffiles[$i] = 'concatena2/'.$filesFs1[$i].'.fasta';
  print ROB ("$ffiles[$i]\n");
  $n++;
}
close (ROB);
print ("Number of files: $n\n");
