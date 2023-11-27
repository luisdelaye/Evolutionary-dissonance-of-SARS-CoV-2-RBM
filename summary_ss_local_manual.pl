#!/usr/bin/perl

# This script summarize the marginalized likelihoods.

# use
# perl summary_ss_local.pl
# out: statistics_ss.tsv


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

use strict;

my $winsi = $ARGV[0]; # longitud de la ventana
my $slide = $ARGV[1]; # numero de columnas que se desplaza la ventana
my $ncols = $ARGV[2]; # numero de columnas en la alineacion multiple

my $n = 1;
my $nfiles = 0;
my %hash;

#---------------------------------------------------------------------------
# Extraigo los valores lnL de los archivos ss

# Este ciclo cicla los archivos ss de cada fragmento individual
for (my $i = 1; ($i + $winsi -1) <= $ncols; $i = $i + $slide){
  # $file es el archivo individual
  my $step = $i + $winsi -1;
  my $file = 'fragment_'.$i.'_to_'.$step.'.ss.out.txt';
  my $key = $i.'_to_'.$step;
  if (-e "marginalL/$file"){
    print ("$file\t");
    $nfiles++;
    my $lnL = 'NA';
    # Extraigo el lnL del archivo
    open (MIA, "marginalL/$file") or die ("Can't open marginalL/$file\n");
    while (my $linea = <MIA>){
      chomp ($linea);
      if ($linea =~ /Mean:\s+([-\d\.]+)/){
        $lnL = $1;
      }
    }
    close (MIA);
    print ("$lnL\n");
    # Guardo el lnL con la llave del archivo individual
    $hash{$key} = $lnL;
    #---------------------------------------------------------------------------
    # Este ciclo cicla los archivos de los concatenados
    my $nfiles2 = 0;
    for (my $j = ($i + $slide -1); ($j + $winsi -1) <= $ncols; $j = $j + $slide){
      my $step2 = $j + $winsi;
      # Todo este desmadre es porque la numeracion final de los archivos no esta ordenada
      # Construyo un nombre parcial del archivo
      my $file2 = 'fragments_'.$i.'_to_'.$step.'_AND_'.($j+1).'_to_'.$step2.'.ss.out.txt';
      my $key2 = $i.'_to_'.$step.'_AND_'.($j+1).'_to_'.$step2;
      #my $nl = 0;
      my $file2Ok = 'NA';
      # Recupero el nombre completo del archivo
      #system ("ls marginalLcon/$file2\* > borrame");
      #open (MIA, "borrame") or die ("Can't open borrame\n");
      #while (my $linea2 = <MIA>){
      #  chomp ($linea2);
      #  if ($linea2 =~ /marginalLcon\/(fragments_\d+_to_\d+_AND_\d+_to_\d+.ss.out.txt)/){
      #    $file2Ok = $1;
      #  }
      #  $nl++;
      #}
      #close (MIA);
      #die ("the borrame file has more than one file\n") if ($nl > 1);
      # Pregunto si el archivo existe
      #if (-e "marginalLcon/$file2Ok"){
      if (-e "marginalLcon/$file2"){
        print ("$file2\t");
        $nfiles2++;
        my $lnL2 = 'NA';
        # Extraigo el lnL del archivo
        open (MIA, "marginalLcon/$file2") or die ("Can't open marginalLcon/$file2\n");
        while (my $linea = <MIA>){
          chomp ($linea);
          if ($linea =~ /Mean:\s+([-\d\.]+)/){
            $lnL2 = $1;
          }
        }
        close (MIA);
        print ("$lnL2\n");
        $hash{$key2} = $lnL2;
      } else {
        # Si no existe el archivo
        print ("$file2\t$file2Ok\n");
        $hash{$key2} = $file2Ok; # NA == $file2Ok
      }
    }
    print ("$nfiles2\n");
    #my $pausa = <STDIN>;
  } else {
    die ("$file <- no existe\n");
  }
  $n++;
}
##print ("numero de archivos analizados: $nfiles\n");

#-------------------------------------------------------------------------------
# Ahora ordeno e imprimo todos los valores

print ("-----\n");
open (ROB, ">statistics_ss.tsv") or die ("Can't open statistics_ss.tsv\n");
print ROB ("frgmt_1s\tlnL1\tfrgmt_2s\tlnL2\tfrgmt_1y2\tlnM0\n");
for (my $i = 1; ($i + $winsi -1) <= $ncols; $i = $i + $slide){
  my $step = $i + $winsi -1;
  my $key1  = $i.'_to_'.$step;
  #my $file = 'outfile_infile_'.$i.'_to_'.$step.'_'.$n.'.ss.log';
  my $choco = $i.'_AND_'.$i;
  print ROB ("$i\t$hash{$key1}\t$i\t$hash{$key1}\t$choco\tNA\n");
  for (my $j = ($i + $slide -1); ($j + $winsi -1) <= $ncols; $j = $j + $slide){
    my $step2 = $j + $winsi;
    #my $file2 = 'outfile_infile_'.$i.'_to_'.$step.'_AND_'.($j+1).'_to_'.$step2.'_';
    my $key1y2 = $i.'_to_'.$step.'_AND_'.($j+1).'_to_'.$step2;
    my $key1y2p = $i.'_AND_'.($j+1);
    my $j_1 = $j+1;
    my $key2 = $j_1.'_to_'.$step2;
    print ROB ("$i\t$hash{$key1}\t$j_1\t$hash{$key2}\t$key1y2p\t$hash{$key1y2}\n");
  }
}
close (ROB);
