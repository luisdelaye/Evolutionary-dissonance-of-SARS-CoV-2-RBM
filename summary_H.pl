#!/usr/bin/perl

# This script summarize the information content analyzed by GALAX.

# use
# perl summary_H.pl
# out: statistics_1.tsv, statistics_2.tsv


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

my $skip = $ARGV[0];

#------------------------------------------------------------------------------
# Read files

my @filesF = sort {$a<=>$b} glob("information/fragment_*.output-galax.txt");

my %filesF;

foreach my $file (@filesF){
  if ($file =~ /(fragment_\d+_to_\d+).output-galax.txt/){
    my $f = $1;
    $file =~ /fragment_(\d+)_to_\d+.output-galax.txt/;
    $filesF{$f} = $1;
  }
}

my @filesF = reverse ordena (\%filesF);

# Create the keys

my @filesFs1;
my %filesFs1;
my %temporal;
my @temporal = sort {$a<=>$b} glob("information/fragments_*.run1.output-galax.txt");
my $menor  = 1000000;
my $mayor  = 0;
my $step   = 1000000;
my $window = 0;
foreach my $file (@temporal){
  if ($file =~ /fragments_(\d+)_to_(\d+)_AND_(\d+)_to_(\d+).run1.output-galax.txt/){
    if ($1 < $menor){
      $menor = $1;
    }
    if ($4 > $mayor){
      $mayor = $4;
    }
    if ($step > ($3 - $1)){
      $step = $3 - $1;
    }
    $window = $2 - $1 + 1;
  }
}
#print ("numero menor: $menor\n");
#print ("numero mayor: $mayor\n");
#print ("intervalo...: $step\n");
#print ("ventana.....: $window\n");
for (my $i = $menor; $i <= $mayor - $window; $i = $i + $step){
  for (my $j = $i+$step; $j <= $mayor - $window +1; $j = $j + $step){
    my $f = 'fragments_'.$i.'_to_'.($i+$window-1).'_AND_'.$j.'_to_'.($j+$window-1);
    #print ("$f\n");
    push (@filesFs1, $f);
    my $k = 'fragments_'.$i.'_to_'.($i+$window-1).'_AND_'.$j.'_to_'.($j+$window-1);
    my $l = $i.'_to_'.($i+$window-1).'_AND_'.$j.'_to_'.($j+$window-1);
    $filesFs1{$k} = $l;
    $temporal{$f} = 0;
  }
}
# Check whether all names correspond to a real file
foreach my $f (@filesFs1){
  my $fi = 'information/'.$f.'.run1.output-galax.txt';
  if (!-e $fi){
    die ("file does not exists in folder information/: $fi\n");
  }
}
# Check whether all files have a name in @filesFs1
foreach my $f (@temporal){
  $f =~ s/information\///;
  $f =~ s/.run1.output-galax.txt//;
  if (!exists $temporal{$f}){
    die ("file does not exists in \@temporal: $f\n");
  }
}

#my @filesFs1 = sort {$a<=>$b} glob("information/fragments_*.run1.output-galax.txt");
#my %filesFs1;
#foreach my $file (@filesFs1){
#  if ($file =~ /(fragments_\d+_to_\d+_AND_\d+_to_\d+).run1.output-galax.txt/){
#    my $f = $1;
#    $file =~ /fragments_(\d+_to_\d+_AND_\d+_to_\d+).run1.output-galax.txt/;
#    $filesFs1{$f} = $1;
#  }
#}
#my @filesFs1 = reverse ordena (\%filesFs1);
#die ("bien!\n");

#-------------------------------------------------------------------------------
# Summarize the information whitin each fragment

my %whitin;

print ("\nInformation within fragments\n");
open (ROB, ">statistics_1.tsv") or die ("Can't open file statistics_1.tsv\n");
print ROB ("fragment\tunique\tcoverage\tH_pr\tH_po\tI\tIpct\tD\tDpct\n");
foreach my $file (@filesF){
  #print ("$file\n");
  my $roww = ();
  my $outputgalax = $file.'.output-galax.txt';
  open (MIA, "information/$outputgalax") or die ("Can't open file information/$outputgalax\n");
  while (my $linea = <MIA>){
    chomp ($linea);
    if ($linea =~ /^\s+merged\s/){
      #print ("\n$file\t$filesF{$file}\t");
      print ROB ("$filesF{$file}\t");
      $roww = "$filesF{$file}\t$filesF{$file}\t";
      #print ("$linea\n");
      my @a = split (/\s+/, $linea);
      print ("@a\n");
      for (my $i = 2; $i <= $#a; $i++){
        if ($i < $#a){
          print ROB ("$a[$i]\t");
          $roww = $roww."$a[$i]\t";
        } else {
          print ROB ("$a[$i]\n");
          $roww = $roww."$a[$i]\n";
        }
      }
    }
  }
  $whitin{$filesF{$file}} =$roww;
  close (MIA);
}
close (ROB);

#-------------------------------------------------------------------------------
# Summarize the information between fragments

my @rows;
my $row;

foreach my $file (@filesFs1){
  my $outputgalax = $file.'.run1.output-galax.txt';
  open (MIA, "information/$outputgalax") or die ("Can't open file information/$outputgalax\n");
  while (my $linea = <MIA>){
    chomp ($linea);
    if ($linea =~ /^\s+merged\s/){
      $filesFs1{$file} =~ /(\d+)_to_\d+_AND_(\d+)_to_\d+/;
      my $frag1 = $1;
      my $frag2 = $2;
      $row = "$frag1\t$frag2\t";
      my @a = split (/\s+/, $linea);
      for (my $i = 2; $i <= $#a; $i++){
        if ($i < $#a){
          $row = $row."$a[$i]\t";
        } else {
          $row = $row."$a[$i]\n";
        }
      }
      push (@rows, $row);
    }
  }
  close (MIA);
}

my $r = 0;
my $s = 0;
my %hashi;
print ("\nInformation between fragments\n");
open (ROB, ">statistics_2.tsv") or die ("Can't open file statistics_2.tsv\n");
print ROB ("frgmt_1\tfrgmt_2\tunique\tcoverage\tH_pr\tH_po\tI\tIpct\tD\tDpct\n");
print ("frgmt_1\tfrgmt_2\tunique\tcoverage\tH_pr\tH_po\tI\tIpct\tD\tDpct\n");
for (my $i = 0; $i <= $#rows; $i++){
  my @a = split (/\t/, $rows[$i]);
  # ititial value
  if ($i == 0){
    $r = $a[0];
  }
  # Order rows
  if ($r != $a[0] && $s == 1){
    my @b = sort {$a<=>$b} keys (\%hashi);
    print ("$whitin{$r}");
    print ROB ("$whitin{$r}");
    for (my $j = 0; $j <= $#b; $j++){
      print ROB ("$hashi{$b[$j]}");
      print ("$hashi{$b[$j]}");
    }
    %hashi = ();
    $r = $a[0];
    $s = 0;
  }
  # Collect rows
  if ($r == $a[0] && $s == 1){
    $hashi{$a[1]} = $rows[$i];
  }
  # Start collecting rows
  if ($r == $a[0] && $s == 0){
    $hashi{$a[1]} = $rows[$i];
    $s = 1;
  }
}
if ($s == 1){
  my @b = sort {$a<=>$b} keys (\%hashi);
  print ("$whitin{$r}");
  print ROB ("$whitin{$r}");
  for (my $j = 0; $j <= $#b; $j++){
    print ROB ("$hashi{$b[$j]}");
    print ("$hashi{$b[$j]}");
  }
  %hashi = ();
  $s = 0;
  my @a = split (/\t/, $rows[-1]);
  print ROB ("$whitin{$a[1]}");
  print ("$whitin{$a[1]}");
}
close (ROB);

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
