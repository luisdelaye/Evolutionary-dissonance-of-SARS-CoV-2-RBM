# Evolutionary-dissonance-of-SARS-CoV-2-RBM
Evolutionary analysis of phylogenetic dissonance of the Recognition Binding Motif (RBM) from SARS-CoV-2

Bitácora spike 14: 250 200

1)  En esta carpeta voy a hacer el análisis de las 14 secuencias del gen de la
    proteína spike que utilicé para la clase de evolución molecular sobre
    inferencia Bayesiana que di junto con Lizbeth en mayo del 2023. Los detalles
    son los siguientes:

    De acuerdo a un análisis que hice, el loop variable está más o menos del
    nucleótido 22814 al 23062 en el genoma MN908947. Y la secuencia es esta:

      attgctgattataattataaattaccagatgattttacaggctgcgttatagcttggaattctaacaatcttgattctaaggttggtggtaattataattacctgtatagattgtttaggaagtctaatctcaaaccttttgagagagatatttcaactgaaatctatcaggccggtagcacaccttgtaatggtgttgaaggttttaattgttactttcctttacaatcatatggtttccaacccact

                              variable            loop en 8GZ5 A                                 variable
      IADYNYKLPDDFTGCVIAWNS->>NNLDSKVGGNYNYLYR->LFRKSNLKPFERDISTEI<-YQAGSTPCNGVEGFNCYFPLQSYGFQPT NGVGY<<-


                     MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFR
                     SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR
                     GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY
                     SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQ
                     GFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFL
                     LKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITN
                     LCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCF
                     TNVYADSFVIRGDEVRQIAPGQTGK IADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN
                     YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPT NGVGYQPY
                     RVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFG
                     RDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAI
                     HADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPR
                     RARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTM
                     YICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFG
                     GFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFN
                     GLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQN
                     VLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA
                     ISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMS
                     ECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAH
                     FPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD
                     SFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELG
                     KYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSE
                     PVLKGVKLHYT

    De acuerdo a ConSurf la región más variable la enmarco entre ->> y <<-
    Esta es justo la región del RBD (Receptor Binding Domain) que contiene los
    hACE2 contact residues (todos menos el primero). Se puede ver en el
    material suplementario de Temmam et al. 2022 en la Figura S3.

    Después me percaté que la estructura 6MOJE, del artículo de Lan et al.
    (2020) Structure of the SARS-CoV-2 spike... Nature... tiene un análisis
    en ConSurf Job. Ahí se ve muy bien cuál es la región hipervariable.


    #Aquí coloco la secuencia que corresponde a la región 1126 a 1376:
    #
    #tct gtc cta tat aat tcc gca tca ttt tcc act ttt aag tgt tat gga gtg tct cct act
    #S   V   L   Y   N   S   A   S   F   S   T   F   K   C   Y   G   V   S   P   T
    #aaa tta aat gat ctc tgc ttt act aat gtc tat gca gat tca ttt gta att aga ggt gat
    #K   L   N   D   L   C   F   T   N   V   Y   A   D   S   F   V   I   R   G   D
    #gaa gtc aga caa atc gct cca ggg caa act gga aag att gct gat tat aat tat aaa tta
    #E   V   R   Q   I   A   P   G   Q   T   G   K   I   A   D   Y   N   Y   K   L
    #cca gat gat ttt aca ggc tgc gtt ata gct tgg aat tct aac aat ctt gat tct aag gtt
    #P   D   D   F   T   G   C   V   I   A   W   N   S->>N   N   L   D   S   K   V
    #ggt ggt aat ta
    #G   G   N
    #
    #Esta región es la que tiene la D más alta de acuerdo al análisis. La región
    #hipervariable comienza en donde está el ->>.


    Después seleccioné las secuencias.

        -> Wuhan-Hu-1/2019 MN908947
        -> BANAL-52        MZ937000.1:21512-25321
        -> BANAL-103       MZ937001.1:21497-25294 (103 y 236 son casi idénticos)
           BANAL-236       MZ937003.1:21537-25334
        -> Guangdong_1     EPI_ISL_410721
        -> bat_RShSTT182   EPI_ISL_852604:21535-25290 Ojo esta secuencia tiene una nota en GISAID
        -> RaTG13          MN996532
        ------------------------
           RpYN06          MZ081381.1:21559-25299 RpYN06 y RmYN02 son muy similares en Lytras et al 2022 Fig.1 (no recombinación)
        -> |RmYN02          EPI_ISL_412977:21544..25227| Según el estudio de Lytras et al 2022 Fig.1 (no recombinación)
        -> PrC31           MW703458.1:21521-25261
        -> |bat_SL_CoVZC45  MG772933| Según el estudio de Lytras et al 2022 Fig.1 (no recombinación)
        -> RacCS203        MW251308.1:21562-25245
           BANAL-116       MZ937002.1:21323-25006 (116 y 247 son casi idénticos)
        -> BANAL-247       MZ937004.1:21518-25201
        ------------------------
        -> Guangxi_P4L     EPI_ISL_410538:21540-25343
        ------------------------
        -> |RsYN04          MZ081380.1:21515..25273| Según el estudio de Lytras et al 2022 Fig.1 (no recombinación)
        ------------------------
        -> Rc_o319         LC556375
        ------------------------

    Hago un archivo que se llama spike14.fasta

    # Antes de seguir tengo que arreglar los nombres de las secuencias...
    #
    #Vamos ahora a alinear las secuencias.
    #
    #    spike14/
    #    spike14.fasta -> MEGA11 -> muscle codon alignment -> varios formatos y
    #    una filogenia con Máxima verosimilitud y bootstrap.
    #    readal -in spike14.fas -out spike14.nex -nexus

    Me traigo las secuencias ya alineadas:

        spike14_250_200/
        cp ../spike14_250_125/spike14* .

    Ahora hago un mapa de la proteína S en la alineación múltiple

                                sin gaps       | con gaps         |
        ------------------------secuencia (aa) | alineación (nuc) | comentario
        péptido señal           | 1 - 12       | 1 - 39           | Temmam et al. 2022; supplementary material
        S1-NTD                  | 13 - 306     | 40 - 918         | Corresponde al dominio NTD en Boni et al. 2020
        RBD                     | 319 - 542    | 985 - 1656       | Corresponde al dominio CTD en Boni et al. 2020; Lan et al. 2020, Fig. 1
        Recombination 10 - 11   | 327          | 1011             | Temmam et al. 2022; supplementary material
        K417                    | 417          | 1279 - 1281      | Temmam et al. 2022;  Lan et al. 2020; supplementary material
        hypervloop (RBM)        | 439 - 505    | 1345 - 1545      | de acuerdo a ConSurf Job 1684375629 (8GZ5_A); Lan et al. 2020, Fig. 1
        G446                    | 446          | 1366 - 1368      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Y449                    | 449          | 1375 - 1377      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Y453                    | 453          | 1387 - 1389      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        L455                    | 455          | 1393 - 1395      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        F456                    | 456          | 1396 - 1398      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        A475                    | 475          | 1453 - 1455      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        F486                    | 486          | 1486 - 1488      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        N487                    | 487          | 1489 - 1491      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Y489                    | 489          | 1495 - 1497      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Q493                    | 493          | 1507 - 1509      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        G496                    | 496          | 1516 - 1518      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Q498                    | 498          | 1522 - 1524      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        T500                    | 500          | 1528 - 1530      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        N501                    | 501          | 1531 - 1533      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        G502                    | 502          | 1534 - 1536      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        Y505                    | 505          | 1543 - 1545      | Temmam et al. 2022; Lan et al. 2020; supplementary material
        S1-(SD1 + SD2)          | 543 - 680    | 1657 - 2070      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        Furin cleavage site     | 681 - 685    | 2071 - 2085      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        S2                      | 686 - 1274   | 2086 - 3855      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        Fusion peptide          | 788 - 806    | 2395 - 2451      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        Furin cleavage site     | 815 - 816    | 2476 - 2481      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        Internal fusion peptide | 817 - 833    | 2482 - 3855      | Temmam et al. 2022; Lan et al. 2020, Fig. 1
        Recombination 11 - 12   | 939          | 2849             | Temmam et al. 2022; Lan et al. 2020, Fig. 1

    # Ahora hago un script de Perl para extraer secciones de las secuencias
    # alineadas fasta
    #
    #     practica_variable_loop/
    #     slices.pl spikeloop.fasta 1 918
    #     out: outfile_1_918.fasta
    #     mv outfile_1_918.fasta spikeloop_1_918_NTD.fasta
    #
    #     practica_variable_loop/
    #     slices.pl spikeloop.fasta 2086 3855
    #     out: outfile_2086_3855.fasta
    #     mv outfile_2086_3855.fasta spikeloop_2086_3855_S2.fasta
    #
    #     practica_variable_loop/
    #     perl joins.pl spikeloop_1_918_NTD.fasta spikeloop_2086_3855_S2.fasta
    #     out: outfile_spikeloop_1_918_NTD.fasta_spikeloop_2086_3855_S2.fasta.fasta
    #     mv outfile_spikeloop_1_918_NTD.fasta_spikeloop_2086_3855_S2.fasta.fasta spikeloop_NTD_and_S2.fasta
    #
    #     practica_variable_loop/
    #     cp spikeloop_1_918_NTD.nex.run?.t NTD/
    #     cp spikeloop_2086_3855_S2.nex.run?.t S2/
    #
    #     practica_variable_loop/NTD/
    #     galax1 --listfile infile --skip 500 --outgroup 5 --details
    #
    #     practica_variable_loop/S2/
    #     galax1 --listfile infile --skip 500 --outgroup 5 --details
    #
    #     practica_variable_loop/NTD_S2/
    #     galax1 --listfile infile --skip 500 --outgroup 5 --details

    Lo que voy a hacer aquí es hacer el análisis de las distintas regiones de
    la proteína Spike por separado. A partir de los anáisis de verosimilitudes
    marginalizadas con el algoritmo steping stone que hice en las carpetas
    spike14_250_125, spike14_250_150 y spike14_250_200, se puede ver que la
    historia evolutiva del gen se puede dividir en tres secciones: 1ra) va de
    la columna 1 a la 751+250 o 801+250; la 2da) va de 876 o 901 o 1001 a la
    2876+250 o 2851+250 o 3001+250. Estas regiones corresponden aproximadamente
    con las recombinaciones detectadas en Temmam et al. (2022) que están en
    la columna 1011 y en 2849.

    Ahora bien, de acuerdo al análisis del algoritmo ss de la carpeta
    spike14_250_150, en el cuadrante m[10,14] el factor Bayes es de -10.44 y
    la disonancia es de 16.5417. Esta región corresponde a la región
    hipervariable que va de la columna 1345 - 1545. El cuadrante m[10,14]
    comprende los fragmentos 1351 a la 1351+250=1601 y 1951 a la 1951+250=2201.

    Si vamos al archivo de GALAX

        spike14_manual/
        cp ../spike14_250_150/information/fragments_1351_to_1600_AND_1951_to_2200.run* m10_14/

    Ahí veo que las especies que entran en conflicto entre las dos regiones son:

        **-*-**--*----
        Order of taxa in split representations:
                   1* WuhanHu1
                   2* RaTG13
                   3  CoVZC45
                   4* Guangdong1
                   5  Rco319
                   6* BANAL52
                   7* BANAL103
                   8  BANAL247
                   9  GuangxiP4L
                  10* RShSTT182
                  11  PrC31
                  12  RacCS203
                  13  RsYN04
                  14  RmYN02


    Entonces, voy a partir la alineación de acuerdo a los puntos de
    recombinación detectados por Temmam et al. (2022) y además voy a hacer
    una bipartición más en las columnas 1345 - 1545 que es en donde parece
    estar la RBM (el loop hipervariable).

    Entonces...

        spike14_manual/
        cp ~/Documents/Scriptoma/slices.pl .
        cp ~/Documents/Scriptoma/joins.pl .
        cp ../spike14_250_150/spike14* .

      # La primer región no recombinante
        spike14_manual/
        perl slices.pl spike14.fas 1 1010
        out: outfile_1_1010.fasta
        mv outfile_1_1010.fasta spike14_1_1010.fas

      # La región antes del RBM de la segunda región no recombinante
        spike14_manual/
        perl slices.pl spike14.fas 1011 1344
        out: outfile_1011_1344.fasta
        mv outfile_1011_1344.fasta spike14_1011_1344.fas

      # El loop hipervariable (RBM)
        spike14_manual/
        perl slices.pl spike14.fas 1345 1545
        out: outfile_1345_1545.fasta
        mv outfile_1345_1545.fasta spike14_1345_1545.fas

      # La región después del RBM de la segunda región no recombinante
        spike14_manual/
        perl slices.pl spike14.fas 1546 2848
        out: outfile_1546_2848.fasta
        mv outfile_1546_2848.fasta spike14_1546_2848.fas

      # La tercer región no recombinante
        spike14_manual/
        perl slices.pl spike14.fas 2849 3855
        out: outfile_2849_3855.fasta
        mv outfile_2849_3855.fasta spike14_2849_3855.fas

      # La segunda región no recombinante incluyendo al RBM
        spike14_manual/
        perl slices.pl spike14.fas 1011 2849
        out: outfile_1011_2849
        mv outfile_1011_2849.fasta spike14_1011_2849.fas

      # La segunda región no recombinante sin incluir al RBM
        spike14_manual/
        perl joins.pl spike14_1011_1344.fas spike14_1546_2848.fas
        out: outfile_spike14_1011_1344.fas_spike14_1546_2848.fas.fasta
        mv outfile_spike14_1011_1344.fas_spike14_1546_2848.fas.fasta spike14_1011_2849_sin_RBM.fas

        Utilizo readal -in filein -out fileout -nexus para hacer los nexfiles
        Ejemplo:
        spike14_manual/
        readal -in spike14_1011_2849_sin_RBM.fas -out spike14_1011_2849_sin_RBM.nex -nexus

    ----------------------------------------------------------------------------
    Nota: faltaría hacer un análisis de GARD para ver que regiones detecta en
    esta alineación múltiple.
    ----------------------------------------------------------------------------

    Corro MrBayes en todos los fragmentos, ejemplo:

        MrBayes/
        mb infile_1345_1545.txt

    Luego hago los análisis de galax1:

        Information/
        #grep -c 'tree gen' ../MrBayes/spike14.nex.run1.t
        #2001
        #galax1 --listfile infile --skip 200 --outgroup 5 --details
        perl goGALAX_manual.pl

    Separo los casos especiales

        Information/
        mv *RBM* Especiales/
        mv *_1011_2849_* Especiales/
        mv spike14_AND_* Especiales/

    Luego hago un resumen de la información

        spike14_manual/
        perl summary_H_manual.pl
        out: statistics_1.tsv, statistics_2.tsv
        mv statistics_* Statistics/

    Como el programa no lo hizo bien, arreglé el archivo statistics_2.tsv
    manualmente. Le hacía falta una arregladida. Y ahora sí visualizo los
    archivos en summary_H.R.

    Ahora sigue hacer el análisis de ss.

        spike14_manual/
        cp MrBayes/*.nex MrBayesML
        rm spike14.nex spike14_1011_2849_sin_RBM.nex spike14_1011_2849.nex

        spike14_manual/
        cp spike14_1_1010.fas spike14_1011_1344.fas spike14_1345_1545.fas spike14_1546_2848.fas spike14_2849_3855.fas fragments/
        perl concatena_manual.pl

    Nota:
        Tengo que tomar en cuenta que cuando las secuencias empalman, sus
        tamaños son más pequeños, ejemplo:
            me falta poner los ejemplos

    Corro todos los ss manualmente, por ejemplo:
        concatena2/
        mb infile_1011_to_1344_AND_1345_to_1545.txt
        MrBayesML/
        mb infile_2849_3855.txt

    Y luego hago el resumen:
        spike14_manual/
        perl minaMLmanual.pl > statistics_ss_manual.tsv
        lo arreglo manualmente y luego...
        mv statistics_ss_manual.tsv Statistics/

    Y luego lo miro con summary_H.R

    Ahora pongo a prueba algunas hipótesis topológicas en el RBM

        RBM/
        cp ../MrBayes/spike14_1345_1545.nex .
        mb
        execute spike14_1345_1545.nex
        lset  nst=6 rates=invgamma ngammacat=4;
        prset brlenspr=unconstrained:GammaDir(1.0, 0.01, 1.0, 1.0) shapepr=exp(1.0) statefreqpr=dirichlet(1.0,1.0,1.0,1.0) revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0);
        constraint RaTG13RsYN04 = RatG13 RsYN04
        constraint noRaTG13RsYN04 negative = RatG13 RsYN04
        prset topologypr=constraints(RaTG13RsYN04)
        ss ngen=2500000 diagnfreq=2500
            Run   Marginal likelihood (ln)
            ------------------------------
              1     -1697.59
              2     -1698.27
            ------------------------------
            Mean:   -1697.87 -> esta es con ngen=250000
            Run   Marginal likelihood (ln)
            ------------------------------
              1     -1698.75
              2     -1699.15
            ------------------------------
            Mean:   -1698.93 -> esta es con ngen=2500000
        prset topologypr=constraints(noRaTG13RsYN04)
        ss
            Run   Marginal likelihood (ln)
            ------------------------------
              1     -1701.06
              2     -1699.67
            ------------------------------
            Mean:   -1700.14 -> esta es con ngen=250000
            Run   Marginal likelihood (ln)
             ------------------------------
               1     -1700.95
               2     -1700.62
             ------------------------------
             Mean:   -1700.77 -> esta es con ngen=2500000
        Resultado: 1700.14 -1697.87 = 2.27  -> esta es con ngen=250000
        Resultado: 1700.77 -1698.93 = 1.84  -> esta es con ngen=2500000

        constraint RaTG13GuangxiP4L = RatG13 GuangxiP4L
        constraint noRaTG13GuangxiP4L negative = RatG13 GuangxiP4L
        prset topologypr=constraints(RaTG13GuangxiP4L)
        ss ngen=250000 diagnfreq=2500
            Run   Marginal likelihood (ln)
            ------------------------------
              1     -1702.16
              2     -1701.56
            ------------------------------
            Mean:   -1701.82
        prset topologypr=constraints(noRaTG13GuangxiP4L)
        ss
            Run   Marginal likelihood (ln)
            ------------------------------
              1     -1700.84
              2     -1700.22
            ------------------------------
            Mean:   -1700.48
        Resultado: 1701.82 -1700.48 = 1.34  -> esta es con ngen=250000

2)  Y bueno, ahora podría hacer todo el análisis que ha sé hacer...

    Me traigo los archivos.

        spike14_250_200/
        cp ../spike14_250_150/*.pl .
        cp ../spike14_250_150/*.R .

    #Calculo el modelo de evolución molecular de ncovsel.aln.fasta y es:
    #
    #  GTR+G+I +I=0.44, +G=0.97, R=2.30, f(A)=0.290, f(T)=0.316, f(C)=0.190, f(G)=0.203
    #
    #Infiero una filogenia por ML en Mega11:
    #
    #  ncov16_rec1/
    #  ncovsel.tree.mtsx, ncovsel.tree.nwk
    #
    #le coloco la raíz en el mismo lugar que en la Figura 1 de Delaune et al.
    #(2021).
    #
    #Ahora preparo el archivo para simular el evento de recombinación:
    #
    #  ncov16_rec1/
    #  cp ncovsel.tree.nwk ncovsel.rec1.nwk
    #
    #Ojo, ncovsel.rec1.nwk es una filogenia sin raíz. Manualmente cambio de
    #clado a la secuencia Wuhan-Hu-1_2019 en el archivo ncovsel.rec1.nwk.
    #
    #Hago el archivo que usaré para simular el evento de recombinación:
    #
    #  ncov16_rec1/
    #  ncovsel.sim.trees
    #
    #Este archivo tiene 4000 nucleótidos en cada extremo 5' y 3' con una
    #topología y 1000 nucleótidos en el centro, con otra topología.
    #
    #Ahora simulo las secuencias:
    #
    #  ncov16_rec1/                      A     C     G     T
    #  seq-gen -m GTR -a 0.97 -i 0.44 -f 0.290 0.190 0.203 0.316 -p 3 -l 9000 -of < ncovsel.sim.trees > ncovsel.rec1.fasta
    #  Random number generator seed: 6392943637274883427

3)  Sigue hacerle el análisis de GALAX.

    Subdivido la alineación múltiple:

        spike14_250_200/
        perl fragment.pl spike14.fas 250 200
        out: fragments/fragment_X_to_Y.fasta
          Statistics:
          Alignment length....: 3855
          Windows size........: 250
          Sliding positions...: 200
          Number of fragments.: 19
          Columns not analyzed: 55

        spike14_250_200/
        perl prepare_list_mcmc.pl
        out: list_for_MrBayes_mcmc.txt
          Statistics:
          numero menor: 1
          numero mayor: 3850
          intervalo...: 200
          ventana.....: 250
          Number of files: 19

        spike14_250_200/
        rm phylogenies/* # si es que tiene archivos

      Hago el análisis de MrBayes mcmc:
        -> Lo estoy haciendo con los siguientes parámetros:
        mcmc ngen=1000000 samplefreq=500 printfreq=10000 burninfrac=0.25 nchains=4 nruns=2

        spike14_250_200/
        perl goMrBayes_local_threads.pl 1 4 1 -> 11:49 hrs. 01/junio/2023
        perl goMrBayes_local_threads.pl 5 8 2
        perl goMrBayes_local_threads.pl 9 12 3
        perl goMrBayes_local_threads.pl 13 16 4
        perl goMrBayes_local_threads.pl 17 19 5
        out:
          fragment_1_to_250.mcmc.ckp
          fragment_1_to_250.mcmc.ckp~
          fragment_1_to_250.mcmc.con.tre
          fragment_1_to_250.mcmc.lstat
          fragment_1_to_250.mcmc.mcmc
          fragment_1_to_250.mcmc.parts
          fragment_1_to_250.mcmc.pstat
          fragment_1_to_250.mcmc.run1.p
          fragment_1_to_250.mcmc.run1.t
          fragment_1_to_250.mcmc.run2.p
          fragment_1_to_250.mcmc.run2.t
          fragment_1_to_250.mcmc.trprobs
          fragment_1_to_250.mcmc.tstat
          fragment_1_to_250.mcmc.vstat
          fragment_1_to_250.mcmc.out.txt
          fragment_1_to_250.nex

    Hago el análisis de GALAX:

        spike14_250_200/
        grep -c 'tree gen' phylogenies/fragment_1_to_2000.mcmc.run1.t
          2001

        spike14_250_200/
        perl goGALAX_local.pl 200
        out:
          information/fragment_1_to_250.output-galax-details.txt
          information/fragment_1_to_250.output-galax-merged.txt
          information/fragment_1_to_250.output-galax.txt
          information/fragments_1_to_250_AND_251_to_500.run1.output-galax-details.txt
          information/fragments_1_to_2000_AND_251_to_500.run1.output-galax-merged.txt
          information/fragments_1_to_2000_AND_251_to_500.run1.output-galax.txt
          information/fragments_1_to_2000_AND_251_to_500.run2.output-galax-details.txt
          information/fragments_1_to_2000_AND_251_to_500.run2.output-galax-merged.txt
          information/fragments_1_to_2000_AND_251_to_500.run2.output-galax.txt

    Hago un resumen de la información

        spike14_250_200/
        perl summary_H.pl
        out: statistics_1.tsv, statistics_2.tsv
        mv statistics_* Statistics/

        spike14_250_150/
        summary_H.R

    Ahora sigue calcular la verosimilitud marginalizada.

        spike14_250_200/
        perl prepare_list_ss.pl
        out: list_for_MrBayes_ss.txt

        spike14_250_150/
        rm marginalL/* # si es que tiene archivos

    Hago el análisis de MrBayes ss:
      -> Lo estoy haciendo con los siguientes parámetros
      ss ngen=1000000 samplefreq=1000 printfreq=10000

        spike14_250_200/
        perl goMrBayesss_local_threads.pl 1 4 1 -> 12:42 hrs. 01/junio/2023
        perl goMrBayesss_local_threads.pl 5 8 2
        perl goMrBayesss_local_threads.pl 9 12 3
        perl goMrBayesss_local_threads.pl 13 16 4
        perl goMrBayesss_local_threads.pl 17 19 5
        out:
          marginalL/fragment_1_to_250.ss.ckp
          marginalL/fragment_1_to_250.ss.ckp~
          marginalL/fragment_1_to_250.ss.mcmc
          marginalL/fragment_1_to_250.ss.out.txt
          marginalL/fragment_1_to_250.ss.run1.p
          marginalL/fragment_1_to_250.ss.run1.t
          marginalL/fragment_1_to_250.ss.run2.p
          marginalL/fragment_1_to_250.ss.run2.t
          marginalL/fragment_1_to_250.ss.ss

        spike14_250_200/
        perl concatena.pl
        out: concatena2/fragments_X_to_Y_AND_Z_to_W.fasta
          Nota:
          Tengo que tomar en cuenta que cuando las secuencias empalman, sus
          tamaños son más pequeños, ejemplo:
          jose$ ls -l concatena2/fragments_1_to_2000_AND_501_to_2500.fasta
          10036 Sep 20 20:48 concatena2/fragments_1_to_2000_AND_501_to_2500.fasta  -> 2500 nt
          jose$ ls -l concatena2/fragments_1_to_2000_AND_1001_to_3000.fasta
          12036 Sep 20 20:48 concatena2/fragments_1_to_2000_AND_1001_to_3000.fasta -> 3000 nt
          jose$ ls -l concatena2/fragments_1_to_2000_AND_1501_to_3500.fasta
          14036 Sep 20 20:48 concatena2/fragments_1_to_2000_AND_1501_to_3500.fasta -> 3500 nt
          jose$ ls -l concatena2/fragments_1_to_2000_AND_2001_to_4000.fasta
          16036 Sep 20 20:48 concatena2/fragments_1_to_2000_AND_2001_to_4000.fasta -> 4000 nt
          jose$ ls -l concatena2/fragments_1_to_2000_AND_2501_to_4500.fasta
          16036 Sep 20 20:48 concatena2/fragments_1_to_2000_AND_2501_to_4500.fasta -> 4000 nt

        spike14_250_200/
        perl prepare_list_ss_con.pl
        out: list_for_MrBayes_ss_con.txt
        Number of files: 171

        spike14_250_150/
        rm marginalLcon/* # si es que tiene archivos

        -> Lo estoy haciendo con los siguientes parámetros
        ss ngen=1000000 samplefreq=500 printfreq=10000

        spike14_250_200/
        perl goMrBayesss_con_local_threads.pl 1 34 1 -> 20:21, 31/mayo/2023
        perl goMrBayesss_con_local_threads.pl 35 69 2
        perl goMrBayesss_con_local_threads.pl 70 104 3
        perl goMrBayesss_con_local_threads.pl 105 139 4
        perl goMrBayesss_con_local_threads.pl 140 171 5
        out:
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.ckp
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.ckp~
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.mcmc
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.out.txt
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.run1.p
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.run1.t
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.run2.p
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.run2.t
          marginalLcon/fragments_10001_to_12000_AND_10501_to_12500.ss.ss

        spike14_250_200/
        perl summary_ss_local.pl 250 200 3855
        out: statistics_ss.tsv
        mv statistics_ss.tsv Statistics/

        spike14_250_200/
        summary_H.R
