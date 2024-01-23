#! /usr/bin/awk -f
# minusK_blob.awk   rough minus K for the output of diamond blast (blobtrools)  , to select one protein per region   Bernardo   25oct2017   v. 27oct
#usage:  minusK_blob.awk  Aqua7shdef_refseqX.m8  > Aqua7shdef_refseqX.minusK1.m8
# diamond blast usage: allow more than one hit to be reported (ideally, 10 or 50)  :
# diamond blastx --block-size 2 --query Aqua7shdef.fasta  --max-target-seqs 100 --sensitive --index-chunks 1  --threads 20 -d /draft1/db/refseq_protein.dmnd  --evalue 1e-25  --outfmt tab --out Aqua7shdef_refseqX.m8

(1==1){if ($7<$8){s_start=$7; s_end=$8}else{s_start=$8; s_end=$7}}

(NR==1)  {
                  gi=$1
                  i ++
                  DNA_start_array[i] = s_start ; DNA_end_array[i] = s_end
                  print $0
                  }

(gi != $1){
                        delete DNA_start_array; delete DNA_end_array
                  i=1
                  gi=$1
                  DNA_start_array[i] = s_start ; DNA_end_array[i] = s_end
                  print $0

          }

(gi == $1){
                   flag_overlap=0
                   for (j=1; j <= i; j++){
                                                if ((s_start >= DNA_start_array[j]) && (s_start <= DNA_end_array[j]) ){flag_overlap ++}   # ref  query ref   query  or  ref query query  ref
                                                if ((s_end   >= DNA_start_array[j]) && (s_end   <= DNA_end_array[j]) ){flag_overlap ++}   # query ref  query ref    or  ref query query  ref
# gi|2    B4MPU6                  83.1    467     45      3       10629   9232    32      465     8.8e-203        720.7
# gi|2    A0A0L0CAM7      76.3    497     75      5       10668   9190    407     864     2.7e-196        699.1
                                                if ((DNA_start_array[j] >= s_start ) && (DNA_start_array[j] <= s_end) ){flag_overlap ++}  # query ref  ref  query
                                                }
                   if (flag_overlap==0){
                                i++
                                                DNA_start_array[i] = s_start ; DNA_end_array[i] = s_end
                                                print $0
                                                }
                        #print "teste", flag_overlap, $0
          }


#program evalutaion:
# wc -l Aqua7shdef_refseqX.minusK1.m8  10703
# wc -l Aqua7shdef_refseqX.m8         487712


# gi|1    XP_553178.3     46.2    2162    737     62      186749  192874  23      1877    0.0e+00 1285.8
# gi|1    XP_319835.4     72.4    1031    176     14      240910  243996  1       924     0.0e+00 1153.7
# gi|1    XP_019533309.1  70.6    1063    197     19      240808  243996  536     1482    0.0e+00 1117.4
# gi|1    XP_019554510.1  70.6    1063    196     19      240808  243996  535     1481    0.0e+00 1117.4
# gi|1    XP_001863809.1  69.0    1063    212     17      240808  243996  57      1002    0.0e+00 1091.3
# gi|1    XP_001652549.2  69.4    1063    208     18      240808  243996  524     1469    0.0e+00 1086.6
# gi|1    XP_311719.5     69.6    771     178     6       3510    1204    4       720     2.2e-301        1050.8
# gi|1    XP_019564547.1  55.8    1256    341     43      116923  120585  245     1321    2.5e-297        1037.3
# gi|1    XP_019558685.1  56.1    1182    330     42      116923  120369  245     1270    1.6e-296        1034.6
# gi|2    XP_314118.4     50.6    2146    616     33      251463  245149  4       1746    0.0e+00 1747.6
# gi|2    XP_001661658.1  72.1    1497    322     22      267322  271755  340     1759    0.0e+00 1726.8
# gi|2    XP_019544369.1  71.9    1493    326     21      267322  271755  408     1822    0.0e+00 1717.2
# gi|2    XP_019544572.1  71.9    1493    326     21      267322  271755  408     1822    0.0e+00 1717.2
# gi|2    XP_314116.4     63.0    1482    42      7       267322  271755  408     1386    0.0e+00 1437.9
# gi|2    XP_310244.7     65.7    1236    278     24      64970   68587   168     1287    0.0e+00 1387.9



# gi|2    B4MPU6                  83.1    467     45      3       10629   9232    32      465     8.8e-203        720.7
# gi|2    A0A0L0CAM7      76.3    497     75      5       10668   9190    407     864     2.7e-196        699.1

