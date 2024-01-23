#!/bin/bash

usage(){
echo '
Written by Isabela Almeida
Last modified January 23, 2024

Description: Runs Blobtools, including the preprocessing steps (BLASTx and
BLASTn) for the desired input fasta file.

Produces temporary files and outputs Y chromosome pure fasta file.

Usage: bash ipda_scriptB_BlobTools.sh [-t target_genome.fasta] [-r reference_genome.fasta] [-c coverage_info.cov] [-p phylum] [-s superkingdom] [-o path/to/output/]

## Input:

-t <target_genome.fasta>    path/to/Target Genome FASTA file
-r <reference_genome.fasta> path/to/Reference Genomes FASTA file
-c <coverage.cov>           path/to/Coverage COV file. If no coverage info is
                            available, create a coverage file with
                            grep ">" ${tar_fasta} | cut -d" " -f1 | cut -d">" -f2 | while read gi; do echo -e "$gi\t100\t300" >> coverage.cov ; done
-p <phylum>                 Phylym to be kept
-s <superkingdom>           Superkingdom to be kept
-o <path/to/output/>        Path to output directory in which files will be
                            saved (Script creates the directory)

## Output:

logfile.txt                 Logfile with commands executed, date, outputs info,
                            and memory/CPU usage for major executions
blastx.hitsfile             DIAMOND BLASTx hitsfile    
blastn.hitsfile             BLASTn hitsfile
blobtools_k3*               BlobTools create blobDB.json files
blast*.png                  BlobTools plot files
*multiplot.stats.txt        BlobTools stats files
blobtools_phylum*           BlobTools view phylum table file
blobtools_superkingdom*     BlobTools view superkingdom table file
blobtools_FinalTable        Combined phylum/superkingdom + BLASTn/BLASTx table
contaminatedGIs             Contaminated IDs identified
contaminationFree_GIs       IDs from the selected phylum and superkingdom
contaminationFreeSeq.fasta  Contamination free target fasta multi line fasta
                            to be used for scriptC (pure sequences)
singleline.fasta            Contamination free target fasta single line fasta
                            to be used for scriptC (pure sequences)
contaminationFree_genomeGIs Contamination free genome IDs
contaminationFree_genome.fasta Contamination free reference fasta

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
'
}

## Display help message
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

#  ____       _   _   _                 
# / ___|  ___| |_| |_(_)_ __   __ _ ___ 
# \___ \ / _ \ __| __| | '_ \ / _` / __|
#  ___) |  __/ |_| |_| | | | | (_| \__ \
# |____/ \___|\__|\__|_|_| |_|\__, |___/
#                             |___/     

## Save execution ID
pid=$BASHPID

## Exit when any command fails
set -e

## Path to databases
blastXdb=/draft2/blast_db/oakdb_prot_01ago19.dmnd
blastNdb=/draft3/dupim_temp/oakdb2/oakdb_genomic_01ago19 #.nal
PTdb_humanContaminants=/home3/blobtools_genomes/Breitwieser_PTN_list.txt # Refseq according to Breitwieser et al 2019 
NTdb_humanContaminants=/home3/blobtools_genomes/Breitwieser_DNA_list.txt # Refseq according to Breitwieser et al 2019 
db_speciesContaminants=/home3/blobtools_genomes/blacklist_26jun19.txt
blobColors=/home3/blobtools_genomes/blobcolors_24ago19.txt

## Path to scripts
fastacmd_loose2=scripts/CarvalhoClark2013/fastacmd_loose2.awk
minusK_blob=scripts/CarvalhoClark2013/minusK_blob.awk

#................................................
#  Input parameters from command line
#................................................

## Get parameters from command line flags
while getopts "t:r:c:p:s:o:" flag
do
    case "${flag}" in
        t) tar_fasta=${OPTARG};;            # path/to/Target Genome FASTA file
        r) ref_fasta=${OPTARG};;            # path/to/Reference Genome FASTA file
        c) cov_file=${OPTARG};;             # path/to/Coverage file file
        p) phylum=${OPTARG};;               # Phylym to be kept
        s) superkingdom=${OPTARG};;         # Superkingdom to be kept
        o) out_path=${OPTARG};;             # Path to output directory
        ?) echo script usage: bash ipda_scriptB_BlobTools.sh [-t target_genome.fasta] [-r reference_genome.fasta] [-c coverage_info.cov] [-p phylum] [-s superkingdom] [-o path/to/output/] >&2
           exit;;
    esac
done

#................................................
#  Extract info from input files
#................................................

## Get input files stem
tar_stem=`echo "$(basename "${tar_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')"`

#................................................
#  Set output path and stem
#................................................

## Create the output directory
mkdir -p ${out_path}

## Get stem for outputs
out_stem=BlobTools_${tar_stem}

#................................................
#  Set Logfile stem
#................................................

## Set Logfile stem
# it contains all the executed commands with date/time;
# and memory/CPU usage for major executions (i.e. DIAMOND BLASTn,
# BLASTx and BlobTools create)
logfile=${out_path}logfile_BlobTools_${out_stem}.txt

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing ipda_scriptB_BlobTools.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs files:"
echo
echo "## Target Genome FASTA file:                  ${tar_fasta}"
echo "## Reference Genome FASTA file:               ${ref_fasta}"
echo "## Coverage info file:                        ${cov_file}"
echo
echo "## Given input parameters/thresholds:"
echo
echo "## Phylum to be kept:                         ${phylum}"
echo "## Superkingdom to be kept:                   ${superkingdom}"
echo
echo "## Output files saved to:                     ${out_path}"
echo
echo "## logfile will be saved as:                  ${logfile}"

# ____  _             _   _                   
#/ ___|| |_ __ _ _ __| |_(_)_ __   __ _       
#\___ \| __/ _` | '__| __| | '_ \ / _` |      
# ___) | || (_| | |  | |_| | | | | (_| |_ _ _ 
#|____/ \__\__,_|_|   \__|_|_| |_|\__, (_|_|_)
#                                 |___/  

#................................................
#  Print Execution info to logfile
#................................................

exec &> "${logfile}"

date
echo "## Executing ipda_scriptB_BlobTools.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs files:"
echo
echo "## Target Genome FASTA file:                  ${tar_fasta}"
echo "## Reference Genome FASTA file:               ${ref_fasta}"
echo "## Coverage info file:                        ${cov_file}"
echo
echo "## Given input parameters/thresholds:"
echo
echo "## Phylum to be kept:                         ${phylum}"
echo "## Superkingdom to be kept:                   ${superkingdom}"
echo
echo "## Output files saved to:                     ${out_path}"
echo
echo "## logfile will be saved as:                  ${logfile}"

set -v

#................................................
#  STEP 1 - TAXONOMY AND ALIGNMENT ANALYSIS
#................................................

date ## Starting Diamond BLASTx & BLASTn at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
diamond blastx --query ${tar_fasta} --db ${blastXdb} --block-size 2 --max-target-seqs 250 \
    --sensitive --index-chunks 1 --evalue 1e-25 --threads 64 \
    --outfmt 6 qseqid staxids bitscore sseqid sstart send qstart qend pident \
    --out ${out_path}blastx_${outPattern}.hitsfile > ${out_path}tmp_blastx_${outPattern}.out &
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
blastn -db ${blastNdb} -query ${tar_fasta} \
    -outfmt '6 qseqid staxids bitscore sseqid sstart send qstart qend pident' \
    -evalue 1e-10 -num_threads 64 > ${out_path}blastn_${outPattern}.hitsfile &
wait
rm ${out_path}tmp_blastx_${outPattern}.out

#................................................
#  STEP 2 - PROCESS ALIGNMENT OUTPUTS
#................................................

date ## Starting to filter BLASTx output [Add Human Contaminants info] at
awk 'BEGIN {OFS="\t"} FNR==1{ f++ } f==1 { a[$1]=$1; next } { if($4 in a); else print $0 }' ${PTdb_humanContaminants} ${out_path}blastx_${outPattern}.hitsfile > ${out_path}tmp_ptHumanInfo_blastx_${outPattern}.txt

date ## Starting to filter BLASTx output [Add Species ID Contaminants info] at
awk 'BEGIN {OFS="\t"} FNR==1{ f++ } f==1 { a[$1]=$1; next } { if($2 in a); else print $0 }' ${db_speciesContaminants} ${out_path}tmp_ptHumanInfo_blastx_${outPattern}.txt > ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2
rm ${out_path}tmp_ptHumanInfo_blastx_${outPattern}.txt

date ## Starting to get BLASTx best hits at
${minusK_blob} k=3 ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2 > ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2_k3
rm ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2

date ## Starting to filter BLASTn output [Add Human Contaminants info] at
awk 'BEGIN {OFS="\t"} FNR==1{ f++ } f==1 { a[$1]=$1; next } { if($4 in a); else print $0 }' ${NTdb_humanContaminants} ${out_path}blastn_${outPattern}.hitsfile > ${out_path}tmp_ntHumanInfo_blastn_${outPattern}.txt

date ## Starting to filter BLASTn output [Add Species ID Contaminants info] at
awk 'BEGIN {OFS="\t"} FNR==1{ f++ } f==1 { a[$1]=$1; next } { if($2 in a); else print $0 }' ${db_speciesContaminants} ${out_path}tmp_ntHumanInfo_blastn_${outPattern}.txt > ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2
rm ${out_path}tmp_ntHumanInfo_blastn_${outPattern}.txt

date ## Starting to get BLASTn best hits at
${minusK_blob} k=3 ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2  > ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2_k3
rm ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2

#................................................
# STEP 3 - RUN BLOBTOOLS
#................................................

date ## Starting BlobTools create with BLASTx processed info at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
blobtools create -i ${tar_fasta} -c  ${cov_file} -t ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2_k3 -o ${out_path}blobtools_blastx_${outPattern}_k3

date ##  Starting BlobTools plot with BLASTx blobDB.json at
blobtools plot -i ${out_path}blobtools_blastx_${outPattern}_k3.blobDB.json -p 10 --noreads --colours ${blobColors} --sort_first "no-hit" --multiplot
mv blastx*.png ${out_path}
mv *multiplot.stats.txt ${out_path}

date ## Starting BlobTools phylum view with BLASTx blobDB.json at
blobtools view -i ${out_path}blobtools_blastx_${outPattern}_k3.blobDB.json
mv blobtools_blastx_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_phylum_blastx_${outPattern}_k3.blobDB.table.txt

date ## Starting BlobTools superkingdom view with BLASTx blobDB.json at
blobtools view -i ${out_path}blobtools_blastx_${outPattern}_k3.blobDB.json -r superkingdom -o superkingdom
mv superkingdom.blobtools_blastx_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_superkingdom_blastx_${outPattern}_k3.blobDB.table.txt

date ## Starting BlobTools create with BLASTn processed info at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
blobtools create -i ${tar_fasta} -c ${cov_file} -t ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2_k3 -o ${out_path}blobtools_blastn_${outPattern}_k3

date ## Starting BlobTools plot with BLASTn blobDB.json at
blobtools plot -i ${out_path}blobtools_blastn_${outPattern}_k3.blobDB.json -p 10 --noreads --colours ${blobColors} --sort_first "no-hit" --multiplot
mv blastn*.png ${out_path}
mv *multiplot.stats.txt ${out_path}

date ## Starting BlobTools phylum view with BLASTn blobDB.json at
blobtools view -i ${out_path}blobtools_blastn_${outPattern}_k3.blobDB.json
mv blobtools_blastn_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_phylum_blastn_${outPattern}_k3.blobDB.table.txt

date ## Starting BlobTools superkingdom view with BLASTn blobDB.json at
blobtools view -i ${out_path}blobtools_blastn_${outPattern}_k3.blobDB.json -r superkingdom -o superkingdom
mv superkingdom.blobtools_blastn_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_superkingdom_blastn_${outPattern}_k3.blobDB.table.txt

#................................................
# STEP 4 - COMBINE BLOBTOOLS OUTPUT TABLES
#................................................

date ## Starting to merge phylum and superkingdom info from BLASTx and BLASTn blobDB.table at
grep -v "^#" ${out_path}blobtools_superkingdom_blastx_${outPattern}_k3.blobDB.table.txt | cut -f6- > ${out_path}tmp_SuperKingBLASTxColumns_${outPattern}.txt
grep -v "^#" ${out_path}blobtools_superkingdom_blastn_${outPattern}_k3.blobDB.table.txt | cut -f6- > ${out_path}tmp_SuperKingBLASTnColumns_${outPattern}.txt
grep -v "^#" ${out_path}blobtools_phylum_blastn_${outPattern}_k3.blobDB.table.txt | cut -f6- > ${out_path}tmp_PhylumBLASTnColumns_${outPattern}.txt
grep -v "^#" ${out_path}blobtools_phylum_blastx_${outPattern}_k3.blobDB.table.txt > ${out_path}tmp_PhylumBLASTxColumns_${outPattern}.txt
paste ${out_path}tmp_PhylumBLASTxColumns_${outPattern}.txt ${out_path}tmp_SuperKingBLASTxColumns_${outPattern}.txt ${out_path}tmp_PhylumBLASTnColumns_${outPattern}.txt ${out_path}tmp_SuperKingBLASTnColumns_${outPattern}.txt > ${out_path}blobtools_FinalTable_${outPattern}_blastx_blastn_k3.txt
rm ${out_path}tmp_PhylumBLASTxColumns_${outPattern}.txt ${out_path}tmp_SuperKingBLASTxColumns_${outPattern}.txt ${out_path}tmp_PhylumBLASTnColumns_${outPattern}.txt ${out_path}tmp_SuperKingBLASTnColumns_${outPattern}.txt
rm ${out_path}blobtools_phylum_blastx_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_phylum_blastn_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_superkingdom_blastx_${outPattern}_k3.blobDB.table.txt ${out_path}blobtools_superkingdom_blastn_${outPattern}_k3.blobDB.table.txt

#................................................
#  STEP 5 - CONTAMINATION INFO AND PURE FILES
#................................................

date ## Starting to get contaminated GIs at
cut -f6,9,12,15 ${out_path}blobtools_FinalTable_${outPattern}_blastx_blastn_k3.txt | tr '\t' '\n' | sort | uniq | grep -v -e "$phylum" -v -e "$superkingdom" -v -e "no-hit" | while read contam ; do grep "$contam" ${out_path}blobtools_FinalTable_${outPattern}_blastx_blastn_k3.txt ; done | sort | uniq | cut -f1 | cut -d"|" -f2 > ${out_path}blobtools_${outPattern}_contaminatedGIs

date ## Starting to check amount of contaminated GIs at
if [ -s ${out_path}blobtools_${outPattern}_contaminatedGIs ]; then
    gis=`wc -l ${out_path}blobtools_${outPattern}_contaminatedGIs | cut -d" " -f1`
    echo $gis # number of contaminated GIs found:

    date ## Starting to create contamination free GIs list at
    grep ">" ${tar_fasta} | cut -d " " -f1 | cut -d"|" -f2 > ${out_path}tmp_mainFastaGIs
    awk -F, 'FNR==NR {f2[$1];next} !($0 in f2)' ${out_path}blobtools_${outPattern}_contaminatedGIs ${out_path}tmp_mainFastaGIs > ${out_path}${outPattern}_contaminationFree_GIs
    rm ${out_path}tmp_mainFastaGIs

    date ## Starting to create contamination free Target fasta at
    ${fastacmd_loose2} gi_type="gi" fasta=${ref_fasta} ${out_path}${outPattern}_contaminationFree_GIs > ${out_path}${outPattern}_contaminationFreeSeq.fasta

    date ## Starting to create single-line contamination free Target fasta at
    cp ${out_path}${outPattern}_contaminationFreeSeq.fasta ${out_path}${outPattern}_tmpfile
    sed -i '/^>/c\NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' ${out_path}${outPattern}_tmpfile
	tr -d '\n' < ${out_path}${outPattern}_tmpfile | sed '$s/ $/\n/' > ${outPath}${target}_singleline.fasta
	sed -i "1 i\>Single-Line_${target}" ${outPath}${target}_singleline.fasta
	echo >> ${outPath}${target}_singleline.fasta
	rm ${out_path}${outPattern}_tmpfile

    date ## Starting to create contamination free Reference Genome GIs list at
    grep ">" ${ref_fasta} | cut -d" " -f1 | cut -d"|" -f2 > ${out_path}tmp_genomeGIs
    awk -F, 'FNR==NR {f2[$1];next} !($0 in f2)' ${out_path}blobtools_${outPattern}_contaminatedGIs ${out_path}tmp_genomeGIs > ${out_path}${outPattern}_contaminationFree_genomeGIs
    rm ${out_path}tmp_genomeGIs

    date ## Starting to create contamination free Reference Genome fasta at
    ${fastacmd_loose2} gi_type="gi" fasta=${ref_fasta} ${out_path}${outPattern}_contaminationFree_genomeGIs > ${out_path}${outPattern}_contaminationFree_genomeSeq.fasta

else:
    echo "# Blobtools analysis did not found any sequences to be from a #" >> ${out_path}blobtools_${outPattern}_contaminatedGIs
    echo "# different phylyum and superkingdom than the given ones      #" >> ${out_path}blobtools_${outPattern}_contaminatedGIs
    echo "# ($phylum and $superkingdom)                                 #" >> ${out_path}blobtools_${outPattern}_contaminatedGIs
    echo "# Your sequence is already trustworthy and pure!              #" >> ${out_path}blobtools_${outPattern}_contaminatedGIs
fi

rm ${out_path}blastx_${outPattern}.hitsfile ${out_path}blastn_${outPattern}.hitsfile
rm ${out_path}tmp_ptHumanAndSpeciesInfo_blastx_${outPattern}.hitsfile2_k3
rm ${out_path}tmp_ntHumanAndSpeciesInfo_blastn_${outPattern}.hitsfile2_k3


# This will remove $VARNAMES from output file with the actual $VARVALUE
# allowing for easily retracing commands
sed -i 's,${tar_fasta},'"${tar_fasta}"',g' "$logfile"
sed -i 's,${ref_fasta},'"${ref_fasta}"',g' "$logfile"
sed -i 's,${cov_file},'"${cov_file}"',g' "$logfile"
sed -i 's,${tar_stem},'"${tar_stem}"',g' "$logfile"
sed -i 's,${phylum},'"${phylum}"',g' "$logfile"
sed -i 's,${superkingdom},'"${superkingdom}"',g' "$logfile"
sed -i 's,${blastXdb},'"${blastXdb}"',g' "$log_file"
sed -i 's,${blastNdb},'"${blastNdb}"',g' "$log_file"
sed -i 's,${PTdb_humanContaminants},'"${PTdb_humanContaminants}"',g' "$log_file"
sed -i 's,${NTdb_humanContaminants},'"${NTdb_humanContaminants}"',g' "$log_file"
sed -i 's,${db_speciesContaminants},'"${db_speciesContaminants}"',g' "$log_file"
sed -i 's,${blobColors},'"${blobColors}"',g' "$log_file"
sed -i 's,${out_path},'"${out_path}"',g' "$logfile"
sed -i 's,${minusK_blob},'"${minusK_blob}"',g' "$logfile"
sed -i 's,${fastacmd_loose2},'"${fastacmd_loose2}"',g' "$logfile"
sed -i 's,${log_file},'"${log_file}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v