#!/bin/bash

usage(){
echo '
Written by Isabela Almeida
Last modified January 23, 2024

Description: Runs YGS (Original Version), including the preprocessing step for
the required input files with Jellyfish.

Produces temporary files and outputs YGS analysis.

Usage: bash ipda_scriptA_YGS.sh [-k 15] [-y "het_reads*.fastq.gz"] [-x "hom_reads*.fastq.gz"] [-t target_genome.fasta] [-q 20] [-m INT] [-f INT] [-g INT] [-v 20] [-p 75] [-o path/to/output/]

## Input:

-k <15>                     k-mer INT size
-y <"het_reads*.fastq.gz">  path/to/Heterozygote FASTQ files (e.g. Male reads)
-x <"hom_reads*.fastq.gz">  path/to/Homozygote FASTQ files (e.g. Female reads)
-t <target_genome.fasta>    path/to/Target Genome FASTA file
-q <20>                     Minimum Phred Score Quality Threshold
                            (it will be converted to ASCII)
-m <INT>                    Heterozygote k-mer Lower Count Threshold
-f <INT>                    Homozygote k-mer Lower Count Threshold
-g <2>                      Genome k-mer Lower Count Threshold
-v <20>                     VSCK Threshold: Number of valid single copy k-mers 
                            used for Y Chromosome GIs selection
-p <75>                     PVSCUK Threshold: Percentage of valid single copy
                            unmatched k-mers used for Y Chromosome GIs selection
-o <path/to/output/>        Path to output directory in which files will be
                            saved (Script creates the directory)

## Output:

logfile.txt                 Logfile with commands executed, date, outputs info,
                            and memory/CPU usage for major executions
k_mer_hash_tables.jelly     Jellyfish count Hash Table files created from the
                            input FASTQ and FASTA files (One for each: 
                            Heterozygote FASTQ, Homozygote FASTQ, and Target
                            Genome FASTA) NOTE: These are large files and will
                            be automatically deleted
input_counts.histo.txt      Jellyfish histo files created from Jellyfish count
                            files (One for each: Heterozygote count,
                            Homozygote count, and Target Genome count)
k_mer_dump.bitvector.gz     Bitvector of Jellyfish dump files created from
                            Jellyfish count files (One for each: Heterozygote
                            count, Homozygote count, and Target Genome count)
                            NOTE: These are large files and will be
                            automatically deleted
output.ygs                  YGS original output analysis file with metrics for
                            each sequence in <target_genome.fasta> and simple
                            execution report
output.tsv                  Tab-delimited file with the metrics present in
                            <output.ygs>
Inferred_ChrY_GIs_Info.txt  Sequence IDs inferred to the Y chromosome based on
                            YGS analysis <output.tsv> and the provided VSCK and
                            PVSCUK thresholds
ChrY.fasta                  FASTA file with the sequences inferred to the Y
                            chromosome based on YGS analysis <output.tsv> and
                            the provided VSCK and PVSCUK thresholds
jellyHisto.png              Histogram plot of Jellyfish histo files
YGS_ChrInferrence.png       Chromosome inferrence plot (Y vs X/Autosomes) based
                            on YGS analysis <output.tsv> and the provided VSCK
                            and PVSCUK thresholds

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

## Path to scripts
YGS=scripts/CarvalhoClark2013/YGS.pl
jfBitVector=scripts/CarvalhoClark2013/jelly_2_bitvector_mxky_1eg.pl
fastacmd_loose2=scripts/CarvalhoClark2013/fastacmd_loose2.awk
plot_YGS=scripts/plots/plot_YGS.py
plot_Jelly=scripts/plots/plotJellyHisto.py

#................................................
#  Input parameters from command line
#................................................

## Get parameters from command line flags
while getopts "k:y:x:t:q:m:f:g:v:p:o:" flag
do
    case "${flag}" in
        k) k_mer="${OPTARG}";;              # k-mer INT size
        y) het_fastq="${OPTARG}";;          # path/to/Heterozygote FASTQ files
        x) hom_fastq="${OPTARG}";;          # path/to/Homozygote FASTQ files
        t) tar_fasta=${OPTARG};;            # path/to/Target Genome FASTA file
        q) quality_threshold=${OPTARG};;    # Minimum Phred Score Quality Threshold
        m) het_lower_threshold=${OPTARG};;  # Heterozygote k-mer Lower Count Threshold
        f) hom_lower_threshold=${OPTARG};;  # Homozygote k-mer Lower Count Threshold
        g) tar_lower_threshold=${OPTARG};;  # Genome k-mer Lower Count Threshold
        v) vsck_threshold=$((${OPTARG}));;  # Number of valid single copy k-mers
        p) pvscuk_threshold=$((${OPTARG}));;# % valid single copy unmatched k-mers
        o) out_path=${OPTARG};;             # Path to output directory
        ?) echo script usage: bash ipda_scriptA_YGS.sh [-k 15] [-y "het_reads*.fastq.gz"] [-x "hom_reads*.fastq.gz"] [-t target_genome.fasta] [-q 20] [-m INT] [-f INT] [-g INT] [-v 20] [-p 75] [-o path/to/output/] >&2
           exit;;
    esac
done

#................................................
#  Extract info from input files
#................................................

## Get input files stem
het_stem=`echo "$(basename "${het_fastq%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`
hom_stem=`echo "$(basename "${hom_fastq%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`
tar_stem=`echo "$(basename "${tar_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')"`

## Get input files size
het_size=`ls -sh ${het_fastq} | cut -d" " -f1`
hom_size=`ls -sh ${hom_fastq} | cut -d" " -f1`
tar_size=`ls -sh ${tar_fastq} | cut -d" " -f1`

## Quality threshold (ASCII)
q=$((quality_threshold+33))
quality=`printf "\x$(printf %x ${q})"`

#................................................
#  Set output path and stem
#................................................

## Create the output directory
mkdir -p ${out_path}

## Get stem for outputs
out_stem=YGS_${tar_stem}_L${tar_lower_threshold}_${hom_stem}_q${quality_threshold}_L${hom_lower_threshold}_${het_stem}_q${quality_threshold}_L${het_lower_threshold}_k${k_mer}

#................................................
#  Set Logfile stem
#................................................

## Set Logfile stem
# it contains all the executed commands with date/time;
# the output files general metrics (such as size);
# and memory/CPU usage for major executions (i.e. Jellyfish count,
# Jellyfish histo, jelly_2_bitvector_mxky_1eg.pl, and YGS.pl)
logfile=${out_path}logfile_YGS_${out_stem}.txt

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing ipda_scriptA_YGS.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs files:"
echo
echo "## Heterozygote FASTQ files:                  ${het_fastq}"
echo "## Homozygote FASTQ files:                    ${hom_fastq}"
echo "## Target Genome FASTA file:                  ${tar_fasta}"
echo
echo "## Given input parameters/thresholds:"
echo
echo "## k-mer INT size:                            ${k_mer}"
echo "## Minimum Phred Score Quality Threshold:     ${quality_threshold} (${quality} in ASCII Table)"
echo "## Heterozygote k-mer Lower Count Threshold:  ${het_lower_threshold}"
echo "## Homozygote k-mer Lower Count Threshold:    ${hom_lower_threshold}"
echo "## Genome k-mer Lower Count Threshold:        ${tar_lower_threshold}"
echo "## Number of valid single copy k-mers:        ${vsck_threshold}"
echo "## % valid single copy unmatched k-mers:      ${pvscuk_threshold}"
echo
echo "## Output files saved to:                     ${out_path}"
echo
echo "## This is logfile:                           ${logfile}"
echo

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
echo "## Executing ipda_scriptA_YGS.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs files:"
echo
echo "## Heterozygote FASTQ files:                  ${het_fastq}"
echo "## Homozygote FASTQ files:                    ${hom_fastq}"
echo "## Target Genome FASTA file:                  ${tar_fasta}"
echo
echo "## Given input parameters/thresholds:"
echo
echo "## k-mer INT size:                            ${k_mer}"
echo "## Minimum Phred Score Quality Threshold:     ${quality_threshold} (${quality} in ASCII Table)"
echo "## Heterozygote k-mer Lower Count Threshold:  ${het_lower_threshold}"
echo "## Homozygote k-mer Lower Count Threshold:    ${hom_lower_threshold}"
echo "## Genome k-mer Lower Count Threshold:        ${tar_lower_threshold}"
echo "## Number of valid single copy k-mers:        ${vsck_threshold}"
echo "## % valid single copy unmatched k-mers:      ${pvscuk_threshold}"
echo
echo "## Output files saved to:                     ${out_path}"
echo
echo "## This is logfile:                           ${logfile}"
echo

set -v

#................................................
#  STEP 1 - GET QUALITY TRIMMED K-MER HASH TABLES
#................................................

date ## Starting Homozygote Hash Table (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
gunzip -c ${hom_fastq} | jellyfish count -t 6 -c 4 -C -s ${hom_size} -m ${k_mer} -o ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.jelly -Q "${quality}" /dev/fd/0

date ## Starting Homozygote Histogram Table (counts) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
jellyfish histo ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.jelly > ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.histo.txt

date ## Starting Homozygote Dump Bitvector (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
perl ${jfBitVector} m_jelly=${k_mer} kmer_size=${k_mer} jelly_file=${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.jelly lower-count=${hom_lower_threshold}
mv ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}m${k_mer}k${k_mer}rep${hom_lower_threshold}.vector_rep.gz ${out_path}${hom_stem}_q${quality_threshold}_L${hom_lower_threshold}_k${k_mer}.bitvector.gz

date ## Starting Heterozygote Hash Table (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
gunzip -c ${het_fastq} | jellyfish count -t 6 -c 4 -C -s ${het_size} -m ${k_mer} -o ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.jelly -Q "${quality}" /dev/fd/0

date ## Starting Heterozygote Histogram Table (counts) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
jellyfish histo ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.jelly > ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.histo.txt

date ## Starting Heterozygote Dump Bitvector (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
perl ${jfBitVector} m_jelly=${k_mer} kmer_size=${k_mer} jelly_file=${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.jelly lower-count=${het_lower_threshold}
mv ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}m${k_mer}k${k_mer}rep${het_lower_threshold}.vector_rep.gz ${out_path}${het_stem}_q${quality_threshold}_L${het_lower_threshold}_k${k_mer}.bitvector.gz

date ## Starting Target Genome Hash Table (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
jellyfish count -t 6 -c 4 -C -s ${tar_size} -m ${k_mer} -o ${out_path}${tar_stem}_k${k_mer}.jelly ${tar_fasta}

date ## Starting Target Genome Histogram Table (counts) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
jellyfish histo ${out_path}${tar_stem}_k${k_mer}.jelly > ${out_path}${tar_stem}_k${k_mer}.histo.txt

date ## Starting Target Genome Dump Bitvector (k-mers and count) at
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
perl ${jfBitVector} m_jelly=${k_mer} kmer_size=${k_mer} jelly_file=${out_path}${tar_stem}_k${k_mer}.jelly lower-count=${tar_lower_threshold}
mv ${out_path}${tar_stem}_k${k_mer}m${k_mer}k${k_mer}rep${tar_lower_threshold}.vector_rep.gz ${out_path}${tar_stem}_L${tar_lower_threshold}_k${k_mer}.bitvector.gz

#................................................
#  STEP 2 - Y GENOME SCAN
#................................................

date ## Starting Y Genome Scan at
ln -s ${tar_fasta} .
ln -s ${out_path}${tar_stem}_L${tar_lower_threshold}_k${k_mer}.bitvector.gz .
ln -s ${out_path}${hom_stem}_q${quality_threshold}_L${hom_lower_threshold}_k${k_mer}.bitvector.gz .
ln -s ${out_path}${het_stem}_q${quality_threshold}_L${het_lower_threshold}_k${k_mer}.bitvector.gz .
/usr/bin/time -f "# Process TIME: %e || MRES: %M || CPUSK: %S || CPUSU: %U" \
perl ${YGS} kmer_size=${k_mer} mode=final_run contig=${tar_fasta} gen_rep=${tar_stem}_k${k_mer}_L${tar_lower_threshold}.bitvector.gz trace=${hom_stem}_q${quality_threshold}_k${k_mer}_L${hom_lower_threshold}.bitvector.gz male_trace=${het_stem}_q${quality_threshold}_k${k_mer}_L${het_lower_threshold}.bitvector.gz
mv ${tar_stem}_${hom_stem}_q${quality_threshold}_k${k_mer}_L${hom_lower_threshold}_${het_stem}_q${quality_threshold}_k${k_mer}_L${het_lower_threshold}.final_result ${out_path}${out_stem}_report.ygs

#................................................
# STEP 3 - CREATE TAB-DELIMITED YGS OUTPUT
#................................................

date ## Starting to get tab separated YGS metrics at
grep ">" ${out_path}${out_stem}_report.ygs | sed 's/ \+/\t/g' > ${out_path}${out_stem}_metrics.tsv
sed -i "1i\GI\tNUM\tMAX_K\tK\tUK\tSC_K\tSC_UK\tP_SC_UK\tVSC_K\tVSC_UK\tP_VSC_UK" ${out_path}${out_stem}_metrics.tsv

#................................................
# STEP 4 - IDENTIFY Y CHR SUSPECTS
#................................................

date ## Starting GIs selection at
cut -f1,9,11 ${out_path}${out_stem}_metrics.tsv | awk -v vsck="${vsck_threshold}" -v pvscuk="${pvscuk_threshold}" -F'\t' '{if($2>=vsck && $3>=pvscuk)print$0}' > ${out_path}tmpfile1
sed '1d' ${out_path}tmpfile1 > tmpfile; mv tmpfile ${out_path}tmpfile1
cut -f1 ${out_path}tmpfile1 | cut -d"|" -f2 | while read gi; do grep -w ">gi|$gi" ${tar_fasta} ; done | cut -d ">" -f2 > ${out_path}tmpfile2
paste ${out_path}tmpfile2 ${out_path}tmpfile1 > ${out_path}${out_stem}_vsck${vsck_threshold}_pvscuk${pvscuk_threshold}_Inferred_ChrY_GIs_Info.txt ; rm ${out_path}tmpfile2

date ## Starting to create FASTA file at
cut -f1 ${out_path}tmpfile1 | cut -d"|" -f2 > ${out_path}tmpfile3
${fastacmd_loose2} gi_type="gi" fasta=${tar_fasta} ${out_path}tmpfile3 > ${out_path}${out_stem}_vsck${vsck_threshold}_pvscuk${pvscuk_threshold}_Inferred_ChrY.fasta
rm ${out_path}tmpfile1 ${out_path}tmpfile3

#file=${out_path}tmpfile3
#awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' ${tar_fasta} | awk -v gisFile="${file}" -F"\t" 'BEGIN{while((getline k < gisFile)>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > ${out_path}${out_stem}_vsck${vsck_threshold}_pvscuk${pvscuk_threshold}_Inferred_ChrY.fasta ; rm ${out_path}tmpfile1 ${out_path}tmpfile3 

#................................................
#  STEP 5 - PLOT RESULTS
#................................................

date ## Starting to plot YGS results at
python3 ${plot_YGS} -f ${out_path}${out_stem}_metrics.tsv -vsck ${vsck_threshold} -pvscuk ${pvscuk_threshold} -o ${out_path}${out_stem}_vsck${vsck_threshold}_pvscuk${pvscuk_threshold}_ChrInferrence.png

date ## Starting to plot Jellyfish Histo at
python3 ${plot_Jelly} -hom ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.histo.txt -het ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.histo.txt -gen ${out_path}${tar_stem}_k${k_mer}.histo.txt -homT "${hom_stem}" -hetT "${het_stem}" -genT "${tar_stem}" -o ${out_path}jellyHisto_k${k_mer}_${tar_stem}_${hom_stem}_q${quality_threshold}_${het_stem}_q${quality_threshold}.png

rm ${out_path}${hom_stem}_q${quality_threshold}_k${k_mer}.jelly ${out_path}${het_stem}_q${quality_threshold}_k${k_mer}.jelly ${out_path}${tar_stem}_k${k_mer}.jelly
rm ${tar_fasta} ${tar_stem}_k${k_mer}_L${tar_lower_threshold}.bitvector.gz ${hom_stem}_q${quality_threshold}_k${k_mer}_L${hom_lower_threshold}.bitvector.gz ${het_stem}_q${quality_threshold}_k${k_mer}_L${het_lower_threshold}.bitvector.gz
rm ${out_path}${tar_stem}_L${tar_lower_threshold}_k${k_mer}.bitvector.gz ${out_path}${hom_stem}_q${quality_threshold}_L${hom_lower_threshold}_k${k_mer}.bitvector.gz ${out_path}${het_stem}_q${quality_threshold}_L${het_lower_threshold}_k${k_mer}.bitvector.gz

# This will remove $VARNAMES from output file with the actual $VARVALUE
# allowing for easily retracing commands
sed -i 's,${hom_fastq},'"${hom_fastq}"',g' "$logfile"
sed -i 's,${hom_size},'"${hom_size}"',g' "$logfile"
sed -i 's,${hom_stem},'"${hom_stem}"',g' "$logfile"
sed -i 's,${hom_lower_threshold},'"${hom_lower_threshold}"',g' "$logfile"
sed -i 's,${het_fastq},'"${het_fastq}"',g' "$logfile"
sed -i 's,${het_size},'"${het_size}"',g' "$logfile"
sed -i 's,${het_stem},'"${het_stem}"',g' "$logfile"
sed -i 's,${het_lower_threshold},'"${het_lower_threshold}"',g' "$logfile"
sed -i 's,${tar_fasta},'"${tar_fasta}"',g' "$logfile"
sed -i 's,${tar_size},'"${tar_size}"',g' "$logfile"
sed -i 's,${tar_stem},'"${tar_stem}"',g' "$logfile"
sed -i 's,${tar_lower_threshold},'"${tar_lower_threshold}"',g' "$logfile"
sed -i 's,${k_mer},'"${k_mer}"',g' "$logfile"
sed -i 's,${vsck_threshold},'"${vsck_threshold}"',g' "$logfile"
sed -i 's,${pvscuk_threshold},'"${pvscuk_threshold}"',g' "$logfile"
sed -i 's,${quality_threshold},'"${quality_threshold}"',g' "$logfile"
sed -i 's,((quality_threshold,'"((${quality_threshold}"',g' "$logfile"
sed -i 's,${quality},'"${quality}"',g' "$logfile"
sed -i 's,${out_path},'"${out_path}"',g' "$logfile"
sed -i 's,${YGS},'"${YGS}"',g' "$logfile"
sed -i 's,${jfBitVector},'"${jfBitVector}"',g' "$logfile"
sed -i 's,${fastacmd_loose2},'"${fastacmd_loose2}"',g' "$logfile"
sed -i 's,${plot_YGS},'"${plot_YGS}"',g' "$logfile"
sed -i 's,${plot_Jelly},'"${plot_Jelly}"',g' "$logfile"
sed -i 's,${log_file},'"${log_file}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v