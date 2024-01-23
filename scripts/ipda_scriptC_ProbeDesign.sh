#!/bin/bash

usage(){
echo '
Written by Isabela Almeida
Last modified January 23, 2024

Description: Runs OligoMiner and remaining OligoY steps for designing and
filtering probes

Produces temporary files and outputs Y chromosome probes file.

Usage: bash ipda_scriptC_ProbeDesign.sh [-t target_chr.fasta] [-r reference_genome.fasta] [-c complete_tar.fasta] [-x "hom_reads*.fastq.gz"] [-o path/to/output/]

## Input:

-t <target_chr.fasta>       path/to/target/singleline.fasta
-r <reference_genome.fasta> path/to/Reference Genomes FASTA file
-x <"hom_reads*.fastq.gz">  path/to/Homozygote FASTQ files (e.g. Female reads)

-c <coverage.cov>           path/to/Coverage COV file. If no coverage info is
                            available, create a coverage file with
                            grep ">" ${tar_fasta} | cut -d" " -f1 | cut -d">" -f2 | while read gi; do echo -e "$gi\t100\t300" >> coverage.cov ; done
-o <path/to/output/>        Path to output directory in which files will be
                            saved (Script creates the directory)

## Output:

logfile.txt                 Logfile with commands executed, date, outputs info,
                            and memory/CPU usage for major executions

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
bP_script=Beliveau2018/blockParse.py # blockParse.py
oC_script=Beliveau2018/outputClean.py # outputClean.py
kF_script=Beliveau2018/kmerFilter.py # kmerFilter.py
sC_script=Beliveau2018/structureCheckpy3.py # structureCheck.py 2to3 to run on python 3.6
nO_script=scripts/noOverlap.py # Isabela Almeida's script to remove probes that overlap for up to n positions to the previous one
nR_script=scripts/nonRedundant.py # Isabela Almeida's script to return the first report of each probe
bedFastq=scripts/bedFastq.py # Isabela Almeida's script to return a .fastq file from .bed

#................................................
#  Input parameters from command line
#................................................

## Get parameters from command line flags
while getopts "t:r:c:x:o:" flag
do
    case "${flag}" in
        t) tar_fasta=${OPTARG};;            # path/to/target/singleline.fasta
        r) ref_fasta=${OPTARG};;            # path/to/Reference Genome FASTA file
        c) complete_tar_fasta=${OPTARG};;   # path/to/target/multiline.fasta
        x) hom_fastq="${OPTARG}";;          # path/to/Homozygote FASTQ files
        o) out_path=${OPTARG};;             # Path to output directory
        ?) echo script usage: bash ipda_scriptC_ProbeDesign.sh [-t target_chr.fasta] [-r reference_genome.fasta] [-c complete_tar.fasta] [-x "hom_reads*.fastq.gz"] [-o path/to/output/] >&2
           exit;;
    esac
done

#................................................
#  Set parameters
#................................................
kmer=18
k_threshold=5
min_len=40
max_len=46
min_TM=47
max_TM=52
hybrid_T=47
bt2parameters="--no-hd -t -k 2 --very-sensitive-local" # change bowtie2 NGS alignment parameters if desired

## Get input files stem
hom_stem=`echo "$(basename "${hom_fastq%%.*}" | sed 's/\(.*\)\..*/\1/')" | sed 's/\*//g'`
ref_stem=`echo "$(basename "${ref_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')"`
tar_stem=`echo "$(basename "${complete_tar_fasta%%.*}" | sed 's/\(.*\)\..*/\1/')"`


#................................................
#  Set output path and stem
#................................................

## Create the output directory
mkdir -p ${out_path}

## Get stem for outputs
out_stem=ProbeDesign_${tar_stem}k${kmer}_kThres${k_threshold}_lengh-${min_len}-${max_len}_TM-${min_TM}-${max_TM}_hybT${hybrid_T}.txt


#................................................
#  Set Logfile stem
#................................................

## Set Logfile stem
# it contains all the executed commands with date/time;
# and memory/CPU usage for major executions (i.e. OligoMiner, Bowtie2)
logfile=${out_path}logfile_ProbeDesign_${out_stem}.txt

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing ipda_scriptC_ProbeDesign.sh"
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
echo "## Executing ipda_scriptC_ProbeDesign.sh"
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
#  STEP 1 - GET HOMOZYGOTE INDEX FILES
#................................................

date ## Starting Homozygote bowtie2 Index at
seqtk seq -a ${hom_fastq} > ${out_path}${hom_stem}.fasta
bowtie2-build ${out_path}${hom_stem}.fasta ${out_path}${hom_stem} -q

#................................................
#  STEP 2 - GET REFERENCE GENOME WITHOUT Y CHROMOSOME (TARGET)
#................................................

date ## Starting to create No Y Chr (Target) Genome GIs list at
grep ">" ${complete_tar_fasta} | cut -d" " -f1 | cut -d"|" -f2 > ${out_path}tmpfile1_${tar_stem}_GIs
sed 's/^/>gi|/' ${out_path}tmpfile1_${tar_stem}_GIs > ${out_path}tmpfile2_${tar_stem}_GIs
sed 's/$/ /' ${out_path}tmpfile2_${tar_stem}_GIs > ${out_path}tmpfile3 ; mv ${out_path}tmpfile3 ${out_path}tmpfile2_${tar_stem}_GIs # To avoid the following e.g.: removes not only >gi|1 but also wrongly removes >gi|12
grep ">" ${ref_fasta} > ${out_path}${ref_stem}_GIs
grep -vf ${out_path}tmpfile2_${tar_stem}_GIs ${out_path}${ref_stem}_GIs | cut -d " " -f1 | cut -d"|" -f2 > ${out_path}${tar_stem}_GIs
rm ${out_path}tmpfile2_${tar_stem}_GIs ${out_path}tmpfile1_${tar_stem}_GIs

date ## Starting to create No Y Chr (Target) Genome fasta & recover its size at
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' ${ref_fasta} | awk -v gisFile="${out_path}${tar_stem}_GIs" -F"\t" 'BEGIN{while((getline k < gisFile)>0)i[k]=1}{gsub("^>","",$0); if(i[$1]){print ">"$1"\n"$2}}' > ${out_path}NoYgenome.fasta
NY_size=`ls -lh ${out_path}NoYgenome.fasta | cut -d " " -f5`

#................................................
# STEP 3 - GET HASH TABLE & INDEX FILES
#................................................

date ## Starting No Y Chr (Target) Genome bowtie2 Index at
bowtie2-build ${out_path}NoYgenome.fasta ${out_path}NoYgenome -q

date ## Starting No Y Chr (Target) Genome Hash Table (kmers and count) at
jellyfish count -s ${NY_size} -m ${kmer} -o ${out_path}NoYgenome_k${kmer}.jelly --out-counter-len 1 -L 2 ${out_path}NoYgenome.fasta
jellyfish histo ${out_path}NoYgenome_k${kmer}.jelly > ${out_path}NoYgenome_k${kmer}.histo.txt

#................................................
# STEP 4 - OLIGOMINER - DESIGN AND FILTER CANDIDATE PROBES
#................................................

date ## Starting to design candidate probes at
python2.7 ${bP_script} -l ${min_len} -L ${max_len} -t ${min_TM} -T ${max_TM} -f ${tar_fasta} -o ${out_path}${tar_stem}_Candidates_bP_overlap_lengh-${min_len}-${max_len}_TM-${min_TM}-${max_TM} -O

date ## Starting to align candidate probes against No Y Genome bowtie2 Index at
bowtie2 -x ${out_path}NoYgenome -U ${out_path}${tar_stem}_Candidates_bP_overlap_lengh-${min_len}-${max_len}_TM-${min_TM}-${max_TM}.fastq ${bt2parameters} -S ${out_path}${tar_stem}_bowtie2_NoYgenomeAlignment_Candidates_bP_overlap_lengh-${min_len}-${max_len}_TM-${min_TM}-${max_TM}.sam

date ## Starting Filter Level 1 - filter candidate probes with zero predicted NGS Alignments at
python2.7 ${oC_script} -0 -f ${out_path}${tar_stem}_bowtie2_NoYgenomeAlignment_Candidates_bP_overlap_lengh-${min_len}-${max_len}_TM-${min_TM}-${max_TM}.sam -o ${out_path}${tar_stem}_F1_oC-ZM

date ## Starting Filter Level 2 - screen for high abundance k-mers at
python2.7 ${kF_script} -f ${out_path}${tar_stem}_F1_oC-ZM.bed -m ${kmer} -j ${out_path}NoYgenome_k${kmer}.jelly -k ${k_threshold} -o ${out_path}${tar_stem}_F2_kF_k${kmer}_kThres${k_threshold}
		
date ## Starting Filter Level 3 - screen for probes with probability of forming secondary structures at
python3.6 ${sC_script} -T ${hybrid_T} -f ${out_path}${tar_stem}_F2_kF_k${kmer}_kThres${k_threshold}.bed -o ${out_path}${tar_stem}_F3_sC_hybT${hybrid_T}

#................................................
#  STEP 5 - REMOVE REDUNDANT PROBES
#................................................

date ## Starting Filter Level 4 - remove overlapping probes at
python ${nO_script} -n ${max_len} ${out_path}${tar_stem}_F3_sC_hybT${hybrid_T}.bed

date ## Starting Filter Level 5 - select first appearence/hit of each probe at
python ${nR_script} ${out_path}${tar_stem}_F3_sC_hybT${hybrid_T}_F4_nO.bed

#................................................
#  STEP 6 - FINAL FILTER AGAINST HOMOZYGOUS READS
#................................................

date ## Starting to convert BED to FASTq file at
python ${bedFastq} ${out_path}${tar_stem}_F3_sC_hybT${hybrid_T}_F4_nO_F5_nR.bed

date ## Starting to align filtered probes against Homozygote bowtie2 Index at
bowtie2 -x ${out_path}${hom_stem} -U ${out_path}${tar_stem}_F3_sC_hybT${hybrid_T}_F4_nO_F5_nR.fastq ${bt2parameters} -S ${out_path}${tar_stem}_bowtie2_${hom_stem}Alignment.sam

date ## Starting Filter Level 6 - select final probes with zero Homozygote predicted NGS Alignments at
python2.7 ${oC_script} -0 -f ${out_path}${tar_stem}_bowtie2_${hom_stem}Alignment.sam -o ${out_path}${tar_stem}_F6_oC-ZM
rm ${out_path}*bt2

#################################################################
#           OLIGOMINER DEFAULT OPTIONS AND EXTRA INFO           #
#################################################################

################################################
##     bP::: blockParse.py Default options    ##
################################################

# bP_file=/PATH/TO/file.fasta #The FASTA file to find probes in
# min_len=36 #The minimum allowed probe length; default is 36
# max_len=41 #The maximum allowed probe length, default is 41
# minGC=20 #The minimum allowed percent G + C, default is 20
# maxGC=80 #The maximum allowed percent G + C, default is 80
# min_TM=42 #The minimum allowed Tm, default is 42
# max_TM=47 #The maximum allowed Tm, default is 47
# prohibitedSeqs='AAAAA,TTTTT,CCCCC,GGGGG' #Prohibited sequence list (separated by commas with no spaces), default is 'AAAAA,TTTTT,CCCCC,GGGGG'
# salt=390 #The mM Na+ concentration, default is 390
# formamide=50 #The percent formamide being used, default is 50
# spacing=0 #The minimum spacing between adjacent probes, default is 0 bases
# dnaC1=25 #Concentration of higher concentration strand [nM] -typically the probe- to use for thermodynamic calculations. Default is 25
# dnaC2=25 #Concentration of lower concentration strand [nM]-typically the target- to use for thermodynamic calculations. Default is 25
# header='' #Allows the use of a custom header in the format chr start-stop. E.g. 'chr2:12500-13500'
# output=/PATH/outputName #Specify the stem of the output filename // include path/to/file

# Usage:
# python blockParse.py -l $min_len -L $max_len -g $minGC -G $maxGC -t $min_TM -T $max_TM -X $prohibitedSeqs -s $salt -F $formamide -S $spacing -c $dnaC1 -C $dnaC2 -f $bP_file -o $output
# python blockParse.py -f $bP_file #same as the datailed one above

################################################
##      bP::: blockParse.py Extra options     ##
################################################

####### For overlapping probes: add -O #Turn on Overlap Mode, which returns all possible candidate probes in a block of sequence including overlaps. Off by default. Note, if selecting this option, the -S/--Spacing value will be ignored
####### For mining progress: add -V #Turn on verbose mode to have probe miningprogress print to Terminal. Off by default
####### For writing a Report file: add -R #Write a Report file detailing the results of each window of sequence considered by the script. The first set of lines give the occurrence of each possible failure mode for quick reference. Off by default. Note, selecting this option will slow the script considerably
####### For writing debung file: add -D #The same as -Report, but prints info to terminal instead of writing a log file. Off by default
####### For writing meta file: add -M #Write a text file containing meta information Off by default. Reports input file <tab> estimated runtime <tab> blockParse version <tab> candidates discovered <tab> span in kb covered by candidate probes <tab> candidate probes per kb

################################################
##    oC::: outputClean.py Default options    ##
################################################

# oc_file=/PATH/TO/file.sam #The sam file to be processed
# prob_threshold=0.5 #The probability threshold for classifying a candidate sequence as likely to have off-target binding using the LDA model. Default=0.5. Selecting smaller values will improve precision (fewer false positives), but at the expense of recall (more false negatives). Selecting larger values will improve recall at the expense of precision.
# hybrid_T=42 #Specify the temperature-specific linear discrimination model to use in LDM. Options are 32, 37, 42, 47, 52, 57. Default=42
# salt=390 #The mM Na+ concentration to be used for Tm calculation, default is 390
# formamide=50 #The percent formamide being used, default is 50
# output=/PATH/outputName #Specify the stem of the output filename // include path/to/file

# Usage:
# python outputClean.py -p $prob_threshold -T $hybrid_T -s $salt -F $formamide -f $oc_file
# python outputClean.py -f $oc_file  #same as the datailed one above

################################################
##     oC::: outputClean.py Extra options     ##
################################################

####### For using LDA Model: add -l #On by default
####### For using Unique Mode: add -u #Only return probes aligning exactly one time. Does not use the LDA model. Off by default
####### For using 'Zero Mode': add -0 # Only return probes aligning zero times. Does not use the LDA model. Off by default.
####### For writing a Report file: add -R #Write a Report file detailing the results of .sam cleaning. Off by default. Note, selecting this option will slow the script.
####### For writing debung file: add -D #The same as -Report, but prints info to terminal instead of writing a log file. Off by default
####### For writing meta file: add -M #Write a text file containing meta information. Off by default. Reports input file <tab> estimated runtime <tab> outputClean version <tab> unique probes identified <tab> number of candidate probes inputted

################################################
##     kF::: kmerFilter.py Default options    ##
################################################

# kf_file=/PATH/TO/file.bed #The probe file to do kmer tallying in
# jf_file=/PATH/TO/dictionary.jf #The Jellyfish .jf file to use
# merLength=18 # The length of kmer used to build the .jf being used, default=18
# k_threshold=5 # Filter probes with kmers occurring => than this many times, default=5
# kf_output='' #The output name prefix
# id='' #Specify an ID to be associated with temporary file names & the default output name. Can be useful for parallelization. Null by default. Will not be automatically included in output file name if -o is flagged

# Usage:
# python kmerFilter.py -f $kf_file -j $jf_file -m $merLength -k $k_threshold -o $kf_output -I $id

################################################
##     kF::: kmerFilter.py Extra options      ##
################################################   

####### For writing a Report file: add -R #Write a Report file detailing the results of the kmer filtering. Off by default. Note, selecting this option will slow the script.
####### For writing debung file: add -D #The same as -Report, but prints info to terminal instead of writing a log file. Off by default
####### For writing meta file: add -M #Write a text file containing meta information. Off by default.Reports input file <tab> estimated runtime <tab> kmerFilter version <tab> kmer occurrence threshold used <tab> kmer length used <tab> number of input probes passing kmer filter <tab> number of input probes

################################################
##   sC:::structureCheck.py Default options   ##
################################################

# sC_file=/PATH/TO/file.bed #The .bed file to check probes in
# formamide=50 #The percent formamide being used, default is 50
# salt=390 #The mM Na+ concentration to be used for Tm calculation, default is 390. NOTE: NUPACK's allowable range is 50-1100 mM
# material=dna1998 #The NUPACK material setting, default is dna1998
# threshold=0.1 #The probability threshold for a probe having a linear structure to use for filtering, default=0.1
# hybrid_T=47 #The temperature at which you want to hybridize your probes. Default=47
# id='' #Specify an ID to be associated with temporary file names & the default output name. Can be useful for parallelization. Null by default. Will not be automatically included in output file name if -o is flagged
# output=/PATH/outputName #Specify the stem of the output filename // include path/to/file
 
################################################
##    sC:::structureCheck.py Extra options    ##
################################################ 

####### For writing a Report file: add -R #Write a Report file detailing the results of the secondary structure filtering. Off by default. Note, selecting this option will slow the script.
####### For writing debung file: add -D #The same as -Report, but prints info to terminal instead of writing a log file. Off by default
####### For writing meta file: add -M #Write a text file containing meta information. Off by default. Reports input file <tab> estimated runtime <tab> structureCheck version ]<tab> hybridization temperature used <tab> number of input probes passing secondary structure filter <tab> number of input probes

# This will remove $VARNAMES from output file with the actual $VARVALUE
# allowing for easily retracing commands
sed -i 's,${tar_fasta},'"${tar_fasta}"',g' "$logfile"
sed -i 's,${ref_fasta},'"${ref_fasta}"',g' "$logfile"
sed -i 's,${complete_tar_fasta},'"${complete_tar_fasta}"',g' "$logfile"
sed -i 's,${hom_fastq},'"${hom_fastq}"',g' "$logfile"
sed -i 's,${out_path},'"${out_path}"',g' "$logfile"
sed -i 's,${bP_script},'"${bP_script}"',g' "$logfile"
sed -i 's,${oC_script},'"${oC_script}"',g' "$logfile"
sed -i 's,${kF_script},'"${kF_script}"',g' "$logfile"
sed -i 's,${sC_script},'"${sC_script}"',g' "$logfile"
sed -i 's,${nO_script},'"${nO_script}"',g' "$logfile"
sed -i 's,${nR_script},'"${nR_script}"',g' "$logfile"
sed -i 's,${bedFastq},'"${bedFastq}"',g' "$logfile"
sed -i 's,${kmer},'"${kmer}"',g' "$logfile"
sed -i 's,${k_threshold},'"${k_threshold}"',g' "$logfile"
sed -i 's,${min_len},'"${min_len}"',g' "$logfile"
sed -i 's,${max_len},'"${max_len}"',g' "$logfile"
sed -i 's,${min_TM},'"${min_TM}"',g' "$logfile"
sed -i 's,${max_TM},'"${max_TM}"',g' "$logfile"
sed -i 's,${hybrid_T},'"${hybrid_T}"',g' "$logfile"
sed -i 's,${bt2parameters},'"${bt2parameters}"',g' "$logfile"
sed -i 's,${hom_stem},'"${hom_stem}"',g' "$logfile"
sed -i 's,${ref_stem},'"${ref_stem}"',g' "$logfile"
sed -i 's,${tar_stem},'"${tar_stem}"',g' "$logfile"
sed -i 's,${log_file},'"${log_file}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v