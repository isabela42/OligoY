#!/usr/bin/awk -f
#fastacmd_loose2.awk   Bernardo  3 Aug 2014 (Cao faz 85 nos)   for Bernardo Lemos project, and general use. v. 29ago2019. 
# modified to add the option gi_type="noregex" to avoid the regex match that may cause errors (e.g., when the gene name starts with a number, but has letters) 
#29ago2019: added PROCINFO["sorted_in"] = "@ind_str_asc" so the fasta files came in a order (useful for interleaved sequnces)
# 4jan2019: added partial_match=PB_well_ID to allow the match between PB well__ID and the full PB defline (each well_ID may contain several subreads)
# based on fastacmd_loose.awk  , but allowing the "mode=reverse" (I need to change the pgm significantly).  
#As fastacmd_loose.awk, the  gi list can be piped, which is usefull for small things
#usage:
#fastacmd_loose2.awk     fasta=CS_CONGO_4361F.dump100c.fasta   mode=reverse  CS_CONGO_4361F.dump100c.with_bls_hit 
# echo "fdy_PB" | fastacmd_loose2.awk  fasta=dmel-all-CDS-r6.03.edited2.cdhit100.fasta mode=reverse gi_type=noregex 



function load_fasta(fasta){
	#print "loading fasta file  " fasta
	while ( (getline  < fasta) > 0 ) {	#loads the fasta file
			if ($1 ~ /^>/){
				  seq ++
				  gi = gensub(/^>/, "" , "g" , $1) ;  #		most general option: simply gets the first word of the defline and remove the ">"   (used with bernardo lemos.)
				  if (gi_type != "noregex"){  #new in 17abr15
						if ( match($1,/>gi\|([0-9]+)/,a) ){gi = a[1]}    #defline of the type  >gi|1226338|
						if ( match($1,/>([0-9_]+[rc]*)/,a) ){gi = a[1]}    #defline of the type  >1226338  or 2310_45_330rc
						}
                  sequence[gi]=$0 "\n"
				  }
			else{sequence[gi] = sequence[gi] $0}
			}
	close(fasta)
return seq}


function usage(dummy){
	print "usage:"
	print "fastacmd_loose2.awk     fasta=CS_CONGO_4361F.dump100c.fasta   mode=reverse  CS_CONGO_4361F.dump100c.with_bls_hit" 
    print "echo \42fdy_PB\42 | fastacmd_loose2.awk  fasta=dmel-all-CDS-r6.03.edited2.cdhit100.fasta mode=reverse gi_type=noregex "   # \42 is "
	print "fastacmd_loose2.awk  gi_type=\42noregex\42   fasta=dmel-all-CDS-r6.03.edited.fasta    dmel-all-CDS-r6.03.longestCDS.gi > dmel-all-CDS-r6.03.longestCDS.fasta"
	print "fastacmd_loose2.awk  partial_match=PB_well_ID"  
	print "partial_match=PB_well_ID to allow the match between PB well__ID and the full PB defline (each well_ID may contain several subreads)"
	print "use gi_type=\42noregex\42 to avoid the regex match that may cause errors (e.g., when the gene name starts with a number, but has letters)"
	flag_exit=1; exit
	}

BEGIN{
	seq=0
	if (ARGC==1){usage(dummy)}
	PROCINFO["sorted_in"] = "@ind_str_asc"
	}

(1==1){array_gi[$1]=""} # I only need the index


END{
	if (flag_exit){exit}
    num_seq=load_fasta(fasta)
	if (partial_match==""){
			if (mode=="reverse"){
			  for (gi in sequence){
				  if (!(gi in array_gi)){print sequence[gi]}
			  }
			}
			if (mode!="reverse"){
			  for (gi in sequence){
				  if (gi in array_gi){print sequence[gi]}
			  }
			}
	}
	if (partial_match=="PB_well_ID"){
			if (mode=="reverse"){
			  for (gi in sequence){
				  gi_well_ID =  gensub(/\/[0-9]+_[0-9]+$/,"","g",gi)
				  if (!(gi_well_ID in array_gi)){print sequence[gi]}
			  }
			}
			if (mode!="reverse"){
			  for (gi in sequence){
				  gi_well_ID =  gensub(/\/[0-9]+_[0-9]+$/,"","g",gi)
				  if (gi_well_ID in array_gi){print sequence[gi]}
			  }
			}
	}	
	
}


