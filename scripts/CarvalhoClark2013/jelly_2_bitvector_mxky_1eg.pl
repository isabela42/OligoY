
    #   jelly_2_bitvector_mxky_1d.pl 
#safe pipe procedures based in http://www.xinotes.org/notes/note/619/    and   http://perldoc.perl.org/perlipc.html#Complete-Dissociation-of-Child-from-Parent
#based on jelly_2_bitvector_mxkx_1b.pl which uses a fixed conversion, and jelly_2_bitvector_k15.pl (variable conversion). This program allows conversion of a jelly file to a smaller k-mer bit-vector, specified by the user. 
#new in 1b:  added --lower-count=2 (to increase speed) , since I am usually saving only the kmers that have at least two copies
#new in 1c: added mode=save_all  This saves all kmers (and not only those with more than one copy). Useful to convert genomic scaffolds to traces (e.g., YNCBI). 
#new in 1d: variable $lower_count_otpion replaced by $lower_count_value; vector file name modified to reflect $lower_count_value (e.g., rep12). Note that --lower-count acts at jellyfish output, not at the save step
#new in 1e: major streamlining of the code, but with full  backward compatibility.   lower-count parameter modified tos "lower-count=" for coherence with other options ( --lower-count prerserved for backward compatibility); convert_mxkx_rep and convert_mxkx_all combined in convert_mxkx. Same for convert_mxky_rep and convert_mxky_all (combined in convert_mxky)
#new in 1eg: done by Isabela Almeida, using 1e code, but changing it to run with jellyfish2 (as in 1g) and also change it to accept the path/to/file rather than just the file names
#usage perl jelly_2_bitvector_mxky_1d.pl m_jelly=18 kmer_size=18 jelly_file=GBR200m_m18_f10.jelly          #or
#      perl jelly_2_bitvector_mxky_1d.pl m_jelly=18 kmer_size=18 jelly_file=GBR200m_m18_f10.jelly --lower-count=2   #Both will save only kmers that occur more than once
#      perl jelly_2_bitvector_mxky_1d.pl m_jelly=15 kmer_size=15 jelly_file=fvirf30m15.jelly --lower-count=12 
#OR    perl jelly_2_bitvector_mxky_1d.pl m_jelly=18 kmer_size=18 jelly_file=GBR200m_m18_f10.jelly mode=save_all


use warnings;
no warnings 'portable';   #avoids a huge number of warnings in the error file (at UW cluster): Binary number > 0b11111111111111111111111111111111 non-portable
use Bit::Vector;     # at UW run in perl-5.14.2 , whose BitVector7.1 was modified (Toolbox.h) 	at DV, runs in perl-5.16.0, whose BitVector7.2 was modified
use Carp::Clan ;
use PerlIO::gzip;

                                                                $program_version = "v_1e "."29 nov 2012 11:40AM";                                         
												

process_command_line("dummy");
print "\n\nstarted:",scalar localtime,    "       program version:  $program_version \n"  ;
#$program =  "/net/gs/vol1/home/abc32/bin/jellyfish";     #for UW cluster: needs full identification of the jellyfish path
$program =  "jellyfish";                                 #for Darthvader
$command = "dump" ; 
if ($lower_count_value eq "not_set"){ $lower_count_value = 2} #changed in mxky_1e  to simplify the code. Keeps the default of saving kmers  with 2 or more copies if lower-count is not specified
if ($mode eq "save_all"){ $lower_count_value = 1}
@options = ("--column","--lower-count=" . $lower_count_value);				
print "command line: perl $0 @ARGV\n";
#print "jellyfish command:  $program  $command  @options  $jelly_full_name \n";

$bits =  4**$kmer_size  ;  #print "\n kmer size used: $kmer_size\n"; # $kmer_size values above 16 require a modified Bit::Vector module, and possible a perl compiled with -Duse64bitint    perlbrew install  perl-5.14.2 -Duse64bitint       
if ($m_jelly  < $kmer_size){die "\n output bit-vector kmer_size should be equal or smaller jelly file. User declared values are  kmer_size:  $kmer_size    m_jelly: $m_jelly\n"}
if ($kmer_size > 16){if ($kmer_size > 19){die "\nexcessive kmer_size: $kmer_size. Normal values between 15 and 19\n"}
                    check_64bits_compatibility($kmer_size);
                    };   
$jelly_kmer_size = get_jelly_kmer_size($jelly_full_name); 
if ($jelly_kmer_size ne $m_jelly ){die "\ndifferent kmer sizes in jelly file ($jelly_kmer_size) and user declared kmer size ($m_jelly).\n"}
else{print "kmer size found in jelly file $jelly_full_name is  $jelly_kmer_size and agrees with user declared kmer size ($m_jelly). \n"} ;

$vector_output = Bit::Vector->new($bits);
if ($lower_count_value == 1) {$out_vector_name = $jelly_prefix . "m" . $m_jelly . "k" . $kmer_size . ".vector"} ;
if ($lower_count_value >  1) {$out_vector_name = $jelly_prefix . "m" . $m_jelly . "k" . $kmer_size . "rep" . $lower_count_value .  ".vector_rep"};
if ($m_jelly == $kmer_size){convert_mxkx("dummy")};
if ($m_jelly  > $kmer_size){convert_mxky("dummy")};
store_vector($vector_output,$out_vector_name);

if ($m_jelly == $kmer_size){            #This is to check if kmers saved by store_vectotr match those present in the jelly file, as selected by lower-count . This is not possible when using different kmer sizes in jelly and vector files (mxky) 
        $distinct_repetitive_kmers = $jelly_Distinct - $jelly_Unique;
        if ($lower_count_value==1){print "jelly file contains  $jelly_Distinct distinct kmers. store_vector should have reported this number as  kmers stored in  file $out_vector_name  \n"};
        if ($lower_count_value==2){print "jelly file contains  $distinct_repetitive_kmers distinct kmers with 2 or more copies. store_vector should have reported this number as  kmers stored in  file $out_vector_name  \n"};   
        if ($lower_count_value> 2){print "Only possible to calculate the number of kmers present in jelly file for lower-count of 1 or 2 , so it is not possible to compare with value reported by store_vector. \n"};        
        }
if ($m_jelly != $kmer_size){print "It is not possible to compare the number of kmers present in jelly file with the value reported by store_vector when using different kmer sizes in jelly and vector files .  \n"};           

print "finished:",scalar localtime,"\n\n";
exit;



# **************************************************************************************************************************************************************

sub convert_mxkx{
$kmer_num=0;
$pid = open(IN, "-|"); # open pipe to read from child
$SIG{PIPE} = sub { die "whoops, $program pipe broken" };
if ($pid) { # parent
    while (<IN>) {              # read line into default varieble $_
        $kmer_num ++ ;
        my @kmer =  split ' ' ;
        $kmer_bin = $kmer[0] ;   
        $kmer_bin =~ s/A/00/g ; $kmer_bin =~ s/T/11/g ; $kmer_bin =~ s/G/10/g ; $kmer_bin =~ s/C/01/g ;
        $kmer_dec = oct("0b" . $kmer_bin); 
        $vector_output->Bit_On($kmer_dec);
        }
    close IN or warn "Child exit code: $?";
}
else {exec($program, $command, @options, $jelly_full_name) or die "Can't exec $program: $!";}  # child
return}


sub convert_mxky{
$kmer_num=0;
my $internal_kmer = $m_jelly - $kmer_size +1 ;  #e.g., a jellyfish 18-mer will be decomposed into 4 15-mers
$pid = open(IN, "-|"); # open pipe to read from child
$SIG{PIPE} = sub { die "whoops, $program pipe broken" };
if ($pid) { # parent
    while (<IN>) {              # read line into default varieble $_
        $kmer_num ++ ;
        my @kmer =  split ' ' ;
        for ($i = 0; $i < $internal_kmer; $i++) {
                $seq_or = substr( $kmer[0] , $i , $kmer_size);  
                $seq_rc = reverse $seq_or ;  $seq_rc =~ tr/ATCG/TAGC/;
                if ($seq_or lt $seq_rc){$kmer_bin = $seq_or}
                else{$kmer_bin = $seq_rc};        
                $kmer_bin =~ s/A/00/g ; $kmer_bin =~ s/T/11/g ; $kmer_bin =~ s/G/10/g ; $kmer_bin =~ s/C/01/g ;
                $kmer_dec = oct("0b" . $kmer_bin); 
                #print "$kmer[0]\t$kmer_dec\n";
                $vector_output->Bit_On($kmer_dec);
                };                    
                            
    }
    close IN or warn "Child exit code: $?";
}
else {exec($program, $command, @options, $jelly_full_name) or die "Can't exec $program: $!";}  # child
return}




sub store_vector{        		#two parameters: vector to be stored, and filename
($vector,$vector_filename) = @_;
$stored_kmer_number = $vector->Norm();
open VECTOR_S, ">:gzip","$vector_filename\.gz" or die "not possible to create the file $vector_filename\n";
$string = $vector->to_Hex();
undef($vector);
print VECTOR_S $string;
undef $string ;								#to really save memory
print "store_vector finished. $stored_kmer_number kmers stored in  file $vector_filename\n";
close (VECTOR_S);
return};


sub check_64bits_compatibility{   #kmer_size bigger than 16 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.14.2 -Duse64bitint)
my ($k_size) = @_;
my($flag_64bits_incompatibility)=0;
my($check_integer) = 2**36 ;                                    #checks perl exponentiation above 2^32
if ($check_integer ne "68719476736"){
        print "ERROR: this perl installation  fails in exponentiation above 2^32\n";
        $flag_64bits_incompatibility ++ ;
        }
else{print "passed  exponentiation above 2^32\n"};
 
my($check_binary) = "10000000000000000000000000000000000" ;     # checks perl bin->dec conversion above 2^32). The binary number is 2^34
$check_integer = oct("0b" . $check_binary);
if ($check_integer ne "17179869184"){
        print "ERROR: this perl installation fails in bin->dec conversion above 2^32 \n";
        $flag_64bits_incompatibility ++ ;
        }
else{print "passed bin->dec conversion above 2^32 \n"};        
        
my($Long_Bits) = Bit::Vector->Long_Bits() ;                     #checks installed Bit::Vector module. 
my($Word_Bits) = Bit::Vector->Word_Bits() ;
if ( ($Long_Bits != 64) || ($Long_Bits != 64)  ){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 16 (Long_Bits / Word_Bits test)\n";
        $flag_64bits_incompatibility ++ ;
        die "Not possible to run the program in this system due to large kmer_size ($k_size).\nkmer_size bigger than 15 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.14.2 -Duse64bitint)\n";
        }
else{print "passed    Long_Bits / Word_Bits=64\n"};      

my($bits) = 4**$k_size  ;
$check_vector = Bit::Vector->new($bits);                   #functional tests of installed Bit::Vector module.
my($vector_size) = $check_vector->Size();
if ($bits != $vector_size){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 15. failed in new() or Size()\n";
        $flag_64bits_incompatibility ++ ;
        }
else{print "passed    functional test of installed Bit::Vector module  new()  Size()\n"};         
        
my($state_before_On,$state_after_On ) = 2;                      #functional tests of installed Bit::Vector module.
$state_before_On = $check_vector->bit_test($vector_size -1);
$check_vector->Bit_On($vector_size -1);
$state_after_On  = $check_vector->bit_test($vector_size -1);
if ( ($state_before_On != 0) || ($state_after_On != 1)  ){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 16. failed in Bit_On() or  bit_test() \n";
        $flag_64bits_incompatibility ++ ;
        }   
else{print "passed    functional test of installed Bit::Vector module  bit_test()  Bit_On()\n"};
 
if ($flag_64bits_incompatibility > 0){die "Not possible to run the program in this system due to large kmer_size ($k_size).\nkmer_size bigger than 16 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.16.0 -Duse64bitint)\n"}
else{print "passed all 64-bits tests\n"};
undef($check_vector);  #to save memory
return}



sub process_command_line{
$jelly_prefix = ""; $jelly_full_name = ""; $kmer_size=0;$m_jelly =0; $mode=""; $lower_count_value="not_set";
if (@ARGV == 0){print "usage: perl jelly_2_bitvector_mxky_1e.pl m_jelly=15 kmer_size=15 jelly_file=fvirf30m15.jelly  lower-count=12  \n"; die}
$num_arg = @ARGV;
for ($i = 0 ; $i < $num_arg ; $i++){
    #print "\n\n num_arg: $num_arg     command argument: @ARGV    \nargument number $i  $ARGV[$i]\n\n"; 
    $arg = $ARGV[$i];
    my($flag_match)=0;
    if ($arg =~ m/jelly_file=([\w\.]+)/ ){$jelly_full_name = (split /\=/, $arg)[1]; $flag_match ++};   
    if ($arg =~ m/jelly_file=(\w+)/ ){ $jelly_prefix = (split /\./, $jelly_full_name)[0]; $flag_match ++};
    if ($arg =~ m/kmer_size=(\d+)/ ){ $kmer_size = $1; $flag_match ++ };
    if ($arg =~ m/m_jelly=(\d+)/ ){ $m_jelly = $1; $flag_match ++ };
    if ($arg =~ m/lower-count=(\d+)/ ){ $lower_count_value = $1; $flag_match ++ };
    if ($arg =~ m/--lower-count=(\d+)/ ){ $lower_count_value = $1; $flag_match ++ };         #old style, preserved for  backward compatibility
    if ($arg =~ m/mode=save_all/ ){$mode="save_all"; $lower_count_value=0 ; $flag_match ++ };   #save_all overrides  lower_count_value   old style, preserved for  backward compatibility
    if ($flag_match==0){die "Unrecognized option $arg . Check spelling. \n"};
    };
return}


sub get_jelly_kmer_size{
my($jelly_file)=@_;
my($jelly_kmer_size)=0;
my(@jellyfish_stats_output) = `$program stats  $jelly_file`;
$jellyfish_stats_output[0] =~ m/Unique:\s*(\d+)/ ;                 $jelly_Unique = $1;
$jellyfish_stats_output[1] =~ m/Distinct:\s*(\d+)/ ;               $jelly_Distinct = $1;
$jellyfish_stats_output[2] =~ m/Total:\s*(\d+)/ ;                  $jelly_Total = $1;print "jellyfish statistics from file $jelly_file : \n";
#getting kmer size info:
my(@jellyfish_info_output) = `$program info  $jelly_file`;
$jellyfish_info_output[0] =~ m/-m\s+(\d+)/ ;   $jelly_kmer_size = $1;
print "jellyfish statistics from file $jelly_file : \n";
print "Unique:     $jelly_Unique   [kmers that occur only once. If mode=single_copy, this number should be reported by store_vector.] \n";
print "Distinct:   $jelly_Distinct    [all different kmers. If all kmers are saved, this number should be reported by store_vector.] \n";
print "Total:      $jelly_Total [total valid kmers (many will be identical) ] \n";
return $jelly_kmer_size}

# [bernardo@localhost GBR200f]$ jellyfish stats -v GBR200f_m15_f10.jelly
# k-mer length (bases): 15
# value length (bytes): 4
# max reprobe         : 1
# hash size           : 1073741824
# Unique:    45288795
# Distinct:  415878572
# Total:     145219366462
# Max_count: 97367652    

