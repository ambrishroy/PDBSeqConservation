#! /usr/bin/perl
use strict;
use common;  ## My module to define paths (e.g. $home, $bindir)

###############################################################################################################################
## This program for calculating residue conservation score (based on Jensen–Shannon divergence score & sum-of-pairs) is      ##
## written by Ambrish Roy                                                                                                    ## 
##                                                                                                                           ##
## This code is inspired from Python code by Capra JA (http://compbio.cs.princeton.edu/conservation/conservation_code.tar.gz)##
##                                                                                                                           ##
## Please address any comments/bug-reports to: ambrish.roy AT gmail.com                                                      ##
###############################################################################################################################

## Scores conservation of each residues as putative binding site residue based on JSD and SOP (Sum of pairs) using
## PsiBlast output
my $myARGC=scalar(@ARGV);
if ($myARGC<1) {&fexit;}


my $s         =$ARGV[0];  # Protein name
my $datadir   =$ARGV[1];  # For finding  PDB file $s.pdb (single chain)
my $record    =$ARGV[2];  # For keeping output

my $rootdir   ="/gpfs1/u";
my $user      ="aroy";
my $nrdb      ="$home/db/nr/nr";     # Non-redundant sequence database 
my $blastdir  ="$rootdir/aroy45/bin/blast";    # Directory containing  BLAST executables

######################################################
sub fexit(){
    print("\nExample\n");
    print("Usage:   ./JSD  1XXXX.pdb  /home/aroy/pdb  /home/aroy/seqconserve/conservation.txt\n");
    print("\n");
    exit(1);
}
my $tag       ="profile_$s"; 
my $random    = int(rand(1000000));
my $workdir   = "/tmp/aroy45/$tag";
if(exists $ENV{'PBS_JOBID'}){
  $workdir = "/tmp/$user/".$ENV{'PBS_JOBID'}."/$tag";
}
if(-d "$workdir"){
    system("rm -rf $workdir") == 0 or die "System call failed: $!";
}
system("mkdir -p $workdir") == 0 or die "System call failed: $!";
chdir "$workdir";

#####################################################
# SETTINGS (Do not change below this line)
my @AA=(
	"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET","PHE",
	"PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASX", "GLX", "UNK");

my %AA2index=(
	      'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12', 
	      'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20', 'B'=>'21', 'Z'=>'22', '-'=>'23');	
my %index2AA=(
	      '1'=>'A','2'=>'R','3'=>'N','4'=>'D','5'=>'C','6'=>'Q','7'=>'E','8'=>'G','9'=>'H','10'=>'I','11'=>'L','12'=>'K', 
	      '13'=>'M','14'=>'F','15'=>'P','16' =>'S','17'=>'T','18'=>'W','19'=>'Y','20'=>'V','21'=>'B','22'=>'Z','23'=>'-');	

my %ts=(
	'ALA'=>'A', 'ARG'=>'R', 'ASN'=>'N', 'ASP'=>'D', 'CYS'=>'C', 'GLN'=>'Q', 'GLU'=>'E', 'GLY'=>'G',
	'HIS'=>'H', 'ILE'=>'I', 'LEU'=>'L', 'LYS'=>'K', 'MET'=>'M', 'PHE'=>'F', 'PRO'=>'P', 'SER'=>'S', 
	'THR'=>'T', 'TRP'=>'W', 'TYR'=>'Y', 'VAL'=>'V', 'ASX'=>'N', 'GLX'=>'Q', 'UNK'=>'G',
	'A'=>'ALA', 'R'=>'ARG', 'N'=>'ASN', 'D'=>'ASP', 'C'=>'CYS', 'Q'=>'GLN', 'E'=>'GLU', 'G'=>'GLY',
	'H'=>'HIS', 'I'=>'ILE', 'L'=>'LEU', 'K'=>'LYS', 'M'=>'MET', 'F'=>'PHE', 'P'=>'PRO', 'S'=>'SER', 
	'T'=>'THR', 'W'=>'TRP', 'Y'=>'TYR', 'V'=>'VAL', 'B'=>'ASX', 'Z'=>'GLX', 'X'=>'UNK');

                # A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
my @BLOSUM62 = (
		[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0],
		[-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, 0, 1,-1],
		[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
		[-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 1, 0,-1],
		[ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2],
		[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
		[-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 0, 2,-1],
		[ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, 0,-2,-1],
		[-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 1, 0,-1],
		[-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1],
		[-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-3,-2,-1],
		[-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1],
		[-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-2, 0,-1],
		[-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1],
		[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2],
		[ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 1, 0, 0],
		[ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0,-1, 0],
		[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2],
		[-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-2,-1,-1],
		[ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1],
		[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
		[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
		[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1]);

# This is the BLOSUM62 background distribution.
my %blsm_bck=(
	      'A'=>0.078, 'R'=>0.051, 'N'=>0.041, 'D'=>0.052, 'C'=>0.024, 'Q'=>0.034, 'E'=>0.059, 'G'=>0.083, 
	      'H'=>0.025, 'I'=>0.062, 'L'=>0.092, 'K'=>0.056, 'M'=>0.024, 'F'=>0.044, 'P'=>0.043, 'S'=>0.059, 
	      'T'=>0.055, 'W'=>0.014, 'Y'=>0.034, 'V'=>0.072, 'B'=>0.041, 'Z'=>0.034, 'X'=>0.083,);
# This is the BLOSUM62 frequency.	 
my @BLOSM_FREQ = ([0.0215,0.0023,0.0019,0.0022,0.0016,0.0019,0.0030,0.0058,0.0011,0.0032,0.0044,0.0033,0.0013,0.0016,0.0022,0.0063,0.0037,0.0004,0.0013,0.0051,0.0019,0.0019,0.0058],
		  [0.0023,0.0178,0.0020,0.0016,0.0004,0.0025,0.0027,0.0017,0.0012,0.0012,0.0024,0.0062,0.0008,0.0009,0.0010,0.0023,0.0018,0.0003,0.0009,0.0016,0.0020,0.0025,0.0017],
		  [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
		  [0.0022,0.0016,0.0037,0.0213,0.0004,0.0016,0.0049,0.0025,0.0010,0.0012,0.0015,0.0024,0.0005,0.0008,0.0012,0.0028,0.0019,0.0002,0.0006,0.0013,0.0037,0.0016,0.0025],
		  [0.0016,0.0004,0.0004,0.0004,0.0119,0.0003,0.0004,0.0008,0.0002,0.0011,0.0016,0.0005,0.0004,0.0005,0.0004,0.0010,0.0009,0.0001,0.0003,0.0014,0.0004,0.0003,0.0008],
		  [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
		  [0.0030,0.0027,0.0022,0.0049,0.0004,0.0035,0.0161,0.0019,0.0014,0.0012,0.0020,0.0041,0.0007,0.0009,0.0014,0.0030,0.0020,0.0003,0.0009,0.0017,0.0022,0.0035,0.0019],
		  [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378],
		  [0.0011,0.0012,0.0014,0.0010,0.0002,0.0010,0.0014,0.0010,0.0093,0.0006,0.0010,0.0012,0.0004,0.0008,0.0005,0.0011,0.0007,0.0002,0.0015,0.0006,0.0014,0.0010,0.0010],
		  [0.0032,0.0012,0.0010,0.0012,0.0011,0.0009,0.0012,0.0014,0.0006,0.0184,0.0114,0.0016,0.0025,0.0030,0.0010,0.0017,0.0027,0.0004,0.0014,0.0120,0.0010,0.0009,0.0014],
		  [0.0044,0.0024,0.0014,0.0015,0.0016,0.0016,0.0020,0.0021,0.0010,0.0114,0.0371,0.0025,0.0049,0.0054,0.0014,0.0024,0.0033,0.0007,0.0022,0.0095,0.0014,0.0016,0.0021],
		  [0.0033,0.0062,0.0024,0.0024,0.0005,0.0031,0.0041,0.0025,0.0012,0.0016,0.0025,0.0161,0.0009,0.0009,0.0016,0.0031,0.0023,0.0003,0.0010,0.0019,0.0024,0.0031,0.0025],
		  [0.0013,0.0008,0.0005,0.0005,0.0004,0.0007,0.0007,0.0007,0.0004,0.0025,0.0049,0.0009,0.0040,0.0012,0.0004,0.0009,0.0010,0.0002,0.0006,0.0023,0.0005,0.0007,0.0007],
		  [0.0016,0.0009,0.0008,0.0008,0.0005,0.0005,0.0009,0.0012,0.0008,0.0030,0.0054,0.0009,0.0012,0.0183,0.0005,0.0012,0.0012,0.0008,0.0042,0.0026,0.0008,0.0005,0.0012],
		  [0.0022,0.0010,0.0009,0.0012,0.0004,0.0008,0.0014,0.0014,0.0005,0.0010,0.0014,0.0016,0.0004,0.0005,0.0191,0.0017,0.0014,0.0001,0.0005,0.0012,0.0009,0.0008,0.0014],
		  [0.0063,0.0023,0.0031,0.0028,0.0010,0.0019,0.0030,0.0038,0.0011,0.0017,0.0024,0.0031,0.0009,0.0012,0.0017,0.0126,0.0047,0.0003,0.0010,0.0024,0.0031,0.0019,0.0038],
		  [0.0037,0.0018,0.0022,0.0019,0.0009,0.0014,0.0020,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0014,0.0047,0.0125,0.0003,0.0009,0.0036,0.0022,0.0014,0.0022],
		  [0.0004,0.0003,0.0002,0.0002,0.0001,0.0002,0.0003,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0008,0.0001,0.0003,0.0003,0.0065,0.0009,0.0004,0.0002,0.0002,0.0004],
		  [0.0013,0.0009,0.0007,0.0006,0.0003,0.0007,0.0009,0.0008,0.0015,0.0014,0.0022,0.0010,0.0006,0.0042,0.0005,0.0010,0.0009,0.0009,0.0102,0.0015,0.0007,0.0007,0.0008],
		  [0.0051,0.0016,0.0012,0.0013,0.0014,0.0012,0.0017,0.0018,0.0006,0.0120,0.0095,0.0019,0.0023,0.0026,0.0012,0.0024,0.0036,0.0004,0.0015,0.0196,0.0012,0.0012,0.0018],
		  [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
		  [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
		  [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378]);


#####################################################
my %AA2sim=(); ## Blosum similarity of amino acids
for(my $i=0;$i<23;$i++){
  for(my $j=0;$j<23;$j++){
    my $x= $i+1;my $y= $j+1;
    #	$AA2sim{$index2AA{$x},$index2AA{$y}}= $BLOSM_FREQ[$i][$j];
    $AA2sim{$index2AA{$x},$index2AA{$y}}= $BLOSUM62[$i][$j];
  }
}

my $PSEUDOCOUNT= 0.0000001;
my $lambda     = 0.5;
my $window     = 3;

###############################################################################
my $pdbid=substr($s,0,4);$pdbid=~s/\s+//g;
my $n=-1000;my $chain='';my $position=0;my %seqQ=();
my $last_res_id='';my $last_res_name='';
####Renumber PDB residues from 1
open(FILE,"<$datadir/$s.pdb") || die "Can't open $datadir/$s.pdb\n";
while(my $line=<FILE>){
  if($line=~/^TER/) {goto pdbend;}  # start over if we find a terminus
  if($line=~/^ATOM/){
    my $alt = substr($line,16,1); #alterate for atom
    if($alt eq " " || $alt eq "A"|| $alt eq "1"){
      $position++;
      if ($n == -1000){
	$n = 1;# start with the starting number in the chain
	$last_res_id   = substr($line,22,5);
	$last_res_name = substr($line,17,3);
	$seqQ{$n}=$ts{$last_res_name};
      }
      my $res_id      = substr($line,22,5);
      my $res_name    = substr($line,17,3);
      if(($res_id ne $last_res_id) || ($res_name ne $last_res_name)){
	$n++;
	$last_res_id   = $res_id;
	$last_res_name = $res_name;
	$seqQ{$n}=$ts{$res_name};
      }
      substr($line,6,5)  = sprintf "%5d",$position;
      substr($line,22,5) = sprintf "%4d ",$n;
    }
  }
}
pdbend:;
close(FILE);
my $Lch=$n;

my $sequence="";
open(FASTA,">$workdir/protein.fasta");
print FASTA ">$s\n";
for(my $i=1;$i<=$Lch;$i++){
  print FASTA "$seqQ{$i}";
  $sequence .=$seqQ{$i};
  if(int($i/60)*60==$i){
    print FASTA "\n";
  }
}
print FASTA "\n";
close(FASTA);


######### Run Blast
if(!-s "$record/blast/$s.gz"){
  print "$blastdir/bin/blastpgp  -b 1000 -j 3 -h 0.001 -d $nrdb -i protein.fasta\n";
  system ("$blastdir/bin/blastpgp  -b 1000 -j 3 -h 0.001 -d $nrdb -i protein.fasta > $s\_blast.out");
  system("gzip -c $s\_blast.out > $s.gz");
  system("/bin/cp $workdir/$s.gz $record/blast/")== 0 or die "System call for file copy failed: $!";	
}
else{
  system("zcat $record/blast/$s.gz > $workdir/$s\_blast.out")== 0 or die "System call for file copy failed: $!";
}

######### Convert BLAST output to Multiple Sequence alignment using alignhits.pl (get it from HHsearch package)
system("$bindir/alignhits.pl -fas $s\_blast.out -q $workdir/protein.fasta $workdir/blast.fasta") == 0 or die "System call failed: $!";

my %ALN=();my $Pcount=0; ## Protein alignment and P-value tracker
open(ALN,"<$workdir/blast.fasta") || die "Cant open blast.fasta: $!";
while(my $line=<ALN>){
  chomp($line);
  if($line=~/^>(\S+)/){
    my $Pname=$1;
    $Pcount++;
    my $Evalue= $1 if($line=~/E=(\S+)/);	
    $ALN{$Pcount,0}=$Pname;
    $ALN{$Pcount,1}=$Evalue;
  }
  else{
    $line=~s/X/-/g;
    $ALN{$Pcount,2}=$line;	
  }
}
close(ALN);

if($Pcount > 10){
  ## Can only use MSA if no of protein is > 10
  ## >1000 protein increases computation without much improvement in performance
  $Pcount=1000 if ($Pcount > 1000);
  my %seq_weights = ();
  &load_sequence_weights(\%ALN,\$Pcount,\%AA2index,\%seq_weights);

  my %jsd_scores = ();my %zscore1=(); my %wscore1=();my %zscore3=();
  &js_divergence_score(\%ALN,\$Pcount,\%blsm_bck,\%AA2index,\%seq_weights,$PSEUDOCOUNT,$lambda,\%jsd_scores);
  &z_scores(\%jsd_scores,\$Lch,\%zscore1);

  #### Window removed because of low performance
  &w_scores(\%jsd_scores,\$Lch,\$window,\$lambda,\%wscore1);
  &z_scores(\%wscore1,\$Lch,\%zscore3);

  my %SOP=();my %zscore2=();my %wscore2=();
  &sum_of_pairs(\%ALN,\$Pcount,\%AA2sim,\%seq_weights,\%SOP);
  &z_scores(\%SOP,\$Lch,\%zscore2);

  #### Window removed because of low performance
  ## JSD:  Jensen–Shannon divergence score
  ## SOPM: Sum-of-pair measure
  ## JSDW: Weighted Jensen–Shannon divergence score
  open(OUT1,">$workdir/$s\_JSD.dat");
  open(OUT2,">$workdir/$s\_SOPM.dat");
  open(OUT3,">$workdir/$s\_JSDW.dat");
  for(my $j=1;$j<=$Lch;$j++){
    printf OUT1 "%-4d  %5s\n",$j,$zscore1{$j};
    printf OUT2 "%-4d  %5s\n",$j,$zscore2{$j};
    printf OUT3 "%-4d  %5s\n",$j,$zscore3{$j};
  }
  close(OUT1);
  close(OUT2);
  close(OUT3);
}
else{
  open(OUT1,">$workdir/$s\_JSD.dat");
  open(OUT2,">$workdir/$s\_SOPM.dat");
  open(OUT3,">$workdir/$s\_JSDW.dat");
  for(my $j=1;$j<=$Lch;$j++){
    printf OUT1 "%-4d  %5s\n",$j,'1.00';
    printf OUT2 "%-4d  %5s\n",$j,'1.00';
  }
  close(OUT1);
  close(OUT2);
  close(OUT3);
}
###### Copy output and clean working directory
system("/bin/cp $workdir/protein.fasta $record/fasta/") == 0 or die "System command failed: $!";
system("/bin/cp $workdir/*.dat $record/conservation") == 0 or die "System command failed: $!";
system("/bin/rm -rf $workdir")== 0 or die "System command failed: $!";
exit();


sub sum_of_pairs
{
    my ($ALN_ref,$Nseq_ref,$blossum_ref,$SW_ref,$SOP)=@_;
    my %align       =%$ALN_ref;
    my $Nseq        =$$Nseq_ref;
    my %sim_matrix  =%$blossum_ref;
    my %weights     =%$SW_ref;  
   
    my @Qres        =split(//,$align{1,2});
    my $Ncol        =$#Qres;    
    my $Qresno=0;my %Qmapping=();
    for(my $j=0;$j<=$#Qres;$j++)
    {
	if($Qres[$j] ne '-'){$Qresno++;$Qmapping{$Qresno}=$j;} 	
    }
    my @ARR=();my $sum_seq_weights = 0;
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	for(my $j=0;$j<=$#res;$j++)
	{
	    $ARR[$i][$j]=$res[$j];
	}
	$sum_seq_weights += $weights{$i};
    }
    my %SOP_scores=();
    for(my $C=0;$C<=$Ncol;$C++)
    {
	my $sum   = 0;
	my $maxsum= 0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	    my $aa1=$ARR[$i][$C];
	    for(my $j=$i+1;$j<=$Nseq;$j++)
	    {
		my $aa2=$ARR[$j][$C];
		if(($ARR[$i][$C] eq '-')||($ARR[$j][$C] eq '-')){next;}
		$maxsum += $weights{$i} * $weights{$j};
		$sum    += $weights{$i} * $weights{$j} * $sim_matrix{$aa1,$aa2};	
	    }	    
	}
	if ($maxsum !=0)
	{
	    $sum   /=$maxsum;
	}
	else
	{
	    $sum=0;
	}
	my $gap_sum=0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	  if($ARR[$i][$C] eq '-')
	  {
	      $gap_sum +=   $weights{$i}; 
	  }
	}
	my $weighted_gap_penalty = (1-  ($gap_sum/$sum_seq_weights));
	$SOP_scores{$C}=($sum*$weighted_gap_penalty);	
    }
    for(my $j=1;$j<=$Qresno;$j++)
    {
	$$SOP{$j}= sprintf("%5.4f",$SOP_scores{$Qmapping{$j}});
    }
    return;
    
}
sub js_divergence_score
{
    my ($ALN_ref,$Nseq_ref,$dist_ref,$AA_ref,$SW_ref,$pseudocount,$lambda,$JSD_ref)=@_;
    my %align  =%$ALN_ref;
    my $Nseq   =$$Nseq_ref;
    my %distr  =%$dist_ref;
    my %weights=%$SW_ref;  
    my %AA2in  =%$AA_ref;    
    my @Qres   =split(//,$align{1,2});
    my $Ncol   =$#Qres;
    
    my $Qresno=0;my %Qmapping=();
    for(my $j=0;$j<=$#Qres;$j++)
    {
	if($Qres[$j] ne '-'){$Qresno++;$Qmapping{$Qresno}=$j;} 	
    }  
    my @ARR=();
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	for(my $j=0;$j<=$#res;$j++)
	{
	    $ARR[$i][$j]=$res[$j];
	}
    }
    my $AAcount = keys %AA2in;
    my %AA_freq=();my %sum_seq_weights=();my %R=();my %JSD_scores=();
    for(my $j=0;$j<=$Ncol;$j++)
    {
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{
	    $AA_freq{$j,$key}=0;
	    $R{$j,$key}      =0;
	}
	$sum_seq_weights{$j} =0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	    $AA_freq{$j,$ARR[$i][$j]} += (1.0 * $weights{$i});
	    $sum_seq_weights{$j}      += $weights{$i}; 
	}
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{	   
	    $AA_freq{$j,$key}=$AA_freq{$j,$key}/($sum_seq_weights{$j} + $AAcount * $pseudocount);
	    if($key ne '-')
	    {
		$R{$j,$key}=($lambda*$AA_freq{$j,$key}) + ((1-$lambda)*$distr{$key});	   
	    }
	}
    }
    for(my $j=0;$j<=$Ncol;$j++)
    {
	my $JSD_score=0;my $weighted_gap_penalty=0;my $gap_sum=0;
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{
	    if($key eq '-'){next;}
	    if($R{$j,$key} != 0.0)
	    {
		if($AA_freq{$j,$key} == 0.0)
		{
		    $JSD_score += ($distr{$key} * log2($distr{$key}/$R{$j,$key}));
		}
		elsif($distr{$key} == 0.0)
		{
		    $JSD_score += ($AA_freq{$key} * log2($AA_freq{$key}/$R{$j,$key}));
		}
		else
		{
		    $JSD_score += ($AA_freq{$j,$key} * log2($AA_freq{$j,$key}/$R{$j,$key})) + ($distr{$key} * log2($distr{$key}/$R{$j,$key}));
		}
	    }
	}	
	$JSD_score = $JSD_score/2;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	  if($ARR[$i][$j] eq '-')
	  {
	      $gap_sum +=   $weights{$i}; 
	  }
	}
	$weighted_gap_penalty = (1-  ($gap_sum/$sum_seq_weights{$j}));
	$JSD_scores{$j}=($JSD_score*$weighted_gap_penalty);	
    }
    for(my $j=1;$j<=$Qresno;$j++)
    {
	$$JSD_ref{$j}= sprintf("%5.4f",$JSD_scores{$Qmapping{$j}});
    }
    return;
}

# Get Seqeuence Weights based on Henikoff-Henikoff schema
sub load_sequence_weights
{
    my ($ALN_ref,$Nseq_ref,$AA_ref,$seq_weights_ref)=@_;

    my %align  =%$ALN_ref;
    my $Nseq   =$$Nseq_ref;
    my %AA2in  =%$AA_ref;

    my %NRC=();my %RC=();my %seen=();
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});

	for(my $j=0;$j<=$#res;$j++)
	{
	    my $AAN=$AA2in{$res[$j]};
	    $RC{$j,$AAN}++;	  
	    if(exists $seen{$j,$AAN}){next;}
	    else{$seen{$j,$AAN}=1;$NRC{$j}++;}	    
	}	
    }   
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	
	for(my $j=0;$j<=$#res;$j++)
	{
	    my $AAN=$AA2in{$res[$j]};
	    $$seq_weights_ref{$i} += 1/($NRC{$j} * $RC{$j,$AAN});    	   
	}
	$$seq_weights_ref{$i} = sprintf("%6.4f",$$seq_weights_ref{$i}/$#res);
    }
    return;  
}
	
sub log2 
{
    my $n = shift;
    return log($n)/log(2);
}

sub w_scores
{# This subroutine takes a list of scores and a length and transforms them so that each position is a weighted average of the surrounding positions. 
     my($score_ref,$L,$W,$lambda,$w_score_ref)=@_;
     my %r_score  =%$score_ref;
     my $Lch   =$$L;
     my $window   =$$W;
     for(my $i=1;$i<=$Lch;$i++)
     {
	 my $sum     = 0;
	 my $n_terms = 0;
	 for(my $j=$i-$window;$j<=$i+$window;$j++)
	 {
	     if(($i != $j) && (exists $r_score{$j}))
	     {
		 $sum +=  $r_score{$j};
		 $n_terms++;
	     }
	 }
	 $$w_score_ref{$i}= sprintf("%5.4f",$sum/$n_terms);
     }
}
    
sub z_scores
{# This subroutine calculates Z-scores for a set of scores (NOTE: hash index starts from 1)   
     my($score_ref,$L,$z_score_ref)=@_;
     my %r_score  =%$score_ref;
     my $Lch   =$$L;
     my $z1_a  = 0; my $z1_a2=0; my %zscore=();
     ### This is a fast way to calculate approx Z-score (not the same as we learn in basic Statistics) 
     for(my $i=1;$i<=$Lch;$i++)
     {	 
	 $z1_a  += $r_score{$i}; 
	 $z1_a2 += $r_score{$i}**2;
     }
     $z1_a/=$Lch;
     $z1_a2/=$Lch;
     my $z1_sd=sqrt($z1_a2-$z1_a**2);
     for(my $i=1;$i<=$Lch;$i++)
     {
         $$z_score_ref{$i}=($r_score{$i}-$z1_a)/$z1_sd;
     }
}

