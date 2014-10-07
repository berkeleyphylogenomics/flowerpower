#!/usr/bin/perl -w
##!/usr/bin/perl -d:DProf # Use this line to profile this script (see dprofpp).

##--- This is a new version of Flowerpower. 

## EDIT: 08/03/09 - DMS added code to cut the number of homologs 
##       back to FPHITS in the case that FPHITS is exceeded.

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/tempdir/;

my $major_version = 2;
my $minor_version = 3;

my $debug1 = 0;

my $seedfile = "";
my $msafile = "";
my $niters = 3;
my $universe = undef;
my $pwid = 0.20;
my $DB = "UniProt/current/protein";
my $MUSCLE_MAXITERS = 2;
my $PSIBLAST_ITERS = 3;
my $PBEVALUE = 10;
my $FPMAXEVALUE = -1;  # Default - depends on seed sequence length.
my $NPSIHITS = 500;
my $FPHITS = 1000;
my $SW_SCORE = 2;
my $COV = 0;
my $HCOV = 0;
my $QCOV = 0;
my $MODE = "global";
my $SPEED = "fast" ;
my $use_queue = 0;
my $USE_QUEUE = "" ;
my $QUEUE = "library";
my $return_value = 0;
my $mult_cmd_tries = 7;
my $ii;
my $temporary_check = 0;
my $verbose = 0;
my $seedmaster = 1;
my $pfam = 0;
my $blast_again = 1;
my $max_insert = 40;
my $max_gap = 40;
my $use_bins = 0;
my $bin_interval = 50;
my $max_length_to_bin = 1200;
my $msa_profile = 0;
my $nooutlog = 0;
my $initial_percent_id = 0.3;
my $initial_max_insert = 30;
my $initial_max_gap = 30;
my $initial_min_coverage = 0;
my $excess_homologs_flag = 0;

# # DMS 08/03/09 - used when the number of homologs exceeds FPHITS
my $uncropped_seed_id = "";

GetOptions (
           "-i=s" => \$seedfile,
           "-a=s" => \$msafile,
           "-n=i" => \$niters,
           "--psiiters=i" => \$PSIBLAST_ITERS,
           "--npsihits=i" => \$NPSIHITS,
           "--psievalue=f" => \$PBEVALUE,
           "--fpmaxevalue=f" => \$FPMAXEVALUE,
           "--fphits=i" => \$FPHITS,
           "-d=s" => \$DB, 
           "-u=s" => \$universe,
           "-s=i" => \$SW_SCORE,
           "-p=f" => \$pwid,
           "--cov=i" => \$COV, 
           "--qcov=f" => \$QCOV, 
           "--hcov=f" => \$HCOV, 
           "--mode=s" => \$MODE,
           "--speed=s" => \$SPEED, 
           "--tempcheck=i" => \$temporary_check,
           "--verbose=i" => \$verbose,
           "--seedmaster=i" => \$seedmaster,
           "--pfam=i" => \$pfam,
           "--blast_again=i" => \$blast_again,
           "--max_insert=i" => \$max_insert,
           "--max_gap=i" => \$max_gap,
           "--use_bins=i" => \$use_bins,
           "--bin_width=i" => \$bin_interval,
           "--longest_bin=i" => \$max_length_to_bin,
           "--use_queue=i" => \$use_queue,
           "--queue=s" => \$QUEUE,
           "--msa_profile=i" => \$msa_profile,
           "--no_output_log=i" => \$nooutlog,
           "--initial_percent_id=f" => \$initial_percent_id,
           "--initial_max_insert=i" => \$initial_max_insert,
           "--initial_max_gap=i" =>\$initial_max_gap,
           "--initial_min_coverage=i" =>\$initial_min_coverage
           );

# -- parse command line options -- #

if ($verbose) {
  $debug1 = 1 ;
}

if ($use_queue) {
  $USE_QUEUE = "--use_queue --queue $QUEUE";
  $temporary_check = 1;
}

print $seedfile;
print $msafile;

unless ( -e $seedfile || -e $msafile ){
   print <<EOF;
Arguments
  -i            Source sequence filename if starting from seed.
  -a            Alignment file if starting from MSA.
  -n            Number of SHMMs iterations (default: 3).
  --psiiters    Number of PSIBLAST iterations (default: 3).
  --npsihits    Number of PSIBLAST hits returned per iteration
  --psievalue   E-value threshold to use in PSIBLAST run (default: 10).
  --fpmaxevalue E-value threshold to use in SHMMs run (default: depends on seed
                length: < 65: 1e-2; <= 100: 1e-3; > 100: 1e-4).
  --fphits      Number of FlowerPower hits to return (default: 1000).
  -d            Database (default: UniProt).
  -u            Universe for flowerpower run (default: will be created by
                running PSIBLAST).
  -s            sw option to use for assign_seqs (default: 2).
  -p            Pwid (default: 0.20).
  --cov         Coverage parameter (default: 0).  If set to 1 less stringent
                parameter will be used; values 2-4 user progressively less stringent values.
  --qcov        Query coverage parameter (default: 0).  Use a fractional setting
                such as 0.65 to override the coverage based on seed length.
  --hcov        Hit coverage parameter (default: 0).  Use a fractional setting
                such as 0.65 to override coverage based on hit length. This
                works only for global option.  Hit coverage is not calculated
                for glocal option.
  --mode        Alignment mode (default: global).  Options: global or glocal.
                If set to glocal: align global to HMM and local to hit.
  --speed       slow or fast (default: fast).  Fast reduces universe size at
                each iteration by removing sequences with e-value > 100).
  --tempcheck   1: use current directory for temporary files (default: 0 - use 
                /tmp)
  --verbose     1: print verbose messages (default: 0).
  --seedmaster  1: keep master-slave alignment with seed as master (default: 1).
  --pfam        1: do not trim the MSA or realign with MUSCLE (default: 0).
  --blast_again 1: run BLAST on supplied universe to get e-values (default: 1).
  --max_insert  Maximum number of contiguous inserts allowed (default: 40).
  --max_gap     Maximum number of contigous gaps allowed (default: 40).
  --use_bins    Whether to use the database binned by length
  --bin_width   Width of the range of sequence lengths to put in one bin
  --longest_bin Length for which longer sequences all go in the last bin
  --no_output_log 1: print logging information to the screen rather than file (default: 0).
  --msa_profile 1: use the input MSA to build an initial profile (default 0).
  --initial_percent_id Percent id cutoff for initial homolog set (default: 0.3)
  --initial_max_insert Largest insert allowable for initial homolog set (default 30)
  --initial_max_gap Largest deletion allowable for initial homolog set (default 30)
  --initial_min_coverage Coverage criterion for initial homolog set (default 0, stringent)
                         use 1-4 for less stringent initial inclusion coverage.

EOF
   exit 0;
}

unless (($MODE eq "global") || $MODE eq "glocal"){
   print "The --mode option has to be either global or glocal\n";
   exit 0;
}

#--- Get path for current working directory 

my $mwd = `pwd`;
chomp $mwd;

# print header for user
print "Using FlowerPower $major_version.$minor_version to retrieve homologs from ";

if ($DB =~ /UniProt/){
    print "UniProt.\n\n";
}
else{
    print "$DB.\n\n";
}

if ($debug1) {
   print "$mwd\n";
}

########################################################################
## Create a temporary working directory;
## Copy seed file or MSA file into it.
## Change to tempdir and do everything within that
########################################################################
if ($temporary_check == 0) {
  my $tempdir = tempdir("fptempXXXX", DIR => "/tmp/", CLEANUP => 0);
  if ($debug1) { 
    print "Copying files to temporary directory: $tempdir\n";
  }

  if ($seedfile) {
   `cp $seedfile $tempdir`;
  }
  if ($msafile) {
   `cp $msafile $tempdir`;
  }
  if ($universe) {
   `cp $universe $tempdir`;
  }

  chdir $tempdir ;
}

my @fa = ();
my $SEED_SEQ = "";
my @seed = ();
my $seedid = "" ;
my $seedseq = "";
my $unaltered_seedid = "";
my $SEQLEN = 0;
my $seed_header = "";

########################################################################
## Check if seedfile contains the "/" character and replace with "_"
########################################################################
my $defline = "";
if ($seedfile) {
    `cp $seedfile submitted-seed.fa`;
    open (NEW, "> seed_modified.fa");
    @seed = fastaFileToArray($seedfile);
    $seed_header = $seed[0];
    $seed_header =~ s/\n.*$//;       # Delete everything after first line.
    $unaltered_seedid = fasta2id($seed[0]);
    $seedid = "seed_" . $unaltered_seedid;
    $seedid =~ s/\//_/;
    $seedseq = fasta2seq($seed[0]);
    open (FILE, "< $seedfile");
    while (my $line = <FILE>) {        
        if ($line =~ /^>/) {
            $defline = $line;
            $line =~ s/>.*/>lcl|$seedid/;
	    
            # # DMS 08/03/09 - used when the number of homologs exceeds FPHITS
	    $uncropped_seed_id = $line;
	    $uncropped_seed_id =~ s/>//;


        }
        if ($line =~ /^>/ && $line =~ /\//) {  
            $line =~ s/\//_/;
            print NEW "$line";
        }
        else {
            print NEW "$line";          
        }
     }
     print NEW "\n";
    `mv seed_modified.fa seed.fa`;
     $seedfile = "seed.fa";

     # print "seedid: $seedid\n";
     # print "defline: $defline\n";
     # print "uncropped seed id: $uncropped_seed_id\n";

} else {
  # Can only keep master-slave alignment to the seed if we have a seed
  $seedmaster = 0
}


########################################################################
## If MSA is used as input find a seed from it to get the universe
########################################################################
if ($msafile) {
   # Read the alignment into an array where each element is a valid fasta entry
   @fa = fastaFileToArray($msafile);

   ########################################################################
   ## Identify the seed as the sequence with highest all-all pwid
   ########################################################################
   my $alignedSeed = identifySeed(@fa);

   ## Write the unaligned seed sequence to seed.fa
   my @seed = split(/\n/, $alignedSeed);
   $seed[1] =~ s/-//g;
   $seed[1] =~ s/\.//g;
   open(F, "> seed.fa");
   print F $seed[0] . "\n" . $seed[1] . "\n";
   close(F);
   $seedfile = "seed.fa";
   $seedid = fasta2id($seed[0]);
   $seedseq = $seed[1];
} else {
  # Can only use the msa to build a profile if we have an msa
  $msa_profile = 0;
}


# BDK open log file for dump
if($nooutlog == 0) { open(BDKLOGFILE, '> flowerpower.screenlog'); }
else { open(BDKLOGFILE, ">&=STDOUT"); }

########################################################################
## Report for seed file
########################################################################
if ($seedfile) {
  #-- Report length of seed ---

  if ($debug1)
  {
    print "\nParameters:\n";
    print "Pwid = $pwid\n";
    print "No. of SHMMs iterations = $niters\n";
  }

  $SEQLEN = length($seedseq);
  if ($FPMAXEVALUE == -1 )
  {
    if ($SEQLEN < 65)
    {
      $FPMAXEVALUE = 1.0e-2;
    }
    elsif ($SEQLEN <= 100)
    {
      $FPMAXEVALUE = 1.0e-3;
    }
    else
    {
      $FPMAXEVALUE = 1.0e-4;
    }
  }
  if ($debug1)
  {
    print "E-value threshold of SHMMs search = $FPMAXEVALUE\n";
    print "No. of PSIBLAST iterations = $PSIBLAST_ITERS\n";
    print "No. of hits returned per PSIBLAST iteration = $NPSIHITS\n";
    print "SW option = $SW_SCORE\n";
    print "Mode = $MODE\n";
    print "Maximum hits from FlowerPower= $FPHITS\n";
    if ($HCOV > 0 || $QCOV > 0)
    {
      print "Query Coverage = $QCOV\n";
      print "Hit Coverage = $HCOV\n";
    }
    else
    {
      print "Coverage parameter = $COV\n";
    }
  }
  print BDKLOGFILE "Seed ID: $seedid\n";
  print BDKLOGFILE "Length of seed : $SEQLEN\n";
  if ($MODE eq "global")
  {
    print BDKLOGFILE "Clustering protocol used: global-global (all sequences must be\n";
    print BDKLOGFILE "roughly the same length and align along their entire lengths)\n";
  }
  elsif ($MODE eq "glocal")
  {
    print BDKLOGFILE "Clustering protocol used: global-local (at least a portion of the\n";
    print BDKLOGFILE "hit sequences match the entire length of the query; hits may\n";
    print BDKLOGFILE "contain additional domains)\n";
  }
  if ($universe)
  {
    print BDKLOGFILE "Sequence database searched = $universe\n";
  }
  else
  {
    if ($DB =~ /UniProt/)
    {
       print BDKLOGFILE "Sequence database searched = UniProt\n\n";	  
    } 
    else
    {
       print BDKLOGFILE "Sequence database searched = $DB\n";
    }
  }
}

##############################################################################
## If MSA file, get fasta file of unaligned sequences and add to universe
##############################################################################
if ($msafile)
{
#   @fa = `cat $msafile` ;
#   my $seq = "";
   open (ALN, "< $msafile");
   open(MSAFA, "> msa.fa");
   my $alnseq = "";
   while (my $line = <ALN>)
   {
     if($line =~ /^>/)
     {
        if ($alnseq)
        {
           print MSAFA "$alnseq\n";
           $alnseq = "";
        }
        print MSAFA "$line";
     }
     else
     {
        $line =~ s/\n//g;
        $line =~ s/\s+//g;
        $line =~ s/\W//g;
        $line =~ tr/a-z/A-Z/;
        $alnseq .= $line;
     }
   }
print MSAFA "$alnseq\n";
close MSAFA;
close ALN;
}

################################################################################
# Check for universe file or create it by running PSIBLAST and using length 
# filter for global option.
################################################################################

if ($universe) {
    unless (-e $universe) {
       print "File for universe missing\n";
       exit 0;
    }
    if ($universe ne "universe.fa") {
       print "Copying universe to universe.fa\n";      
       `cp $universe universe.fa`;
    }
    if ($msafile) {
       `cat msa.fa >> universe.fa`;
    } 
    #print "Removing duplicate IDs, if any, from the universe\n";
    makeunique("universe.fa");

    # # DMS 08/03/09 - force seed into universe.
    `cat uniq-universe.fa $seedfile > universe.fa`;

    print "Formatting universe.fa\n";
    $return_value = 1;
    $ii = 0;
    while ( $return_value != 0 && $ii < $mult_cmd_tries ) {
       if ( $ii > 0 ) {
          print BDKLOGFILE "Re-running...\n";
       }
       $return_value = system("formatdb -o T -i universe.fa");
       $ii++;
    }
    if ($return_value != 0) {
       print "ERROR: formatdb exited with non-zero value $return_value\n";
       exit 1;
    }
    if($msafile) {
      if ($blast_again == 1) {
        blast("msa.fa", "pb");
      } else {
        `touch pb`;
      }
    }
    else {
        blast($seedfile, "pb");  
    }
} else {
    print "  Starting PSI-BLAST ($PSIBLAST_ITERS iterations) to retrieve candidate set... ";  
    if ($msa_profile) {
      blastpgp($seedfile, "pb", $msafile); 
    } else {
      blastpgp($seedfile, "pb", ""); 
    }
    print "done.\n";
    if ($debug1) {
      print "Creating Universe from PSIBLAST file\n";
    }
    createUniverse("pb", $seedfile); 
    if ($msafile) {
      `cat msa.fa >> universe.fa`;
       print "Making unique\n";
       makeunique("universe.fa");
    }
    `formatdb -o T -i universe.fa`;
 } 

########################################################################
# Read in the lengths of each of the sequences in the universe
########################################################################
my %protein_lengths = ();
my %protein_seqs = ();
my %protein_headers = ();
my @universe_fa=fastaFileToArray("universe.fa");
for (my $i=0; $i <= $#universe_fa; $i++)
{
  my $protein_header = $universe_fa[$i];
  $protein_header =~ s/\n.*$//;       # Delete everything after first line.
  my $protein_id = fasta2id($universe_fa[$i]);
  my $protein_seq = fasta2seq($universe_fa[$i]);
  $protein_lengths{$protein_id}=length($protein_seq);
  $protein_seqs{$protein_id}=$protein_seq;
  $protein_headers{$protein_id}=$protein_header;
}
$protein_headers{$seedid}=$seed_header;
$protein_seqs{$seedid}=$seedseq;
$protein_lengths{$seedid}=length($seedseq);

########################################################################
# Create first round of homologs
########################################################################
print "\n  Iteration 0: Selecting the first set of sequences...";
print BDKLOGFILE "***STARTING Iteration 0 ***\n";
if($msafile) {
   # Remove columns with > 50% gaps.
   trimMSA($msafile, "initial-sel.mus", 0.50);
}
else {
  print BDKLOGFILE "Selecting sequences for MUSCLE alignment\n";  
  if (createHomologs("pb", "universe.fa")) {
     if ($debug1) { 
      print "Creating HMM from seed\n";
     } 
      $return_value = 1;
      $ii = 0;
      while ( $return_value != 0 && $ii < $mult_cmd_tries ) {
         if ( $ii > 0 ) {
            print BDKLOGFILE "Re-running...\n";
         }
         $return_value = system("w0.5 $seedfile initial.mod >& w0.5.out");
         $ii++;
      }
      if ($return_value != 0) {
         print "ERROR: w0.5 exited with non-zero value $return_value\n";
         exit 1;
      }
     if ($debug1) {
      print "Aligning initial set of homologs to seed HMM\n"; 
     }
      $return_value = 1;
      $ii = 0;
      while ( $return_value != 0 && $ii < $mult_cmd_tries ) {
         if ( $ii > 0 ) {
            print BDKLOGFILE "Re-running...\n";
         }
         $return_value = system("align2model initial-homologs -i initial.mod -db blast-homologs.fa -sw $SW_SCORE -adpstyle 5 >& align2model.out");   
         $ii++;
      }
      if ($return_value != 0) {
         print "ERROR: align2model exited with non-zero value $return_value\n";
         exit 1;
      }

  } 
  else
  {
     print BDKLOGFILE "\n\nThis FlowerPower run is unable to find any homologs to your query that meet the \ncriteria you selected. You may want to resubmit the sequence with modified \nparameters to reduce stringency, or select a subregion of your sequence for a \nnew FlowerPower run. See the online help page for recommended parameter settings.\n";
    `cp submitted-seed.fa $mwd/final.a2m`;
    `cp psiblast-hits.fa $mwd/psiblast-hits.fa`;
    `cp universe.fa $mwd/initial-universe.fa`;
    print "FlowerPower was unable to find any homologs that meet the selected criteria.";
    exit 0;
  } 
   ##---Get the seedseq from initial-homologs.a2m to compare for pwid ----

     @fa = (); 
     @fa = fastaFileToArray("initial-homologs.a2m");
     my $aligned_seed = "";                  #--seed sequence in alignment
     for (my $i = 0; $i <= $#fa; $i++) {
        my $hitid = fasta2id($fa[$i]);
        if ($hitid eq $seedid || $hitid eq $unaltered_seedid) {
           $aligned_seed = $fa[$i];
           $SEED_SEQ = fasta2seq($aligned_seed);
           last;
        }
     }
    if ($debug1) {  
     print "Running QA on initial-homologs.a2m\n";
    } 
    # The length of the initial HMM is the same as that of the seed
    my $hmm_length = 0;
    if ($MODE eq "glocal") {
	$hmm_length = $SEQLEN;
    }
    else {
	$hmm_length = $protein_lengths{$seedid};
    }
    my $nseqs = 0;
     open (FILE, "> blast-sel.fa");
     for (my $i = 0; $i <= $#fa ; $i++) {
          if (($MODE eq "global") 
                 && testCoverageConservative($fa[$i], $hmm_length) 
                 && (pwid($aligned_seed, $fa[$i]) >= $initial_percent_id)) {
               #testInsert($seedseq, $fa[$i], 85))
               my $id = fasta2id($fa[$i]);
               my $idseq = $protein_seqs{$id};
               if($idseq) {
               } else {
                if($debug1) {
                  print "Could not find $id in universe\n";
                }
                 exit 0;
               }  
               print FILE "$protein_headers{$id}\n";
               print FILE "$idseq\n";
               #print FILE `fastacmd -s $id -d universe.fa`;
               $nseqs++;
          } elsif (($MODE eq "glocal") 
                      && testCoverageConservative($fa[$i], $hmm_length) 
                      && (pwid($aligned_seed, $fa[$i]) >= $initial_percent_id)) {
                 my ($def, $seq) = split(/\n/, $fa[$i]);
                 my $id = fasta2id($def);
                 if ($pfam == 0) {
                   $seq =~ s/\.//g;
                   if ($seedmaster == 0 || $id ne $seedid) {
                     # Remove initial (N-terminal) lowercase characters
                     $seq =~ s/^[a-z]+//g;
                     # Remove final (C-terminal) lowercase characters
                     $seq =~ s/[a-z]+$//g;
                     # Uppercase the remaining characters
                     $seq =~ tr/a-z/A-Z/;
                   }
                   $seq =~ s/\-//g;
                 }
                 $fa[$i] = $def . "\n" . $seq;
                 print FILE "$fa[$i]\n";
                 $nseqs++;
          } else {
            splice(@fa, $i, 1);
            $i--;
          }
     }
     close FILE;
     if ($debug1)
     {
       print  " ". ($#fa + 1) . " passed criteria\n";
     }

    print " done. " . ($#fa + 1) . " sequences were selected.\n";

    if ($nseqs <= 1)
    {
      print BDKLOGFILE "\n\nThis FlowerPower run is unable to find any homologs to your query that meet the criteria you selected. You may want to resubmit the sequence with modified parameters to reduce stringency, or select a subregion of your sequence for a new FlowerPower run. See the online HELP page for recommended parameter settings.\n";
      `cp submitted-seed.fa $mwd/final.a2m`;
      `cp psiblast-hits.fa $mwd/psiblast-hits.fa`;
      `cp universe.fa $mwd/initial-universe.fa`;
      print "FlowerPower was unable to find any homologs that meet the selected criteria.";
      exit 0;
    } 
    ########################################################################
    ## Create initial alignment file
    ########################################################################
    open(FILE, "> initial-sel.a2m");
    for (my $i = 0; $i <= $#fa ; $i++)
    {   
       print FILE  "$fa[$i]\n";
    }
    close FILE;
    if($nseqs > 1 && $pfam == 0)
    { 
       print BDKLOGFILE "Aligning initial set with MUSCLE\n";
       muscle("blast-sel.fa", "initial.mus");
       print BDKLOGFILE "Masking the alignment to remove columns with >50% gap characters\n";
       trimMSA("initial.mus", "initial-sel.mus", 0.50);
    } elsif ($pfam == 1) {
      `prettyalign blast-sel.fa -f > initial-sel.mus`;
    }
}


if($debug1)
{
  print "Creating fp-workspace\n\n";
}

 unless (-e 'fp-workspace') {
   `mkdir fp-workspace`;
 }
 `cp $seedfile fp-workspace`;
 `cp universe.fa fp-workspace`;
 `cp initial-sel.mus fp-workspace`;
  chdir "fp-workspace";
  `cp universe.fa initial-universe.fa`; 


########################################################################
## Start first iteration
########################################################################

  print "\n  Subsequent iterations will use SCI-PHY subfamily identification\n";
  print "  and subfamily HMM scoring to retrieve additional sequences.\n\n";

  my $count = 0;
  $count++;
  my $dirname = runIteration($count, "initial-sel.mus", $seedfile);
  print "  Iteration 1:";
  print BDKLOGFILE "\n\n*** STARTING $dirname *****\n"; 
  chdir $dirname;                         ##--- We are in "iter1" directory" ---
#  print "RUNNING SHMMS\n";
  runShmms("initial-sel.mus", $FPMAXEVALUE, $count);
  if ($debug1)
  {
    print "RUNNING QA on big.fa\n";
  }
  my $filename = "$dirname" .  "." . "a2m";
  
########################################################################
## Check best e-value from shmms.score and run QA based on that
#######################################################################

  my $best_evalue = testEvalue("shmms.score");
  if($debug1)
  {
    print "Best e-value was $best_evalue\n";
  }

  if (runQA($best_evalue, $pwid, "initial-sel.mus", %protein_lengths))
  {
     `cp big-final.a2m $filename`;
  }
  elsif (runQA($FPMAXEVALUE, $pwid, "initial-sel.mus", %protein_lengths))
  {
     `cp big-final.a2m $filename`;
  }
  else
  {
     print BDKLOGFILE "No more sequences added in $dirname. Finishing FlowerPower.\n";
     replacedefline("initial-sel.mus", $seedid, $defline);
     `cp initial-sel.mus $mwd/final.a2m`;
     `cp ../../psiblast-hits.fa $mwd/psiblast-hits.fa`; 
     `cp ../../universe.fa $mwd/initial-universe.fa`;
      print "No new sequences.\n";
      exit 0;
  }

  ## XX BDK - count number of new sequences added in iteration 1, and print that
  my $infacount = `grep ">" big-final.a2m | wc -l`;
  my $newseqcount = $infacount - ($#fa + 1);

  print " $newseqcount additional sequences retrieved.\n";

  `cp $filename ..`;


########################################################################
## Reduce universe size : remove hits with e-value >  100
########################################################################
 if ($SPEED eq "fast")  
 {
    my %id_evalue = ();
    open (FILE, "< shmms.score") || die "Could not open file shmms.score: $!\n" ;
    while (my $line = <FILE>){
        last if ( $line =~ /^align2model/ || $line =~/^\[/ );
         if ($line !~ /^\#/ && $line !~ /^\s+/)
         {
             chomp $line;
             (my $id, my $evalue) = (split (/\s+/, $line))[0,7];
             $id_evalue{ $id } = $evalue;
         }
     }

    close FILE;

   for(my $i=0; $i<=$#universe_fa; $i++)
   {
      my $id = fasta2id($universe_fa[$i]);

      # assign_seqs_to_shmms truncates IDs at 15 characters, so do so here, too.

      $id = substr( $id, 0, 15 );
      if ($id_evalue{$id} > 100) 
      {
         splice(@universe_fa, $i, 1);
         $i--; 
      }
   }
 
   open (FILE, "> universe.fa");
   for(my $i=0; $i<=$#universe_fa; $i++)
   {
      print FILE  "$universe_fa[$i]\n";
   }
   close FILE;
   if ($debug1)
   {
      print "Reduced universe contains  " . ($#universe_fa + 1) . " sequences\n";
   }
 } 
  `cp universe.fa ..`;	 

########################################################################
## Start next iteration
########################################################################
#  print "Finished ITER1 and starting next iteration\n";

  for (my $i = 1; $i < $niters; $i++)
   {
        $count++;
        chdir "..";
        $dirname = runIteration($count, $filename, $seedfile);
        print "  Iteration $count: ";
        print BDKLOGFILE "\n\n*** STARTING $dirname ****\n";  
        chdir $dirname;             ##--- We are in next iter directory --
        my $nseqs_univ = createNewUniverse($filename, "universe.fa");
        if ($nseqs_univ == 0) 
        {
          print BDKLOGFILE "Exhausted all sequences in PSIBLAST set. No more sequences to search; Stopping here. Bye!\n";
          print "No remaining candidate sequences.\n";
          last;
          #replacedefline($filename, $seedid, $defline);
          #`cp $filename $mwd/final.a2m`;
          #`cp ../../psiblast-hits.fa $mwd/psiblast-hits.fa`;
          #`cp ../../universe.fa $mwd/initial-universe.fa`;
          # exit 0;
        }

        print `formatdb -o T -i universe.fa` ;
        ###--- Determine number of sequences in the alignment.
         my $nseqs_iter = 0 ;
         #open (FILE, "< $filename");
         $return_value = 1;
         $ii = 0;
         while ( $return_value != 0 && $ii < $mult_cmd_tries )
         {
           if ( $ii > 0 )
           { 
             print BDKLOGFILE "Re-running...\n";
           }
           $return_value = system("make_nr_at_100_with_dict.py ${filename}-uniq $filename >& make_nr_at_100_with_dict.out");
           $ii++;
        }
        if ($return_value != 0)
        {
           print "ERROR: make_nr_at_100_with_dict.py exited with non-zero value $return_value\n";
           exit 1;
        }
        system("ln -s ${filename}-uniq_nr100.fa ${filename}-uniq.a2m");
        open (FILE, "< ${filename}-uniq.a2m");
        while (my $line = <FILE>)
        {
          if ($line =~ /^>/)
	  {
	    $nseqs_iter++ ;
          } 
        }
	if($debug1)
	{
          print "Number of sequences in alignment is $nseqs_iter\n";
          print "Number of FP hits to be returned is $FPHITS\n"; 
	}

        if ($nseqs_iter >= $FPHITS)
        {
          print BDKLOGFILE "Number of homologs exceeds $FPHITS; Stopping here. Bye!\n";

	  # # DMS 07/30/09 -- set a flag when the number of homologs exceeds FPHITS
	  $excess_homologs_flag = 1;

          last;
          #replacedefline($filename, $seedid, $defline);
          #`cp $filename $mwd/final.a2m`;
          #`cp ../../psiblast-hits.fa $mwd/psiblast-hits.fa`;
          #`cp ../../universe.fa $mwd/initial-universe.fa`;
          # exit 0;
        }
        my $alnfile = $filename;
        if ($debug1)
        {
           print "Running SHMMS score on $filename\n";
        } 
        runShmms($filename, $FPMAXEVALUE,$count);
        if ($debug1)
        {
          print "Running QA on big.fa\n";
        } 
        $best_evalue = testEvalue("shmms.score");
        my $pass = 0;  ###--- Counter to check if the runQA worked
        if ($best_evalue < $FPMAXEVALUE)
        {
           $best_evalue = ($best_evalue/0.1); 
          if($debug1)
          {
            print "Best e-value for $dirname is $best_evalue\n";
          } 
           while ($best_evalue <= $FPMAXEVALUE) ###--- Iterate till some evalue below evalue cut-off returns sequences. Else quit.
           {
             if (runQA($best_evalue, $pwid, $filename, 
                        %protein_lengths))
             {
                $filename = "$dirname" .  "." . "a2m";
                `cat $alnfile big-final.a2m > $filename` ;
                 print `prettyalign $filename -f > pretty.pa 2>prettyalign.out`;
                 print `mv pretty.pa $filename`;
                 $pass++;
                 last;
             }
             else 
             {   
                $best_evalue = ($best_evalue/0.1);
                if($debug1)
                {
                   print "Current Best evalue is $best_evalue\n";
                }
             }
           }
        } 
        if ($pass == 0) 
        {
            print BDKLOGFILE "No new sequences added in $dirname\n";
            print "No new sequences.\n";
            last;
        }

        else
        {
            ## XX BDK - count number of new sequences added in iteration 1, and print that
            my $infacount = `grep ">" big-final.a2m | wc -l`;
            my $newseqcount = $infacount / 1.0;
            #$newseqcount=~s/^\s+//;
            #$newseqcount=~s/\s+$//;
            print "$newseqcount additional sequences retrieved.\n";

        }

        `cp $filename ..`;
        `cp universe.fa ..`;
        if ( $debug1 ) 
        {
            print "DEBUG: Filename is $filename\n"; 
        }
   }  
 
########################################################################
## Final realignment of selected sequences
########################################################################
if($debug1)
{ 
  print "Executing the final realignment\n";
}
  # # DMS 07/30/09 - Remove homologs, if they number more than FPHITS.
  if (1 == $excess_homologs_flag) 
  {
      # print "Trimming $filename to $FPHITS sequences...\n";
      # print BDKLOGFILE "uncropped_seed_id: $uncropped_seed_id\n";
      # print BDKLOGFILE "filename: $filename\n;

      my $return_value = system("fp_top_by_pwid_with_seed.py $filename trimmed_last.a2m $FPHITS '$uncropped_seed_id'");
      if (0 != $return_value)
      {
	  print "ERROR: Trimming $filename failed";
	  if (512 == $return_value)
	  {
	      print " ( missing seed: $defline - $seedid )";
	  }
	  print "\n";
	  exit 1;
      }
      `cp trimmed_last.a2m ../last.a2m`;
  }
  else
  {
      `cp $filename ../last.a2m` ;
  }
  `cp $filename ../../untrimmed_last.a2m` ;
  # # END DMS EDIT

  chdir "..";
  unless (-e 'final') {
    `mkdir final`;
  }
  `cp last.a2m final`;
  `cp ../universe.fa final`;
   chdir "final";
   print `formatdb -o T -i universe.fa`;

  ########################################################################
  ## Get unaligned sequences from last.a2m
  ########################################################################
  open (FFILE, "> final.fa");
  if ($seedmaster) {
    open(SDFILE, "> seed_dups.fa");
  }
  @fa = fastaFileToArray("last.a2m");
  for (my $i=0; $i<=$#fa ; $i++)
  {
    my $id = fasta2id($fa[$i]);
    my $idheader = $protein_headers{$id};
    my $idseq = $protein_seqs{$id};
    # For master-slave alignment, the seed must be aligned 1-1 with the HMM
    if ($seedmaster == 0 || $idseq ne $seedseq) {
      print FFILE "$idheader\n";
      print FFILE "$idseq\n";
    } elsif ($seedmaster && $idseq eq $seedseq) {
      print SDFILE "$idheader\n";
      print SDFILE "$idseq\n";
    }
    #`fastacmd -s $id -d universe.fa >> final.fa`;
  }
  if ($seedmaster) {
    close SDFILE;
  }
  close FFILE;
  ########################################################################
  ## Running SCI-PHY one last time
  ########################################################################  
   $return_value = 1;
   $ii = 0;
   if (-e 'last.mod') {
     `rm last.mod`;
   }
   while ( $return_value != 0 && $ii < $mult_cmd_tries )
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("SCI-PHY last -i last.a2m >& SCI-PHY.out");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: SCI-PHY exited with non-zero value $return_value\n";
      exit 1;
   }

  ########################################################################
  ## Running Assign_seqs one last time
  ########################################################################  

   $return_value = 1;
   $ii = 0;
   while ( $return_value != 0 && $ii < $mult_cmd_tries )
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("assign_seqs_to_shmms -f final.fa --reuse --sw $SW_SCORE -d 5 --dbsize 100000 -b $USE_QUEUE --models *.mod > shmms.score");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: Assign_seqs exited with non-zero value $return_value\n";
      exit 1;
   }

  ########################################################################
  ## Copying the final files required 
  ########################################################################  
  if ($debug1)
  { 
    print "Copying final files to working directory\n";
  }
  if ($seedmaster) {
    `cat seed_dups.fa final.ass/big.fa > big_with_seed_dups.a2m`;
    `prettyalign big_with_seed_dups.a2m -f > ../final.a2m 2> prettyalign.big_with_seed_dups.out`;
  } else {
    `cp final.ass/big.fa ../final.a2m`;
  }
   replacedefline("../final.a2m", $seedid, $defline);
   my $finalnumber = `grep '>' ../final.a2m | wc -l`;
   chomp $finalnumber;
  `cp ../final.a2m $mwd/final.a2m`;
  if($temporary_check == 0)
  {
    `cp ../../psiblast-hits.fa $mwd/psiblast-hits.fa`;
  }
   my $fnum = `grep '>' $mwd/final.a2m | wc -l`;
   chomp $fnum;
   if ($fnum == 0) {
     `cp last.a2m $mwd/final.a2m`;
   }
  `cp ../../universe.fa $mwd/initial-universe.fa`;
   chdir $mwd ;  

  print "\nFlowerPower retrieved $finalnumber sequences.\n";

exit 0;

sub fasta2id
{
    my $fa = shift;
    my $id = $fa;

    $id =~ s/^> />/;        # Get rid of blank after ">".
    $id =~ s/\n.*$//;       # Delete everything after first line.
    $id =~ s/ .*$//g;       # Delete everything after first blank.
    $id =~ s/^[^|]+\|//;    # Delete non-pipe ("|") chars up through first pipe.
    $id =~ s/^>//;          # Delete the ">".
    $id =~ s/\|.*//;        # Delete pipe and everything after.
    
    if($MODE eq "global")
    {
        $id =~ s/\/.*//;        # Delete "/" and after.
    }

    $id =~ s/[\(\)\[\]]//g; # Delete parens and square brackets.
    $id =~ s/\r//g;         # Delete CRs.
    $id =~ s/,//g;          # Delete commas.

    return $id;
}

sub fasta2seq
{
    my $fa = shift;
    my $seq = $fa;

    $seq =~ s/[^\n]*\n//g;
    $seq =~ s/\n//g;

    return $seq;
}

sub testLen
{
    my $query   = shift;
    my $subject = shift;

    my $querySeq   = fasta2seq($query);
    my $subjectSeq = fasta2seq($subject);
    my $querylen = length ($querySeq);
    my $subjectlen = length($subjectSeq);
    my $covmin  = covmin($querylen, $initial_min_coverage);
   
    #if($QCOV > 0)
    #{ # user-supplied value overrides
    #  $covmin = $QCOV;
    #}

    if ($subjectlen >= ($covmin * $querylen))
    { # subject is long enough to cover query sequence
      if (($MODE eq "global") && (($subjectlen * $covmin) > $querylen))
      { # query is too short to cover subject, no good
        return 0; 
      }  
      else
      { # subject and query are okay
        return 1;
      }
    }

    return 0;
}

# ------------------------------------------------------------------------------
sub testCoverage {

   # See if this sequence meets coverage criteria.

   my $query    = shift;
   my $hmm_length = shift;

   my $queryseq = fasta2seq($query);
   my $queryid  = fasta2id($query);
   my $queryLen = $protein_lengths{$queryid};

   my $upper      = countUpper( $queryseq );
   # We require that the hit align to at least $qcov_check fraction of the HMM 
   my $qcov_check = covmin( $hmm_length, $COV ); 
   # In global-global mode, we require that at least $hcov_check fraction of the
   # hit align to the HMM
   my $hcov_check = covmin( $queryLen, $COV );
    
   if ($QCOV) {
      $qcov_check = $QCOV;  
   } 
   if ($HCOV) {
      $hcov_check = $HCOV;
   }

   if ( $debug1 ) {
      print ("DEBUG: $queryid\t$upper\t$queryLen\t$qcov_check\t$hcov_check\t" , (($hmm_length) *$qcov_check), "\t", (($queryLen)/$hcov_check),"\n"); 
   }

   # Do query coverage for both global and glocal mode.  
   # (aka "HMM coverage" = aligned chars / HMM length)

   if ( $upper < $hmm_length*$qcov_check ) {   

      # Fail: aligned chars / HMM length < qcov_check

      if ( $debug1 ) {
            print "DEBUG: $queryid failed query coverage test\n";
      }
      return 0;

   }
   # Check that no deleted area of the HMM is greater than $max_gap residues.
   if ( $max_gap > 0 ) {
     my $longest_deleted_region_len = longestGap( $queryseq );
     if ( $longest_deleted_region_len > $max_gap ) {
       if ( $debug1 ) {
         print "DEBUG: $queryid failed long deleted region test\n";
       }
       return 0;
     }
     if ( $debug1 ) {
       print "DEBUG: $queryid passed long deleted region test\n";
       print "DEBUG: longest deleted region: $longest_deleted_region_len\n";
     }
   }

   # Do hit-coverage test only for the "global" mode.  
   # (??? aka "Sequence coverage = aligned chars / sequence length)

   my $longest_inserted_region_len = 0;

   if ($MODE eq "global") {
      if ( $upper < $queryLen*$hcov_check ) {

         # Fail: aligned chars / sequence length < qcov_check

         if ( $debug1 ) {
            print "DEBUG: $queryid failed hit coverage test\n";
         }
      return 0;
      }

      # Also in global mode, if flag set, check that no inserted area of the
      # query is greater than $max_insert residues.

      if ( $max_insert > 0 ) { 
         $longest_inserted_region_len = longestLowercase( $queryseq );
         if ( $longest_inserted_region_len > $max_insert ) {
            if ( $debug1 ) {
              print "DEBUG: $queryid failed long inserted region test\n";
            }
            return 0;
         }
      }
   } else {
      # Check that no inserted area of the query that lies within the aligned
      # region is greater than $max_insert residues.
      if ( $max_insert > 0 ) {
         $longest_inserted_region_len = longestInternalLowercase($queryseq);
         if ( $longest_inserted_region_len > $max_insert ) {
            if ( $debug1 ) {
              print "DEBUG: $queryid failed long internal inserted region test\n";
            }
            return 0;
         }
      }
   }
    
   if ( $debug1 ) {
      print "DEBUG: $queryid passed coverage test\n";  
      if ($MODE eq "global") {
        print "DEBUG: longest inserted region: $longest_inserted_region_len\n";
      } else {
        print "DEBUG: longest internal inserted region: ";
        print "$longest_inserted_region_len\n";
      }
   }
   return 1;
}

sub testCoverageConservative
{
   # See if this sequence meets coverage criteria.

   my $query    = shift;
   my $hmm_length = shift;

   my $queryseq = fasta2seq($query);
   my $queryid  = fasta2id($query);
   my $queryLen = $protein_lengths{$queryid};

   my $upper      = countUpper( $queryseq );
   # We require that the hit align to at least $qcov_check fraction of the HMM 
   my $qcov_check = covmin( $hmm_length, $initial_min_coverage ); 
   # In global-global mode, we require that at least $hcov_check fraction of the
   # hit align to the HMM
   my $hcov_check = covmin( $queryLen, $initial_min_coverage );
    
   if ($QCOV) {
      $qcov_check = $QCOV;  
   } 
   if ($HCOV) {
      $hcov_check = $HCOV;
   }

   if ( $debug1 ) {
      print ("DEBUG: $queryid\t$upper\t$queryLen\t$qcov_check\t$hcov_check\t" , (($hmm_length) *$qcov_check), "\t", (($queryLen)/$hcov_check),"\n"); 
   }

   # Do query coverage for both global and glocal mode.  
   # (aka "HMM coverage" = aligned chars / HMM length)

   if ( $upper < $hmm_length*$qcov_check ) {   

      # Fail: aligned chars / HMM length < qcov_check

      if ( $debug1 ) {
            print "DEBUG: $queryid failed query coverage test\n";
      }
      return 0;

   }
   # Check that no deleted area of the HMM is greater than $max_gap residues.
   if ( $initial_max_gap > 0 ) {
     my $longest_deleted_region_len = longestGap( $queryseq );
     if ( $longest_deleted_region_len > $initial_max_gap ) {
       if ( $debug1 ) {
         print "DEBUG: $queryid failed long deleted region test\n";
       }
       return 0;
     }
     if ( $debug1 ) {
       print "DEBUG: $queryid passed long deleted region test\n";
       print "DEBUG: longest deleted region: $longest_deleted_region_len\n";
     }
   }

   # Do hit-coverage test only for the "global" mode.  
   # (??? aka "Sequence coverage = aligned chars / sequence length)

   my $longest_inserted_region_len = 0;

   if ($MODE eq "global") {
      if ( $upper < $queryLen*$hcov_check ) {

         # Fail: aligned chars / sequence length < qcov_check

         if ( $debug1 ) {
            print "DEBUG: $queryid failed hit coverage test\n";
         }
      return 0;
      }

      # Also in global mode, if flag set, check that no inserted area of the
      # query is greater than $max_insert residues.

      if ( $initial_max_insert > 0 ) { 
         $longest_inserted_region_len = longestLowercase( $queryseq );
         if ( $longest_inserted_region_len > $initial_max_insert ) {
            if ( $debug1 ) {
              print "DEBUG: $queryid failed long inserted region test\n";
            }
            return 0;
         }
      }
   } else {
      # Check that no inserted area of the query that lies within the aligned
      # region is greater than $max_insert residues.
      if ( $initial_max_insert > 0 ) {
         $longest_inserted_region_len = longestInternalLowercase($queryseq);
         if ( $longest_inserted_region_len > $initial_max_insert ) {
            if ( $debug1 ) {
              print "DEBUG: $queryid failed long internal inserted region test\n";
            }
            return 0;
         }
      }
   }
    
   if ( $debug1 ) {
      print "DEBUG: $queryid passed coverage test\n";  
      if ($MODE eq "global") {
        print "DEBUG: longest inserted region: $longest_inserted_region_len\n";
      } else {
        print "DEBUG: longest internal inserted region: ";
        print "$longest_inserted_region_len\n";
      }
   }
   return 1;

}


# ------------------------------------------------------------------------------
sub longestLowercase {
  # Find length of longest region of continuous lowercase.
  # A gap or dot does not end the continuous lowercase region.
  # Example: seed (=HMM): AC......EFGH.KLMN
  #          query:       ACpqrstv----a-LMN
  # The length of this lowercase region is 7.
  my $querystr = shift;
  # Remove dashes and dots
  $querystr =~ s/[\.\-]//g;
  # Replace consecutive runs of uppercase characters by '#'
  $querystr =~ s/[A-Z]+/#/g;
  # Split along '#'
  # The split regions will consist of consecutive runs of lowercase
  my @lowercase_regions = split(/#/, $querystr);
  my $longest_lowercase_region_length = 0;
  my $lowercase_region_length = 0;
  foreach (@lowercase_regions) {
    $lowercase_region_length = length($_);
    if ($lowercase_region_length > $longest_lowercase_region_length) {
      $longest_lowercase_region_length = $lowercase_region_length;
    }
  }
  return $longest_lowercase_region_length;
}


# ------------------------------------------------------------------------------
sub longestInternalLowercase {
  # Find length of longest region of continuous lowercase within the aligned
  # region.  A gap or dot does not end the continuous lowercase region.
  # Example: seed (=HMM): ..........AC......EFGH.KLMN...
  #          query:       alvmlvstlvACpqrstv----a-LM-haw
  # The length of the longest continuous internal lowercase region is 7.
  my $querystr = shift;
  # Remove dots
  $querystr =~ s/\.//g;
  # Remove initial (N-terminal) lowercase characters
  $querystr =~ s/^[a-z]+//g;
  # Remove final (C-terminal) lowercase characters
  $querystr =~ s/[a-z]+$//g;
  # Remove dashes
  $querystr =~ s/\-//g;
  # Replace consecutive runs of uppercase characters by '#'
  $querystr =~ s/[A-Z]+/#/g;
  # Split along '#'
  # The split regions will consist of consecutive runs of lowercase
  my @lowercase_regions = split(/#/, $querystr);
  my $longest_lowercase_region_length = 0;
  my $lowercase_region_length = 0;
  foreach (@lowercase_regions) {
    $lowercase_region_length = length($_);
    if ($lowercase_region_length > $longest_lowercase_region_length) {
      $longest_lowercase_region_length = $lowercase_region_length;
    }
  }
  return $longest_lowercase_region_length;
}

# ------------------------------------------------------------------------------
sub longestGap {
  # Find length of longest region of continuous gaps.
  # A dot or lowercase character not end the continuous region of gaps.
  # Example: seed (=HMM): AC......EFGH.KLMN
  #          query:       ACpqrstv----a-LMN
  # The length of this region of gaps is 5.
  my $querystr = shift;
  # Remove dots and lowercase characters
  $querystr =~ s/[\.a-z]//g;
  # Replace consecutive runs of uppercase characters by '#'
  $querystr =~ s/[A-Z]+/#/g;
  # Split along '#' 
  # The split regions will consist of consecutive runs of dashes
  my @gappy_regions = split(/#/, $querystr);
  my $longest_gapped_region_length = 0;
  my $gapped_region_length = 0;
  foreach (@gappy_regions) {
    $gapped_region_length = length($_);
    if ($gapped_region_length > $longest_gapped_region_length) {
      $longest_gapped_region_length = $gapped_region_length;
    }
  }
  return $longest_gapped_region_length;
}

# ------------------------------------------------------------------------------
sub testEvalue
{
  my $infile = shift;
  my $cnt10 = 0;
  my $cnt8 = 0;
  my $cnt6 = 0;
  my $evalue = 0; 
 
 open (FILE, "< $infile") || die "Could not open file $infile: $!\n" ;
  while (my $line = <FILE>)
  {
      last if ( $line =~ /^align2model/ || $line =~/^\[/ );
      if ($line !~ /^\#/ && $line !~ /^\s+/)
      {
          chomp $line;
          ($evalue) = (split (/\s+/, $line))[7];
         if ($evalue <= 1e-10)
         {
            $cnt10++;
         } 
         if ($evalue <= 1e-08)
          {
             $cnt8++;
          }
         if ($evalue <= 1e-06)
          {
             $cnt6++;
          }
      }
  }

 close FILE;
  
  if ($cnt10 >= 20)
  {
    $evalue = 1e-10;
  }
  elsif ($cnt8 >= 20)                 
   {
     $evalue = 1e-08;
   }
  elsif ($cnt6 >= 20)                 
   {
     $evalue = 1e-06;
   } 
  else                 
   {
     $evalue = 1e-05;
   }
   return $evalue;
}


sub countUpper
{
    my $seq = shift;

    return ($seq =~ s/([A-Z])/$1/g);
}

sub countLower
{
    my $seq = shift;

    return ($seq =~ s/([a-z])/$1/g);
}

sub countUpperLower
{
    my $seq = shift;
    return ($seq =~ s/([A-Za-z])/$1/g);
}

sub countDashes
{
    my $seq = shift;
    return ($seq =~ s/(\-)/$1/g);
}


sub testPwid
{
    my $query   = shift;
    my $subject = shift;
    my $cutoff_pwid = shift;
    my $subject_id = fasta2id($subject);
    my $pwid = pwid($query, $subject);
    if ($pwid < $cutoff_pwid)
    {
        if ( $debug1 )
        {
           print "DEBUG: $subject_id failed pwid test\n"; 
        }
        return 0;
    }
    if ( $debug1 )
    {
        print "DEBUG: $subject_id passed pwid test. Pwid: $pwid\n"; 
    }
    return 1;
}

sub covmin {   
# This subroutine returns a real value, which is the minimum required
# fractional overlap, given the sequence length. [.55 - .85]

    my $len  = shift; 
    my $criterion = shift;

    my $baseline_prop = 0.60;

    if ($criterion == 1) {
        $baseline_prop = 0.55;
    }

    elsif ($criterion == 2) {
        $baseline_prop = 0.50;
    }

    elsif ($criterion == 3) {
        $baseline_prop = 0.45;
    }

    elsif ($criterion == 4) {
        $baseline_prop = 0.40;
    }

    if    ($len < 100) { return $baseline_prop; }
    elsif ($len < 200) { return $baseline_prop + 0.05; }
    elsif ($len < 250) { return $baseline_prop + 0.10; }
    elsif ($len < 300) { return $baseline_prop + 0.13; }
    elsif ($len < 350) { return $baseline_prop + 0.15; }
    elsif ($len < 400) { return $baseline_prop + 0.18; }
    elsif ($len < 450) { return $baseline_prop + 0.20; }
    elsif ($len < 500) { return $baseline_prop + 0.23; }
    return $baseline_prop + 0.25;
}

sub checkUsage {

    unless ($ARGV[0] && -e $ARGV[0])
    {
        print "Usage: klusterexpand KLUSTERDIR\n";
        exit;
    }
}

sub fastaFileToArray
{
    my $fn = shift;

    my @s;

    open(F, $fn);

    while (<F>)
    {
        chomp;
        s/\n//g;
        s/\r//g;

        if (/^>/)
        {
            push(@s, $_ . "\n");
        }
        elsif ($_ !~ //)
        {
            $s[$#s] .= $_;
        }
    }

    close(F);

    return @s;
}


sub averageAllAllPwidPerSequence
{
    ## Only arg is a fasta array
    my @fa = @_;

    ## Store each pwid comparison here (keys are row #'s)
    my %pwids;

    if ($#fa == 0)
    {
        $pwids{$#fa} = 1.0;
        return %pwids;
    }

    ## Initialize %pwids to zero for each row
    for (my $row = 0; $row <= $#fa; $row++)
    {
        $pwids{$row} = 0;
    }

    ## Compute pwid for lower diagonal
    for (my $row = 0; $row <= $#fa; $row++)
    {
        for (my $sub = $row + 1; $sub <= $#fa; $sub++)
        {
            my $pwid = pwid($fa[$row], $fa[$sub]);
            $pwids{$row} += $pwid;
            $pwids{$sub} += $pwid;
        }
    }

    ## Average the pwid totals
    foreach my $key (keys %pwids)
    {
        $pwids{$key} /= $#fa;
    }

    return %pwids;
}

sub pwid
{
    my $astr = shift;
    my $bstr = shift;
    my $pwid = 0;

    $astr = fasta2seq($astr);
    $bstr = fasta2seq($bstr);
    
    my $aln = $astr;
       $aln =~ s/\.//g;
       $aln =~ s/[a-z]//g;
    my $aln_len  = length($aln);
    
    my $bln = $bstr;
       $bln =~ s/\.//g;
       $bln =~ s/[a-z]//g;
   
    my @a = split(//, $aln);
    my @b = split(//, $bln);

    if ($#a != $#b)
    {
        if ( $debug1 )
        {
            print "DEBUG: Lengths do not match\n";
        }
        return 0.0;
    }

    # $aligns counts the number of columns of the alignment where the sequences
    # both have uppercase characters
    my $aligns  = 0;
    # $matches counts the number of columns of the alignment where the
    # sequences have identical uppercase characters
    my $matches = 0;

    for (my $col = 0; $col <= $#a; $col++)
    {
        my $aaa = $a[$col];
        my $aab = $b[$col];
        if ($a[$col] ne '-' && $b[$col] ne '-')
        {
            $aligns++;
            if ($a[$col] eq $b[$col])
            {
                $matches++;
            }
        }
    }
   
   if ($MODE eq "global")
   {
      $pwid = $matches/$aln_len;
      return $pwid;
   }
   if ($MODE eq "glocal")
   {
     if ($aligns > 0)
     {
        $pwid = $matches/$aligns ;
     }
     else
     {
       $pwid = 0;
     } 
     return $pwid;
   }
}

sub isUpper
{
    my $c = shift;

    return ($c =~ /[A-Z]/);
}

sub isLower
{
    my $c = shift;

    return ($c =~ /[a-z]/);
}

sub muscle
{
    my $ifn = shift;
    my $ofn = shift;

    my $cmd = "muscle -in $ifn -out $ofn ";
    $cmd .= "-maxiters $MUSCLE_MAXITERS >/dev/null 2>&1";
    print `$cmd`;
}

sub identifySeed
{
    my @fa = @_;

    ## Get the all-all average pwid per seq
    my %pwids = averageAllAllPwidPerSequence(@fa);

    ## Identify the seq with the highest pwid (call this the seed)
    my $pwidSum     = 0.0;
    my $highestPwid = 0.0;
    my $highestKey  = 0;

    foreach my $key (keys %pwids)
    {
        $pwidSum += $pwids{$key};

        if ($pwids{$key} > $highestPwid)
        {
            $highestPwid = $pwids{$key};
            $highestKey  = $key;
        }
    }

    return $fa[$highestKey];
}

# ------------------------------------------------------------------------------
sub binspec
{
  my $seed_length = shift;
  my $remainder = $max_length_to_bin % $bin_interval;
  my $floor_num_bins = ($max_length_to_bin - $remainder) / $bin_interval;
  my $num_narrow_bins = $floor_num_bins;
  if ($remainder != 0) {
    $num_narrow_bins += 1;
  }
  my $first_length_in_wide_bin = $num_narrow_bins * $bin_interval + 1;
  if ( $seed_length >= $first_length_in_wide_bin ) {
    return sprintf("_%04dandmore", $first_length_in_wide_bin);
  }
  my $rem = $seed_length % $bin_interval;
  my $i = ($seed_length - $rem) / $bin_interval;
  my $first_seed_length;
  my $last_seed_length;
  if ($rem == 0) {
    $first_seed_length = ($i-1) * $bin_interval + 1;
    $last_seed_length = $seed_length;
  } else {
    $first_seed_length = $i * $bin_interval + 1;
    $last_seed_length = ($i+1) * $bin_interval;
  }
  return sprintf("_%04dthrough%04d", $first_seed_length, $last_seed_length);
}

# ------------------------------------------------------------------------------
sub blastpgp
{
    my $ifn = shift;
    my $ofn = shift;
    my $input_msa = shift;
    my $nb = 500;
    my $nv = 500;

    if (($NPSIHITS > 500) && ($NPSIHITS <= 3000))
    {
      $nb = $NPSIHITS;
      $nv = $NPSIHITS;
    }      

    my $cmd = "blastpgp";
    if ($input_msa ne "") {
     $return_value = 1;
     $ii = 0;
     while ( $return_value != 0 && $ii < $mult_cmd_tries ) 
     {
        if ( $ii > 0 )
        {
           print BDKLOGFILE "Re-running...\n";
        }
        $return_value = system("w0.5 $input_msa msa_profile.mod");
        $ii++;
     }
     if ($return_value != 0)
     {
        print "ERROR: w0.5 in blastpgp exited with nonzero value $return_value\n";
        exit 1;
     }
     system("model_convert.pl msa_profile.mod msa_profile.psiblast");
      
     $cmd .= " -R msa_profile.psiblast ";
    }
    
    $cmd .= " -i $ifn -o $ofn -m 9 -j $PSIBLAST_ITERS -I T -b $nb "
                                                       . "-v $nv -e $PBEVALUE ";
    $cmd .= "-d $DB -F F";
    if ($use_bins) {
      $cmd .= binspec($SEQLEN);
    }

   $return_value = 1;
   $ii = 0;
   while ( $return_value != 0 && $ii < $mult_cmd_tries ) 
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("$cmd");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: PSIBLAST exited with non-zero value $return_value\n";
      exit 1;
   }
}


# ------------------------------------------------------------------------------
sub blast
{
    my $ifn = shift;
    my $ofn = shift;
    my $cmd = "blastall -p blastp -i $ifn -o $ofn -m 9 -e 100 -b 1000 -v 1000 " 
                                                       . "-F F -d universe.fa";

    print `$cmd`;
}

sub testQueryCoverage
{
    my $query   = shift;
    my $subject = shift;

    my $querySeq   = fasta2seq($query);
    my $subjectSeq = fasta2seq($subject);

    my $qcov = countUpper($subjectSeq) / countUpperLower($querySeq);

    if ($qcov < covmin($querySeq, $COV))
    {
        return 0;
    }

    return 1;
}

sub testHitCoverage
{
    my $subject = shift;

    #print "Subject: $subject\n";

    my $subjectId = fasta2id($subject);

    print `fastacmd -d $DB -s $subjectId > $subjectId.fa 2>/dev/null`;

    unless (`grep '>' $subjectId.fa`)
    {
        print "Could not find $subjectId in UniProt...checking local\n";
       my $cmd = "formatdb -i acceptedseqs-orig.fa -o T -l /dev/null";
        print `$cmd`;
        print `fastacmd -d acceptedseqs-orig.fa -s $subjectId > $subjectId.fa`;
    }

    my @subjectFa = fastaFileToArray("$subjectId.fa");
    `rm $subjectId.fa`;
    my $fullSubjectSeq = fasta2seq($subjectFa[0]);

    my $subjectSeq = fasta2seq($subject);

    my $hcov = countUpper($subjectSeq) / countUpperLower($fullSubjectSeq);

    if ($hcov < covmin($fullSubjectSeq, $COV))
    {
        #print "Failed hit coverage: " . $hcov . " < " . covmin($subjectSeq) . "\n";
        return 0;
    }

    #print "Passed hit coverage: " . $hcov . " > " . covmin($subjectSeq) . "\n";
    return 1;
}


sub createUniverse
{
   my $psiblast_file = shift;
   my $seed_file = shift;

######################################################################################
## Perform E-value QA on the PSI-BLAST results : remove any sequence with
## e-value >= $PBEVALUE
######################################################################################

  my %id_evalue = ();
  my %id_unique = ();
  open (FILE, "< pb");
  open (PB_ID_FILE, "> psiblast-hits.id");
  while (my $line = <FILE>)
  {
     if ($line !~ /^#/)
     {
        chomp $line;
        (my $id, my $evalue) = (split(/\s+/, $line))[1,10];
        my $id_munge = $id ; 
        if ($id =~ /\|/)
        {
          $id_munge = (split(/\|/, $id))[1];
        }
        if(!exists($id_unique{$id_munge}))
        {
          $id_unique{$id_munge} = 1;
          print PB_ID_FILE "$id\n";
        }
     }
  }
  close PB_ID_FILE;
  close FILE;
  my $nseqs = keys(%id_unique);
  print BDKLOGFILE "$nseqs sequences retrieved\n";
  print `fastacmd -i psiblast-hits.id -o psiblast-hits.fa -d $DB`;
  print `cp psiblast-hits.fa blast.fa`;

########################################################################
## Perform Length QA on the PSI-BLAST results
########################################################################

 my @fa = fastaFileToArray("psiblast-hits.fa");
   for (my $i = 0; $i <= $#fa; $i++)
   {
     if (testLen($seedseq, $fa[$i]))
     {
     }
     else
     {
        splice(@fa, $i, 1);
        $i--;
     }
   }
   print BDKLOGFILE "" . ($#fa + 1) . " sequences from PSI-BLAST run passed initial criteria\n";

   open(F, "> blastlen.fa");
   foreach my $seq (@fa)
   {
     my @seq = split(/\n/, $seq);
     $seq[1] =~ s/-//g;
     $seq[1] =~ s/\.//g;
     print F "$seq[0]\n$seq[1]\n";
  }
  close(F);
 

 ##--- Check if seed is already present in blast results --

 my $seed_id = fasta2id($seed[0]);
 my $count = 0;
 @fa = (); 
 @fa = fastaFileToArray("blastlen.fa");
 for (my $i = 0; $i <= $#fa; $i++)
 {
   my $check_id = fasta2id($fa[$i]);
   if ($check_id eq $seed_id)
   {
     $count++;
     if($debug1)
     {
      print "Seed Id is $seed_id\nCheck Id is $check_id\n";
     }
     print BDKLOGFILE "Seed is present in PSI-BLAST set\n\n";
    `mv blastlen.fa universe.fa`;
   }
 }
  if($count == 0)
  {   
     print BDKLOGFILE "Seed is not present in PSI-BLAST set. Adding seed\n";
     `cat seed.fa blast.fa > universe.fa`;
      print `formatdb -o T -i universe.fa`;
  }
  if ($debug1)
  {
   print "Created universe with ". ($#fa + (1-$count)) . " sequences\n";	 
  }
}

# ------------------------------------------------------------------------------
sub createHomologs {
  ##--- Read in first iteration of psiblast : IDs and e-values to a hash

 my $psiblast_file = shift; #psiblast output file
 my $nseqs = 0; 
 my @loose_initial_homologs_fa = ();
 my @strict_initial_homologs_fa = ();
 my $fasta_record = "";
 my %id_evalue = ();
 open (FILE, "< $psiblast_file");
 while (my $line = <FILE>)
 { 
    last if ($line =~ /^# Iteration: 2/);
    if ($line !~ /^#/)
    {  
      (my $id, my $evalue) = (split(/\s+/, $line))[1,10];
      if ($evalue <= 1e-05)
      {
        if($id =~ /\|/)
        { 
          $id  = (split(/\|/, $id))[1]; 
        }
        # Check that the id exists in the universe
        # In particular, it may not have passed testLen
        if (exists ($protein_seqs{$id})){
          $nseqs++;
          if ($debug1)
          {
             print "$id : $evalue\n";
          } 
          $id_evalue{$id} = $evalue;
          
          $fasta_record = "$protein_headers{$id}\n$protein_seqs{$id}\n";
          push(@loose_initial_homologs_fa, $fasta_record);
          if ($evalue <= 1e-5) {
            push(@strict_initial_homologs_fa, $fasta_record);
          }
        }
      }
    }
    last if ($nseqs >= 200);
 }
 close FILE; 

  if ($debug1) {
    print "DEBUG: No. of seqs from first iteration is $nseqs\n";
  }
  if (@loose_initial_homologs_fa > 0) {
    open (FILE, "> blast-homologs.fa");
    foreach my $fasta_record (@loose_initial_homologs_fa) {
      print FILE "$fasta_record";
    }
    close(FILE); 
    if (!exists($id_evalue{$seedid})) {
      `cat seed.fa >> blast-homologs.fa`;
    }
    if ($debug1) {
      print "First set of " . ($#loose_initial_homologs_fa + 1) . 
        " homologs found.\n"; 
    }
    return 1;
  } else {
    if ($debug1) { 
       print "No blast homologs retrieved at 1e-05.\n";
    }
    return 0;
  }
}

# ------------------------------------------------------------------------------
sub runIteration {
  my $count = shift;
  my $aln = shift;
  my $seedfile = shift;

  my $dirname = "iter" . "$count";

  unless (-e $dirname) {
    `mkdir $dirname`;
  }
  `cp $aln $dirname`;
  `cp $seedfile $dirname`;
  `cp universe.fa $dirname`;
   return $dirname;

}

# ------------------------------------------------------------------------------
sub runShmms {
  my $aln = shift;  ##--- Input alignment for SCI-PHY
  my $evalue = shift;  ##--- e-value cutoff for assignseqs
  my $count = shift;   ##--- counter to check for which iteration in order to remove singletons

  ################################################################################
  ## Edit input alignment to remove unique sequences
  ################################################################################
if ($debug1)
{
   print "Making alignment non-redundant at 100%\n";
}
   $return_value = 1;
   $ii = 0;
   while ( $return_value != 0 && $ii < $mult_cmd_tries ) 
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("make_nr_at_100_with_dict.py trimmed-uniq $aln >& make_nr_at_100_with_dict.out");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: make_nr_at_100_with_dict exited with non-zero value $return_value\n";
      exit 1;
   }
   system("ln -s trimmed-uniq_nr100.fa trimmed-uniq.a2m");

#   system("cp $aln trimmed-uniq.a2m");
   
   #######################################################################
   ## Run SCI-PHY on trimmed-uniq.a2m
   ## Try several times.
   ########################################################################
   print BDKLOGFILE "Running SCI-PHY on MSA from previous iteration to identify subfamilies\nConstructing subfamily HMMs\n";
   #`cp $aln trimmed-uniq.a2m`;  ### Just doing this to keep rest of the names happy!!
   if (-e 'trimmed-uniq.mod') {
     `rm trimmed-uniq.mod`;
   }
   $return_value = 1;
   $ii = 0;
   while ( $return_value != 0 && $ii < $mult_cmd_tries ) 
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("SCI-PHY trimmed-uniq  -i trimmed-uniq.a2m > SCI-PHY.out 2>&1");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: SCI-PHY exited with non-zero value $return_value\n";
      exit 1;
   }

   my $subfam_no = `grep "%" trimmed-uniq.subfam |wc -l`;
   chomp $subfam_no;
   if ($debug1)
   {
     print "   $subfam_no subfamilies found\n";
   }

   ######################################################################################
   ## From Iteration 2 :  Check number of sequences per subfamily and remove singletons
   ######################################################################################
   if ($count >= 2)
   {
      my $subfam = "";
      my $subfam_count = 0;
      my $seq_count = 0;
      my %subfam_seq = ();

      open(FILE, "< trimmed-uniq.subfam");
      while (my $line = <FILE> )
      {
         if ( $line =~ /^%/ )
           {   
              if ($seq_count > 0)
              {
                  $subfam_seq{$subfam} = $seq_count; 
              } 
              chomp $line;
              ($subfam) = (split(" ", $line))[1];

               #-- Increment count
               $subfam_count++;
               $seq_count = 0;
           } 
         elsif ($line =~ /^>/) 
          {
              $seq_count++;
          }
      }

     $subfam_seq{$subfam} = $seq_count;

     close(FILE);
    if ($debug1)
    {
      print "Found $subfam_count subfamilies in specified file\n";
    } 

     while ((my $first, my $last) = each(%subfam_seq))
      {
           if ($last == 1)
            {
              if ($debug1)
              { 
                 print "Removing $first from list of models\n";
              } 
               print `mv trimmed-uniq.$first.mod trimmed-uniq.$first.hmm`;
            }
      }  
   }

  run_assign_seqs_to_shmms($evalue);
}

sub run_assign_seqs_to_shmms
{
   my $evalue = shift;

   ########################################################################
   ## Run assignseqsshmms  
   ## Try several times in hope of getting around "Stale file handle" errors.
   ## The -b ("bigalign") option to assign_seqs_to_shmms causes it to align
   ## each sequence to its closest SHMM using align2model --adpstyle 5
   ## All the aligned sequences meeting the e-value cutoff (-c $evalue)
   ## are written to big.a2m, which is run through prettyalign to create big.fa
   ########################################################################

   print BDKLOGFILE "Scoring the PSI-BLAST selected set with subfamily HMMs\nAligning these sequences to their closest subfamily HMM\n";
   $return_value = 1;
   $ii = 0;
   while ( $return_value != 0 && $ii < $mult_cmd_tries ) 
   {
      if ( $ii > 0 )
      {
         print BDKLOGFILE "Re-running...\n";
      }
      $return_value = system("assign_seqs_to_shmms -f universe.fa --reuse --sw $SW_SCORE -d 5 -c $evalue --dbsize 100000 -b $USE_QUEUE --models *.mod > shmms.score 2>assignseqs.log");
      $ii++;
   }
   if ($return_value != 0)
   {
      print "ERROR: assign_seqs_to_shmms exited with non-zero value $return_value\n";
      exit 1;
   }
}


sub runQA
{
########################################################################
## Run QA on big.fa
########################################################################
    
###---- Parse shmms.scores to a hash with ID as key and evalue as reference
if ($debug1)
{
    print "    Parsing shmms.score file\n";
}

    my $cutoff_evalue = shift;   ##--- Evalue cut off for comparing to hit e-value
    my $cutoff_pwid = shift;     ##--- Pwid cutoff for comparing to hit pwid
    my $aln_file = shift;        ##--- This is the alignment file used as input to SCI-PHY. 
    my $id = " ";                ##--- This is id of hit to check for e-value and reduce universe
    my $hmm_len = 0;             ##--- Length of SHMM to which hit is assigned
    my $evalue = 0;              ##--- Evalue of hit
    my %id_evalue = ();          ##--- hash that stores hitid to evalue map
    my %hmm_lengths = ();        ##--- hash that stores hitid to HMM length map

    open (FILE, "< shmms.score") || die "Could not open file shmms.score: $!\n" ;
    while (my $line = <FILE>){
         last if ($line =~ /^align2model/);
         if ($line !~ /^\#/ && $line !~ /^\s+/){
             chomp $line;
             ($id, $hmm_len, $evalue) = (split (/\s+/, $line))[0,2,7];
             $id_evalue{ $id } = $evalue;
             $hmm_lengths{ $id } = $hmm_len;
         }
     }
    close FILE;

###---- Create a hash with IDs from input alignment to SCI-PHY. These IDs do not go through Coverage and other tests as they are already accepted.

  my %check_aln_id = ();
  my @fa = ();
  @fa = fastaFileToArray($aln_file);
  for (my $i = 0; $i <= $#fa; $i++)
  {
     my $aln_id = fasta2id($fa[$i]);
     $check_aln_id{$aln_id} = 1 ;
  } 
   
  my $bigaln = "universe.ass/big.fa";
  my $hitid = " ";  ###--- Hit IDs in big.fa
  my %idhash= ();   ###--- This hash stores the ids of sequences that passed QA. This will be used to lookup the seqs from big.fa
if ($debug1) 
{
  print "    Doing QA on big.fa\n";
}
  @fa = ();
  @fa = fastaFileToArray($bigaln);
  for (my $i = 0; $i <= $#fa; $i++)
  {
     $hitid = fasta2id($fa[$i]);
###--- Do QA only on new sequences in big.fa 
     if (exists($check_aln_id{$hitid}))
     {
       if ($debug1)
       {
          print "$hitid already present in $aln_file\n";
       } 
       $idhash{$hitid} = 1;
     }
     elsif ((testCoverage($fa[$i],$hmm_lengths{$hitid})) 
             && ($id_evalue{ $hitid } <= $cutoff_evalue))
     {
       if ( $debug1 )
       {
           print "DEBUG: $hitid passed all tests\n";
       }
       $idhash{$hitid} = 1;
     }
  }

  my @hits_passed = keys(%idhash);
  if($debug1)
  {
    print "        " . ($#hits_passed + 1) . " passed quality tests\n";
  }

  if (@hits_passed == 0)
  {
     print BDKLOGFILE "No new sequences passed cut-offs.\n";
     return 0;

  }

##----Before doing PWid test, write out the sequence from original big.fa without trimming ----
  
  if ($debug1)
  {
    print "Writing big-sel.a2m\n";
  }
  open (FILE, "> big-sel.a2m");
  for (my $i = 0; $i <= $#fa; $i++)
  {
     my $bigid = fasta2id($fa[$i]);
     if (exists($idhash{$bigid}))
     {
         print FILE  "$fa[$i]\n";
     } 
  }
  `cp big-sel.a2m test-sel.a2m`;
  removeDottyColumns ("big-sel.a2m", "big-sel.edit");
  `mv big-sel.edit big-sel.a2m`;
   if (testPwidAll($aln_file, "big-sel.a2m", $cutoff_pwid))
   {
       return 1;
   }
  return 0;  
} 


sub testInsert
{
   my $seed = shift;
   my $hit = shift;
   my $lc_value = shift;

   my $seedseq = fasta2seq($seed);
   my $hitseq = fasta2seq($hit); 
   my $hitid = fasta2id($hit);

   my $seed_lc = countLower($seedseq);
   my $hit_lc = countLower($hitseq);
  
   if ( $debug1 )
   {
       print "DEBUG : Lowercase : $seed_lc  :  $hit_lc\n";
   }
  
   my $seed_len = countUpperLower($seedseq); 
   my $inscutoff = (0.25 * $seed_len);
  
   if($lc_value > $inscutoff)
   {
     $lc_value = $inscutoff;
   }

  if ($hit_lc > $lc_value)
  {
    if ( $debug1 )
    {
      print "DEBUG : $hitid failed check insert test\n";
    }
    return 0;
  }
  if ( $debug1 )
  {
    print "DEBUG: $hitid passed check insert test\n";
  }
  return 1;
}

sub testConsInsert
{
  ###--- Check the number of consequent lowercase letters in hit with respect to seed ---
  ###--- Cutoff is 80 ---

   my $hit = shift;

   my $hitseq = fasta2seq($hit);
   my $hitid = fasta2id($hit);

   my $tmp1 = $hitseq;
   my $tmp3 = $hitseq;
   
   $tmp1 =~ s/[A-Z].*//g;
   my $tmp1_len = countLower($tmp1);
   if ($tmp1_len >= 100)
   {
     if ($debug1)
     {
       print "$hitid : $tmp1_len ; failed Cons insert test\n";
     }
     return 0; 
   } 
   
   $tmp3 =~ s/[a-z][^A-Z]*.//;
   my @array = split(//, $tmp3); 
   my $count = 0;
   for (my $i=0; $i<=$#array ; $i++)
   {
     if ($array[$i] =~ /[a-z]/)
     {
       $count++;
     }
    elsif($array[$i] =~ /[A-Z]/)
    {
      $count = 0;
    }
    if($count >= 100)
    {
      if ($debug1)
      {
        print "$hitid failed testConsinsert\n";
      }     
      return 0;
     }
   }
   
   if ($debug1)
   {
        print "$hitid passed Cons insert test\n";
   }
   return 1;
}

# ------------------------------------------------------------------------------
sub createNewUniverse {
  # Remove the homologs that were already accepted from the universe.
  my $alnfile = shift;
  my $univfile = shift;
  my %alnid    = ();
  
  my @fa = fastaFileToArray($alnfile);
  for (my $i = 0; $i<=$#fa; $i++)
  {
    my $id = fasta2id($fa[$i]);
    $alnid{$id} = 1 ;
  }

  @fa = ();
  @fa = fastaFileToArray($univfile);
  for (my $i = 0; $i<=$#fa; $i++)
  {
      my $id = fasta2id($fa[$i]);
      if (exists ($alnid{$id}))
      {
         splice(@fa, $i, 1);
         $i--;
      }
  }
 
  open (FILE, "> universe.fa");
  for (my $i = 0; $i<=$#fa; $i++)
  {
     print FILE "$fa[$i]\n";
  }
  close FILE;
  if($debug1)
  { 
    print "Universe has " . ($#fa + 1) . " sequences\n";
  }
  my $nseqs_newuniv = ($#fa + 1);
  return $nseqs_newuniv;
}


# ------------------------------------------------------------------------------
sub trimMSA {
  # Removes columns from the alignment that have more than $gappycol fraction of
  # delections, provided that either we are not keeping master-slave alignment
  # to the seed, or, the seed itself has a deletion in that column
  my $input = shift;
  my $output = shift;
  my $gappycol = shift;

  #print "Reading input file into array\n";
  my %h_fasta = ();
  my $id = "";
  my $seq = "";
  my @fa = fastaFileToArray($input);
  foreach my $fa (@fa)
  {
    $id = fasta2id($fa);
    $seq = fasta2seq($fa);
    $h_fasta{$id} = $seq;
  }

  ## Remove gappy columns

  if ($debug1) {
     print "Removing gappy cols\n";
  }
  my $length = length($seq);
  #print "$length\n";
  my $gid = "";
  my $seed_has_gap = 0;
  my $seed_fasta;
  if ($seedmaster) {
    if (exists($h_fasta{$unaltered_seedid}) && 
        defined($h_fasta{$unaltered_seedid})) {
      $seed_fasta = $h_fasta{$unaltered_seedid};
    } elsif (exists($h_fasta{$seedid}) &&
            defined($h_fasta{$seedid})) {
      $seed_fasta = $h_fasta{$seedid};
    } else {
      print "Error in trimMSA: unable to find seed in $input.";
    }
  }
  for (my $col = 0; $col < $length; $col++) {
    my $gaps = 0;
    #if ($debug1) {
    #   print "Processing column $col out of $length\n";
    #} 
    # if keeping master-slave alignment to seed, only remove gappy columns
    # which are also gapped in the aligned seed
    if ($seedmaster) {
      $seed_has_gap = (substr($seed_fasta, $col, 1) eq "-");
    }
    if ($seedmaster == 0 || $seed_has_gap) {
      foreach $gid ( keys %h_fasta ) {
          $gaps++ if ( substr($h_fasta{$gid}, $col, 1) eq "-" );
      }
     
      #if ($debug1) {
      #   printf "%5d %5d %5d %8.2f\n",  $col, $gaps, ($#fa + 1), $gaps/($#fa + 1);
      #} 
      if (($gaps / ($#fa + 1)) > $gappycol || $seedmaster && $seed_has_gap) {
         foreach $gid ( keys %h_fasta )
         {
             substr($h_fasta{$gid}, $col, 1) =~ tr/A-Z/a-z/;
             substr($h_fasta{$gid}, $col, 1) =~ tr/-/./;
         }
      }
    }
 }

 ## Output the results
 if($debug1) {
  print "Writing output file\n";
 }
 open (FILE, "> $output");

 foreach $gid (keys %h_fasta) {
    print FILE ">$gid\n";

    print FILE "$h_fasta{$gid}\n";
 }
 close FILE;
}

# ------------------------------------------------------------------------------
sub testPwidAll
{
  my $refaln = shift;
  my $testaln = shift;
  my $cutoff_pwid = shift;

  my %refseqs = ();  ###--- This is a hash to map accid to sequence from reference alignment

  my @fa = fastaFileToArray($refaln);
  for(my $i=0; $i<=$#fa; $i++)
  {
    my $accid = fasta2id($fa[$i]);
    $refseqs{$accid} = $fa[$i];
  }
 
  @fa = ();
  @fa = fastaFileToArray($testaln);
  for(my $i=0; $i<=$#fa; $i++)
  { 
    my $count = 0;
    my $hitid = fasta2id($fa[$i]);
    if (exists $refseqs{$hitid})
    {
      $count++;
      if ($debug1)
      {
          print "$hitid already exists in $refaln\n";
      }
      next ;
   } 
    my $accid = "" ;
    foreach $accid (sort keys (%refseqs))
    {
       my $accseq = $refseqs{$accid};
       my $checkpwid = pwid($accseq, $fa[$i]);
#       if ((pwid($accseq, $fa[$i])) > $cutoff_pwid)
       if ($checkpwid > $cutoff_pwid)
       {
          $count++;
          if ($debug1)
          {
             print "$hitid passed pwid-all test with $accid\n";
          } 
          last; 
       }
    }
    if ($count == 0)
    {
      if ($debug1)
      {
       print "$hitid failed pwid-all test\n";
      } 
      splice(@fa, $i, 1);
      $i--; 
    }
  }   

###--- Output results
 
  if (@fa == 0)
  {
   if($debug1)
   {
      print "No new sequences passed pwid-all test\n";
   }
    return 0;
  }
  else
  {
      my $sequence_sequences = "sequence";
      if ( $#fa > 0 )
      {
         $sequence_sequences .= "s";
      }

      my $additional = "";
      if ( $count >= 2 ) 
      {
          $additional = " additional";
      }
      print BDKLOGFILE "" . ($#fa+1)
         . "$additional $sequence_sequences passed user-specified criteria.\n";
      open (FILE, "> big-final.a2m");
      for(my $i=0; $i<=$#fa; $i++)
      {
          print FILE "$fa[$i]\n";
      }      
      `cp big-final.a2m test-final.a2m`;
      removeDottyColumns("big-final.a2m", "edit.a2m");
      `mv edit.a2m big-final.a2m`;
  }
  
  return 1; 
}

sub removeDottyColumns
{
   my $input = shift;
   my $output = shift;
   my %h_fasta = ();

   my $def = "";
   my $seq = "";
   my $cnt = 0;
   open (INFILE, "< $input");
   while (my $line = <INFILE>)
   {
     chomp $line;
     if ($line =~ /^>/)
     {
       if ($cnt > 0)
       {
         $h_fasta{$def} = $seq;
         $seq = "";
       }
       $cnt++;
       $def = $line;
       $h_fasta{$def} = "";
     }
     else
     {
       $seq .= $line;
     }
   }
    $h_fasta{$def} = $seq;
   close INFILE;
   #print "SEQ is \n$seq\n";
   my $length = length($seq);
   #print "Length of alignment = $length\n";

   for (my $col = 0; $col < $length; $col++)
   {
      my $dot = 0;
      #if ($debug1)
      #{
      #   print "Processing column $col\n";
      #}
      foreach $def (keys %h_fasta)
      {
        if (substr($h_fasta{$def}, $col, 1) eq ".")
        {
            $dot++;
        }
    }
    #if ($debug1)
    #{
    #  print "$col ; $dot ; ". scalar(keys %h_fasta) . "\n";
    #}
    if ($dot == (scalar keys %h_fasta))
    {
        foreach $def (keys %h_fasta)
        {
           substr($h_fasta{$def}, $col, 1) = "";
        }

        $col--;
        $length--;
    }
  }
  open (FILE, "> $output");
  foreach my $def (keys %h_fasta)
  {
    print FILE "$def\n";
    print FILE "$h_fasta{$def}\n";
  }
}
sub removeInserts
{
  my $fasta = shift;
  my ($def, $seq) = split(/\n/, $fasta);
  $seq =~ s/\.//g;
  $seq =~ s/[a-z]//g;
  $fasta = $def . "\n" . $seq;

  return $fasta;
}

sub replacedefline
{
  my $file = shift;
  my $seedid = shift;
  my $def = shift;
 
  #print "Checking $seedid\n"; 
  open (NEW, "> replace-defline");
  open (FILE, "< $file");
  while (my $line = <FILE>)
  {

     # "[", "]", "(", or ")" in seed ID would give trouble in regular
     # expression, but deleted in fasta2id().

     my $escaped_seedid = $seedid;
     ####$escaped_seedid =~ s/([\[\]\(\)])/\\$1/g;

     if (($line =~ /^>lcl\|$escaped_seedid/))
     {
        $line = $defline;
        print NEW "$line";
     }
     else
     {
        print NEW "$line";
     }
  }
  `cp replace-defline $file`;
}

sub makeunique
{
 ##check for duplicate ids and write out only the representative one.
 my $infile = shift;
 my %uniqidhash = ();
 my @seq = fastaFileToArray($infile);
 for (my $i=0; $i<=$#seq; $i++)
 {
   my $seqid = fasta2id($seq[$i]);
   if (exists $uniqidhash{$seqid})
   {
   }
   else
   {
     $uniqidhash{$seqid} = $seq[$i];
   }
 } 
 open (OFILE, "> uniq-universe.fa");
 foreach my $key (keys (%uniqidhash))
 {
  print OFILE "$uniqidhash{$key}\n";
 } 
 `cp uniq-universe.fa universe.fa`;
}

sub getseq
{ 
  my $id = shift;
  my $infile = shift;
 
  my @seq = fastaFileToArray($infile);
  for (my $i=0; $i<=$#seq; $i++)
  {
    my $seqid = fasta2id($seq[$i]);
    if ($id eq $seqid)
    {
      return $seq[$i]; 
    } 
  } 
  return 0;
}
