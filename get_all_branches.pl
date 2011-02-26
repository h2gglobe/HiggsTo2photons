#!/usr/bin/perl -w
use strict;

sub getlines_infile {
# input parameter file name
# returns list of lines
    my ($file)=$_[0];
    my (@lines);

#    print "FILE $file\n";
    open(INFO,$file);
    @lines=<INFO>;
    close (INFO);
    return @lines;
}

sub getfiles_indir {
# input parameter directory name
# returns list of files
    my ($dir)=$_[0];
    my ($name);
    my (@files);

#    print "DIRECTORY $dir\n";
    opendir(DIR, "$dir");
    while($name = readdir(DIR)) {
	unshift(@files, $name);
    }
    closedir(DIR);
    return @files;
}

sub check_node {
#returns 0 for lxplus, 1 for pcmssd04, -1 unknown
    my ($returnvalue);
    $returnvalue=-1;
    my ($myuname) = qx! uname -a !;
    if ($myuname =~ /pcmssd04/) {
	$returnvalue=1;
    }
    if ($myuname =~ /lxplus/) {
	$returnvalue=0;
    }
    return $returnvalue;
}

my $DIR="src";
my $DIRH="interface";

my @files = (
#"GlobeAnalyzer",
	     "GlobeCommon",
	     "GlobeCaloTowers",
	     "GlobeEcalClusters",
	     "GlobeEcalHits",
	     "GlobeElectrons",
	     "GlobeGenerator",
	     "GlobeGenJets",
	     "GlobeGenJets",
	     "GlobeGenJets",
	     "GlobeHcal",
	     "GlobeHLT",
	     "GlobeHT",
	     "GlobeJets",
	     "GlobeJets",
	     "GlobeJets",
	     "GlobeJets",
	     "GlobeJets",
	     "GlobeJets",
	     "GlobeL1",
	     "GlobeLeptons",
	     "GlobeMET",
	     "GlobeMuons",
	     "GlobePhotons",
	     "GlobeReducedGen",
	     "GlobeSimTracks",
	     "GlobeSimHits",
	     "GlobeTracks",
	     "GlobeTrackingParticles",
	     "GlobeVertex",
	     "GlobeVertex",
	     );

my @names = (
#"_",
	     "_",
	     "_",
	     "_",
	     "_",
	     "_std_",
	     "_",
	     "_algo1_",
	     "_algo2_",
	     "_algo3_",
	     "_",
	     "_",
	     "_",
	     "_algo1_",
	     "_algo2_",
	     "_algo3_",
	     "_algoPF1_",
	     "_algoPF2_",
	     "_algoPF3_",
	     "_",
	     "_",
	     "_",
	     "_glo_",
	     "_",
	     "_",
	     "_",
	     "_",
	     "_",
	     "_",
	     "_",
	     "_std_",
	     "_pix_",
	     );
my $filename;

#if(!$ARGV[0]) {
#    print "need input 0\n";
#    die;
#}

#my @sentlines=getlines_infile("ru.txt");
#my $sentline;
#my @lines = qx! grep "rubuilder" ru.txt  !;
	
my $output="h2gglobe/branchdef/treebranch.h";
system "rm $output";
system "touch $output";

my $output1="h2gglobe/branchdef/branchdef.h";
system "rm $output1";
system "touch $output1";

my $output2="h2gglobe/branchdef/setbranchaddress.h";
system "rm $output2";
system "touch $output2";

my $output3="h2gglobe/branchdef/getbranch.h";
system "rm $output3";
system "touch $output3";

my $output4="h2gglobe/branchdef/getentry.h";
system "rm $output4";
system "touch $output4";

my $output5="h2gglobe/branchdef/treedef.h";
system "rm $output5";
system "touch $output5";

my $output6="h2gglobe/branchdef/newclonesarray.h";
system "rm $output6";
system "touch $output6";

#my @files = getfiles_indir("$DIR");

my $line;

my $i=0;

#header files
foreach $filename (@files) {
    print "filename $DIRH/$filename.h $names[$i] \n";
    my @lines = qx! cat $DIRH/$filename.h  !;
    my $var1="EMPTY";
    my $var2="EMPTY";
    my $var3="EMPTY";
    my $var4="EMPTY";

    foreach $line (@lines) {
	my $process=1;
	if(($line=~ /\/\/.*[a-z].*[a-z].*[a-z]/)) {
	   $process=0;
	}
	if(!($line=~ /[a-z].*[a-z].*[a-z]/)) {
	   $process=0;
	}
	if(($line=~ /\(.*\)/)) {
	   $process=0;
	}
	if(($line=~ /\#include/)) {
	   $process=0;
	}
	if(($line=~ /\#if/)) {
	   $process=0;
	}
	if(($line=~ /\#endif/)) {
	   $process=0;
	}
	if(($line=~ /\#define/)) {
	   $process=0;
	}
	if(($line=~ /protected:/)) {
	   $process=0;
	}
	if(($line=~ /private:/)) {
	   $process=0;
	}
	if(($line=~ /public:/)) {
	   $process=0;
	}
	if(($line=~ /edm::/)) {
	   $process=0;
	}
	if(($line=~ /void /)) {
	    $process=0;
	}
	if(($line=~ /TClonesArray/)) {
	   #$process=0;
	}
	if(($line=~ /debug_level/)) {
	   $process=0;
	}
	if(($line=~ /char\*/)) {
	   $process=0;
	}
	if(($line=~ /GlobeCuts/)) {
	   $process=0;
	}
	if(($line=~ /namespace/)) {
	   $process=0;
	}
	#if(($line=~ /std::/)) {
	if(($line=~ /std::/ && $line!~/std::vector/ && $line!~/std::map/)) {
	   $process=0;
	}
	if(($line=~ /class /)) {
	   $process=0;
	}
	if(($line=~ /typedef /)) {
	   $process=0;
	}
	if(($line=~ /bool do/)) {
	   $process=0;
	}
	if(($line=~ /HepMC::/)) {
	   $process=0;
	}
	if(($line=~ /struct/)) {
	   $process=0;
	}
	if(($line=~ /TrackerHitAssociator\*/)) {
	   $process=0;
	}
	if(($line=~ /enum CorrectionType/)) {
	   $process=0;
	}
	if(($line=~ /[a-z]\*/)) {
	   #$process=0;
	}
	if(($line=~ /float leptons_cut/)) {
	   $process=0;
	}
	if($process) {
      if(($line=~ /Float_t/ || $line=~ /Int_t/ || $line=~ /Bool_t/ || $line=~ /Double_t/ )) {
		#print "$line";
		if(($line=~ /(.+_t )(.+)_(.+)_(.+)_(.+)_(.+)_(.+)/)) {
        system "echo '$1$2$names[$i]$3_$4_$5_$6_$7' >> $output5";
		}
		elsif(($line=~ /(.+_t )(.+)_(.+)_(.+)_(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3_$4_$5_$6' >> $output5";
		}
		elsif(($line=~ /(.+_t )(.+)_(.+)_(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3_$4_$5' >> $output5";
		}
		elsif(($line=~ /(.+_t )(.+)_(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3_$4' >> $output5";
		}
		elsif(($line=~ /(.+_t )(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3' >> $output5";
		}
		elsif(($line=~ /(.+_t )(.+)/)) {
		    system "echo '$1$2' >> $output5";
		}
		else {
		    print "ERROR _t $line\n";
		}
	    }
	    elsif(($line=~ /TClonesArray/)) {
		if(($line=~ /(.*TClonesArray.*\* *)(.+)_(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3_$4' >> $output5";
		}
		elsif(($line=~ /(.*TClonesArray.*\* *)(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3' >> $output5";
		}
		elsif(($line=~ /(.*TClonesArray.*\* *)(.+)/)) {
		    system "echo '$1$2' >> $output5";
		}
		else {
		    print "ERROR TClonesArray $line\n";
		}
	    }
	    elsif(($line=~ /TVector3/)) {
		if(($line=~ /(.*TVector3.*\* *)(.+)_(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3_$4' >> $output5";
		}
		if(($line=~ /(.*TVector3.*\* *)(.+)_(.+)/)) {
		    system "echo '$1$2$names[$i]$3' >> $output5";
		}
		elsif(($line=~ /(.*TVector3.*\* *)(.+)/)) {
		    system "echo '$1$2' >> $output5";
		}
		else {
		    print "ERROR TVector3 $line\n";
		}
	    }
		  elsif(($line=~ /std::vector<std::short>/)) {
          system "echo '$line' >> $output5";
		  }
		  elsif(($line=~ /std::vector<unsigned short>/)) {
          system "echo '$line' >> $output5";
		  }
		  elsif(($line=~ /std::vector<int>/)) {
          system "echo '$line' >> $output5";
		  }
		  elsif(($line=~ /std::map/)) {
          system "echo '$line' >> $output5";
		  }
	    else {
		print "****UNACCOUNTED**** $line";
	    }
	}

    }
    $i++;
}

$i=0;
foreach $filename (@files) {

    #if ($filename =~ /\.cc$/) 
    #{
#remember %d
#remember TcloneArray
	
    print "filename $DIR/$filename.cc $names[$i] \n";
    my @lines = qx! cat $DIR/$filename.cc  !;
	
    my $var1="EMPTY";
    my $var2="EMPTY";
    my $var3="EMPTY";
    my $var4="EMPTY";
	
    foreach $line (@lines) {


#get new TClonesArray
	if(!($line=~ /\/\/.+new.+TClones/)) {

	    if(($line=~ /\s*(.+)\s*=\s*new\s*TClonesArray(.*)/)) {

		#print "$1 $2 --- $line";
		my $varname = $1;
		my $rest = $2;
		
		if($varname=~ /(.+)_(.+)_(.+)_(.+)/) {
		    system "echo '  loops->$1$names[$i]$2_$3_$4 = new TClonesArray$rest' >> $output6";
		}
		elsif($varname=~ /(.+)_(.+)_(.+)/) {
		    system "echo '  loops->$1$names[$i]$2_$3 = new TClonesArray$rest' >> $output6";
		}
		elsif($varname=~ /(.+)_(.+)/) {
		    system "echo '  loops->$1$names[$i]$2 = new TClonesArray$rest' >> $output6";
		}
		else {
		    #print "bbb $varname = new TClonesArray$rest bbb\n";		    }
		    system "echo '  loops->$varname = new TClonesArray$rest' >> $output6";
		}
	    }
	}

	if(!($line=~ /\/\/.+Branch/)) {
	    if($line=~/sprintf.+a1.+\"(.+)_.s_(.+)\"/) {
		#print "A1 $1_$2 $line";
#		if(!($names[$i] =~/^$/)) {
#		    $var1="${1}_$2";
#		}
#		else {
		$var1="$1$names[$i]$2";
#		}
	    }
	    if($line=~/sprintf.+a2/) {
		#if($line=~/sprintf.+a2.+\"(.+)_.s_(.+)\[(.+)_.s_(.+)\].+[.d](.+)\".+\(MAX.+\);/) {
		if($line=~/sprintf.+a2.+\"(.+)_.s_(.+)\[(.+)_.s_(.+)\].+[.d](.+)\".+(MAX.+)\)/) {
		   
        #print "$line\n";
		    #print "A2 $1 $2 $3 $4 $5 $6\n";
		    # sprintf(a2, "vtx_%s_tkind[vtx_%s_n][%d]/I", nome, nome, MAX_VERTEX_TRACKS);
        # A2 vtx tkind vtx n ]/I MAX_VERTEX_TRACKS

		    $var2="${1}$names[$i]${2}[$3$names[$i]$4][$6$5";
		    $var2="${1}$names[$i]${2}[$3$names[$i]$4][";
		    $var4=$6;
		    my $vartemp5=$5;
		    
		    my $maxlin=qx! grep -h $var4 $DIR/../interface/Limits.h  !;
		    $maxlin=~/$var4 (\d+)/;
		    $var4=$1;
		    $var2="$var2$var4$vartemp5";
		    #print "HERE A2 $var4 $var2  $line";
		}
		elsif($line=~/sprintf.+a2.+\"(.+)_.s_(.+)\[(.+)_.s_(.+)\](.+)\"/) {
		    #print "HERE A2 array $1_${2}[$3_$4]$5  $line";
#		    if(!($names[$i] =~/^$/)) {
#			$var2="${1}_${2}[${3}_$4]$5";
#		    }
#		    else {
		    $var2="${1}$names[$i]${2}[$3$names[$i]$4]$5";
#		    }
		}
		elsif ($line=~/sprintf.+a2.+\"(.+)_.s_(.+)\"/) {
#       print "A2 $1_$2 scalar $line";
#		    if(!($names[$i] =~/^$/)) {
#			$var2="${1}_${2}";
#		    }
#		    else {
		    $var2="$1$names[$i]${2}";
#		    }
		}
		elsif ($line=~/sprintf.+a2.+\"(.+)\[.d\](.+)\".+(MAX.+)\)/) {
		    $var2=$1;
		    $var4=$3;
		    my $vartemp5=$2;
		    my $maxlin=qx! grep -h $var4 $DIR/../interface/Limits.h  !;
		    $maxlin=~/$var4 (\d+)/;
		    $var4=$1;
		    $var2="$var2\[$var4\]$vartemp5";
		}
	    }
      
      if(($line =~/std::map/ || $line =~/std::vector/) && $line =~/tree/) {
		    if ($line=~/tree.+a1.+"(.+).+,.+&(.+)_(.+).+;/){
          #print "$1 $2 $3\n";
          $var2 = $1;
          $var3 = "$2$names[$i]$3";
        }
	      elsif ($line=~/tree.+Branch+."(.+)_(.+)".+"(.+).+,(.+)/){
          #print "$1 $2 $3\n";
          $var3 = "$1$names[$i]$2";
          $var2 = $3;
        }
		    system "echo '   outputTree-\>Branch(\"$var3\", \"$var2\", &loops->$var3);' >> $output";
		    system "echo '   TBranch \*b_$var3;' >> $output1";		    
		    system "echo '   fChain->SetBranchAddress(\"$var3\", \&$var3, \&b_$var3); ' >> $output2";		    
		    system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
		    system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
        #print "var1 = $var1\n";
        #print "var2 = $var2\n";
        #print "var3 = $var3\n";
        #print "var4 = $var4\n\n";
      }
	    elsif($line=~/tree.+Branch.+TClonesArray.+,(.+),.+,.+\)/) {
		#print "HERE1 $1 $line";

		$var3="$1"; #was this
		#$var3="$1";

		$var3 =~ s/\_/$names[$i]/;
		if($var1=~/EMPTY/) {
		    #print "HERE2 $line";
		    $line =~ s/tree/outputTree/;
		    $line =~ s/\s//g;
		    
		    $line =~ s/Array\",&/AMPER/;
		    $line =~ s/Array\",/COMMA/;
		    $line =~ s/COMMA/Array\",loops->/;
		    $line =~ s/AMPER/Array\",&loops->/;
		    		    
#		    $line =~ s/,&/AMPER/;
#		    $line =~ s/,/COMMA/;
#		    $line =~ s/COMMA/,loops->/;
#		    $line =~ s/AMPER/,&loops->/;
			
		    chomp($line);
		    #system "echo     '      $line' >> $output";
		    $var3=~s/\&//;
		    $var3=~s/\s//g;
		    system "echo '      outputTree-\>Branch(\"$var3\", \"TClonesArray\",&loops->$var3, 32000, 0);' >> $output";
		    system "echo '    TBranch \*b_$var3;' >> $output1";		    
		    system "echo '   fChain->SetBranchAddress(\"$var3\", \&$var3, \&b_$var3); ' >> $output2";		    
		    system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
		    system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		}
		
		
		else {
		    #print "HERE3 $line";
		    my $temp=$var3;
		    $temp =~ s/\s//g;
		    if($temp =~/\&/) {
			$temp =~ s/\&/\&loops->/g;
		    }
		    else {
			$temp= "loops->$temp";
		    }
		    
		    system "echo     'outputTree-\>Branch(\"$var1\", \"TClonesArray\",$temp, 32000, 0);' >> $output";
		    $var3=~s/\&//;
		    $var3=~s/\s//g;
		    system "echo '    TBranch \*b_$var3;' >> $output1";		    
		    system "echo '   fChain->SetBranchAddress(\"$var3\", \&$var3, \&b_$var3); ' >> $output2";
		    system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
		    system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		}
	    }
	    elsif ($line=~/tree.+Branch.+,.+"TVector3",(.+),(.+),.+\)\;/) {
		$var3="$1";
		$var3 =~ s/\_/$names[$i]/;
		
		if($var1=~/EMPTY/) {
		    $line =~ s/tree/outputTree/;
		    $line =~ s/\s//g;
		    $line =~ s/\"\,\&/AMPER/;
		    $line =~ s/\"\,/COMMA/;
		    $line =~ s/COMMA/\"\,/;
		    $line =~ s/AMPER/\"\,\&loops->/;
		    
		    if ($line =~/a2\)/) {
			$line =~ s/a2\)/\"$var2\"\)/g;
		    }
		    chomp($line);
		    system "echo '      $line' >> $output";
		    $var3=~s/\&//;
		    $var3=~s/\s//g;
		    system "echo '    TBranch \*b_$var3;' >> $output1";
		    system "echo '   fChain->SetBranchAddress(\"$var3\", $var3, \&b_$var3); ' >> $output2";
		    system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
		    system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		}
	    }
	    elsif ($line=~/tree.+Branch.+,(.+),(.+)\)\;/) {
		#print "poss1 $line";
		$var3="$1";
		my $lastvar=$2;
		$var3 =~ s/\_/$names[$i]/;
		
		if($var1=~/EMPTY/) {
		    #print "poss1EMPTY1 $line $lastvar\n";
		    $line =~ s/tree/outputTree/;
		    $line =~ s/\s//g;
		    $line =~ s/\"\,\&/AMPER/;
		    $line =~ s/\"\,/COMMA/;
		    $line =~ s/COMMA/\"\,loops->/;
		    $line =~ s/AMPER/\"\,\&loops->/;
		    
		    if ($line =~/a2\)/) {
			$line =~ s/a2\)/\"$var2\"\)/g;
		    }

		    if ($lastvar =~/a2/) {
			$lastvar =~ s/a2/\"$var2\"/g;
		    }
		    #print "poss1EMPTY2 $line $lastvar\n";

		    chomp($line);
		    #system "echo '      $line' >> $output";
		    $var3=~s/\&//;
		    $var3=~s/\s//g;
		    system "echo '    TBranch \*b_$var3;' >> $output1";
		    system "echo '      outputTree->Branch(\"$var3\",\&loops->$var3,$lastvar); ' >> $output";
		    if($line =~ /\[.+\]/) {


			system "echo '   fChain->SetBranchAddress(\"$var3\", $var3, \&b_$var3); ' >> $output2";
			system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
			system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		    }
		    else {
			system "echo '   fChain->SetBranchAddress(\"$var3\", \&$var3, \&b_$var3); ' >> $output2";
			system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
			system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		    }
		}
		else {
		    #print "poss1FULL  $line\n";
		    #print "FINAL $var1 $var2 $var3 \n";
		    my $temp=$var3;
		    $temp =~ s/\s//g;
		    if($temp =~/\&/) {
			$temp =~ s/\&/\&loops->/g;
		    }
		    else {
			$temp= "loops->$temp";
		    }
		    system "echo '    outputTree-\>Branch(\"$var1\",$temp, \"$var2\");' >> $output";
		    $var3=~s/\&//;
		    $var3=~s/\s//g;
		    system "echo '    TBranch \*b_$var3;' >> $output1";
		    
		    if($var2 =~ /\[.+\]/) {
			system "echo '   fChain->SetBranchAddress(\"$var3\", $var3, \&b_$var3); ' >> $output2";
			system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
			system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		    }
		    else {
			system "echo '   fChain->SetBranchAddress(\"$var3\", \&$var3, \&b_$var3); ' >> $output2";
			system "echo '   b_$var3 = fChain->GetBranch(\"$var3\"); ' >> $output3";
			system "echo '   b_$var3->GetEntry(jentry); ' >> $output4";
		    }
		    
		}
	

		$var1="EMPTY";
		$var2="EMPTY";
		$var3="EMPTY";
		$var4="EMPTY";
		
	    }
	}
    }

    #system "echo ' ' >> $output";
    #system "echo ' ' >> $output1";
    #system "echo ' ' >> $output2";
    #system "echo ' ' >> $output3";
    #system "echo ' ' >> $output4";

    $i++;
}

print "Done\n";
