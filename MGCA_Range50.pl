use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Copy;

## CONSTANTS ##
my $TRUE              = 1;
my $FALSE             = 0;
my $DEBUG             = $FALSE;
my $EXITSTATUS        = 0;

# Default umask
#umask 027;

# Define variables


##############################
#   M A I N   P R O G R A M  #
##############################

if(scalar(@ARGV) != 3)
{
	print STDERR "Incorrect number of arguments\n";
	exit(1);
}

my $fail = $FALSE;

my $file1 = $ARGV[0];
if(!-e $file1)
{
	print STDERR "File $file1 didn't exist\n";
	$fail = $TRUE;
}

my $file2 = $ARGV[1];
if(!-e $file2)
{
	print STDERR "File $file2 didn't exist\n";
	$fail = $TRUE;
}

my $file3 = $ARGV[2];
if(!-e $file3)
{
	print STDERR "File $file3 didn't exist\n";
	$fail = $TRUE;
}

if($fail)
{
	exit(1);
}

my $AllGene = $ARGV[0];
my $InputGene = $ARGV[1];
my $dataset = $ARGV[2];

#读取总的基因名称和基因数。
my %Position;
my @WGM = ();
my $TotalGeneNumb = 0;
open(F0IN, $AllGene) or die("Could not open $AllGene\n");#读取总基因列表。
while(my $line = <F0IN>)#while循环，逐行读取。
{
	chomp($line);#去除换行，方便直接读取。
	$WGM[$TotalGeneNumb] = $line;
	$Position{$line} = $TotalGeneNumb;
	$TotalGeneNumb+=1;
}	
close(F0IN);
my %rhash = reverse %Position;

#读取导入基因名称和基因数。
my @Gene = ();
my @GenePosition = ();
my $ASDGeneNumb = 0;
open(F1IN, $InputGene) or die("Could not open $InputGene\n");#读取导入基因列表。
while(my $line = <F1IN>)#while循环，逐行读取。
{
	chomp($line);#去除换行，方便直接读取。
	if ($line ~~ @WGM)#if语句，如果导入基因存在于总基因列表中。
	{
		$Gene[$ASDGeneNumb] = $line;
		$GenePosition[$ASDGeneNumb] = $Position{$line};		
		$ASDGeneNumb+=1;
	}
}	
close(F1IN);

my $Data_Line = 0;
my @row = ();
my @Ave_ASD = ();
my @Numb = ();
open(F2IN, $dataset) or die("Could not open $dataset\n");#读取数据矩阵。
while(my $line = <F2IN>)#while循环，逐行读取。
{
	chomp($line);#去除换行，方便直接读取。
	@row = split /\s+/, $line;#以空格拆分成列。
	#$a = scalar(@row);#统计列数。
	
    my $Sum_ASD = 0;
	for (my $a=0; $a<$ASDGeneNumb; $a++)
	{
		my $b = $GenePosition[$a];
		$Sum_ASD = $Sum_ASD + $row[$b];
	}
	my $tmp_ASD = $Sum_ASD/$ASDGeneNumb;
	$Ave_ASD[$Data_Line] = eval sprintf('%.4f', $tmp_ASD);
		
    my @RandGene = ();
    my @mRand = ();
    my $range = 50;
    my $Count = 0;
    my $Ave_mRand = 0;
    for (my $Cycle_mRand=0; $Cycle_mRand<100000; $Cycle_mRand++)#for循环，从0循环到999999,循环1000000次.
    {
        my $Sum_mRand = 0;
	    print $AllGene . "\t" .  "Calculate_No:" . $Data_Line .  "\t" . "Select_mRand:" . $Cycle_mRand . "\n";
	    for (my $Cycle_ASD=0; $Cycle_ASD<$ASDGeneNumb; $Cycle_ASD++)#for循环，从0循环到导入基因数减一。
	    {
	    	my $MidPosition = $GenePosition[$Cycle_ASD];
		    Start:
		    my $random_number = (int(rand($range)) + 1) * (int(rand(2))-0.5)/0.5;#在正负50之内选取随机整数。
		    my $RandPos = $MidPosition + $random_number;#匹配随机基因的位置。
		    if ($RandPos < 0)#如果选到0之前的位置，则重选。
		     {
		    	goto Start;
		     }
		    if ($RandPos > ($TotalGeneNumb-1))#如果选到基因总数减一之后的位置，则重选。
		     {
		    	goto Start;
		     }
             #if ($RandPos ~~ @GenePosition)
		 	#{
				#goto Start;
			 #}
		    my $Selected_gene = $rhash{$RandPos};#根据位置信息，读取匹配随机基因的基因名。
		    if ($Selected_gene ~~ @RandGene)#如果选取的基因已经在匹配随机基因列表中，则重选。
		     {
			    goto Start;
		     }
		     $RandGene[$Cycle_ASD] = $Selected_gene;
            $Sum_mRand = $Sum_mRand + $row[$RandPos];
	    }
        my $tmp_mRand = $Sum_mRand/$ASDGeneNumb;
		$Ave_mRand = eval sprintf('%.4f',$tmp_mRand);
        if ($Ave_mRand > $Ave_ASD[$Data_Line])
		{
			$Count+=1;
		}		
	}
    $Numb[$Data_Line] = $Count;
    $Data_Line+=1;
 }

my $Output1 = "$AllGene-$InputGene-MatchedRange50-0-1M-Results";
open(F1OUT, "> $Output1") or die("Could not write to file: $Output1\n");
for (my $Cycle_TOTAL=0; $Cycle_TOTAL<$Data_Line; $Cycle_TOTAL++)
{
	print F1OUT $Numb[$Cycle_TOTAL] . "\n";
#	print "Result_writing:" . $Cycle_TOTAL . "\n";
}
close (F1OUT);

exit(1); 







