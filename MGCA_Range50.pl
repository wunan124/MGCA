#!/usr/bin/perl

##############################################################################
#
#                                 MGCA_Range50
#
# Author:  Xiaobing Yuan (xbyuan@brain.ecnu.edu.cn)
# Written: 25th Sept 2019
# 
# DESCRIPTION:
# This script is designed to calculate P-value of gene by MGCA under gene expression level, gene length or GC content ranked. The number of groups of mRand gene is 100000, and the matched range of mRand genes is 50. The input files of this script are AllGene list file, the input gene list file, and the correlation matrix dataset. AllGene list needs to be ranked based on on gene expression levels, gene size or GC content. Matrix dataset are correlation coefficiencies of ranked all genes.
#
# Copyright 2021 Xiaobing Yuan
#                Institute of Brain Functional Genomics
#                ECNU, Shanghai
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see:
# <http://www.gnu.org/licenses/>.
###########################################################################

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

# Check to see if we received 3 inputs in the arguments
if(scalar(@ARGV) != 3)
{
	print STDERR "Incorrect number of arguments\n";
	exit(1);
}

my $fail = $FALSE;

# Check to see if we received 1 file in the first argument
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

# If the file didn't exist, let's kill it
if($fail)
{
	exit(1);
}

# Read the file names in from the command line. The first file is the total gene list. The second file is the gene of interest. The third file is the matrix of genome-wide gene pairwise correlations.
my $AllGene = $ARGV[0];
my $InputGene = $ARGV[1];
my $dataset = $ARGV[2];

#get gene names of the WGM gene data and total gene number
my %Position;
my @WGM = ();
my $TotalGeneNumb = 0;

open(F0IN, $AllGene) or die("Could not open $AllGene\n");
while(my $line = <F0IN>)
{
	chomp($line);
	$WGM[$TotalGeneNumb] = $line;
	$Position{$line} = $TotalGeneNumb;
	$TotalGeneNumb+=1;
}	
close(F0IN);
#make the reverse harsh of %Position in order to get the key (gene name) from the value (rank) quickly in later steps
my %rhash = reverse %Position;

#get gene names and gene number of the input geneset
my @Gene = ();
my @GenePosition = ();
my $InputGeneNumb = 0;
open(F1IN, $InputGene) or die("Could not open $InputGene\n");
while(my $line = <F1IN>)
{
	chomp($line);
	if ($line ~~ @WGM)
	{
		$Gene[$InputGeneNumb] = $line;
		$GenePosition[$InputGeneNumb] = $Position{$line};		
		$InputGeneNumb+=1;
	}
}	
close(F1IN);

my $Data_Line = 0;
my @row = ();
my @Ave_Input = ();
my @Numb = ();
#Read the dataset and calculate the P-value
open(F2IN, $dataset) or die("Could not open $dataset\n");
while(my $line = <F2IN>)
{
	chomp($line);
	@row = split /\s+/, $line;
    my $Sum_Input = 0;
	for (my $a=0; $a<$InputGeneNumb; $a++)
	{
		my $b = $GenePosition[$a];
		$Sum_Input = $Sum_Input + $row[$b];
	}
	my $tmp_Input = $Sum_Input/$InputGeneNumb;
	$Ave_Input[$Data_Line] = eval sprintf('%.4f', $tmp_Input);
		
	# set genes for the 100000 random genesets of matched expression, size or GC content (50 upper or lower in the rank for each gene in the input gene list)
    my @RandGene = ();
    my @mRand = ();
    my $range = 50;
    my $Count = 0;
    my $Ave_mRand = 0;
    for (my $Cycle_mRand=0; $Cycle_mRand<100000; $Cycle_mRand++)
    {
        my $Sum_mRand = 0;
	    print $AllGene . "\t" .  "Calculate_No:" . $Data_Line .  "\t" . "Select_mRand:" . $Cycle_mRand . "\n";
	    for (my $Cycle_Input=0; $Cycle_Input<$InputGeneNumb; $Cycle_Input++)
	    {
			#Random selection of matched gene from witin ±50 range of each input gene. Avoid potential overlap with input genes.
	    	my $MidPosition = $GenePosition[$Cycle_Input];
		    Start:
		    my $random_number = (int(rand($range)) + 1) * (int(rand(2))-0.5)/0.5;
		    my $RandPos = $MidPosition + $random_number;
		    if ($RandPos < 0)
		     {
		    	goto Start;
		     }
		    if ($RandPos > ($TotalGeneNumb-1))
		     {
		    	goto Start;
		     }
            #if ($RandPos ~~ @GenePosition)
		 	#{
				#goto Start;
			#}
		    my $Selected_gene = $rhash{$RandPos};
		    if ($Selected_gene ~~ @RandGene)
		     {
			    goto Start;
		     }
		     $RandGene[$Cycle_Input] = $Selected_gene;
            $Sum_mRand = $Sum_mRand + $row[$RandPos];
	    }
        my $tmp_mRand = $Sum_mRand/$InputGeneNumb;
		$Ave_mRand = eval sprintf('%.4f',$tmp_mRand);
        if ($Ave_mRand > $Ave_Input[$Data_Line])
		{
			$Count+=1;
		}		
	}
    $Numb[$Data_Line] = $Count/100000;
    $Data_Line+=1;
 }

my $Output1 = "$AllGene-$InputGene-MGCA_Range50-Results";
open(F1OUT, "> $Output1") or die("Could not write to file: $Output1\n");
for (my $Cycle_TOTAL=0; $Cycle_TOTAL<$Data_Line; $Cycle_TOTAL++)
{
	print F1OUT $WGM[$Cycle_TOTAL-1] . "\t";
	print F1OUT $Numb[$Cycle_TOTAL] . "\n";
}
close (F1OUT);

exit(1); 







