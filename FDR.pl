#!/usr/bin/perl

##############################################################################
#
#                                 FDR
#
# Author:  Xiaobing Yuan (xbyuan@brain.ecnu.edu.cn)
# Written: 25th Sept 2018
# 
# DESCRIPTION:
# This script is designed to calculate FDR of whole genes by MGCA. The input files of this script are AllGene list file, the input gene list file, and the correlation matrix dataset. AllGene list needs to be ranked based on on gene expression levels, gene size or GC content. Matrix dataset are correlation coefficiencies of ranked all genes.
#
# Copyright 2018 Xiaobing Yuan
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

if($fail)
{
	exit(1);
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

my @Matrix2;
my @Matrix3;
my @Matrix4;
my @InputGene_Ranked;
my @Row = ();
my @WGM = ();

#get gene names of the WGM gene data and total gene number
my $TotalGeneNumb = 0;
my %Position;
my $Matrix_Row=0;
my $Matrix_Column=0;

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
my $InputGeneNumb=0;
my @Gene = ();
my @GenePosition = ();
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

my @mRandGene = ();
my $range = 50;

# set genes for the 5000 random genesets of matched expression, size or GC content (50 upper or lower in the rank for each gene in the input gene list), calculate the CC matrix
for (my $Cycle_mRand=1; $Cycle_mRand<5001; $Cycle_mRand++)	      
{	
	for (my $Cycle_Input=0; $Cycle_Input<$InputGeneNumb; $Cycle_Input++)
	{
#Random selection of matched gene from witin ±50 range of each input gene. Avoid potential overlap with input genes.
		 my $MidPosition = $GenePosition[$Cycle_Input];
		 Start:
		 my $mRand_number = (int(rand($range)) + 1) * (int(rand(2))-0.5)/0.5;
		 my $mRandPos = $MidPosition + $mRand_number;
		 if ($mRandPos < 0)
		 {
			goto Start;
		 }
		 if ($mRandPos > ($TotalGeneNumb-1))
		 {
			goto Start;
		 }
		 my $Selected_gene = $rhash{$mRandPos};
#Avoid duplicated selection of matched random gene
		 if ($Selected_gene ~~ @mRandGene)
		 {
			goto Start;
		 }
		 $mRandGene[$Cycle_Input] =  $Selected_gene;	
	}

#Write selected random genes into the mRand-$Cycle file 
	my $tmp ="$AllGene-mRandomGene5000-$Cycle_mRand";	
	open(F1OUT, "> $tmp") or die("Could not write to file: $tmp\n");	
	for (my $Cycle_Input=0; $Cycle_Input<$InputGeneNumb; $Cycle_Input++)
	{
		print F1OUT $mRandGene [$Cycle_Input] . "\n";
	}
	close (F1OUT);

    #Calculate the P-value of whole genome genes based on mRand-$Cycle geneset.
	my $status = `perl MGCA_Range50.pl $AllGene $tmp $dataset`;

    #set the P-value as the value of the summary matrix
    my $input2 = "$AllGene-$tmp-MGCA_Range50-Results";
    open(F2IN, '<', $input2) or die("Could not open: $input2\n");
    $Matrix_Row = 0;
    while(my $line = <F2IN>)
    {
	    chomp($line);		
        @Row = split /\s+/, $line;
	    $Matrix_Column = scalar(@Row);
	    $Matrix2[$Matrix_Row-1][$Cycle_mRand-1] = $Row [$Matrix_Column-1];
	    $Matrix_Row+=1;
	}
	close(F2IN);

	unlink $tmp;
}

my @FDR = ();
my $Count = 0;
#Calculate FDR and write into the output file
my $Output2 = "$AllGene-FDR-Sum";

open(F2OUT, "> $Output2") or die("Could not write to file: $Output2\n");	
for (my $Input=1; $Input<$TotalGeneNumb+1; $Input++)
{
	print F2OUT $WGM[$Input-1] . "\t";
	for (my $Column_Numb=1; $Column_Numb<5001; $Column_Numb++)
	{
        if($Matrix2[$Input-1][$Column_Numb-1] < 0.001)
        {
            $Count+=1;
        }	
	}
    $FDR[$Input] = $Count/5000;
    print F2OUT $FDR[$Input] . "\n";
}
close (F2OUT);

#delete the matched random genesets and matrix derived from each matched random geneset
for (my $Cycle_mRand=1; $Cycle_mRand<5001; $Cycle_mRand++)   
{
	my $tmp ="$AllGene-mRandomGene5000-$Cycle_mRand";
	my $input1 = "$AllGene-$tmp-MGCA_Range50-Results"; 
	unlink $input1;
}

exit(1); 