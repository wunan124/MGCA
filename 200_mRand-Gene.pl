#!/usr/bin/perl

##############################################################################
#
#                                 200_mRand-Gene
#
# Author:  Xiaobing Yuan (xbyuan@brain.ecnu.edu.cn)
# Written: 25th Sept 2018
# 
# DESCRIPTION:
# This script is designed to get 200 random genesets of matched gene expression level, gene length or GC content relative to the geneset of interest and calculate their internal CC, CC to the input geneset and CC to whole genes. The input files of this script are AllGene list file, the input gene list file, and the correlation matrix dataset. AllGene list needs to be ranked based on on gene expression levels, gene size or GC content. Matrix dataset are correlation coefficiencies of ranked all genes.
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

# set genes for the 200 random genesets of matched expression, size or GC content (50 upper or lower in the rank for each gene in the input gene list), calculate the CC matrix
for (my $Cycle_mRand=1; $Cycle_mRand<6; $Cycle_mRand++)	      
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
		 if ($mRandPos ~~ @GenePosition)
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
	my $tmp ="$AllGene-mRandomGene200-$Cycle_mRand";	
	open(F1OUT, "> $tmp") or die("Could not write to file: $tmp\n");	
	for (my $Cycle_Input=0; $Cycle_Input<$InputGeneNumb; $Cycle_Input++)
	{
		print F1OUT $mRandGene [$Cycle_Input] . "\n";
	}
	close (F1OUT);

#Calculate the internal CC of the mRand-$Cycle geneset, its CC to the input geneset, and the CC of whole genome genes to mRand-$Cycle geneset.
	my $status = `perl Matrix.pl $AllGene $tmp $tmp $dataset`;
	   $status = `perl Matrix.pl $AllGene $InputGene $tmp $dataset`;
	   $status = `perl Matrix.pl $AllGene $AllGene $tmp $dataset`;


#set the average CC of matched random dataset-derived matrix as the value of the summary matrix
	my $input2 = "$tmp-$tmp-$dataset-Matrix";
 	my $input3 = "$InputGene-$tmp-$dataset-Matrix";
	my $input4 = "$AllGene-$tmp-$dataset-Matrix";

	open(F2IN, '<', $input2) or die("Could not open: $input2\n");
	$Matrix_Row = 0;
	while(<F2IN>)
	{
		chomp($_);		
		if ($Matrix_Row > 0)
		{			
			@Row = split /\s+/, $_;
			$Matrix_Column = scalar(@Row);
			$Matrix2[$Matrix_Row-1][$Cycle_mRand-1] = $Row [$Matrix_Column-1];
		}
		$Matrix_Row+=1;
	}
	close(F2IN);

	open(F3IN, '<', $input3) or die("Could not open: $input3\n");
	$Matrix_Row = 0;
	while(<F3IN>)
	{
		chomp($_);		
		if ($Matrix_Row > 0)
		{			
			@Row = split /\s+/, $_;
			$Matrix_Column = scalar(@Row);
			$Matrix3[$Matrix_Row-1][$Cycle_mRand-1] = $Row [$Matrix_Column-1];
			$InputGene_Ranked[$Matrix_Row-1] = $Row [0];
		}
		$Matrix_Row+=1;
	}
	close(F3IN);

	open(F4IN, '<', $input4) or die("Could not open: $input4\n");
	my $Matrix_Row1 = 0;
	while(<F4IN>)
	{
		chomp($_);		
		if ($Matrix_Row1 > 0)
		{			
			@Row = split /\s+/, $_;
			$Matrix_Column = scalar(@Row);
			$Matrix4[$Matrix_Row1-1][$Cycle_mRand-1] = $Row [$Matrix_Column-1];
		}
		$Matrix_Row1+=1;
	}
	close(F4IN);

	unlink $tmp;
}

#write the matrix value into the output file
my $Output2 = "$AllGene-MatchedRand200-MatchedRand200-Sum";
my $Output3 = "$AllGene-$InputGene-MatchedRand200-Sum";
my $Output4 = "$AllGene-MatchedRand200-Sum";

open(F2OUT, "> $Output2") or die("Could not write to file: $Output2\n");
open(F3OUT, "> $Output3") or die("Could not write to file: $Output3\n");	
for (my $Input=1; $Input<$Matrix_Row; $Input++)
{
	my $tmp = "$AllGene-MatchedRandomGene200-$Input";
	print F2OUT $tmp . "\t";
	print F3OUT $InputGene_Ranked[$Input-1] . "\t";
	for (my $Column_Numb=1; $Column_Numb<6; $Column_Numb++)
	{
		print F2OUT $Matrix2[$Input-1][$Column_Numb-1] . "\t";
		print F3OUT $Matrix3[$Input-1][$Column_Numb-1] . "\t";

	}
	print F2OUT "\n";
	print F3OUT "\n";

}
close (F2OUT);
close (F3OUT);

open(F4OUT, "> $Output4") or die("Could not write to file: $Output4\n");	
for (my $Input=1; $Input<$TotalGeneNumb+1; $Input++)
{
	print F4OUT $WGM[$Input-1] . "\t";
	for (my $Column_Numb=1; $Column_Numb<6; $Column_Numb++)
	{
		print F4OUT $Matrix4[$Input-1][$Column_Numb-1] . "\t";
	}
	print F4OUT "\n";
}
close (F4OUT);

#delete the matched random genesets and matrix derived from each matched random geneset
for (my $Cycle_Rand=1; $Cycle_Rand<6; $Cycle_Rand++)   
{
	my $tmp ="$AllGene-mRandomGene200-$Cycle_Rand";
	my $input1 = "$tmp-$tmp-$dataset-Matrix"; 
	unlink $input1;
	my $input2 = "$InputGene-$tmp-$dataset-Matrix"; 
	unlink $input2;
	my $input3 = "$AllGene-$tmp-$dataset-Matrix"; 
	unlink $input3;
}

exit(1); 