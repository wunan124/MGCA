#!/usr/bin/perl

##############################################################################
#
#                                 Matched Gene Analysis
#
# Author:  Xiaobing Yuan (xbyuan@brain.ecnu.edu.cn)
# Written: 25th Sept 2018
# 
# DESCRIPTION:
# This script is designed to get 200 random genesets of matched gene length # and/or expression level relative to the geneset of interest and calculate # their internal CC and CC to the input geneset. The input files of this   # script are allGene list file, the input gene list file, and the          # correlation matrix dataset. AllGene list needs to be ranked based on     # either gene size or gene expression levels. Matrix dataset are           # correlation coefficiencies of ranked all genes. There should be a        # /Matrix/ folder under the working directory.
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
my @row = ();
my @WGM = ();

#get gene names of the WGM gene data and total gene number
my $TotalGeneNumb = 0;
my %Position;
my $i=0;
my $k=0;
my $A=0;

open(F0IN, $AllGene) or die("Could not open $AllGene\n");
while(my $line = <F0IN>)
{
	chomp($line);
	$WGM[$i] = $line;
	$Position{$line} = $i;
	$i+=1;
}	
close(F0IN);
$TotalGeneNumb = $i;
#make the reverse harsh of %Position in order to get the key (gene name) from the value (rank) quickly in later steps
my %rhash = reverse %Position;

#get gene names of the input geneset and its position in the all gene list
$i=0;
my @Gene = ();
my @GenePosition = ();
open(F1IN, $InputGene) or die("Could not open $InputGene\n");
while(my $line = <F1IN>)
{
	chomp($line);
	if ($line ~~ @WGM)
	{
		$Gene[$i] = $line;
		$GenePosition[$i] = $Position{$line};
		$i+=1;
	}
}	
close(F1IN);

my @RandGene = ();
my $range = 50;

# set genes for the 100000 random genesets of matched size or expression (50 upper or lower in the rank for each gene in the input gene list), calculate the CC matrix
for (my $j=1; $j<100001; $j++)	      
{	
	for (my $a=0; $a<$i; $a++)
	{
#Random selection of matched gene from witin -50-50 range of each input gene. Avoid potential overlap with input genes.
		 my $MidPosition = $GenePosition[$a];
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
		# if ($RandPos ~~ @GenePosition)
		# {
		#	goto Start;
		# }
		 my $Selected_gene = $rhash{$RandPos};
#Avoid duplicated selection of random gene
		 if ($Selected_gene ~~ @RandGene)
		 {
			goto Start;
		 }
		 $RandGene[$a] =  $Selected_gene;	
	}

#Write selected random genes into the MatchedRand-$j file 
	my $tmp ="$AllGene-MatchedRand100000-$j";	
	open(F1OUT, "> $tmp") or die("Could not write to file: $tmp\n");	
	for ($a=0; $a<$i; $a++)
	{
		print F1OUT $RandGene [$a] . "\n";
	}
	close (F1OUT);

#Calculate the internal CC of the MatchedRand-n geneset, its CC to the input geneset, and the CC of whole genome genes to MatchedRand-n geneset.
	my $status = `perl Matrix.pl $AllGene $tmp $tmp $dataset`;
	   $status = `perl Matrix.pl $AllGene $InputGene $tmp $dataset`;


#set the average CC of random dataset-derived matrix as the value of the summary matrix
	my $input2 = "$tmp-$tmp-$dataset-Matrix";
 	my $input3 = "$InputGene-$tmp-$dataset-Matrix";

	open(F2IN, '<', $input2) or die("Could not open: $input2\n");
	$k = 0;
	while(<F2IN>)
	{
		chomp($_);		
		if ($k > 0)
		{			
			@row = split /\s+/, $_;
			$A = scalar(@row);
			$Matrix2[$k-1][$j-1] = $row [$A-1];
		}
		$k+=1;
	}
	close(F2IN);

	open(F3IN, '<', $input3) or die("Could not open: $input3\n");
	$k = 0;
	while(<F3IN>)
	{
		chomp($_);		
		if ($k > 0)
		{			
			@row = split /\s+/, $_;
			$A = scalar(@row);
			$Matrix3[$k-1][$j-1] = $row [$A-1];
		}
		$k+=1;
	}
	close(F3IN);

	unlink $tmp;
}

#write the matrix value into the output file
my $Output2 = "$AllGene-MatchedRand100000-MatchedRand100000-Sum";
my $Output3 = "$AllGene-$InputGene-MatchedRand100000-Sum";

open(F2OUT, "> $Output2") or die("Could not write to file: $Output2\n");
open(F3OUT, "> $Output3") or die("Could not write to file: $Output3\n");	
for ($a=1; $a<$k; $a++)
{
	my $tmp = "$AllGene-MatchedRand100000-$a";
	print F2OUT $tmp . "\t";
	print F3OUT $Gene[$a-1] . "\t";
	for (my $b=1; $b<100001; $b++)
	{
		print F2OUT $Matrix2[$a-1][$b-1] . "\t";
		print F3OUT $Matrix3[$a-1][$b-1] . "\t";

	}
	print F2OUT "\n";
	print F3OUT "\n";

}
close (F2OUT);
close (F3OUT);

#delete the matched random genesets and matrix derived from each random geneset
for ($a=1; $a<100001; $a++)   
{
	my $tmp ="$AllGene-MatchedRand100000-$a";
	my $input1 = "$tmp-$tmp-$dataset-Matrix"; 
	unlink $input1;
	my $input2 = "$InputGene-$tmp-$dataset-Matrix"; 
	unlink $input2;
}

exit(1); 





