#!/usr/bin/perl

##############################################################################
#
#                                 Matrix
#
# Author:  Xiaobing Yuan (xyuan@@hussmanautism.org)
# Written: 10th Aug 2015
# 
# DESCRIPTION:
# This script is designed to get the matrix of correlation coefficient
# between two interested gene sets. 
# 
# Copyright 2017 Xiaobing Yuan
#                Hussman Institute for Autism
#                Baltimore, USA
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
##############################################################################

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

# Check to see if we received 4 files in the arguments
if(scalar(@ARGV) != 4)
{
	print STDERR "Incorrect number of arguments\n";
	exit(1);
}

my $fail = $FALSE;
foreach my $file (@ARGV)
{
	if(!-e $file)
	{
		print STDERR "File $file didn't exist\n";
		$fail = $TRUE;
	}
}


# If the file didn't exist, let's kill it
if($fail)
{
	exit(1);
}

# Read the file names in from the command line

my $GeneList = $ARGV[0];

my $Geneset1 = $ARGV[1];

my $Geneset2 = $ARGV[2];

my $dataset = $ARGV[3];


my %oldpair = ();

open(F1IN, $GeneList) or die("Could not open $GeneList\n");

my @Gene = ();
my @Gene_expression = ();
my @Interested_Gene1 = ();
my @Interested_Gene2 = ();
my @data = ();
my $i = 0;
my $Sum = 0;
my $Ave = 0;
#my $status = `mkdir Matrix`;
#my $Todir = "Matrix";

#get the gene names of the WGM gene data 
while(my $line = <F1IN>)
{
	chomp($line);
	$Gene [$i] = $line;
	$i+=1;
}	
close(F1IN);


$i = 0;
#get the gene names of the interested geneset1 
open(F2IN, $Geneset1) or die("Could not open $Geneset1\n");
while(my $line = <F2IN>)
{
	chomp($line);
	$Interested_Gene1 [$i] = $line;
	$i+=1;
}	
close(F2IN);

$i = 0;
#get the gene names of the interested geneset2 
open(F2IN, $Geneset2) or die("Could not open $Geneset2\n");
while(my $line = <F2IN>)
{
	chomp($line);
	$Interested_Gene2 [$i] = $line;
	$i+=1;
}	
close(F2IN);

my $Output = "$Geneset1-$Geneset2-$dataset-Matrix"; 
open(F1OUT, "> $Output") or die("Could not write to file: $Output\n");	

#Print the first row, the gene names in order
print F1OUT "Gene" . "\t";
for ($i=0; $i<scalar(@Gene);$i++)	      
{
	if ($Gene[$i] ~~ @Interested_Gene2)
	{
		print F1OUT $Gene [$i] . "\t";
	}
}
print F1OUT "\n";

# start to print the matrix data after the first row of gene name
open(F3IN, $dataset) or die("Could not open $dataset\n");
my $j = 0;
while (<F3IN>)
{
	print "GeneNo:\t" . $j . "\n";
	chomp($_);
	my $temp = $Gene[$j];
	if ($temp ~~ @Interested_Gene1)
	{
		print F1OUT $temp . "\t";
		$i=0;
		$Sum = 0;
	   	@data = split /\s+/, $_;
		for ($i=0; $i<scalar(@Gene);$i++)	      
		{

			if ($Gene[$i] ~~ @Interested_Gene2)
			{
				print F1OUT $data [$i] . "\t";
				$Sum = $Sum + $data [$i]; 
			}
		}
		$Ave = $Sum/scalar(@Interested_Gene2);
		my $num = eval sprintf('%.4f', $Ave);
		print F1OUT $num;
		print F1OUT "\n";
	}
	$j+=1;
}

close (F1OUT);
#my $status = `mv $Output $Todir`;

close(F3IN);

exit(1); 




