#!/usr/bin/perl

##############################################################################
#
#                                 Cross-Correlation
#
# Author:  Xiaobing Yuan (xbyuan@brain.ecnu.edu.cn)
# Written: 10th Aug 2015
# 
# DESCRIPTION:
# This script is designed to calculate the cross-correlation of gene expression 
# 
# 
# Copyright 2021 Xiaobing Yuan
#                Institute of Brain Functional Genomics
#                ECNU, Shanghai

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


## CONSTANTS ##
my $TRUE              = 1;
my $FALSE             = 0;
my $DEBUG             = $FALSE;
my $EXITSTATUS        = 0;

# Default umask
umask 027;

# Define variables


##############################
#   M A I N   P R O G R A M  #
##############################

# Check to see we received 1 files in the arguments
if(scalar(@ARGV) != 1)
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

my $Input_Profile = shift(@ARGV);

my $Output = "Correlation_Coefficient_$Input_Profile"; 


# First count how many lines are there in the input file
my $lineCounter = 0;

open(F1IN, $Input_Profile) or die("Could not open $Input_Profile\n");
open(F1OUT, "> $Output") or die("Could not write to file: $Output\n");
while(my $line = <F1IN>)
{
		chomp($line);
		$lineCounter++;
}	
close(F1IN);
print $lineCounter . "\n";

my $x = [];
my @TestGene = ();
			
for(my $i=0;$i<$lineCounter;$i++)
{
# Get the test gene expression data 

		open(F1IN, $Input_Profile) or die("Could not open $Input_Profile\n");
		print "Line NO.: \t" . $i . "\n";
		my $j = 0;
		while(my $line = <F1IN>)
		{
			chomp($line);
			if ($j == $i)
			{
				@TestGene = split /\s+/, $line;

				last;
			}
 
			$j++;
		}	
		close(F1IN);

my $DataNum = @TestGene;
#print $DataNum . "\n";
$DataNum += 1;
	
# Compare test gene data with every line to get the correlation coefficient  

		open(F1IN, $Input_Profile) or die("Could not open $Input_Profile\n");



		while(my $line = <F1IN>)
		{
			#chomp($line);
			my @EachGene = split /\s+/, $line;

				# convert data into column 
				for (my $n=1;$n<$DataNum;$n++)
				{
					  $x->[$n][1]= $TestGene[$n-1];
										   
					  $x->[$n][2]= $EachGene[$n-1];

						#print "data is :" . $x->[$n][1] . "\t" . $x->[$n][2] . "\n";
				}

			# Get the correlation coefficient and print data into the output file column by column seperated by \t 

			my $correlation = correlation($x);


 

			#generate an anonymous 2D array where $x->[1] is the row
			#$x->[1][1] is the value in row 1 column 1 and $x->[1][2] is the value of row 1 column 2
			#once you build the entire array, pass it to the correlation subroutine as above
			#my $corrleation = correlation($x)
 
			#if you want to see what's inside $x use the code below
			#for (my $i = 1; $i < scalar(@{$x}); ++$i){
			#   my $line_to_print = '';
			#   for (my $j = 1; $j < scalar(@{$x->[$i]}); ++$j){
			#      $line_to_print .= "$x->[$i]->[$j]\t";
			#   }
			#   $line_to_print =~ s/\t$//;
			#   print "$line_to_print\n";
			#}
						
			print F1OUT $correlation . "\t";
		}
		# At the end of a data line, print \n to jump to the next line  	
		
		close(F1IN);
                print F1OUT "\n";

}		

close (F1OUT);

exit(1); 


sub mean 
{
   my ($x)= @_;
   my $num = @{$x};
#	print $num . "\n";
   my $sum_x = 0;
   my $sum_y = 0;
   for (my $i = 1; $i < $num; $i++)
   {
      $sum_x = $sum_x + ($x->[$i][1]);
      $sum_y = $sum_y + ($x->[$i][2]);
   }


   my $mu_x = $sum_x / ($num-1);
   my $mu_y = $sum_y / ($num-1);


   return($mu_x,$mu_y);
}
 
### ss = sum of squared deviations to the mean
sub ss 
{
   my ($x,$mean_x,$mean_y,$one,$two)= @_;
   my $num1 = scalar(@{$x});
   my $sum = 0;
   for (my $d = 1;$d < $num1; $d++)
   {
     $sum += ($x->[$d][$one]-$mean_x)*($x->[$d][$two]-$mean_y);
   }
   return $sum;
}
 
sub correlation 
{
   my ($x) = @_;
   my ($mean_x,$mean_y) = mean($x);
   my $ssxx=ss($x,$mean_x,$mean_y,1,1);
   my $ssyy=ss($x,$mean_x,$mean_y,2,2);
   my $ssxy=ss($x,$mean_x,$mean_y,1,2);
   my $correl=correl($ssxx,$ssyy,$ssxy);
   my $xcorrel=sprintf("%.4f",$correl);
   return($xcorrel);
 
}
 
sub correl 
{
   my($ssxx,$ssyy,$ssxy)= @_;
   my $sign=$ssxy/abs($ssxy);
   my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
   return $correl;
} 


