#!/usr/bin/perl

use strict;

open(D,$ARGV[0]) || die "$!\n";
#Protein 888. Allele HLA-C*08:02. Number of high binders 0. Number of weak binders 0. Number of peptides 9

my %prot2binding=();

while(<D>)
{
	chomp($_);
	#Protein 14. Allele HLA-A*11:01. Number of high binders 0. Number of weak binders 0. Number of peptides 9
	if($_=~/^Protein\s+(\d+).+(HLA\-\S+)\.\s+Number\ of\ high\ binders\ (\d+)\.\ Number\ of\ weak\ binders\ (\d+)/)
	{
		my $prot = $1;
		my $all = $2;
		my $hb = $3;
		my $wb = $4;

		print STDERR "$prot\t$all\t$hb\t$wb\n";

		if($hb>0)
		{
			$prot2binding{$prot} = 2;
		}
		elsif($wb>0 & ((!defined $prot2binding{$prot}) | $prot2binding{$prot}==0))
		{
			$prot2binding{$prot} = 1;
		}
		elsif($hb==0 & $wb ==0 & ((!defined $prot2binding{$prot}) | $prot2binding{$prot}==0))
		{
			$prot2binding{$prot} = 0;
		}

	}
}	
close(D);

for my $key (sort keys %prot2binding)
{
	print "$key\t$prot2binding{$key}\n";
}
