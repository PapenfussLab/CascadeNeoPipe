#!/usr/bin/perl

use strict;

#open(D,$ARGV[0]) || die "$!\n";
my $folder = $ARGV[0];

my @dirs = `ls $folder/\*\_epitopes\.out`;

my $id;

my %epi2id=();
my %epi2sample=();

for(my $i=0;$i<@dirs;$i++)
{
	chomp($dirs[$i]);
	open(D,$dirs[$i]) || die "$!\n";

	$dirs[$i]=~/$folder\/(\S+)\_epitopes.out/;

	print STDERR "Working with $dirs[$i]\n";

	my $sample_name = $1;

	while(<D>)
	{
		chomp($_);
		if($_=~/\>(ENST00\S+\ \S+)/)
		{
			$id = $1;
		}
		else
		{
			if($epi2id{$_}!~/$id\;/)
			{
				#my $s_epi=substr($_,0,15);
				$epi2id{$_}.=$id . ';';
			}
			if($epi2sample{$_}!~/$sample_name\;/)
			{
				$epi2sample{$_}.=$sample_name . ';';
			}	
		}
	}
	close(D);
}	

open(ET,">$folder/epitope2transcripts.list");
open(EX,">$folder/epitope2transcripts.excluded.log");
open(IN_NHC,">$folder/input_netMHC_noRedundancy.fsa");

my $count_id=0;

my $epitope_length = (17 +1)/2;

for my $key (sort keys %epi2id)
{
	$count_id++;
	#exclude cases where the stop codon occurs within the first window OR when the extracted kmer is shorter than the window
	if (($key=~/\*/ & (index($key,'*') +1 <= $epitope_length)) | length($key)<$epitope_length)
	{
		print EX "$count_id\t$key\t$epi2id{$key}\t$epi2sample{$key}\n";
		next;
	}
	print ET "$count_id\t$key\t$epi2id{$key}\t$epi2sample{$key}\n";
	print IN_NHC '>',$count_id,"\n";
	print IN_NHC $key,"\n";

	
}
close(ET);
close(EX);
close(IN_NHC);
