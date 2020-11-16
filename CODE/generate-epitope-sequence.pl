#!/usr/bin/perl

use strict;

sub min
{
	my $a=shift;
	my $b=shift;
	if($a<$b)
	{
		return $a;
	}
	else
	{
		return $b;
	}
}


sub max
{
        my $a=shift;
        my $b=shift;
        if($a>$b)
        {
                return $a;
        }
        else
        {
                return $b;
        }
}

open(GVF,$ARGV[0]) || die "$!\n";
#17	CooVar	SNV	32598166	32598166	.	+	.	ID=snp_812;Variant_seq=G;Reference_seq=T;Variant_type=non_conservative_missense_codon;Variant_effect=non_conservative_missense_codon 0 mRNA ENST00000394627;Note=ttg>tGT_L>C_aa44_codon_loc2_RADICAL(198)


open(VARPEP,$ARGV[1]) || die "$!\n";

open(TSV,$ARGV[2]) || die "$!\n";

my $epitope_length = $ARGV[3];

my $window = ($epitope_length+1)/2;

print STDERR "For length $epitope_length generate $window\-mers\n";


my %trans2var=();
my %trans2coord=();

#Need to get deletion information, as coovar doesn't provide it

my $index_codonchange=0;
my $index_varianttype=0;
my $index_feature=0;

my %som_coord=();
my %del_trans2var=();
my %del_trans2coord=();

my %seen=();

while(<TSV>)
{
	chomp($_);
	
	my @line = split(/\t/,$_);
	if($_=~/^\#/)
	{
		for(my $i=0;$i<@line;$i++)
		{
			if($line[$i] eq 'HGVSp')
			{
				$index_codonchange = $i;
			}
			if($line[$i] eq 'Variant_Type')
			{
				$index_varianttype= $i;
			}
			if($line[$i] eq 'Feature')
			{
				$index_feature = $i;
			}
		}
	}
	next if ($_=~/^\#/);
	next if ($_!~/frame|missense/);

	next if ($line[$index_varianttype]!~/SNV|deletion/);
	my $pos = $line[2];

	if($line[$index_varianttype]=~/deletion/)
	{
		$pos+=1;
	}
	next if (defined $seen{"$line[$index_feature]\t$line[1]\-$pos"});
	$seen{"$line[$index_feature]\t$line[1]\-$pos"}++;


	$line[$index_codonchange]=~/p\.\S\S\S(\d+)/;
	my $codon = $1;

	if($line[$index_varianttype]=~/deletion/)
	{
		$del_trans2coord{$line[$index_feature]}.= $line[1] . "\-" . $pos . ";" ;
		$del_trans2var{$line[$index_feature]}.=$codon . ";";
	}
	$som_coord{"$line[$index_feature]\t$line[1]\-$pos"}++;

	print STDERR "$line[$index_feature]\t$line[1]\t$pos\t$codon\n";

}
close(TSV);

while(<GVF>)
{
	chomp($_);

	next if ($_!~/missense|frame/);
	my @line = split(/\t/,$_);
	my @vars = split(/\;/,$line[8]);
	

	if($_=~/missense/)
	{
		$_=~/Variant_effect.+mRNA\ (\S+)\;Note\=(\S+)/;

		my $coord = $line[0] . "\-" . $line[3]; 
		my $transcripts = $1;
		my $changes = $2;

		my @trans = split(/\,/,$transcripts);
		my @chang = split(/\,/,$changes);

		if(scalar(@trans)!=scalar(@chang))
		{
			print STDERR "Warn: $line[0] $line[3] vector of trans and vector of changes not same length";
		}
		else
		{
			for(my $i=0;$i<@trans;$i++)
			{
				$chang[$i]=~s/\S+aa(\d+)\_codon\S+/$1/;
				next if(!defined $som_coord{"$trans[$i]\t$coord"});
				$trans2var{$trans[$i]}.=$chang[$i] . ";";
				$trans2coord{$trans[$i]}.=$coord . ";";
			}
		}
	}
	elsif ($_=~/frame/)
	{
		$_=~/Variant_effect.+mRNA\ (\S+)\s*/;

		my $coord = $line[0] . "\-" . $line[3]; 
		my $transcripts = $1;

		my @trans = split(/\,/,$transcripts);
		for(my $i=0;$i<@trans;$i++)
		{
			next if(!defined $som_coord{"$trans[$i]\t$coord"});
			$trans2var{$trans[$i]}.=$del_trans2var{$trans[$i]};
			$trans2coord{$trans[$i]}.=$del_trans2coord{$trans[$i]};
		}
	}
	else
	{
		print STDERR "Not identified as missense or frameshift deletion: $_\n";
		last;
	}
	
}
close(GVF);

print STDERR keys(%trans2var) . " transcripts to consider\n";

my $transcript;

while(<VARPEP>)
{
	chomp($_);
	if($_=~/^\>(\S+)/)
	{
		$transcript = $1;
	}
	else
	{
		next if (!defined $trans2var{$transcript});
		my @sequence = split(//,$_);
		my @vars = split(/\;/,$trans2var{$transcript});
		my @coords = split(/\;/,$trans2coord{$transcript});

		for(my $i=0;$i<@vars;$i++)
		{
			
			my $final_seq="";
			for(my $j=max($vars[$i]-$window+1,0);$j<=min($vars[$i]+$window-1,$#sequence+1);$j++)
			{
				$final_seq.=$sequence[$j-1];
			}
			next if (($final_seq=~/\*/ & (index($final_seq,'*') +1 <= $window)) | length($final_seq)<$window);

			if($final_seq=~/\*/ & (index($final_seq,'*') +1 > $window))
			{
				$final_seq = substr($final_seq,0,index($final_seq,'*'));
			}
			print "\>$transcript\_$vars[$i]\ $coords[$i]\n";
			print "$final_seq\n";
		}
	}
}
close(VARPEP);

