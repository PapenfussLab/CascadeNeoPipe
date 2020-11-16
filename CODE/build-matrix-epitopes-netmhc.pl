#!/usr/bin/perl

use strict;

open(EXP,"DATA/TCGA_expression_medians.out") || die "$!\n";

my $case = $ARGV[0];
my $len = $ARGV[1];

my %gene2exp=();

while(<EXP>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	$gene2exp{$line[0]}=$line[1];
}
close(EXP);

my @snvs = `ls DATA/CALLS/$case\_\*.coovar`;
#1-4_filtered_10X.coovar
my %snv2sample=();

my %samples=();
my %snvs=();


for(my $i=0;$i<@snvs;$i++)
{
	chomp($snvs[$i]);
	#next if ($case eq 'E' & $snvs[$i]=~/E_SWISS/);
	#next if ($case eq 'F' & $snvs[$i]=~/F_SWISS/);
	#next if ($case!~/F_SWISS/);

	open(D,$snvs[$i]) || die "$!\n";

	$snvs[$i]=~/DATA\/CALLS\/(\S+)\.coovar/;

	my $sample = $1;

	$samples{$sample}++;

	while(<D>)
	{
		chomp($_);
		my @line = split(/\t/,$_);
		$snv2sample{"$sample\t$line[1]\-$line[2]"}++;
		$snvs{"$line[1]\-$line[2]"}++;
	}
	close(D);

}

my %snv2epitope=();
my %snv2epitopeid=();

open(D,"EPITOPES\_$len/$case/epitope2transcripts.list") || die "$!\n";

my %seen=();

my %snv2transcript=();

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	my @snv_data =split(/\;/,$line[2]);
	
	for(my $i=0;$i<@snv_data;$i++)
	{
		if($snv_data[$i]=~/(\S+)\_\d+\s+(\S+)$/)
		{
			my $trans = $1;
			my $var = $2;
			$snv2transcript{"$var\t$trans"}++;
			next if (defined $seen{"$var\t$line[1]\t$trans"});
			$seen{"$var\t$line[1]\t$trans"}++;
			$snv2epitope{"$var\t$trans"}.= $line[1] . " ";
			$snv2epitopeid{"$var\t$trans"}.= $line[0] . " ";
		}
	}
}
close(D);

my @tsvs = `ls DATA/TSV/\*\.tsv`;

my %snv2gene=();

for(my $i=0;$i<@tsvs;$i++)
{
	chomp($tsvs[$i]);
	open(D,$tsvs[$i]) || die "$!\n";

	my $index_symbol="-1";
	my $index_feature = "-1";

	while(<D>)
	{
		chomp($_);
		if($_=~/^\#CHROM/)
		{
			my @line  =split(/\t/,$_);
			for(my $j=0;$j<@line;$j++)
			{
				if($line[$j] eq 'SYMBOL')
				{
					$index_symbol = $j;
				}
				if($line[$j] eq 'Feature')
				{
					$index_feature = $j;
				}
			}
		}
		if($_!~/^#/)
		{
			my @line  =split(/\t/,$_);
			if(defined $snv2transcript{"$line[0]\-$line[1]\t$line[$index_feature]"})
			{
				$snv2gene{"$line[0]\-$line[1]"} = $line[$index_symbol];
			}
		}
	}
}

my %snv2epitope_ref=();
my %snv2epitopeid_ref=();

open(D,"EPITOPES_REF\_$len/$case/epitope2transcripts.list") || die "$!\n";

my %seen_ref=();

my %snv2transcript_ref=();

while(<D>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	my @snv_data =split(/\;/,$line[2]);
	
	for(my $i=0;$i<@snv_data;$i++)
	{
		if($snv_data[$i]=~/(\S+)\_\d+\s+(\S+)$/)
		{
			my $trans = $1;
			my $var = $2;
			$snv2transcript_ref{"$var\t$trans"}++;
			next if (defined $seen_ref{"$var\t$line[1]\t$trans"});
			$seen_ref{"$var\t$line[1]\t$trans"}++;
			$snv2epitope_ref{"$var\t$trans"}.= $line[1] . " ";
			$snv2epitopeid_ref{"$var\t$trans"}.= $line[0] . " ";
		}
	}
}
close(D);


open(BIN,"EPITOPES\_$len/$case/I_input_netMHC_noRedundancy.tab") || die "$!\n";

my %epitopeid2score=();

while(<BIN>)
{
	chomp($_);
	my @line = split(/\t/,$_);

	$epitopeid2score{$line[0]} = $line[1];	
}
close(BIN);

open(BIN,"EPITOPES_REF\_$len/$case/I_input_netMHC_noRedundancy.tab") || die "$!\n";

my %epitopeid2score_ref=();

while(<BIN>)
{
	chomp($_);
	my @line = split(/\t/,$_);

	$epitopeid2score_ref{$line[0]} = $line[1];	
}
close(BIN);

print "Gene\;Coord\;Mutated_Epitope";

for my $key (sort keys %samples)
{
	print "\t",$key;
}
print "\tReference_Epitope\tControl\tAvgExpression\n";

my %seen_print=();


for my $key (sort keys %snv2transcript)
{
	next if (!defined $snv2epitope{$key});
	#print $key;

	my @ids = split(/\ /,$snv2epitopeid{$key});
	my @epis = split(/\ /,$snv2epitope{$key});

	my @ids_ref = split(/\ /,$snv2epitopeid_ref{$key});
	my @epis_ref = split(/\ /,$snv2epitope_ref{$key});


	$key=~/(\S+)\t\S+/;

	my $var = $1;



	my @to_print;

	for(my $i=0;$i<@epis;$i++)
	{
		next if (defined $seen_print{"$snv2gene{$var}\;$var\;$epis[$i]"});
		$seen_print{"$snv2gene{$var}\;$var\;$epis[$i]"}++;

		print "$snv2gene{$var}\;$var\;$epis[$i]";

		for my $key2 (sort keys %samples)
		{
			if(defined $snv2sample{"$key2\t$var"} & defined $epitopeid2score{$ids[$i]})
			{
				print "\t",$epitopeid2score{$ids[$i]};
			}
			else
			{
				print "\t-1";
			}
		}
		print "\t$epis_ref[0]\t$epitopeid2score_ref{$ids_ref[0]}";
		if(defined $gene2exp{$snv2gene{$var}})
		{
			print "\t",$gene2exp{$snv2gene{$var}};
		}
		else
		{
			print "\t-1";
		}
		print "\n";
	}
}
