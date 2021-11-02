
#This script is used to extract SNPs in CDS regions from hapmap files.

use strict;
use warnings;

my $in=shift;
my $in2=shift;

die "perl $0 <gff> <hmp>\n" unless ($in2);

open IN,$in or die;
my (%info,%g);
while (<IN>)
{	next if($_=~/#/);
	chomp;
	my @a=split /\t/,$_;
#	next unless ($a[2] eq "gene");
	$a[0]=~s/Chr//g;
	my ($gene,$size);
	if ($a[2] eq "mRNA")
	{
	$gene=(split /=|.v/,$a[-1])[1];
	$size=($a[4]-$a[3]+1)/1000;
	$info{$gene}="$a[0]\t$a[3]\t$a[4]\t$size";
	}
	next unless ($a[2] eq "CDS" );#or $a[2]=~/UTR/);
	my $k="$a[0]\t$a[3]\t$a[4]\t$size";
#	my $gene=(split /=|.v/,$a[-1])[1];
	push @{$g{$gene}},$k;
	#print "$gene\t$k\n";
}

print "BLOCK\tCHR\tBP1\tBP2\tSIZE\tNSNPS\tSNPS\n";

open IN2,$in2 or die;
#open IN3,"head -n 1000 $in2|" or die;
my %chr;
while (<IN2>)
{
	next if ($.==1);#=~/#/);
	my @hmp=(split /\t/,$_,5);
	my $snp="$hmp[3]\t$hmp[0]";
	push @{$chr{$hmp[2]}},$snp;
	#print "$snp a \n";
}
close IN2;

foreach my $gene(sort {$a<=>$b} keys %g)
{
	
	my @k=split /\t/,$info{$gene};
	my %snp;
	foreach (@{$g{$gene}})
	{
	my @k=split /\t/,$_;
#	my %snp; #print "$k[0]\n";
	foreach my $hmp(@{$chr{$k[0]}})
	{#	print "$hmp\ta";
		my @hmp=split /\t/,$hmp;
#		print "$hmp[0]\t$k[1]\t$k[2]\t$hmp[1]\n";
		next unless ($hmp[0] >=$k[1] && $hmp[0] <=$k[2]);
	#	push @snp,$hmp[1];
#		print "$hmp[0]\t$k[1]\t$k[2]\t$hmp[1]\n";
		$snp{$hmp[1]}=1;
	}
	}
#	my $size=($k[2]-$k[1]+1)/1000;
	my $size=$k[-1];
	my @snp= keys %snp;
	my $NSNPS=@snp;
	my $SNPS=join "|",@snp;
	next unless (@snp>0);
#	next if (@snp>30);
	print "$gene\t$k[0]\t$k[1]\t$k[2]\t$size\t$NSNPS\t$SNPS\n";
}



