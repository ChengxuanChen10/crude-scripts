
#This script is used to filter low quality SNPs and convert to maker-type table.

use strict;
use warnings;

my $in=shift;
my $miss_rate=shift;
my $chi_p=shift;

die "perl $0 <SNP.table> <miss_rate:0.25> <chi_p:0.05> > <output> \n" unless ($in);

my $miss_rate=0.25 unless ($miss_rate);
my $chi_p=0.05 unless ($chi_p);


#my @output;
open IN ,$in or die;
while (<IN>)
{
	chomp;
	
	if ($.==1)
	{
	my @head=split(/\s+/,$_);
	s/\s+//g foreach @head;
	@head=@head[6..$#head];
	print "Maker\tType\t";
	print  join "\t",@head;
	print "\n";
	next;
	}
	
	my @a=split /\t/,$_;
	$a[0]=~s/Chr//g;
	$a[0]="$a[0]_$a[1]";
#	next if ($a[4] eq $a[5]);
	my ($type,@g);
	if ( ($a[4] eq "H") && ($a[5] eq "H") )	
	{
		$type="<hkxhk>";my ($ml,%count);
		for $i (6..$#a)
		{	
			my $g;
			if ($a[$i] eq "H"){$g="hk";}
			elsif ($a[$i] eq $a[2]){$g="hh";}
			elsif($a[$i] eq $a[3]){$g="kk";}
			else {$g="--";}
			#print $g;exit;
			$count{$g}+=1 unless ($g eq "--");
			$ml+=1 unless ($g eq "--");
			push @g,$g;	
		}
		
		next if ($ml < (1-$miss_rate)*@g);  # miss rate check
		my (@obs,@rat,@chi);
		@obs=($count{"hh"},$count{"hk"},$count{"kk"});
        	@rat=(1,2,1);
       		@chi=&chisq(\@obs,\@rat);
		next if ($chi[-1] < $chi_p); # chi-test pvalue check

		my $geno=join "\t",@g;
		my $out= "$a[0]\t$type\t$geno\n";
		print "$out";next;
	#	push @output,$out;
	}
	
    elsif( ($a[4] eq "H")&& (($a[5] eq $a[2]) || ($a[5] eq $a[3])))
	{
		$type="<lmxll>";my ($ml,%count);
		for $i (6..$#a)
		{	
			my $g;
			if ($a[$i] eq "H"){$g="lm";}
			elsif ($a[$i] eq $a[2]){$g="ll";}
			elsif($a[$i] eq $a[3]){$g="ll";}
			else {$g="--";}
			$count{$g}+=1 unless ($g eq "--");
			$ml+=1 unless ($g eq "--");
			push @g,$g;	
		}
		
		next if ($ml < (1-$miss_rate)*@g);  # miss rate check
		my (@obs,@rat,@chi);
		@obs=($count{"lm"},$count{"ll"});
        	@rat=(1,1);
        	@chi=&chisq(\@obs,\@rat);
		next if ($chi[-1] < $chi_p); # chi-test pvalue check

		my $geno=join "\t",@g;
		my $out= "$a[0]\t$type\t$geno\n";
		print "$out";next;
	#	push @output,$out;
	}
	
	elsif(($a[5] eq "H") &&  (($a[4] eq $a[2]) || ($a[4] eq $a[3])))
 	{
		$type="<nnxnp>";my ($ml,%count);
		for $i (6..$#a)
		{	
			my $g;
			if ($a[$i] eq "H"){$g="np";}
			elsif ($a[$i] eq $a[2]){$g="nn";}
			elsif($a[$i] eq $a[3]){$g="nn";}
			else {$g="--";}
			$count{$g}+=1 unless ($g eq "--");
			$ml+=1 unless  ($g eq "--");
			push @g,$g;	
		}
		
		next if ($ml < (1-$miss_rate)*@g);  # miss rate check
		my (@obs,@rat,@chi);
		@obs=($count{"np"},$count{"nn"});
        	@rat=(1,1);
        	@chi=&chisq(\@obs,\@rat);
		next if ($chi[-1] < $chi_p); # chi-test pvalue check

		my $geno=join "\t",@g;
		my $out= "$a[0]\t$type\t$geno\n";
		print "$out";next;
	#	push @output,$out;
	}
}

sub chisq
{
    use Statistics::Distributions qw< chisqrprob >;
    use List::Util qw< sum >;
    my $ref_array_obs=shift;
    my $ref_array_rat=shift;
    my $chisq=0;
    my $pval=1;
    my $total=sum(@$ref_array_obs);
    my $totalexp=sum(@$ref_array_rat);
    my @ratios=map{sprintf("%.2f",$_/$totalexp)} @$ref_array_rat;
    my $degree= scalar @$ref_array_obs -1;

    for(my$i=0;$i<scalar@$ref_array_obs;$i++)
	{
       $chisq += ( ( $$ref_array_obs[$i] - $ratios[$i]*$total )**2 ) / ( $ratios[$i]+0.00000000000000000000001);
    }
    return ($chisq,$degree,chisqrprob($degree,$chisq));
}                       
