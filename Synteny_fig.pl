
#This script is used to display the synteny map between 3 genomes using the perl SVG module.

use SVG ;
$in=shift;
$in2=shift;
$len=shift;
$len2=shift;
$len3=shift;

die "perl $0 <syn_result1.txt> <syn_result2.txt> <chr_len1.txt> <chr_len2.txt> <chr_len3.txt>\n" unless ($len3);

my $svg = SVG->new(width=>1500, height=>1000);

open L,$len or die;
while (<L>)
{
	chomp;
	my @a=split /\t/,$_;
	$l{$a[0]}=$a[1];
}
open L2,$len2 or die;
while (<L2>)
{
        chomp;
        my @a=split /\t/,$_;
        $l2{$a[0]}=$a[1];
}

open L3,$len3 or die;
while (<L3>)
{
        chomp;
        my @a=split /\t/,$_;
        $l3{$a[0]}=$a[1];
}


@chr=("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10");

$x=20;
foreach (@chr)
{
my $lenB=$l{$_}/100000;
my $lenE=$l2{$_}/100000;
my $lenR=$l3{$_}/100000;
$loc{$_}=$x;
$loc2{$_}=$x+50;
$loc3{$_}=$x+50+50;
$svg->rect(x => $x, y => 50, width =>20 , height => $lenB,fill=>"blue",stroke=>"black");
$svg->text(x => $x+50, y => 40,width => 20, height => 5,  "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>16,"-cdata" => "$_" );
$svg->rect(x => $x+50, y => 50, width =>20 , height => $lenE,fill=>"green",stroke=>"black");
$svg->rect(x => $x+50+50, y => 50, width =>20 , height => $lenR,fill=>'rgb(93,15,249)',stroke=>"black");
#$svg->text(x => $x+75+5, y => 40, "font-family"=>"Arial", "text-anchor"=>"end","font-size"=>16,"-cdata" => "$_" );
$x=$x+50+50+75;
}

open IN,$in or die;
while(<IN>)
{	next if ($.==1);
	chomp;
	my @a=split /\t/;
	next unless ($a[0] eq $a[3]);
	my $xv = [$loc{$a[0]}+22,$loc2{$a[0]}-2,$loc{$a[0]}+22,$loc2{$a[3]}-2];
	my $yv = [50+$a[4]/100000,50+$a[1]/100000,50+$a[5]/100000,50+$a[2]/100000];
	my $points = $svg->get_path(x=>$xv,y=>$yv,-type =>'polygon');
	$svg->polygon(%$points,style => {'stroke'=>'grey'});
	#$svg->line(x1 =>$loc{$a[0]}+22,y1=> 50+$a[1]/100000,x2=>$loc2{$a[0]}-2,y2=>50+$a[4]/100000,style=>{'stroke-width'=>'0.1','stroke'=>'grey'} );
	#$svg->line(x1 =>$loc{$a[0]}+22,y1=> 50+$a[2]/100000,x2=>$loc2{$a[3]}-2,y2=>50+$a[5]/100000,style=>{'stroke-width'=>'0.1','stroke'=>'grey'});
}


open IN2,$in2 or die;
while(<IN2>)
{       next if ($.==1);
        chomp;
        my @a=split /\t/;
        next unless ($a[0] eq $a[3]);
        my $xv = [$loc2{$a[0]}+22,$loc3{$a[0]}-2,$loc2{$a[0]}+22,$loc3{$a[3]}-2];
        my $yv = [50+$a[1]/100000,50+$a[4]/100000,50+$a[2]/100000,50+$a[5]/100000];
        my $points = $svg->get_path(x=>$xv,y=>$yv,-type =>'polygon');
        $svg->polygon(%$points,style => {'stroke'=>'grey'});
        #$svg->line(x1 =>$loc{$a[0]}+22,y1=> 50+$a[1]/100000,x2=>$loc2{$a[0]}-2,y2=>50+$a[4]/100000,style=>{'stroke-width'=>'0.1','stroke'=>'grey'} );
        ##$svg->line(x1 =>$loc{$a[0]}+22,y1=> 50+$a[2]/100000,x2=>$loc2{$a[3]}-2,y2=>50+$a[5]/100000,style=>{'stroke-width'=>'0.1','stroke'=>'grey'});
      }
        


my $out = $svg->xmlify;
open SVGFILE, ">BER.svg";
print SVGFILE $out;

