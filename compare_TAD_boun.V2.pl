#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long qw(:config no_ignore_case);
my $USAGE = qq{
Name:
	$0
Fuction:
	This is used to compare two sample boundry;
Useage:
	$0 options
Options:
	-h|help		Print this help
	-o|out 	<str>   dir  for output file;
	-in	<str>	TAD file list and DI file for samples;
	-res <int>   resolution for TAD;
	-num <int>   number for random distribution;
	-lev <float>  level for select;
	-Len <int> threshold for  TAD length,
	-fai <str> fa.fai file,
	-gff <str> gff file,
	-fun <str> functional file,
Author:
	luhongfeng\@novogene.cn
Version:
	1.0; Wed Jan 13 16:33:46 CST 2016
 ;
};
my ($help,$in,$out,$res,$num,$lev,$len,$fai,$gff,$fun);
GetOptions(
	"h|help" =>\$help,
	"in=s" =>\$in,
	"o|out=s" =>\$out,
	"res=i" =>\$res,
	"num=i" =>\$num,
	"lev=f" =>\$lev, 
	"Len=i" =>\$len,
	"fai=s" =>\$fai,
	"gff=s" =>\$gff,
	"fun=s" =>\$fun,
);
$out ||="Result";

die "$USAGE" if ($help);

`mkdir -p $out/00.bin`;
`mkdir -p $out/01.result`;
`mkdir -p $out/02.plot`;


my $bin="/lustre/liyan/scripts/Pipeline/TAD.boundary/bin/";
my $Rscript="/lustre/liyan/01.software/03.R/bin/Rscript";
my $perl="/lustre/liyan/01.software/applink/perl";

my @sam;
my @TAD;
my @DI;
open IN,"$in";
while(<IN>){
	chomp;
	my @arr=split(/\t/,$_);
	push @sam,$arr[0];
	push @TAD,$arr[1];
	push @DI,$arr[2];
}
close IN;
my $n=@sam;

if ($n==2){
	open OUT, ">$out/00.bin/ran.sh";
	my $s1=$sam[0];
	my $s2=$sam[1];
	my $name=join("_",@sam);
	print OUT "$Rscript $bin/random.R $DI[0] $DI[1] $num $out/01.result/$s1.$s2.random\n";
	print OUT "$perl $bin/getThr.pl $out/01.result/$s1.$s2.random $lev $out/01.result/$s1.$s2.random.th\n";
	print OUT "$perl $bin/getBoundry.pl $TAD[0] $out/01.result  $len $s1\n";
	print OUT "$perl $bin/getBoundry.pl $TAD[1]  $out/01.result  $len $s2\n";
	print OUT "$perl $bin/CombineMidV2.pl $out/01.result $s1,$s2 TAD.boun.txt $out/01.result/$s1.$s2.combine.boun.txt\n";
	print OUT "$perl $bin/getMid.pl $out/01.result/$s1.$s2.combine.boun.txt $res $out/01.result/$s1.$s2.combine.boun.Mid.txt\n";
	print OUT "$perl $bin/getDI.pl $out/01.result/$s1.$s2.combine.boun.Mid.txt $DI[0] $DI[1] $s1 $s2 $out/01.result/$s1.$s2.Mid.DI\n";
	print OUT "$Rscript $bin/stat.R $out/01.result/$s1.$s2.Mid.DI $out/01.result/$s1.$s2.Mid.DI.cor\n";
	print OUT "$perl $bin/stat_TAD_boun.pl $out/01.result/$s1.$s2.random.th $out/01.result/$s1.$s2.Mid.DI.cor $out/01.result  $s1,$s2\n";
	print OUT "mkdir -p $out/01.result/$s1\nmkdir -p $out/01.result/$s2\n";
	print OUT "$perl $bin/getGene.pl $gff $fun $out/01.result/$s1.spe.TAD.boun.txt $out/01.result/$s1\n";
	print OUT "$perl $bin/getGene.pl $gff $fun $out/01.result/$s2.spe.TAD.boun.txt $out/01.result/$s2\n";
	print OUT "rm $out/01.result/$s1.TAD.boun.txt $out/01.result/$s2.TAD.boun.txt $out/01.result/$s1.$s2.random  $out/01.result/$s1.$s2.combine.boun.Mid.txt  $out/01.result/$s1.$s2.Mid.DI\n";
	close OUT;

	open OUT1, ">$out/02.plot/list";
	print OUT1 "$s1\t$TAD[0]\n$s2\t$TAD[1]\n";
	close OUT1;

	open OUT, ">$out/00.bin/plot.sh";
	print OUT "$Rscript $bin/TAD.whole.genome.R $out/02.plot/list $fai $out/02.plot 1 9 10 $name\n";
	close OUT;

	open IN, ">$out/00.bin/do.sh";
	print IN "qsub -q tangqianzi -l nodes=1:ppn=4 -l mem=1gb $out/00.bin/ran.sh\n";
	print IN "qsub -q tangqianzi -l nodes=1:ppn=4 -l mem=1gb $out/00.bin/plot.sh\n";
	close IN;
}elsif($n==3){
	open OUT, ">$out/00.bin/ran.sh"; 
	my $s1=$sam[0];
	my $s2=$sam[1];
	my $s3=$sam[2];
	my $name=join("_",@sam); 
	print OUT "$Rscript $bin/random.R $DI[0] $DI[1] $num $out/01.result/$s1.$s2.random\n";
	print OUT "$Rscript $bin/random.R $DI[0] $DI[2] $num $out/01.result/$s1.$s3.random\n";
	print OUT "$Rscript $bin/random.R $DI[1] $DI[2] $num $out/01.result/$s2.$s3.random\n";
	print OUT "$perl $bin/getThr.pl $out/01.result/$s1.$s2.random $lev $out/01.result/$s1.$s2.random.th\n";
	print OUT "$perl $bin/getThr.pl $out/01.result/$s1.$s3.random $lev $out/01.result/$s1.$s3.random.th\n";
	print OUT "$perl $bin/getThr.pl $out/01.result/$s2.$s3.random $lev $out/01.result/$s2.$s3.random.th\n";
	print OUT "cat $out/01.result/*.random.th >$out/01.result/random.all.th\n";
	print OUT "$perl $bin/getBoundry.pl $TAD[0] $out/01.result  $len $s1\n$perl $bin/getBoundry.pl $TAD[1] $out/01.result  $len $s2\n$perl $bin/getBoundry.pl $TAD[2] $out/01.result  $len $s3\n";
	print OUT "$perl $bin/CombineMidV2.pl $out/01.result $s1,$s2,$s3  TAD.boun.txt $out/01.result/$s1.$s2.$s3.combine.boun.txt\n";
	print OUT "$perl $bin/getMid.pl $out/01.result/$s1.$s2.$s3.combine.boun.txt $res $out/01.result/$s1.$s2.$s3.combine.boun.Mid.txt\n";
	print OUT "$perl $bin/getDI_for_three.pl $out/01.result/$s1.$s2.$s3.combine.boun.Mid.txt $DI[0] $DI[1] $DI[2] $s1 $s2 $s3 $out/01.result/$s1.$s2.$s3.Mid.DI\n";
	print OUT "$Rscript $bin/stat_for_three.R $out/01.result/$s1.$s2.$s3.Mid.DI $out/01.result/$s1.$s2.$s3.Mid.DI.cor\n";
	print OUT "$perl  $bin/stat_TAD_boun_for_three.pl $out/01.result $out/01.result/$s1.$s2.$s3.Mid.DI.cor $out/01.result $s1,$s2,$s3\n";
	print OUT "rm $out/01.result/*.random $out/01.result/*.random.th  $out/01.result/*.TAD.boun.txt  $out/01.result/$s1.$s2.$s3.combine.boun.Mid.txt $out/01.result/$s1.$s2.$s3.Mid.DI\n";
	close OUT;

	open OUT1, ">$out/02.plot/list";
	print OUT1 "$s1\t$TAD[0]\n$s2\t$TAD[1]\n$s3\t$TAD[2]\n";
	close OUT1;

	open OUT, ">$out/00.bin/plot.sh";
	print OUT "$Rscript $bin/TAD.whole.genome.R $out/02.plot/list $fai $out/02.plot 1 9 10 $name\n";
	close OUT;
	
	open IN, ">$out/00.bin/do.sh";
	print IN "qsub -q tangqianzi -l nodes=1:ppn=4 -l mem=1gb $out/00.bin/ran.sh\n";
	print IN "qsub -q tangqianzi -l nodes=1:ppn=4 -l mem=1gb $out/00.bin/plot.sh\n";	
	close IN;
}
