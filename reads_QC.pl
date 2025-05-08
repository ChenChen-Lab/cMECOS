#!/usr/bin/perl 
=head1 Description
	sequence quality control and nohost
=head1 Usage
	perl reads_QC.pl [options]
	general input and output:
	-F 		    F_adapter sequence
	-R 		    R_adapter sequence
        -bl  <a file>       barcode_list.txt(reads barcodes)
	-sl  <a file>       list of sample name and fastq file of read(file suffix can be ".fq.gz"), one strain per line, PE read seperated by ",", and different by "\n"
	-host               host_type,'human' or 'mice' or 'none'
	-o   <a directory>  output dir, default current directory [./01.QC]
	-thd <num>          thread for dsub
	-h                  show help

=head1 Example
	perl reads_QC.pl -F -R -bl barcode_list.txt -sl sample.list -host human -o 

=head1 Version
        Author: guchaoyang0826@163.com
        Date: 2022-08-17
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw/abs_path/;

my ($F_adapter,$R_adapter,$barcodelist,$samplelist,$hosttype,$outdir,$thd,$help);

GetOptions (
	"F:s"=>\$F_adapter,
	"R:s"=>\$R_adapter,
        "bl:s"=>\$barcodelist,
	"sl:s"=>\$samplelist,
	"host:s"=>\$hosttype,
	"o:s"=>\$outdir,
	"h"=>\$help
);
die `pod2text $0` if ($help || !$samplelist);

unless ($outdir) {
	$outdir="./01.QC";
}

unless (-e $outdir) {
	mkdir $outdir;
}

unless ($thd){
	$thd=8;
}

$outdir=abs_path($outdir);

#============================================================================================================================================================

if (!$samplelist){
	print "No input read files\n";
	exit;
}

if (!$hosttype){
        print "No input read type: 'human' or 'mice' or 'none'\n";
        exit;
}

my $split_barcode="/datapool/bioinfo/guchaoyang/pipeline/stLFR/bin/split_barcode3.pl";	
my $SOAPfilter="/datapool/bioinfo/lijr/project/2022_culture/scripts/SOAPfilter_v2.2";
my $bowtie2="/datapool/software/anaconda3/bin/bowtie2";
my $cutadapt="/datapool/bioinfo/guchaoyang/software/cutadapt-4.1/bin/cutadapt";
my $samtools="/datapool/software/anaconda3/envs/qiime2/bin/samtools";
#my $splits = '\n\n';
open (OUT, ">$outdir/SOAPfilter.sh")||die;
if ($samplelist && $hosttype eq "human"){
	open (IN,$samplelist)||die;
	open (LIST,">$outdir/fq.list") ||die;
	while (<IN>){
	chomp;
	my $sample=(split/\t/,$_)[0];	
	my $dir="$outdir/$sample";
	(-d $dir) || `mkdir $dir`;
        my $sample_fq1=(split/\,/,((split/\t/,$_)[-1]))[0];
	my $sample_fq2=(split/\,/,((split/\t/,$_)[-1]))[1];
	my $cmd.="cd $dir\nperl $split_barcode $barcodelist $sample_fq1 $sample_fq2 split_reads  > cutoff_barcode.log\ncp /datapool/bioinfo/guchaoyang/pipeline/stLFR/bin/lane.lst $dir\n$SOAPfilter -t 20  -y  -F $F_adapter  -R $R_adapter  -p -M 2 -f -1 -Q 10 lane.lst stat_SOAP.txt > SOAP.log\nmv split_reads.1.fq.gz.clean.gz split_reads.1.clean.fq.gz\nmv split_reads.2.fq.gz.clean.gz split_reads.2.clean.fq.gz\n$bowtie2  --very-sensitive  -p 20 -x \/datapool\/db\/hg38\/hg38  -1 split_reads.1.clean.fq.gz  -2 split_reads.2.clean.fq.gz  --al-conc-gz  $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -@ 10 -b -o $sample.bam $sample.sam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz\n####rm lane.lst\n";
 	print OUT $cmd."\n";
	print LIST "$sample\t$dir/$sample.unmap.1.fq.gz\,$dir/$sample.unmap.2.fq.gz\n";
	}

}

if ($samplelist && $hosttype eq "mice"){
        open (IN,$samplelist)||die;
 	open (LIST,">$outdir/fq.list") ||die;
	while (<IN>){
        chomp;
	my $sample=(split/\t/,$_)[0];
        my $dir="$outdir/$sample";
	(-d $dir) || `mkdir $dir`;
	my $sample_fq1=(split/\,/,((split/\t/,$_)[-1]))[0];
        my $sample_fq2=(split/\,/,((split/\t/,$_)[-1]))[1];
        my $cmd.="cd $dir\nperl $split_barcode $barcodelist $sample_fq1 $sample_fq2 split_reads > cutoff_barcode.log\ncp /datapool/bioinfo/guchaoyang/pipeline/stLFR/bin/lane.lst $dir\n$SOAPfilter -t 20  -y -F $F_adapter  -R $R_adapter  -p -M 2 -f -1 -Q 10 lane.lst stat_SOAP.txt > SOAP.log\nmv split_reads.1.fq.gz.clean.gz split_reads.1.clean.fq.gz\nmv split_reads.2.fq.gz.clean.gz split_reads.2.clean.fq.gz\n$bowtie2 --very-sensitive  -p 20  -x \/datapool\/db\/mm39\/mm39  -1  split_reads.1.clean.fq.gz -2 split_reads.2.clean.fq.gz  --al-conc-gz  $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -@ 10 -b -o $sample.bam $sample.sam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz\n####rm lane.lst\n";
          
	print OUT $cmd."\n";
	print LIST "$sample\t$dir/$sample.unmap.1.fq.gz\,$dir/$sample.unmap.2.fq.gz\n";
        }
#	 my $num=`less SOAPfilter.sh  |grep \"\^$\" |wc -l`;
#         chomp($num);

}

if ($samplelist && $hosttype eq "none"){
        open (IN,$samplelist)||die;
        open (LIST,">$outdir/fq.list") ||die;
	while (<IN>){
        chomp;
 	my $sample=(split/\t/,$_)[0];
	my $dir="$outdir/$sample";
        (-d $dir) || `mkdir $dir`;
	my $sample_fq1=(split/\,/,((split/\t/,$_)[-1]))[0];
        my $sample_fq2=(split/\,/,((split/\t/,$_)[-1]))[1];
        my $cmd.="cd $dir\nperl $split_barcode $barcodelist $sample_fq1 $sample_fq2 split_reads  > cutoff_barcode.log\ncp /datapool/bioinfo/guchaoyang/pipeline/stLFR/bin/lane.lst $dir\n$SOAPfilter -t 20  -y -F $F_adapter  -R $R_adapter  -p -M 2 -f -1 -Q 10 lane.lst stat_SOAP.txt > SOAP.log\nmv split_reads.1.fq.gz.clean.gz $sample.unmap.1.fq.gz\nmv split_reads.2.fq.gz.clean.gz $sample.unmap.2.fq.gz\n####rm lane.lst\n";
                
	print OUT $cmd."\n";
	print LIST "$sample\t$dir/$sample.unmap.1.fq.gz\,$dir/$sample.unmap.2.fq.gz\n";
        }
}

	close IN;
	close OUT;
	close LIST;
	`perl $Bin/dsub_batch.pl -thd 10 -mem 4 $outdir/SOAPfilter.sh`;

