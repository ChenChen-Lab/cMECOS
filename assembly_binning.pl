#!/usr/bin/perl 
=head1 Description
	assembly binning
=head1 Usage
	perl assembly_binning.pl [options]
	general input and output:
	-assembly_list	    list of athena assembly fasta
	-o   <a directory>  output dir, default current directory [./02.Assembly]
	-thd <num>          thread for dsub
	-h                  show help
	-notrun             only write the shell, but not run
=head1 Example
	perl assembly_binning.pl  -fl fq.list  -o 

=head1 Version
        Author: guchaoyang0826@163.com
        Date: 2023-03-23
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw/abs_path/;

my ($assembly_list,$outdir,$thd,$help,$notrun);

GetOptions (
	"assembly_list:s"=>\$assembly_list,
	"o:s"=>\$outdir,
	"h"=>\$help,
	"notrun"=>\$notrun
);
die `pod2text $0` if ($help || !$assembly_list);

unless ($outdir) {
	$outdir="./03.Binning_new";
}

unless (-e $outdir) {
	mkdir $outdir;
}

unless ($thd){
	$thd=8;
}

$outdir=abs_path($outdir);

#============================================================================================================================================================

if (!$assembly_list){
	print "No input read files\n";
	exit;
}

my $bowtie2_build="/datapool/software/anaconda3/bin/bowtie2-build";
my $bowtie2="/datapool/software/anaconda3/bin/bowtie2";
my $samtools="/datapool/software/anaconda3/bin/samtools";
my $jgi_summarize_bam_contig_depths="/datapool/software/anaconda3/bin/jgi_summarize_bam_contig_depths";
my $metabat2="/datapool/software/anaconda3/bin/metabat2";
my $checkm="/datapool/software/anaconda3/envs/qiime2/bin/checkm";
my $checkm_data="/datapool/software/anaconda3/envs/qiime2/checkm_data";
my $gunc="/datapool/software/anaconda3/envs/qiime2/bin/gunc";
my $gunc_db="/datapool/bioinfo/guchaoyang/test/athena/binning/gunc/db/progenomes/gunc_db_progenomes2.1.dmnd";
#my $splits = '\n\n';

open (OUT, ">$outdir/assembly_binning.sh")||die;
if ($assembly_list){
	open (IN,$assembly_list)||die;
	while(<IN>){
		chomp;
		my $sample=(split/\t/,$_)[0];
		my $dir="$outdir/$sample";
        	(-d $dir) || `mkdir  $dir`;
		(-d "$dir/metabat2/bin_input") || `mkdir -p  "$dir/metabat2/bin_input"`;
		(-d "$dir/checkm") || `mkdir "$dir/checkm"`;
		(-d "$dir/gunc") || `mkdir "$dir/gunc"`;
		(-d "$dir/gunc/output") || `mkdir -p "$dir/gunc/output"`;
		my $sample_fq1=(split/\,/,((split/\t/,$_)[1]))[0];
		my $sample_fq2=(split/\,/,((split/\t/,$_)[1]))[1];
		my $sample_athena=(split/\t/,$_)[-1];
		my $cmd="export PATH=\"/datapool/software/anaconda3/envs/qiime2/bin/\:\$PATH\"\ncd $dir/metabat2\nln -s $sample_athena  $sample.fasta\n$bowtie2_build  -f $sample.fasta $sample.fasta\n$bowtie2 -x $sample.fasta  -1 $sample_fq1  -2 $sample_fq2  -p 20 -S $sample.sam 2> bowtie2.log\n$samtools view -F 4 -Sb $sample.sam  > $sample.bam\n$samtools  sort $sample.bam  -o $sample.sort.bam\n$samtools index  $sample.sort.bam\n$jgi_summarize_bam_contig_depths  --outputDepth $sample.depth.txt $sample.sort.bam\n$metabat2 -i $sample.fasta -a $sample.depth.txt -o $sample --sensitive -t 20 -v > $sample.log.txt\ncd $dir/metabat2/bin_input\nln -s $dir/metabat2/$sample*.fa . \n";
		$cmd.="cd $dir/checkm\n$checkm data setRoot $checkm_data\n$checkm lineage_wf $dir/metabat2/bin_input $dir/checkm/output -x fa -t 20  --nt --tab_table  -f bins_qa.txt\n$checkm qa -f qa.tsv -o 2 --tab_table  output/lineage.ms  $dir/checkm/output\nless qa.tsv |awk -F \"\\t\" '{print \$1\"\\t\"\$6\"\\t\"\$7\"\\t\"\$9\"\\t\"\$11\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14\"\\t\"\$15\"\\t\"\$16\"\\t\"\$17\"\\t\"\$18}' |sed '1d' |awk '\$2>90&&\$3<5' >$sample.checkm.HG.list\nless qa.tsv |awk -F \"\\t\" '{print \$1\"\\t\"\$6\"\\t\"\$7\"\\t\"\$9\"\\t\"\$11\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14\"\\t\"\$15\"\\t\"\$16\"\\t\"\$17\"\\t\"\$18}' |sed '1d' |awk '\$2>50&&\$3<10' >$sample.checkm.MG.list\n";
		$cmd.="cd $dir/gunc\n$gunc run  -d $dir/metabat2/bin_input/ -o  $dir/gunc/output -r $gunc_db  -e .fa -t 20\n$gunc merge_checkm -g output/GUNC.progenomes_2.1.maxCSS_level.tsv -c $dir/checkm/qa.tsv -o $dir/gunc\nless GUNC_checkM.merged.tsv |awk -F \"\\t\" '{print \$1\"\\t\"\$6\"\\t\"\$15\"\\t\"\$16\"\\t\"\$20}' |awk '\$2<5&&\$3>90&&\$4<5&&\$5==\"True\"' >$sample.HQ.list\n####rm\n";

 		print OUT $cmd."\n";
		
	}

}

	close IN;
	close OUT;
	`perl $Bin/dsub_batch2.pl -thd 10 -mem 4 $outdir/assembly_binning.sh`;

$notrun && exit;
