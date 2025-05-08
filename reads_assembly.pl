#!/usr/bin/perl 
=head1 Description
	sequence assembly
=head1 Usage
	perl reads_assembly.pl [options]
	general input and output:
	-fl  <a file>       list of fastq from 01.QC(file suffix can be ".fq.gz"), one strain per line, PE read seperated by ",", and different by "\n"
	-type               assembly_type, default 'metaspades_athena' [metaspades and athene-meta]
	-o   <a directory>  output dir, default current directory [./02.Assembly]
	-thd <num>          thread for dsub
	-h                  show help

=head1 Example
	perl reads_assembly.pl  -fl fq.list -type  -o 

=head1 Version
        Author: guchaoyang0826@163.com
        Date: 2022-08-18
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw/abs_path/;

my ($fqlist,$assemblytype,$outdir,$thd,$help);

GetOptions (
	"fl:s"=>\$fqlist,
	"type:s"=>\$assemblytype,
	"o:s"=>\$outdir,
	"h"=>\$help
);
die `pod2text $0` if ($help || !$fqlist);

unless ($outdir) {
	$outdir="./02.Assembly";
}

unless (-e $outdir) {
	mkdir $outdir;
}

unless ($thd){
	$thd=8;
}

$outdir=abs_path($outdir);

#============================================================================================================================================================

if (!$fqlist){
	print "No input read files\n";
	exit;
}

if (!$assemblytype){
        print "No input read type: 'human' or 'mice' or 'none'\n";
        exit;
}

my $cutadapt="/datapool/bioinfo/guchaoyang/software/cutadapt-4.1/bin/cutadapt";
my $fastq_pair="/datapool/bioinfo/guchaoyang/software/fastq-pair-1.0/build/fastq_pair";
my $metaspades="/datapool/software/anaconda3/bin/metaspades.py";
my $bwa="/datapool/bioinfo/guchaoyang/software/bwa-0.7.17/bwa";
my $samtools="/datapool/bioinfo/guchaoyang/software/samtools-1.10/samtools";
my $athena="/datapool/software/anaconda3/envs/qiime1/bin/athena-meta";

#my $splits = '\n\n';

open (OUT, ">$outdir/assembly.sh")||die;
if ($fqlist && $assemblytype eq "metaspades_athena"){
	
	open (IN,$fqlist)||die;
	while(<IN>){
		chomp;
		my $sample=(split/\t/,$_)[0];
		my $dir="$outdir/$sample";
        	(-d $dir) || `mkdir $dir`;
		(-d "$dir/athena") || `mkdir "$dir/athena"`;
		my $sample_fq1=(split/\,/,((split/\t/,$_)[-1]))[0];
		my $sample_fq2=(split/\,/,((split/\t/,$_)[-1]))[1];
		open (Config,">$dir/athena/$sample.config") || die;
		print Config "\{\n    \"input_fqs\"\: \"$dir/$sample.interleaved.clean.fq\"\,\n    \"ctgfasta_path\"\: \"$dir/metaspades/contigs.fasta\"\,\n    \"reads_ctg_bam_path\"\: \"$dir/metaspades/align-reads.metaspades-contigs.bam\"\n\}\n";
		my $cmd.="cd $dir\nzcat $sample_fq1 |sed \'N\;N\;N \;s\/\\n\/\\t\_\|\_\/g\'|grep -v \'0\_0\_0\' > $sample.unmap.clean.tmp.1.fq\nsort -k 2 $sample.unmap.clean.tmp.1.fq |sed \'s\/\\t\\_\|\\_\/\\n\/g\' >$sample.unmap.clean.1.fq\nzcat $sample_fq2 |sed \'N\;N\;N \;s\/\\n\/\\t\_\|\_\/g\'|grep -v \'0\_0\_0\' > $sample.unmap.clean.tmp.2.fq\nsort -k 2 $sample.unmap.clean.tmp.2.fq |sed \'s\/\\t\\_\|\\_\/\\n\/g\' >$sample.unmap.clean.2.fq\n$fastq_pair $sample.unmap.clean.1.fq $sample.unmap.clean.2.fq\n$cutadapt  -o $sample.interleaved.clean.fq  --interleaved $sample.unmap.clean.1.fq.paired.fq $sample.unmap.clean.2.fq.paired.fq\n$metaspades -1 $sample.unmap.clean.1.fq.paired.fq -2 $sample.unmap.clean.2.fq.paired.fq  -o $dir/metaspades\ncd $dir/metaspades\n$bwa  index  contigs.fasta\n$bwa  mem -C  contigs.fasta $dir/$sample.unmap.clean.1.fq.paired.fq  $dir/$sample.unmap.clean.2.fq.paired.fq \| $samtools  sort -o align-reads.metaspades-contigs.bam \-\n$samtools index align-reads.metaspades-contigs.bam\n";

#print OUT_config "\{\n    \"input_fqs\"\: \"$dir/$sample.interleaved.clean.fq\"\,\n\"ctgfasta_path\"\: \"$dir/contigs.fasta\"\,\n\"reads_ctg_bam_path\"\: \"$dir/align-reads.metaspades-contigs.bam\"\n\}\n";
 		print OUT $cmd."\#\$ \-S /bin/bash\nconda activate qiime1\nexport PATH=\"/datapool/bioinfo/guchaoyang/software/samtools-1.10/bin/\:\$PATH\"\ncd $dir/athena\n$athena --config $dir/athena/$sample.config  --threads 128 > athena.log\n#rm $sample.unmap.clean.tmp.1.fq\n#rm $sample.unmap.clean.1.fq\n#rm $sample.unmap.clean.tmp.2.fq\n####rm $sample.unmap.clean.2.fq\n\n";
	}

}

if ($fqlist && $assemblytype eq "metaspades"){
        open (IN,$fqlist)||die;
        while (<IN>){
	        chomp;
		my $sample=(split/\t/,$_)[0];
	        my $dir="$outdir/$sample";
		(-d $dir) || `mkdir $dir`;
		my $sample_fq1=(split/\,/,((split/\t/,$_)[-1]))[0];
        	my $sample_fq2=(split/\,/,((split/\t/,$_)[-1]))[1];
        	my $cmd.="cd $dir\nzcat $sample_fq1 |sed \'N\;N\;N \;s\/\\n\/\\t\_\|\_\/g\'|grep -v \'0\_0\_0\' > $sample.unmap.clean.tmp.1.fq\nsort -k 2 $sample.unmap.clean.tmp.1.fq |sed \'s\/\\t\\_\|\\_\/\\n\/g\' >$sample.unmap.clean.1.fq\nzcat $sample_fq2 |sed \'N\;N\;N \;s\/\\n\/\\t\_\|\_\/g\'|grep -v \'0\_0\_0\' > $sample.unmap.clean.tmp.2.fq\nsort -k 2 $sample.unmap.clean.tmp.2.fq |sed \'s\/\\t\\_\|\\_\/\\n\/g\' >$sample.unmap.clean.2.fq\n$fastq_pair $sample.unmap.clean.1.fq $sample.unmap.clean.2.fq\n$cutadapt  -o $sample.interleaved.clean.fq  --interleaved $sample.unmap.clean.1.fq.paired.fq $sample.unmap.clean.2.fq.paired.fq\n$metaspades   -1 $sample.unmap.clean.1.fq.paired.fq -2 $sample.unmap.clean.2.fq.paired.fq    -o $dir/metaspades\ncd $dir/metaspades\n$bwa  index  contigs.fasta\n$bwa  mem -C  contigs.fasta $dir/$sample.unmap.clean.1.fq.paired.fq  $dir/$sample.unmap.clean.2.fq.paired.fq \| $samtools  sort -o align-reads.metaspades-contigs.bam \-\n$samtools index align-reads.metaspades-contigs.bam\n#rm $sample.unmap.clean.tmp.1.fq\n#rm $sample.unmap.clean.1.fq\n#rm $sample.unmap.clean.tmp.2.fq\n####rm $sample.unmap.clean.2.fq\n";
		print OUT $cmd."\n";
        }
#	 my $num=`less SOAPfilter.sh  |grep \"\^$\" |wc -l`;
#         chomp($num);

}


	close IN;
	close OUT;
	close Config;
	`perl $Bin/dsub_batch2.pl -thd 10 -mem 4 $outdir/assembly.sh`;

