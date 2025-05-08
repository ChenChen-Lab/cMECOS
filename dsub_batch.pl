#!/usr/bin/perl


=head1 Description
	parallelly run by dsub
=head1 Usage
	perl dsub_batch.pl [options] <shell file>
		--mem maximum memory to use [Gb], default 2
		--thd maximum thread to run, default 2 in each node;
		--split default "\n"
		--ppn CPU cores used, default 1
		--sub default: 1==dsub; 2: multiple threads at current node and "thd" is required
		--tim maximum run time [hours], default 10000 hours
		--env the environment of software, current default base
		--h show help
=head1 Example

=head1 Version
	Author: Du Pengcheng dupengcheng@icdc.cn
	Date: 2015-12-15

=cut

use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ($mem, $th_num,$split, $ppn, $run_time, $sub, $env, $help);
GetOptions (
	"mem:s"=>\$mem,
	"thd:s"=>\$th_num,
	"split:s"=>\$split,
	"ppn:s"=>\$ppn,
	"tim:s"=>\$run_time,
	"sub:s"=>\$sub,
	"env:s"=>\$env,
	"h"=>\$help,
);
die `pod2text $0` if ($help || @ARGV==0);

my $sh = $ARGV[0];
$sh = abs_path($sh);
my $cur_dir = abs_path("./");

$mem = 2 unless ($mem);
$ppn = 1 unless ($ppn);
$sub = 1 unless ($sub);
$run_time = 10000 unless ($run_time);
$env = "base" unless ($env);
$split = "\n" unless ($split);

my $time = sub_format_datetime(localtime(time()));
my $sh_dir = $cur_dir."/".basename($sh)."\.".$time;
mkdir $sh_dir;

my $job_num = `wc -l $sh` || die "Input shell does not exist!\n";
$job_num = (split /\s/, $job_num)[0];

my $sample_N = `grep \"####\" $sh \|wc -l ` || die "Sample Number does not exist!\n";
$sample_N = (split /\s/, $sample_N)[0];

my $per_N;
$per_N = int($job_num/$sample_N);

##if ($th_num) {
##	$per_N = $job_num%$th_num==0 ? $job_num/$th_num : int($job_num/$th_num)+1;
##}
##else {
##	my $node_num = `cluster_free | grep -c batch`;
##	chomp $node_num;
##	$per_N = int($job_num/($node_num*2))+1;
##}


my ($i, $j) = (1, 0);
my %shell;
open (IN, $sh) || die;
while (<IN>) {
	my $sub_sh;

	if ($i == 1 or $i > $per_N) {
		$i=1;
		$j++;
		$sub_sh = sprintf "%03.0f", $j;
		$sub_sh = $sh_dir."/work_".$sub_sh.".sh";
		$shell{basename($sub_sh)} = 1;

		if ($sub==1) {
			open (OUT, ">$sub_sh") || die;
		#	print OUT "#!/bin/sh\n";
		#	print OUT "#PBS -q batch\n";
		#	print OUT "#PBS -l mem=".$mem."gb,nodes=1:ppn=".$ppn.",walltime=".$run_time.":00:00\n";
		#	print OUT "#HSCHED -s hschedd\n";
		#	print OUT "cd $cur_dir\n";
			#print OUT "source /home/dupc/.bashrc\n";
			#print OUT 'export PATH=/datapool/software/anaconda3/bin:$PATH';
			#my $cmd_1="#"."\$"." -S /bin/bash"."\n"."conda activate R4.1"."\n";
		#	print OUT "#"."\$"." -S /bin/bash"."\n";
		###	print OUT "conda activate $env";
		#	print OUT "\n";
		}
	}
	else {
		$sub_sh = sprintf "%03.0f", $j;
		$sub_sh = $sh_dir."/work_".$sub_sh.".sh";
	}

	if (/\s\.\.\//) {
		s/\s\.\.\// $cur_dir\/\.\.\//g;
	}
	
	if (/\s\.\//) {
		s/\s\.\// $cur_dir\//g;
	}
	
	open (OUT, ">>$sub_sh") || die;
	print OUT $_;
	$i++;
}
close IN;

if ($sub==1) {
	my %jobs;
	chdir $sh_dir;
	foreach my $h (sort keys %shell) {
		my $tmp = `qsub -cwd $h`;
		$jobs{$1}=1 if ($tmp=~/(\d+)\.manager/);
	}
	chdir $cur_dir;

	OUT: while (1) {
		my $qstat = `sleep 10; qstat`;
		my %run;
		foreach my $line (split /\n/, $qstat) {
			if ($line=~/(\d+)\.manager/ and exists $jobs{$1}) {
				$run{$1}=1;
			}
		}
		
		if (0+keys %run == 0) {
			print "\nThis work is completed!\n$sh\n";
			last OUT;
		}
	}
}
else {
	chdir $sh_dir;
	open (OUT, ">tmp.sh") || die;
	foreach my $h (sort keys %shell) {
		print OUT "sh $h &\n";
	}
	close OUT;
	`sh tmp.sh`;
	chdir $cur_dir;
	print "\nThis work is completed!\n$sh\n";
}

sub sub_format_datetime
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d.%02d.%02d.%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
