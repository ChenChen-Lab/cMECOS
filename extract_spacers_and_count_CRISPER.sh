#!/bin/env bash
#AUTHER:LUCAS
#DATE:20230328
#GOAL:获得输入基因组中的spacers序列，并输出包含的满足层级CRISPER的数目
#USAGE:script 输入序列 输出目录名 最小evidence层级 线程数
##注意：使用之前必须激活CRISPRCasFinder的conda环境！！！

set -e
set -u
IN=$1
OUT=$2
level=$3
THRE=$4

CRI_FIND='/disk_RAID5_27T/software/CRISPRCasFinder-master/CRISPRCasFinder.pl'
SO_FILE='/disk_RAID5_27T/software/CRISPRCasFinder-master/src/sel392v2.so'
SEQK='/home/train/anaconda3/bin/seqkit'

#获得输入文件的绝对路径
ABS_IN=`readlink -f ${IN}`

#由于CRISPRCasFinder的原因，会在运行的地方生成影响下一个CRISPRCasFinder执行的中间文件。因此先将执行的输入fasta文件在新的目录中创建软连接

mkdir -p ${OUT} &> /dev/null
cd ${OUT} 
ln -s ${ABS_IN}
#使用seqkit处理“.”
IN2=`basename ${ABS_IN}`
cat ${IN2}|~/anaconda3/bin/seqkit replace -p '\.' -r '' > REPLACE_${IN2}

perl ${CRI_FIND} -in REPLACE_${IN2} -so ${SO_FILE} -lMin ${level} -keep -cpuM ${THRE} -out CRISPRCasFinder_result
#RAN1
RAN1=`head -n 20 /dev/urandom | cksum |cut -f1 -d " "`

cat CRISPRCasFinder_result/TSV/Crisprs_REPORT.tsv|sed '1d'| awk -F '\t' '$27 >= 1 {print $2}' > tmp_${RAN1}
#CRISPER的数目
NUM_CRISP=`cat tmp_${RAN1}|wc -l`

for i in `cat tmp_${RAN1}`;do cat CRISPRCasFinder_result/CRISPRFinderProperties/${i}_properties/Spacers/*|${SEQK} rmdup -s|${SEQK} replace -p '^' -r "${i}:_:";done > CRISPRCasFinder_result/All_spacers.fa
trap 'rm -fr tmp_${RAN1}' ERR EXIT INT TERM

echo $NUM_CRISP > CRISPRCasFinder_result/CRISPR_NUM
