#!/bin/env python3
#AUTHER:LUCAS
#DATE:20230702
#GOAL:根据提供的MAG列表和首次聚类结果使用MASH进行二次聚类，并输出聚类结果；可以用来确定株水平聚类(99%)
#USAGE:SCRIPT 自己MAG文件名文件  dRep的Cdb文件 存放包括自己MAG文件和UHGG所有MAG文件的目录 UHGG的metadata 输出目录 阈值(通常0.99为株，0.95为种) 线程数
#warning: 注意所有除了自己的其他的序列结尾应该为.fna结尾在dREP的Cdb文件中。
#        2. 存放UHGG的目录中UHGG的基因组结尾为.fa
import sys
import pandas as pd
import warnings
import os
import time
from multiprocessing.dummy import Pool as ThreadPool
warnings.filterwarnings("ignore")
argin=[i for i in sys.argv]

#读取自己的MAG文件（每一行为MAG的文件名，是Cdb中的自己的MAG名）
My_MAG_path=argin[1]
with open(My_MAG_path,"r") as f:
    My_MAG_path_list=f.read().splitlines()

Cdb_path=argin[2]

ALL_MAG_DIR=argin[3]
UHGG_metadata_path=argin[4]
OUT_DIR_main=argin[5]
cut_off=argin[6]
threads=argin[7]
work_path=str("/".join(sys.argv[0].split("/")[:-1])).strip()

print("{0}\nStart  accordding to dRep result from |{1}|".format(time.strftime('%Y.%m.%d.%H:%M:%S',time.localtime(time.time())),Cdb_path))
print("Your own MAG file is |{}|".format(My_MAG_path))
print("UHGG metadata is |{}|".format(UHGG_metadata_path))
print("All MAG fasta file is |{}|".format(ALL_MAG_DIR))
print("The cluster cut-off is |{}|".format(cut_off))
print("The threads used is |{}|".format(threads))

##创建输出目录
os.environ['OUT_DIR_main']=str(OUT_DIR_main)
os.environ['ALL_MAG_DIR']=str(ALL_MAG_DIR)
os.system('mkdir -p ${OUT_DIR_main} > /dev/null')

#step1，取出含有自己MAG的dREP簇
Cdb_raw=pd.read_csv(Cdb_path)
Cdb_my_secondary_cluster=(Cdb_raw[Cdb_raw.genome.isin(My_MAG_path_list)])["secondary_cluster"].tolist()  #仅包含自己MAG的dREP簇的列表
#step2，判断每一个簇是否包含UHGG的MAG，若包含则对每一个簇其执行株水平聚类，若不包含则将株水平分类全部赋予0
UHGG_metadata=pd.read_csv(UHGG_metadata_path,sep="\t") #读取UHGG metadata
#输出输入的所有MAG包含的文件名
os.system('ls ${ALL_MAG_DIR} > ${OUT_DIR_main}/ALL_MAG_list.txt')
UHGG_IN_SP_ID=[ z.split(".fna")[0].split(".fa")[0] for z in Cdb_raw[Cdb_raw.secondary_cluster.isin(Cdb_my_secondary_cluster)].genome.tolist()]
UHGG_IN_MAG=(UHGG_metadata[UHGG_metadata.Species_rep.isin(UHGG_IN_SP_ID)]).Genome.tolist()
#需要的所有
IN_genome_list_all=set([z+".fa" for z in UHGG_IN_MAG+UHGG_IN_SP_ID] )
#输出IN_genome_list_raw.txt中为需要的所有基因组
with open("{0}/Need_genome_list.txt".format(OUT_DIR_main),"w") as f:
    for x in IN_genome_list_all:
        f.write(x+"\n")

with open("{}/ALL_MAG_list.txt".format(OUT_DIR_main)) as f:
    INPUT_genome=f.read().splitlines()

if len( set(IN_genome_list_all) - set(INPUT_genome)) > 0:
    print("Error! Input all Genomes can not include needed query genomes in {0}! Exit.\nPlease check the {1}".format(ALL_MAG_DIR,"Need_genome_list.txt"))
    sys.exit(1)
para_go_list_1=[]
para_go_list_2=[]
for i in set(Cdb_my_secondary_cluster):
    Cdb_i=Cdb_raw[Cdb_raw.secondary_cluster==i]
    #判断是否该簇的所有集合均为自己的MAG
    if Cdb_i.shape[0]==sum(Cdb_i.genome.isin(My_MAG_path_list)):
        i_pd=pd.DataFrame()
        Cluster_sklearn=0
        i_pd["ID"]=Cdb_i.genome
        i_pd["Cluster_dRep"]=i
        i_pd["Cluster_sklearn"]=Cluster_sklearn
        os.environ['OUT_DIR']=str(OUT_DIR_main)+"/"+str(i)
        os.system('mkdir -p ${OUT_DIR} > /dev/null')
        i_pd.to_csv(str(OUT_DIR_main)+"/"+str(i)+"/novel_cluster.tsv",sep="\t",index=None,header=None)
    else:
        #赋予IN_path为该文件的名字，并赋予shell该变量
        ##设置linux中的变量
        os.environ['OUT_DIR']=str(OUT_DIR_main)+"/"+str(i)
        os.system('mkdir -p ${OUT_DIR} > /dev/null')
        #os.environ["sec_threshold"]=str(cut_off)
        #根据匹配的种水平ID获得该种中的所有基因组
        UHGG_genome=(UHGG_metadata[UHGG_metadata.Species_rep.isin([ z.split(".fna")[0] for z in Cdb_i.genome])]).Genome.tolist()
        IN_genome_list=set([f+".fa" for f in UHGG_genome + [z.split(".fna")[0].split(".fa")[0] for z in Cdb_i.genome.tolist()]])
        with open("{0}/{1}/IN_genome_list".format(OUT_DIR_main,i),"w") as f:
            for x in IN_genome_list:
                f.write(x+"\n")
        IN_MAG_file="{0}/{1}/IN_genome_list".format(OUT_DIR_main,i)
        #os.environ['IN_MAG_file']=str(IN_MAG_file)
        #os.environ['Cluster_dRep']=str(i)
        #os.environ['work_path']=str(work_path)
        #执行MASH的操作生成上三角矩阵
        para_go_list_1.append("{0}/MASH_sketch_and_triangle.sh {1} {2} 1 {3} &> {3}/MASH_sketch_triangle.log".format(work_path,IN_MAG_file,ALL_MAG_DIR,str(OUT_DIR_main)+"/"+str(i)))
        para_go_list_2.append("{0}/obtain_cluster_from_sklearn.py {1}/up_triangle.txt {2} {3} {1} &> {1}/sklearn_cluster.log ".format(work_path,str(OUT_DIR_main)+"/"+str(i),cut_off,i))
        #os.system('${work_path}/MASH_sketch_and_triangle.sh ${IN_MAG_file} ${ALL_MAG_DIR} 1 ${OUT_DIR}')
        #os.system('${work_path}/obtain_cluster_from_sklearn.py ${OUT_DIR}/up_triangle.txt ${sec_threshold} ${Cluster_dRep} ${OUT_DIR} &> ${OUT_DIR}/sklearn_cluster.log')
def para_go(i):
    os.system(i)

with ThreadPool(int(threads)) as p: # 指定n个线程
    p.map(para_go, para_go_list_1)
with ThreadPool(int(threads)) as p: # 指定n个线程
    p.map(para_go,para_go_list_2)
##合并株水平分类结果
os.system('find ${OUT_DIR_main} -name "*.tsv"|xargs cat|sed "1i Genome\tdRep_cluster\tMash_cluster" > ${OUT_DIR_main}/Result_MASH_secend_cluster.txt')

print("{0}\nCluster end!\nThe result is stored in |{1}| named |Result_MASH_secend_cluster.txt|".format(time.strftime('%Y.%m.%d.%H:%M:%S',time.localtime(time.time())),OUT_DIR_main))
