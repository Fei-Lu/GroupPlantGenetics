#!/bin/bash
#允许的进程数
THREAD_NUM=16
#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#预先写入指定数量的换行符，一个换行符代表一个进程
for ((i=0;i<$THREAD_NUM;i++))
do
    echo -ne "\n" 1>&9
done
if [ $# != 1 ] ;then
        echo "The parameters you enter is not correct !";
        exit -1;
fi

while read line
do
{
    #进程控制
    read -u 9
    {
        #
        prefix=${line/.all_snpgene_pairs.txt.gz/}
        prefix1=${prefix/"/data2/xiaohan/metasoft/v6_allpairs/GTEx_Analysis_v6p_all-associations/"/}
        python /data2/xiaohan/metasoft/shell/extract_pairs.py ${line} /data2/xiaohan/metasoft/meta-sig.combined_signifpairs.txt.gz $prefix1 -o /data2/xiaohan/metasoft/v6_extract/
        echo -ne "\n" 1>&9
    }&
}
done < $1
wait
echo "执行结束"
rm tmp
