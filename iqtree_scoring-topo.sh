for i in $(ls ${1})
do
phylo=$(echo $i | cut -d "." -f 1 | cut -d "_" -f 3)
iqtree -s ${1}/${i} -m LG+F -z ../topologies/${phylo}.tre -pre ${1}/${i} -n 0
rm -f ${1}/${i}.log
rm -f ${1}/${i}.ckp.gz
rm -f ${1}/${i}.trees
rm -f ${1}/${i}.treefile
done
