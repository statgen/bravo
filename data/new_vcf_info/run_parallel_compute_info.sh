#!/bin/bash

compute_script=compute_info.py
merge_script=merge_gzvcfs_by_pos.py

vcf=$1
samples=$2
chunk=$3
output=$4

starts=()
ends=()

chromosomes=(`tabix --list-chroms ${vcf}`)

if [ ${#chromosomes[@]} -ne 1 ]
then
   echo "Multiple chromosomes detected. Only one is allowed."
   exit
fi

chromosome=${chromosomes[0]}

echo "Chromosome ${chromosome} detected."

start_position=`gzip -dc ${vcf} | grep -v "^#" -m1 | cut -f2`
end_position=$((start_position + chunk))

more=`tabix ${vcf} ${chromosome}:${star_position} | head -n1 | wc -l` 

while (( ${more} > 0 ))
do
   starts+=(${start_position})
   ends+=(${end_position})
   start_position=$((end_position + 1))
   end_position=$((start_position + chunk))
   more=`tabix ${vcf} ${chromosome}:${start_position} | head -n1 | wc -l`
done

touch ${output}.list.temp
for i in ${!starts[@]}
do
   python ${compute_script} --in-gzvcf ${vcf} --in-samples ${samples} --chrom ${chromosome} --start ${starts[${i}]} --end ${ends[${i}]} --out-gzvcf ${output}.${i}.temp.gz &
   pid=$! 
   echo "Computing region ${i} (${chromosome}:${starts[${i}]}-${ends[${i}]}) PID = ${pid}..."
   echo ${output}.${i}.temp.gz >> ${output}.list.temp
done

original_lc=`pigz -dc ${vcf} | grep -v "^#" | wc -l`

wait

echo "All computations finished."
echo "Merging regions..."

python ${merge_script} --in ${output}.list.temp --out ${output}

merged_lc=`pigz -dc ${output} | grep -v "^#" | wc -l`

if [ $original_lc -ne $merged_lc ]
then
   echo "Region merging failed."
   exit
fi

rm ${output}.*.temp.gz ${output}.list.temp

echo "Region merging finished."
