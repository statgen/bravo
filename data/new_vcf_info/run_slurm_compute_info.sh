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

echo "#!/bin/bash" > ${output}.compute.sbatch
echo "#SBATCH --array=0-$(( ${#starts[@]} - 1 ))" >> ${output}.compute.sbatch
echo "#SBATCH --job-name=compute_${output}" >> ${output}.compute.sbatch
echo "#SBATCH --output=${output}-compute-sbatch-%a.log" >> ${output}.compute.sbatch
echo "declare -a jobs" >> ${output}.compute.sbatch
for i in ${!starts[@]}
do
    echo -e "jobs[${i}]=\"python ${compute_script} --in-gzvcf ${vcf} --in-samples ${samples} --chrom ${chromosome} --start ${starts[${i}]} --end ${ends[${i}]} --out-gzvcf ${output}.${i}.temp.gz\"" >> ${output}.compute.sbatch
    echo ${output}.${i}.temp.gz >> ${output}.list.temp
done
echo -e "\${jobs[\${SLURM_ARRAY_TASK_ID}]}" >> ${output}.compute.sbatch

echo "SBATCH compute jobs file ${output}.compute.sbatch is ready."

echo "#!/bin/bash" > ${output}.merge.sbatch
echo "#SBATCH --job-name=merge_${output}" >> ${output}.merge.sbatch
echo "#SBATCH --output=${output}-merge-sbatch.log" >> ${output}.merge.sbatch
echo -e "original_lc=\`pigz -dc ${vcf} | grep -v \"^#\" | wc -l\`" >> ${output}.merge.sbatch
echo "python ${merge_script} --in ${output}.list.temp --out ${output}" >> ${output}.merge.sbatch
echo -e "merged_lc=\`pigz -dc ${output} | grep -v \"^#\" | wc -l\`" >> ${output}.merge.sbatch
echo -e "if [ \$original_lc -ne \$merged_lc ]; then echo \"Region merging failed.\"; exit 1; fi" >> ${output}.merge.sbatch
echo "rm ${output}.*.temp.gz ${output}.list.temp" >> ${output}.merge.sbatch

echo "SBATCH merge job file ${output}.merge.sbatch is ready."

echo "#!/bin/bash" > ${output}.run.sh
echo -e "message=\$(sbatch -p nomosix --mem=8000 --time=14-0 ${output}.compute.sbatch)" >> ${output}.run.sh
echo -e "if ! echo \${message} | grep -q \"[1-9][0-9]*\$\"; then echo \"Compute job(s) submission failed.\"; exit 1; fi" >> ${output}.run.sh
echo -e "job=\$(echo \${message} | grep -oh  \"[1-9][0-9]*\$\")" >> ${output}.run.sh
echo -e "echo \"$(( ${#starts[@]} - 1 )) compute job(s) submitted in \${job} job array.\"" >> ${output}.run.sh

echo -e "message=\$(sbatch -p nomosix --mem=8000 --time=14-0 --depend=afterok:\${job} ${output}.merge.sbatch)" >> ${output}.run.sh
echo -e "if ! echo \${message} | grep -q \"[1-9][0-9]*\$\"; then echo \"Merge job submission failed.\"; exit 1; fi" >> ${output}.run.sh
echo -e "job=\$(echo \${message} | grep -oh  \"[1-9][0-9]*\$\")" >> ${output}.run.sh
echo -e "echo \"Merge job \${job} submitted.\"" >> ${output}.run.sh

echo "Submission script ${output}.run.sh is ready. Execute it!"
