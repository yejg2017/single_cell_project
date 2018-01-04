save_path=$1
printf "%-10s %-6s %-6s\n" File Rows Cols 
for gz in `ls *.gz`;
do
rows=`zcat $gz | wc -l`
cols=`zcat $gz | awk '{ print NF}'| awk 'NR==2{ print}' `
printf "%-6s %-6s\n" $gz $rows $cols >> $save_path
done

