cd `pwd`	
data_path=$1

for file in `cat $data_path`;
do
name=`basename $file`
if [ -e $name ];
then 
echo "$name exists!"
else
wget $file
fi
done

