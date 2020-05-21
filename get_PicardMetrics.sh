######################################################################################################
# this script processes metrics from Picard (i.e. 11037_S7_Aligned.sortedByCoord.out.bam.metrics.txt)
######################################################################################################

if [ $# != 2 ] ; then
    echo -e "ERROR: 2 arguments are required: \
	(1) path to data folder [with bam.metrics.txt files] \
	(2) and path to folder for output files ... Exiting"
    exit 1
fi    

echo "NOTICE: MAKE SURE THAT THE OUTPUT FILES DO NOT ALREADY EXIST IN THE WORKING DIRECTORY"
echo "(out_l2.txt, out_l7.txt, out_l8.txt, out_metrics, out_filenames.txt and out_filenames_FINAL.txt)"
# assign variables
data_dir=$1   						# | 11026_S12_Aligned.sortedByCoord.out.bam.metrics.txt
out_dir=$2

# get rows of interest
for f in $data_dir/*metrics*; do  head -n 2 $f | tail -n 1 >> $out_dir/out_l2.txt ; done   # => parameters
for f in $data_dir/*metrics*; do  head -n 7 $f | tail -n 1 >> $out_dir/out_l7.txt ; done   # => heading
for f in $data_dir/*metrics*; do  head -n 8 $f | tail -n 1 >> $out_dir/out_l8.txt ; done   # => actual data

# get unique heading
uniq $out_dir/out_l7.txt >> $out_dir/out_metrics.txt

# concatenate the heading with the actual data
cat $out_dir/out_l8.txt >> $out_dir/out_metrics.txt  # => output

# get file names
for f in $data_dir/*metrics.txt*; do echo $f >> $out_dir/out_filenames.txt ; done

# drop ‘_Aligned.sortedByCoord.out.bam.metrics.txt’
cat $out_dir/out_filenames.txt | rev | cut -d'/' -f 1 | rev | cut -d "_" -f-2 > $out_dir/out_filenames_FINAL.txt

# add 'ID' heading
echo 'ID'$'\n'"$(cat $out_dir/out_filenames_FINAL.txt)" > $out_dir/out_filenames_FINAL_BIS.txt

# combine in excel: out_metrics.txt and out_filenames_FINAL_BIS.txt
paste -d , $out_dir/out_filenames_FINAL_BIS.txt $out_dir/out_metrics.txt > $out_dir/final.csv

# remove unnecessary files
cd $out_dir
rm out_filenames.txt out_l2.txt out_l7.txt out_l8.txt out_filenames_FINAL.txt out_filenames_FINAL_BIS.txt out_metrics.txt

echo FINISHED !

