######################################################################################################
# this script processes metrics from Trimmomatic (i.e. submit-2-trim.sh.e226300.1)
######################################################################################################

if [ $# != 2 ] ; then
    echo -e "ERROR: 2 arguments are required: \
	(1) path to data folder [with submit-2-trim.sh.e226300.1 files] \
	(2) and path to folder for output files ... Exiting"
    exit 1
fi    

echo "NOTICE: MAKE SURE THAT THE OUTPUT FILES DO NOT ALREADY EXIST IN THE WORKING DIRECTORY"
echo "(out_l9.txt)"

# assign variables
data_dir=$1
out_dir=$2

# get rows of interest
for f in $data_dir/*trim.sh.e*; do head -n 6 $f | tail -n 1 >> $out_dir/out_l6.txt ; done

# save only one line
head -n 1 $out_dir/out_l6.txt > $out_dir/out_l6_bis.txt 

# get headings:
line_1=$(echo `cut -d" " -f 1-2 $out_dir/out_l6_bis.txt` , `cut -d" " -f 4 $out_dir/out_l6_bis.txt` , `cut -d" " -f 7 $out_dir/out_l6_bis.txt`)

echo $line_1 > $out_dir/line_1.txt

# get column entries
cut -d" " -f 3 $out_dir/out_l6.txt > $out_dir/col_1.txt		# I column
cut -d" " -f 5-6 $out_dir/out_l6.txt > $out_dir/col_2.txt	# II column
cut -d" " -f 8-9 $out_dir/out_l6.txt > $out_dir/col_3.txt	# III column

# concatenate all columns
paste -d "," $out_dir/col_1.txt $out_dir/col_2.txt $out_dir/col_3.txt > $out_dir/cols_all.txt

# concatenate the heading with the actual data
cat $out_dir/line_1.txt $out_dir/cols_all.txt > $out_dir/out_stats.txt  # => output

# extract line 2 to get file names
for f in $data_dir/*trim.sh.e*; do head -n 2 $f | tail -n 1 >> $out_dir/out_l2.txt ; done

cat $out_dir/out_l2.txt | cut -d" " -f 5 | rev | cut -d"/" -f 1 | rev | cut -d"." -f 1 > $out_dir/out_filenames.txt

# add 'ID' heading
echo 'ID'$'\n'"$(cat $out_dir/out_filenames.txt)" > $out_dir/out_filenames_FINAL.txt

# combine: out_metrics.txt and out_filenames_FINAL.txt
paste -d "," $out_dir/out_filenames_FINAL.txt $out_dir/out_stats.txt > $out_dir/final.csv

# remove unnecessary files
cd $out_dir
rm cols_all.txt col_1.txt col_2.txt col_3.txt line_1.txt out_l6.txt out_l6_bis.txt out_filenames.txt out_filenames_FINAL.txt out_l2.txt out_stats.txt

echo FINISHED !

