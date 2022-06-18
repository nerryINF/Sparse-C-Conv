#/bin/bash
num_vox=$(wc -l < features_numpoints.csv) 
file=../data.h
# num
printf "int num_sparse = $num_vox; \n" > $file
# features
printf "\nfloat features[$num_vox][5]={\n" >> $file
sed 's/^/{/' features_numpoints.csv | sed 's/$/},/' | sed '$ s/.$/};/'>> $file
# indices
printf "\nint indices[$num_vox][3]={\n" >> $file
sed 's/^/{/' indices.csv | sed 's/$/},/' | sed '$ s/.$/};/'>> $file
# shape
printf "\nint sparse_shape[3]=" >> $file
sed 's/^/{/' sparse_shape.csv | sed 's/$/},/' | sed '$ s/.$/;/'>> $file
# kernel
# generate
echo "27"> kernel.csv; for i in {1..27}; do echo "1," >> kernel.csv; done
# add to .h
k_vol=$(head -n 1 kernel.csv)
printf "int KL_VOL = $k_vol; \n" >> $file
printf "\nfloat kernel[$k_vol]={\n" >> $file
tail -n +2 kernel.csv >> $file
printf "};" >> $file
cat $file

