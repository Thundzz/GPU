i='100'
values='100 200 500 750 1000 2500 5000 7500 10000 15000 20000 35000'
a=4

if [ "$#" -ne "$a" ]
	then 
		echo "Format : 'bin_filename' '-s/-O/-d' 'output_device' 'dat_output_file'"
	else
		echo -e "nb_atoms  \t  time/iteration (usec)" > $4.dat
	for j in $values 
	do
		echo "$j"
		echo -e "$j\t" > out
		../ocl/fichiers/bin/$1 -v $2 $3 -i $i -n $j 2>&1 | grep -E '\[PERF\][ \t]*[0-9]+' | grep -o -E "$j[[:space:]]+[0-9]+[[:space:]]" | awk '{print $2}' >> out
		tr '\n' ' ' < out > tmp
		echo " " >> tmp
		cat tmp >> $4.data
	done
fi
rm out
rm tmp
