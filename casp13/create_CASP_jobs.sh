repeat=20
filename="run"
copyname="run1"
infile=$1
rounds=$2
starttemp=$3

starttemp=${starttemp:=600}

echo "The program has started!"
echo "1"
cd $copyname
python3 CASP_jobs_sub.py $rounds
cd ..
for ((i=2; i<= $repeat; i++))
do
        echo "$i"
        mkdir $filename$i
        cp ./$copyname/*   ./$filename$i
        cd ./$filename$i
        sed "s/create $starttemp\.0 20000/create $starttemp\.0 $RANDOM/"  $infile > temp.in
        cp temp.in $infile
        rm -rf temp**
        python3 CASP_jobs_sub.py $rounds
	cd ..
done
