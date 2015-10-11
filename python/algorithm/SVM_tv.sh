#! /bin/bash
# Author: Hom 2015.8.18
# Use "python SVM_trainvalid.py data_file modelnum trainradio" to generate training and validation sets first.

if [ -z $1 ];then
	echo "No input data file is given!"
	exit;
fi
fpre=${1%.*}
fext=${1##*.}

setnum=$2
if [ -z $2 ];then
	setnum="1";
fi

cvallist=$3

if [ -z $3 ];then
	cvallist="1 50 100 150 200 300 400 500 700 1000";
fi

echo $1 > ${fpre}_out.log


for cval in $cvallist
do
for i in `seq $setnum`
do

./svm_multiclass_learn.exe -c $cval ${fpre}_train_${i}.${fext} ${fpre}_model_${i}_c${cval}.txt >/dev/null
echo -n ${fpre}_train_${i} $cval " : ">> ${fpre}_out.log
./svm_multiclass_classify.exe ${fpre}_train_${i}.${fext} ${fpre}_model_${i}_c${cval}.txt ${fpre}_train_${i}_out.txt |grep "Zero/one-error" >> ${fpre}_out.log
echo -n ${fpre}_valid_${i} $cval " : ">> ${fpre}_out.log
./svm_multiclass_classify.exe ${fpre}_valid_${i}.${fext} ${fpre}_model_${i}_c${cval}.txt ${fpre}_valid_${i}_out.txt |grep "Zero/one-error" >> ${fpre}_out.log

done
done

if [ -z $3 ];then
	rm -rf ${fpre}_model_*_c*.txt
	rm -rf ${fpre}_train_*_out.txt ${fpre}_valid_*_out.txt
fi