#! /bin/bash
# Author: Hom 2015.8.18
# Use "python SVM_trainvalid.py data_file modelnum trainradio" to generate training and validation sets first.
# Use "./svm_test.sh data_file test_file modelnum" to test with many models.

if [ -z $1 ];then
	echo "No input data file is given!"
	exit;
fi
fpre=${1%.*}
fext=${1##*.}

setnum=$3
if [ -z $3 ];then
	setnum="1";
fi

testf=$2
if [ -z $2 ];then
	echo "No input test file is given!"
	exit;
fi
tpre=${2%.*}
text=${2##*.}

echo $1 $2 > ${tpre}_testout.txt


for i in `seq $setnum`
do
echo -n ${fpre}_model_${i} " : ">> ${tpre}_testout.txt
./svm_multiclass_classify.exe $testf ${fpre}_model_${i}.txt ${tpre}_model_${i}_out.txt |grep "Zero/one-error" >> ${tpre}_testout.txt
done
