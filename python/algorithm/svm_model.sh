#! /bin/bash
# Author: Hom 2015.8.18
# Use "python SVM_trainvalid.py data_file modelnum trainradio" to generate training and validation sets first.
# Use "./svm_model.sh data_file modelnum cval" to generate many model.

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

cval=$3
if [ -z $3 ];then
	cval="1";
fi

for i in `seq $setnum`
do
./svm_multiclass_learn.exe -c $cval ${fpre}_train_${i}.${fext} ${fpre}_model_${i}.txt > /dev/null
done