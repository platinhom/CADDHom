#! /bin/bash

for f in `cat $1`
	do
		mv $f HH${f}
	done
