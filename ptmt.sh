#!/bin/bash
###################################
# Author: Xiaoyu Zhou
# Email:  xyzh@biomed.au.dk
###################################

function readfile() { 
for file in `ls $1/*.fna`
do 
		if [ -f $1"/"$file ];then
			#echo $1"/"$file 
			name=`basename $file`
			tag=${name%.fna}
			specie=${name%_*.fna}
			typeset -l specie
			#`perl cut.pl $tag`
			#echo "$tag"
			db="blast/"$tag"_nu.txt"
			q="blast/"$tag"_mt.txt"
			p="blast/"$tag"_pt.txt"
			#echo "$db"
			#echo "$q"
			#if [ -s $db ];then
			#	if [ -s $q ];then
					#`mv $db $q /home/zhouxy/CCDB/NUMT/`
					#`/home/zhouxy/RepeatMasker/RepeatMasker -parallel 10 -species $specie -html -gff -dir repeat $name > log 2>> err`
					#echo "$tag RepeatMasker Done"
			#		blastn="blast/"$tag"_blastn.txt"
			#		reblastn="blast/"$tag"_reblastn.txt"
					#`cd /home/zhouxy/CCDB/fna/blast/`
			#		`~/ncbi-blast-2.6.0+/bin/makeblastdb -in $db -dbtype nucl`
					#`cd /home/zhouxy/CCDB/fna/blast/`
			#		`~/ncbi-blast-2.6.0+/bin/blastn -db $db -query $q -out $blastn -word_size 11 -evalue 0.0001 -outfmt 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 `
			#		echo "$tag BLAST Done"
					#`cd /home/zhouxy/CCDB/fna/blast/`	
			#		`~/ncbi-blast-2.6.0+/bin/makeblastdb -in $q -dbtype nucl`
					#`cd /home/zhouxy/CCDB/fna/blast/`
			#		`~/ncbi-blast-2.6.0+/bin/blastn -db $q -query $db -out $reblastn -word_size 11 -evalue 0.0001 -outfmt 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 `
					#`cd /home/zhouxy/CCDB/fna/blast/`
			#		`mv $blastn $reblastn /home/zhouxy/CCDB/NUMT/blast_result/`
			#		echo "$tag REBLAST and Move Done"
					#`cd /home/zhouxy/CCDB|perl data.pl $tag >log 2>>err`
					#echo "$tag ALL DONE"
			#		`perl numt.pl $tag`
			#		echo "perl numt.pl $tag Done"
			#	fi
			#fi
			if [ -s $q ];then
                                if [ -s $p ];then
                                        #`mv $db $q /home/zhouxy/CCDB/NUMT/`
                                        #`/home/zhouxy/RepeatMasker/RepeatMasker -parallel 10 -species $specie -html -gff -dir repeat $name > log 2>> err`
                                        #echo "$tag RepeatMasker Done"
                                        blastn_ptmt="blast/"$tag"_blastn_ptmt.txt"
                                       reblastn_ptmt="blast/"$tag"_reblastn_ptmt.txt"
                                        #`cd /home/zhouxy/CCDB/fna/blast/`
                        #                `~/ncbi-blast-2.6.0+/bin/makeblastdb -in $q -dbtype nucl`
                                        #`cd /home/zhouxy/CCDB/fna/blast/`
                        #                `~/ncbi-blast-2.6.0+/bin/blastn -db $q -query $p -out $blastn_ptmt -word_size 11 -evalue 0.0001 -outfmt 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 `
                        #                echo "ptmt $tag BLAST Done"
                                        #`cd /home/zhouxy/CCDB/fna/blast/`
                        #                `~/ncbi-blast-2.6.0+/bin/makeblastdb -in $p -dbtype nucl`
                                        #`cd /home/zhouxy/CCDB/fna/blast/`
                        #                `~/ncbi-blast-2.6.0+/bin/blastn -db $p -query $q -out $reblastn_ptmt -word_size 11 -evalue 0.0001 -outfmt 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 `
                                        #`cd /home/zhouxy/CCDB/fna/blast/`
                        #                `mv $blastn_ptmt $reblastn_ptmt /home/zhouxy/CCDB/NUMT/blast_result/`
                        #                echo "ptmt $tag REBLAST and Move Done"
                                        #`cd /home/zhouxy/CCDB|perl data.pl $tag >log 2>>err`
                                        #echo "$tag ALL DONE"
                                        `perl ptmt.pl $tag`
                                        echo "perl ptmt.pl $tag Done"
                                fi
                        fi
	
		fi
done
}
#read -p "Please input the directory or file you want to their raed name:" n;
[ $# -eq 0 ] && { echo "Usage: $0 file path"; exit 1; }
#folder="./test"
readfile $1
