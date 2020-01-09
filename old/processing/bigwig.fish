set folder "/hps/nobackup/stegle/users/ricard/Guo_2017/acc/not_used/ESC"
cd $folder

for i in *.bw
	echo $i
	# bigWigToWig $i $i.wig | wig2bed --input=fmt - > $i.bed
	set j (echo $i | sed 's/\.[^.]*$//')
	# bigWigToWig $i $j.wig
	sed '/^#\|^lambda/ d' $j.wig | sortBed | cut -f 1,2,4 > $j.tsv
	# wig2bed < $j.wig > $j.tsv 
end

# rm *.wig
pigz -p 4 *.tsv


