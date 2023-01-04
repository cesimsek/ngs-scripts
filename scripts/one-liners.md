# NGS/fasta/fastq manipulations


*  Find a file in a directory: find . -iname "bbyeniheader" -print

*  Number of sequences in fasta file: grep -c "^>" d.fasta

*  Number of sequences in a fastq file: zcat ---.fastq.gz | echo $((`wc -l`/4))
 //OR
			source activate qiime1
  			count_seqs.py –i R1.fastq (no fastq.gz) 
 		    conda deactivate

*  Get the average sequence length of a fastq file:
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' foo.fastq

* bash variable with an underscore:  "${line}"_ORF.fas

*  Extract the first “N” sequences from a fasta file:
awk "/^>/ {n++} n>N {exit} {print}" foo.fasta

*  Sorting numerically according to col 2: 
sort  -k2n output.txt > output_sorted.txt 

*  Move all except Tux.png: mv ~/!(Tux.png) ~/New/

*  Get the common strings for two files: awk 'NR==FNR{A[$1];next}$1 in A' file1 file2

//OR
*  Get 1’s for the matches and 0’s for the mismatches:
while read line; do grep -c "$line" VirFinder/node_names_VF.txt; done < MetaPhinder/node_names_MF.txt > same.txt
Number of 1’s and 0’s: sort < same.txt | uniq -c

*  Number of sequences shorter than a threshold (1000 here): 
// requires bioawk
bioawk -cfastx 'BEGIN{ shorter = 0} {if (length($seq) < 1000) shorter += 1} END {print "shorter sequences", shorter}' Gorilla_merged_contigs.fasta

*  Get the sequences above a threshold: 
bioawk -c fastx '{ if(length($seq) > 600) { print ">"$name; print $seq }}' some_seq.fastq

*  Print a certain column only: awk ‘{print $5}’ text.txt

*  Add something to end of all header lines (loop not possible):
sed 's/>.*/&WHATEVERYOUWANT/' file.fa > outfile.fa

*  Prepend the fasta header names (possible with a loop ☺). Prefix only removes sampleName_1,2,3 etc. Part of BBMap package. 
rename.sh in=x.scaffolds.fasta out=header.x.scaffolds.fasta prefix=sampleName addprefix=t
rename.sh in=x.scaffolds.fasta out=header.x.scaffolds.fasta prefixonly

*  To extract FASTA ids:
grep "^>" 1000-4200-picobir.fasta > headernodes2.txt

*  Remove duplicated sequences:
sed -e '/^>/s/$/@/' -e 's/^>/#/' file.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t $'\t' -f -k 2,2  | sed -e 's/^/>/' -e 's/\t/\n/'

*  Get the sequences of the provided header from a FASTA file:
1.	xargs samtools faidx Gorilla.scaffolds.fasta < nodenames.txt > seqs.fasta (the names have to exactly match)

2.	samtools faidx fff.fasta header_name

//OR: input fasta, headerIDs you want to extract, output file:
3.	getsequencesfromHeader.py picobir.fasta picobirnaBLASThits.txt out.fasta

//OR: give the substring of names to be extracted if the main file has longer headers:
4.	filterbyname.sh in=nn.fasta out=ff.fasta names=kk.txt substring=t include=t overwrite=t


* Check the header names (if they are node names) and cut from the last _ to only see the sample names:
grep '>' all_scaffolds.fasta | awk -F_ '{print $NF}' | sort -u

* Cut the last 2nd ‘_’  from the file
awk 'BEGIN{OFS=FS="_"}{if(/^>/){NF-=2}}{print $0}' vp4_prodigal.fasta > vp4_prodigal-2.fasta

*  Sort according to a column (here 2nd): sort -k2 -n yourfile

*  Get only the query IDs (here contigs) with a hit from the blastn output (if the outfmt was 7):
cat blastn.out |awk '/hits found/{getline;print}' | grep -v "#" > top_hits.txt

*  From a blast output where the bitscore is at col12, the e-value is at col11 and %id match is at col3 – only get the best hits for query IDs (contig IDs):
sort -k1,1 -k12,12gr -k11,11g -k3,3gr blastout.txt | sort -u -k1,1 --merge > bestHits

*  Remove the content of file2 from file1 and save as file3:
awk 'NR==FNR{a[$0];next} !($0 in a)' file2 file1 > file3
//OR:
grep -v -f file2 file1 > file3 (apparently a bit slow)

*  Filter based on a column and print the lines suitable:
awk '$2>300 {print $0; }' allb.cont.magnitudes > 300-allb.cont.magnitudes

*  Filter out contig names smaller than 1000: 
// add to bashrc: 
alias FASTAgrep="awk '{gsub(\"_\",\"\\\_\",\$0);\$0=\"(?s)^>\"\$0\".*?(?=\\\n(\\\z|>))\"}1' | pcregrep -oM -f -"
echo "NODE_\d+_length_(\d){4,}_" | FASTAgrep --buffer-size=100000000 contigs.fa | grep -i ">" > output.txt

*  Filter out contigs with coverage less than 20:
echo "NODE_\d+_length_\d+_cov_((\d){3,}|[2-9]\d)\." | FASTAgrep --buffer-size=100000000 contigs.fa | grep -i ">" > out

*  ETE TOOLS: get the taxonomy from the tax IDs (we already have python scripts to do this in gitlab)

cut -f2 CR10994.tab | ete3 ncbiquery --info > CR10994_taxpath.txt

*  Change text in a file
sed -i 's/old-text/new-text/g' input.txt

*  Get the lines that contain a specific string (here 12058) in the a specific column (2 here, change the $2 accordingly if another column. Also specify a field separator if it is different from space and/or tab)
awk '{if($2~/12058/){print}}' AllSamples_RV+.tab.taxidpath.txt > whatever.txt

*  Delete lines with specific word: 
awk ‘!/query/' file  
grep -v "pattern" file > temp && mv temp file
sed -n '/pattern/!p' file   # worked fine

*  Count number of blank lines in a file
grep -cvP '\S'

*  Remove blank lines: sed -i '/^$/d'

*  Remove specific word from a file (but leaves an extra space?)
sed -e 's/\<foo\>//g' file.txt
sed -e 's/>//'    (to remove > )

*  Get the mean of column two:
awk '{ total += $2 } END { print total/NR }' file.txt

*  Get the sum of second column
 awk '{sum+=$2} END {print sum}' magnitudes_virus_all.txt

*  Get the sum of column (2nd here) of a specified row (NODE_9 is the row entry here) by defining a separator (-F option) and print the row entry and the sum
awk -F ',' '$1 ~ /NODE_9/{sum += $2} END {print $1 "," sum}' abund.bbmap.txt

*  Split multi-fasta file into individual fasta files:
split_multifasta.pl –input_file –output_dir --seqs_per_file=10 (without this, each contig is a .fsa file)

*  Divide the second column of the magnitudes file with the total number of trimmed reads :
awk '{print $2/474763769}' Gorilla.magnitudes > relativeCountsGorilla.magnitudes

*  Merge several fasta sequences into one (in one header) :
awk '
   	 /^>/ { 
       	 # print the first header
        	if (c++ == 0) {print} 
        	next
    	} 
   	 /^$/ {next} 
   	 {printf "%s", $0} 
   	 END {print ""}
' a.fasta > b.fasta

* Prepend string to a column in a file with a shell variable : //not tried
awk -v line="$line" 'BEGIN { OFS = "\t" } { $1 = "line" $1; print }' line.m8 > hey

* Delete until the nth _ from the beginning:
awk 'BEGIN{OFS=FS="_"}{if(/^>/){NF-=5}}{print $0}'  gg.fasta > shit.fasta

* Change the ‘_’ in the filenames in a directory to ‘.’ (they have ‘ORF’ in them – specify that to find them easier)
find . -name '*ORF*' -exec sh -c 'f="{}"; mv -- "$f" "${f%_ORF.fas}.ORF.fas"' \;

* Replace the lines of vp6.ORF1.fas with the lines in ‘headers.txt’ (has to be ordered!)
(starting with > here, you can change it in the second line)
awk '/^>/{getline <"headers.txt"}1' vp6.ORF1.fas > file.fasta

* Print the filenames in a directory and put their names (after cutting extensions etc) in a names folder for future loops
find . -name '*ORF*' | cut -d O -f1 | sed s'/.\///g' | tr -d . | sort -r > names.txt

* Change new line (and F) to other character in OS X terminal:
sed -e ':a' -e 'N' -e '$!ba' -e 's/\nF/,F/g'