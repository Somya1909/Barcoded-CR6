
## Barcoded CR6

This repository contains code to reproduce analyses in the following manuscript: **"Interferons and tuft cell numbers are bottlenecks for persistent murine norovirus infection"**. 

**Authors**: Somya Aggarwal<sup>1</sup>, Forrest C. Walker<sup>1</sup>, James S. Weagley<sup>1</sup>, Dylan Lawrence<sup>1</sup>, Broc T. McCune<sup>2</sup>, Xiaofen Wu<sup>1</sup>, Lawrence A. Schriefer<sup>1</sup>, Heyde Makimaa<sup>1</sup>, Pratyush Sridhar<sup>1</sup>, Megan T. Baldridge<sup>1,3</sup>.  <sup>1</sup> Division of Infectious Diseases, Department of Medicine, Edison Family Center for Genome Sciences & Systems Biology, Washington University School of Medicine, St. Louis, MO, USA; <sup>2</sup> Department of Microbiology, University of Texas Southwestern Medical Center, Dallas, Texas, USA; <sup>3</sup> Department of Molecular Microbiology, Washington University School of Medicine, St. Louis, MO, USA.

### 1. Unzip fastq file

This is to unzip fastq files by converting fastq.gz (compressed/zipped) files to fast (unzipped) files

```bash
gunzip *.gz
```

### 2. Create file list

This is to create a file that lists all of the file names to feed into bbtools.

```bash
ls *.fastq | sed 's/_[IR]*[12]*_001.fastq//g' | uniq > files.txt
```

### 3. Merge paired reads

Feed merged file list into [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/). Depending on where you have it saved, the path to the script will need to change. 

```bash
for file in `cat files.txt`; do 
/Users/baldridgelab/Desktop/SOMYA/Data_BarcodeVirus/220927_Baldridge_BC001_TissuesNo1_bc_10_8/bbmap/bbmerge.sh in1=$file'_R1_001.fastq' in2=$file'_R2_001.fastq' out=$file'_merged.fastq'
done
```

### 4. Search fastq and barcode

This is the command for having it search each .fastq for each barcode, including 6 nt upstream and downstream. You may need to modify this file if using a different set of barcodes. Additionally, it has the WT, un-inserted CR6 sequence for background. But it should run pretty well, it just sometimes gets confused on the names of files, so edit the sed commands as necessary to make sure it works.

```bash
for f in `cat context.txt`; 
do for f2 in *_merged.fastq; 
do name=$(echo ${f2} | sed s/\.fastq//g); 
count=$((grep -H ${f} ${f2}) | cut -d ':' -f 1 | uniq -c| sed -e 's/^ *//;s/ /,/'|cut -d ',' -f 1); 
echo -e "${f}\t ${name}\t ${count}" >> all_counts.txt; done; 
echo "${f} complete"; 
done
```

### 5. Clean up file names

This next line isn’t strictly necessary, but it does help reduce the obnoxiously long file names in your output. This may also need to be modified slightly if your file names are wonky, but this all_counts_fixed.txt should pretty much be your final output file.

```bash
sed 's/Baldridge_//g' all_counts.txt | sed 's/_[ATCG]*_[ATCG]*_S[1234567890]*_merged//g' > all_counts_fixed.txt
```

### 6. Count the total number of reads

This last is to count up the total number of reads for each one to see how much trash is in there. It is unnecessary to run.  This gives 4*the total # of reads. It needs to be divided by 4 to get the total number of reads.

```bash
for f in *_merged.fastq;
do (wc -l ${f}) >> total_reads.txt;
done

sed 's/Baldridge_//g' total_reads.txt | sed 's/_[ATCG]*_[ATCG]*_S[1234567890]*_merged.fastq//g' > total_reads_fixed.txt
```
