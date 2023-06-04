# GenomicsCourse
Code used in Genomics in Biodiversity course, University of Oulu 2023

## Log into Puhti:

    ssh akiebler@puhti.csc.fi

    cd /scratch/project_2007265/


## 0) Getting our nanopore file into the correct directory

created a folder for me (ange) in the project, and a subfolder named 'fast5'

    cd /scratch/project_2007265
    mkdir ange
    mkdir ange/fast5
    cd ange/fast5 

*link* the raw fast5 nanopore run file into the new fast5 folder 

    ln -s ../fast5/FAV95239_b49beda4_bbfa58a9_67.fast5

Load all module:

    module load biopythontools
    module load biokit
    module load medaka



## 1a) Basecalling + Demultiplexing.	tool: guppy

You can first allocate some GPU power to your project: "one GPU for 30 minutes"

    sinteractive --account project_2007265 -g 1 --time 00:30:00
    n

running of 'guppy' basecaller with the 'fast' protocol:

     /scratch/project_2007265/ont-guppy/bin/guppy_basecaller -i ../fast5/ -s ../fastq/ -r -c dna_r10.4.1_e8.2_400bps_fast.cfg --barcode_kits "SQK-RBK114-24" --num_callers 1 -x auto --disable_pings


>In the created "fastq" folder there should be some data now.
The most important folder there is "pass", where the fastq files are sorted by barcodes.


## 1b) Demultiplexing.			tool: guppy

... already. Barcode02.


## 2) Quality Control

Creatinges quality control plots and displays other parameters, like N50 and median read length and read quality.

    cd fastq/pass/barcode02
    NanoPlot --fastq fastq_runid_bbfa58a95bda5f200f2528e784b729a78fdb668d_0_0.fastq -o barcode02_nanoplot -p barcode02

>downloading the summary report barcode02_NanoPlot-report.html file in a new command control window using "scp". Example:

    scp akiebler@puhti.csc.fi:/scratch/project_2007265/ange/fastq3/pass/barcode02/barcode02_nanoplot/barcode02_NanoPlot-report.html "C:Downloads"

> Mean read quality Q10.7, 1300 reads; 1,653,973 total bases. 
> Saved in results: barcode02_fastguppy_NanoPlot-report.html

Command to check number of reads in a fastq file:

    awk '{s++}END{print s/4}' *.fq


## Redo 1a) and 2) with superior basecalling [extra]

Redo basecalling with superior guppy model:

    sinteractive --account project_2007265 -g 1 --time 00:30:00

    /scratch/project_2007265/ont-guppy/bin/guppy_basecaller -i fast5/ -s fastq_sup/ -r -c dna_r10.4.1_e8.2_400bps_sup.cfg --barcode_kits "SQK-RBK114-24" --num_callers 1 -x auto --disable_pings

    cd fastq_sup/pass/barcode02/

    NanoPlot --fastq *.fastq -o barcode02_nanoplot -p barcode02_supguppy

    scp akiebler@puhti.csc.fi:/scratch/project_2007265/ange/fastq_sup/pass/barcode02/barcode02_nanoplot/*NanoPlot-report.html "C:Downloads"

> Mean read quality 17.8, 1605 reads; 2,004,736 total bases
> Saved in results: barcode02_supguppy_NanoPlot-report.html

#### INFO:

From here on, we got complete basecalled [ZMUO.062935.fq] file from Marko.

> Compare Nanoplot reports: 
barcode02_supguppy_NanoPlot-report & 
ZMUO.062935_nanoplotNanoPlot-report 
(in C:\Users\kiebl\Desktop\GenomicsCourse\Results)


## 4) Reference based mapping. tool: bwa, samtools

> working in folder /Group_4 now!

Map all reads ZMUO.062935.fq against reference from NCBI MZ726800.fas:

    bwa index MZ726800.fas index_bwa

    bwa mem index_bwa ZMUO.062935.fq > ZMUO.062935.fq.sam

    samtools sort ZMUO.062935.fq.sam > ZMUO.062935.fq_sorted.bam

    samtools index ZMUO.062935.fq_sorted.bam 

> output: .bam and .bai file. View on Tablet [ZMUO.062935.fq_sorted.bam(.bai)]


### Extracting reads and convert files:

    samtools view -b -F 0x904 ZMUO.062935.fq_sorted.bam > ZMUO.062935.mapped.sorted.bam

    bam2fastx -M ZMUO.062935.mapped.sorted.bam -o ZMUO.062935.sorted.extracted.fq

> output: ZMUO.062935.sorted.extracted.fq


## 5) De novo Asssembly

### Extract reads and Assembly with flye:

First assembly, with extracts [ZMUO.062935.sorted.extracted.fq]

Assemble extracted reads:

    flye --nano-hq ZMUO.062935.sorted.extracted.fq -g 40K -o ZMUO.062935.flye -t 2


Copy the assembly.fasta file into our working directory, here [Group_4]:

    cp ZMUO.062935.flye/assembly.fasta Group_4/

### mapping with bwa:

    bwa index assembly.fasta

    bwa mem assembly.fasta ZMUO.062935.fq | samtools sort -o ZMUO.062935_consenus_sorted.bam

    samtools index ZMUO.062935_consenus_sorted.bam

> Look at the .bam and .bai files with Tablet.
> 23.8k bp assembly length. Coverage still uneven. 6182 reads.

=> Now we have more reads that are mapped to the new reference.
=> We can extact those reads again and redo all steps from "extracting reads", create anotherde novo assembly with more reads!



### Second assembly pipeline: extract, assemble (flye/medaka), map

    samtools view -b -F 0x904 ZMUO.062935_consenus_sorted.bam > extract2.bwa.bam
    bam2fastx -M extract2.bwa.bam -o extract2.bwa.fq

    flye --nano-hq extract2.bwa.fq -g 40K -o flye.extract2 -t 3

    medaka_consensus -i ZMUO.062935.fq -d flye.extract2/assembly.fasta -m r1041_e82_400bps_sup_g615 -o flye.extract2.medaka.consensus

> medaka_consensus -i "file with the reads that you used to create the assembly" -d "the assembly.fasta file that was created with flye/wtdbg2" -m "model name of software" -o "output directory"

    bwa index flye.extract2.medaka.consensus/consensus.fasta
    bwa mem flye.extract2.medaka.consensus/consensus.fasta ZMUO.062935.fq | samtools sort -o extract2.flye.bam
    samtools index extract2.flye.bam

> Download .bam and .bai files. View on Tablet.
> 24.8k length = slightly longer



#### EXTRA:

We could try to get larger assemblies by defining a minimum read overlap. The default is 1000bp, but could define it as 4000bp too!

example: flye --nano-hq extracted_reads.fq -g 40K -m 4000 -o extracted_reads.flye -t 3


## wtdbg2

> As I remember and understand, here we use the ZMUO.062935.fq reads to assemble the mt genome with a different program (wtdbg2) and then use medaka again for polishing. 

Assemble extracted reads:

    wtdbg2 -x ont -g 50k -i ZMUO.062935.sorted.extracted.fq -t 2 -fo ZMUO.062935_wtdbg2

Derive consensus:

    wtpoa-cns -t 2 -i ZMUO.062935_wtdbg2.ctg.lay.gz  -fo ZMUO.062935_wtdbg2.raw.fa

And run medaka (version 1.7.2):

    medaka_consensus -i ZMUO.062935.fq -d ZMUO.062935_wtdbg2.raw.fa -m r1041_e82_400bps_sup_g615 -o wtdbg2

> Downloaded .bam and .bai files in wtdbg2 folder.
> 22k bp length assembly, coverage also still uneven.


## 6) Annotation

...with MITOS online tool.

For annotation we've used the medaka consensus assembly with a length of 24.8k bp. 
