---
title: Variants in Alignment
author: Xiaoju (Ju) Zhang
---


AVAlign - Worknote Build whole genome alignment
----
Edited / Authored by Xiaoju (Ju) Zhang
Updated: 8/11/2015

###SECTION 1: **BACKGROUND KNOWLEDGES**

**`Terms:`**
`chains and nets`:
http://genome.ucsc.edu/goldenpath/help/chain.html
http://genomewiki.cse.ucsc.edu/index.php/Chains_Nets
http://www.personal.psu.edu/zuz17/blogs/psu_life/2011/02/understand-ucsc-netchain-alignment-1.html
`liftover`:
http://genome.sph.umich.edu/wiki/LiftOver

**`How to align?`**
**Pages**
[Vertebrate Multiz Alignment & Conservation (46 Species)](https://cgwb.nci.nih.gov/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=cons46way&hgta_table=multiz46way&hgta_doSchema=describe+table+schema)
[HowTo: Syntenic Net or Reciprocal Best][3]
[3]: <http://genomewiki.ucsc.edu/index.php/HowTo:_Syntenic_Net_or_Reciprocal_Best>
[**Whole genome alignment howto**][5]
[5]: <http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto>
[LASTZ manual](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html) 
	- BLASTZ parameters
	- scoring matrix
	- options 
[UCSC lastz parameter documentation for hg19 100way alignment](http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters)
[Table Browser User's Guide ](https://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html),  inluding:
	- Displaying CDS FASTA alignments.
	- Explanation of CDS FASTA header format
[UCSC Genome Bioinformatics FAQ for Data and Downloads](http://genome.ucsc.edu/FAQ/FAQdownloads.html#TOP)
	- Repeat-masking data
	- Availability of repeat-masked data
[Minimal Browser Installation](http://genomewiki.ucsc.edu/index.php/Minimal_Browser_Installation)
[Preparing whole genome alignments](http://genomeview.org/manual/Preparing_whole_genome_alignments)
**Mailing list**
[Mailing list:: Multiple alignments of 99 vertebrate genomes with Human, Valya Burskaya][6]
[6]:<https://groups.google.com/a/soe.ucsc.edu/forum/#!msg/genome/S1Z8EQb-hlk/jc2oP3cX-vQJ>
[Mailing list:: Multiple alignment, UCSC mailing list, Hong Lu](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/pvc0O8Civpg/XTnJv4ME6-wJ)
[Mailing list:: UCSC human multiple alignment (protein), UCSC mailing list, Hong Lu](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/0cqDoYRtG0Y/Qk4l-f5oH1wJ)
[Mailing list:: **Building your own multiple sequence alignment**, Nimrod Rubinstein](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/xHJLq1yKMyI/YLfeRWQ1ldwJ)
[Mailing list:: Retrieve the multiple seq alignment with defined regions.](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/F_YjGiYMcDY/xKans0pSV-UJ)
[Mailing list:: How to generate pairwise alignment data](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/zcHuWtmw-LE/QxI5cx_6UVMJ)
[Mailing list:: alignments of primate genomes, Qianfeng Cliff Zhang](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/2Q8gSkwuQkA/CbIVBLGayL0J)
> If I have to generate it by myself, is "doBlastzChainNet.pl" the correct script that I should use? What parameters would you recommend? 

[Mailing list:: How to get multiple alignment with Multiz?, Jingting Guan](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/iPhOIWpf3iA/7JT8gDuESN4J)
>I have gotten local alignments by using Blastz, then how to get a multiple alignment with Multiz and what the command is.

[Mailing list:: Stand alone copy of Multiz and 4 species alignments with TBA and Multiz and the Gmaj viewer](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/7d4Pe4BJd30/Al1QYt59yQAJ)
>Hi I am trying to do a 4 species, man, mouse, rat, dog (or chicken) of  the DMP1 and MEPE genes, where I include 10kb of 5'flanking and the transcription unit.  I have all the sequences. I am not quite sure where to go from here.

[Mailing list:: Application of MULTIZ, Christian Rödelsperger](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/9PIkDnGfVRg/N599WO99qUAJ)
>I am trying to use the multiz program together with pairwise alignments from UCSC. 


[**Mailing list:: Whole genome alignment Howto Question, Xiaoju Zhang**](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/BZjUyNgRKnM/2YR0HceVcuMJ)

> Hi Ju,
> 
> Thank you for your questions about the process of creating a Multiz
> multiple alignment.
> 
> >Chains and nets introduction page (http://genomewiki.ucsc.edu/index.php/Chains_Nets) explains these two
> terms. "Net is hierarchical collection of chains". If I want to go to
> multiz step to generate multiple alignment, do I have to do netting,
> or I can start from all.chain file? The LiftOver files download from
> UCSC are chain files but not netted, right?
> 
> Yes, you will need to generate net files, and the type of net files
> will depend on the quality of the assemblies you are using. You can
> see some examples of the different types of net files in the
> "alignment type" column here:
> http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics
> 
> >the Whole genome alignment HowTo tutorial actually only shows how to generate pairwise alignment, and then generate .maf format of pairwise
> alignment, is that right?
> 
> Yes, if you want more than just two species in your alignment you will
> need to use something like Multiz to create your N-way multiple
> alignment.
> 
> >With results from pairwise alignment (Liftover), do you have any suggestion for me to learn how to create multiple alignment?
> 
> Our "make docs" contain all of the steps we used when generating a
> N-way multiple alignment and the associated phastCons and phyloP data.
> You can the see the steps used to create the 20-way multiple alignment
> on hg38 here:
> http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob_plain;f=src/hg/makeDb/doc/hg38/multiz20way.txt;hb=HEAD.
> 
> You can use the steps described in the make doc as a guide for
> creating your own multiple alignments. 
> 
> >however for two species alignment, how do I define the tree?
> 
> As noted in my response to your last question, you can find this
> information in the make doc. You can use the "tree_doctor" program
> included in the PHAST package, http://compgen.cshl.edu/phast/, to
> extract the tree for a subset of species from a larger tree. 
> 
> For example, we have a tree for a 10-way multiple alignment on ce9
> here:
> http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob_plain;f=src/hg/utils/phyloTrees/ce9.10way.nh
> 
> This file has the contents:
> 
> ((((((ce9:0.003,((cb3:0.005,caeRem3:0.003):0.004,caePb2:0.013):0.002):0.001,
> caeJap3:0.023):0.04,haeCon1:0.06):0.06,priPac1:0.12):0.06,  
> (melHap1:0.14,melInc1:0.15):0.03):0.06,bruMal1:0.24);
> 
> Using 'tree_doctor' from the PHAST package:
> 
> $ tree_doctor --prune-all-but ce9,cb3 ce9.10way.nh
> 
> The resulting tree is:
> 
> (ce9:0.003000,cb3:0.011000);
> 
> You can then substitute ce10 for ce9.
> 
> We have other phyloTrees that include other species available here: 
> http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=tree;f=src/hg/utils/phyloTrees
> 
> I hope this is helpful. If you have any further questions, please
> reply to gen...@soe.ucsc.edu. All messages sent to that address are
> archived on a publicly-accessible Google Groups forum. If your
> question includes sensitive data, you may send it instead to
> genom...@soe.ucsc.edu.
> 
> Matthew Speir UCSC Genome Bioinformatics Group
> 






<br>
###**SECTION 2: AVAlign - Pairwise alignment**
**Dataset** 
The following datasets for *C. elegans* and  *C. briggsae* are downloaded from UCSC genome browser website.

*C. elegans-ce10:*
In the original script, ce9 was used. However, it was replaced with ce10.
ce10.2bit, ce10.chrom.sizes, cd10.chromFaMasked.tar.gz:
http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/
http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit
http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.chrom.sizes
http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/chromFaMasked.tar.gz


*C. briggsae- cb3:*
cb3.2bit, cb3.chrom.sizes, cb3.chromFaMasked.tar.gz
http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/
http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/cb3.2bit
http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/cb3.chrom.sizes
http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/chromFaMasked.tar.gz


<br>
#### **Build Pairwise Whole Genome Alignment Approach 1**
--- modified version of tutorial based on the original one for _C. elegans-ce10_ and  _C. briggsae- cb3_
> **Date: 7/16/2015** 
> **work log: sandbox/align_howto5.1, 6.1 analysis/work_3a/align_ce10.cb3_step_by_step, align_ce10.caePb1_step_by_step** @usav1svBIOd1
align_howto5.1: aligning for only chrI.
align_howto6.1: aligning for whole genomes.
align_ce10.caePb1_step_by_step: align whole genome for ce10 and cb3
align_ce10.caePb1_step_by_step: align whole genome for ce10 and caePb1

**OBJECTIVE:**
Run the whole genome alignment methods and tools, following the tutorial: [Whole genome alignment howto](http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto). Do a whole genome alignment of _C. elegans-ce10_ and  _C. briggsae- cb3_ (_C. intestinalis_ V2 and _C.savignyi_ V2 are used in the original tutorial, which some of the vital data is missing, e.g. the database to fetch tandem repeat mark of _C.savignyi_ .)


>Due to the time consuming of alignment, and the fact that C. elegans_ and  _C. briggsae_ show an almost complete conservation of synteny, with 1:1 orthologs present on a single chromosome in one species also found on a single chromosome in the other, herein for this test, chromosome 1 is chosen to represent the whole genome for both species.

reference: {Comparison of C. elegans and C. briggsae Genome Sequences Reveals Extensive Conservation of Chromosome Organization and Synteny, LaDeana Hillier, et. al PLOS 2007} http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050167

*Test chrI alignment hit a dead end at the step calling database, which was originally generated in the whole genome manner. It needs extra testing to figure out how to solve this issue. 

<br>
**MATERIAL AND METHODS:**
**Methods:** 
[Whole genome alignment howto](http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto)


**Tools:** 
**- Mandatory tools**
>**kentUtils**: 
> kentUtil: https://github.com/ENCODE-DCC/kentUtils
> This package includes the key tools for the alignment. The tools from the tutorial  are almost all from kenUtils.
>```bash
$ cd kentUtils
$ make # install general tools
$ cd src/hg #make hg utils
$ make
$ cd src/utils #make other utils
$ make
>```
*if during the aligning steps, some command is missing, firstly to check kenUtils has it or not, and it might not compiled.

>**LASTZ**: 
>``` bash
$ wget http://www.bx.psu.edu/miller_lab/dist/lastz-1.02.00.tar.gz
$ tar -zxvf lastz-1.02.00.tar.gz
>```
Change folder from lastz-xxx.xx to lastz; editing the make‑include.mak file to set the definition of installDir. Also, be sure to add the directory installDir that you choose to your $PATH in `~/.bash_profile`.
>```bash
$ cd lastz
$ make
$ make install
>```

>**Tandem Repeat Finder**:http://tandem.bu.edu/trf/trf.download.html
Download Tandem Repeat Finder and put it into executable path (for example /usr/local/bin), and name it as "trf", make it executable `chmod +x trf`, this tool is used in the second step mask the tandem repeat.

>**Scoring Matrix**: 
>HoxD55.q is a scoring matrix and it can get from the following website. There are other scoring matrices available.
http://genomewiki.ucsc.edu/index.php/Ce10_conservation_lastz_parameters

>**UCSC Genome Browser**: 
http://genomewiki.ucsc.edu/index.php/Source_tree_compilation_on_Debian/Ubuntu
Details read: https://github.com/maximilianh/browserInstall
>This shell script makes installation much easier, and it also offers command line to download and import database you are referring. The details instruction can be found in the repo's README.
>```bash
$ sudo -i
$ wget https://raw.githubusercontent.com/maximilianh/browserInstall/master/browserInstall.sh
$ sudo bash browserInstall.sh
	#use the same script to download and import assemblies database.
$ sudo bash  browserInstall.sh hg19 ce10 cb3
	#this may complains about the size of disk if the partition of the UCSC genome browser and mysql folder are on root mount. Do following.
	#create those two folders in my working folder.
$ cd /home/xzhang/Tools/browser
$ mkdir gbdb
$ mkdir mysql
	# create symlink 
$ ln -s /gbdb    /home/xzhang/Tools/browser/gbdb
$ ln -s /var/lib/mysql      /home/xzhang/Tools/browser/mysql
>```
**Phast**:
This package contains phastCons to compute Phastcon scores.
Install instruction
    Part 1 - Installing Clapack - (If you already have Clapack installed, skip to Part 2)
    1. Download Clapack from the following URL http://www.netlib.org/clapack/clapack.tgz
    2. Unzip clapack.tgz with the command 'tar -xvzf clapack.tgz'
    3. Go into the newly created Clapack directory (i.e. 'cd CLAPACK-3.2.1')
        and type 'cp make.inc.example make.inc && make f2clib && make blaslib && make lib'
    Part 2 - Installing Phast
    4. Download a copy of Phast from http://compgen.bscb.cornell.edu/phast/
        and extract the contents of phast*.tgz using 'tar -xvzf phast*.tgz'
    5. Change directory to 'phast/src/' and run 'make CLAPACKPATH=/usr/local/software/clapack'
        replacing '/usr/local/software/clapack' with the path of your
        Clapack install (e.g., CLAPACKPATH=/home/username/CLAPACK-3.2.1)
    6. The Phast binaries should be created in the '../bin/' directory
>```bash
$ wget http://www.netlib.org/clapack/clapack.tgz
$ tar -xvzf clapack.tgz
$ mv CLAPACK.3.2.1 CLAPACK && cd CLAPACK
$ cp make.inc.example make.inc && make f2clib && make blaslib && make lib
$ wget http://compgen.bscb.cornell.edu/phast/downloads/phast.v1_3.tgz
$ tar -xvzf phast.v1_3.tgz
$ mv  phast-1.3  phast && cd phast/src
$ make CLAPACKPATH=/home/xzhang/Tools/CLAPACK
	# then add /home/xzhang/phast/bin to $PATH in ~/.bash_profile
>```

______
**- Optional Tools**
Samtools library is required for UCSC genome browser.
>**Samtools:** https://github.com/samtools/samtools
>``` bash
$ wget https://github.com/samtools/samtools
$ unzip develop.zip #unzip the package
$ cd samtools && sudo make install #install
>```
>**Samtools library**:  https://github.com/samtools/htslib 
>``` bash
$ wget https://github.com/samtools/htslib 
$ unzip develop.zip #unzip the package
$ cd htslib 
$ autoconf && ./configure && make  && make install #install
>```
>**MySQL Server**: http://askubuntu.com/questions/489815/cannot-install-mysql-server-5-5-the-following-packages-have-unmet-dependicies
>http://www.rackspace.com/knowledge_center/article/installing-mysql-server-on-ubuntu
> downgrade MySQL from 5.6 to 5.5
>```bash
sudo apt-get purge mysql-client-core-5.6
sudo apt-get autoremove
sudo apt-get autoclean
sudo rm  -r /var/lib/mysql # this step is added based on an error message
sudo apt-get install mysql-client-core-5.5
sudo apt-get install mysql-server  
>```


<br>
**RESULTS:**
Some of the detailed command lines are skipped, which can be found in the original Approach 1 note or refer to original tutorial.
1. Download dataset into `/sandbox/align_howto5.1/genome`,  `/sandbox/align_howto6.1/genome`. For this approach, the sequence files are in FASTA format. cb3 has two sets of FASTA files, random and regular.  Here all the .fa files are used, since alignment tools will take care of overlaps. And all of the FASTA files have been already masked, which allow us to skip Mask tandem repeat step.
Explanation from [UCSC genome FAQ: chrN_random tables](http://genome.ucsc.edu/FAQ/FAQdownloads.html) about random FASTA files:
>In the past, these tables contained data related to sequence that is known to be in a particular chromosome, but could not be reliably ordered within the current sequence.

2. Decompress the sequence files. For some FASTA file is in one large file, we need to split them into smaller ones for better parallel computing, and Lastz is good for short sequence alignment. For FASTA file in this test, all the sequences are already splitted in chromosomes. One can further split them to optimize the alignment efficiency and accuracy using  `faSplit` if more cpu clusters are available. Check its usage, by type `faSplit`.  This step is equivalent to the scripted approach, using perl script `partitionSequence.pl` to make a list of files suitable for creating a parasol job list. **If the genome FASTA file is already fragmented in chromosomes, this step can be skipped.**

	
	```bash
	#creat cb3 and ce10 sub folders, put a copy of the  genome FASTA files of each species into their subfolders individually. Size of the fragment can refer to script method, and it uses the following commands. 
	$ cd cb3
	$ faSplit byName ../cb3.fa ./
	$ cd ../ce10
	$ faSplit base ../ce10.fa ./

	```
The script method split the target and query sequence with different schema.
>`echo "1. partitioning the two genomes into:"`
`echo "   a. 10,000,000 overlapping 10,000 chunks for the target sequence"`
`echo "   b. 20,000,000 no overlap chunks for the query sequence"`


3. **For this test, since FASTA files are already masked, this step can be SKIPPED.**  Mask tandem repeat. In the tutorial, the cs2.fa is not masked. For unmasked sequences do following. 
Note that `trfBig` needs `trf` available. `trf` does not come with kenUtils, it can be downloaded from http://tandem.bu.edu/trf/trf.download.html, change the tool name to `trf`, details see above **Tandem Repeat Finder** installing instruction.

	```bash
	$ mkdir trf
	$ for i in cs2/*.fa; do trfBig $i trf/`basename $i`; done
	```

9. **Since the genome size file are available in database and downloaded, this step could be skipped.** However, here in this test,  chr I is assumed as the whole genome, it is necessary to calculate the genome size for chr I. In addition, the hashtag names have been changed due to the split, the genome need to be merged together, instead of using the original chrI file. 
Determine the genome sizes of the two species
	- merge all the fragmented .fa for each species
	- measure the detailed size of each chromosome. 
	```bash
	$ #mkdir genome
	$ for spec in ce10 cb3; do cat $spec/*.fa > ../${spec}.fa; faSize -detailed ../${spec}.fa > ${spec}.chrom.sizes; done  

	```

4. Convert fasta format sequence to nib format, since LASTZ ignore lower case characters which might be seen in fasta file. Script approach uses 2bit format, by using tool `faToTwoBit` Since 2bit file is more compact and efficient, in this test, 2bit format will be used. Instead, the original tutorial uses .nib format. The description of these format can be found [**at UCSC genome FAQ**](http://genome.ucsc.edu/FAQ/FAQformat.html).
	
	```
	$ for i in ce10/*.fa cb3/*.fa; do faToTwoBit $i `echo $i | sed -e s/.fa/.2bit/`; done
	# nib file will be used for chaining.
	$ for i in ce10/*.fa cb3/*.fa; do faToNib $i `echo $i | sed -e s/.fa/.nib/ `; done 
	```
5. Alignment with LASTZ. Align every sequence of one species (cb3) onto every sequence of another species (ce10) with LASTZ. This step will be computationally expensive. HoxD55.q is a scoring matrix file for alignment computing.  Two species with different evolutionary distance may need to use different scoring matrix. [UCSC ref](http://genomewiki.ucsc.edu/index.php/Ce10_conservation_lastz_parameters), [Blastz ref](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#options_scoring)
>HoxD55
>|-|**A**|**C**|**G**|**T**|
|:--:|:--:|:--:|:--:|:--:|
|**A**|91|-90|-25|-100|
|**C**|-90|100|-100|-25|
|**G**|-25|-100|100|-90|
|**T**|-100|-25|-90|91|

	Bash scripts to generate job scripts and run the alignment. At this step actually the order of input species matters. Target should be at first position and Query should be the following one. For instance: `lastz Target.nib Query.nib H=2000 Y=3400 L=6000 K=2200 Q=HoxD55.q > Target-Query.lav`. The original tutorial did the reversed alignment, and an extra step was used to swap the order. Options meaning look at [this page](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#options_scoring); referece [UCSC lastz parameter documentation for hg19 100way alignment](http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters)

	```bash
	#Divide the job into separate parts, and run parallelly. 
	#!/bin/bash
    mkdir -p lav
    batch=batch_lastz.bat
    echo ''> $batch
    for i in $( ls -v ce10/*.2bit ); do
        id=`basename $i .2bit`
    cat <<EOF > lastz_$id.sh
    for j in cb3/*.2bit; do
        echo "lastzing $i onto \$j ..."
        lastz $i \$j  H=2000 Y=3400 L=6000 K=2200 Q=HoxD55.q > lav/\`basename $i .2bit\`-\`basename \$j .2bit\`.lav;
    done
    EOF
    chmod +x lastz_$id.sh
    echo lastz_$id.sh '&' >> $batch
    done
    chmod +x $batch
	```
	```
	$ batch
	```

6. Convert the lav files into the more compact psl format.

	```bash
	 $ cd lav
	 $ for i in *.lav; do echo $i; lavToPsl $i `basename $i .lav`.psl; done; 
	 $ rm -f *.lav
	```
7. Merge and split all the psl file. pslSplitOnTarget may need to compile. It is in //kentUtils/src/hg/pslSplitOnTarget 

	``` bash
	$ cd lav
	$ cat *.psl > ../all.psl.1
	#remove the line starts with '#', since pslSplitOnTarget does not skip it and complains. pslSwap does the '#' removal. 
	$ cd ..
	$ grep -v '^#' all.psl.1 > all.psl 
	#same result with `pslSwap all.psl.1 all.psl.1-swap & pslSwap all.psl.1-swap all.psl`
	$ mkdir -p psl 
	$ pslSplitOnTarget all.psl psl/ 
	```
option -lump, is useful with scaffolds, hashes on targ name to lump together. In other words, it lump all the target chromosome fractions, and generate one single psl, and name merged file with systematic indexing number. Otherwise, without -lump option, it will leave target fractions, only merge same query chromosomes. Try `pslSplitOnTarget all.psl psl/ -lump`, since here we assume chromosome 1 as the whole genome, and it should have multiple .psl output, the -lump is not chose. The default method should use `-lump`

	> [7]. Swap the target and query order, if in previous step [5], the order is not kept correctly.
	
	>``` bash
	$ cd psl
	$ cat *.psl > ../all.psl
	$ pslSwap ../all.psl ../all-swap.psl
	$ cd .. & mkdir psl
	$ pslSplitOnTarget all-swap.psl psl/ -lump
	```
8. Now the chaining 	
	
	```bash
	$ mkdir -p chain
	$ for i in psl/*.psl; do echo $i; axtChain $i ce10 cb3 chain/`basename $i .psl`.chain -linearGap=loose -psl; done
	
	
	```

10. Sort and filter the chains with the chromosome size information provided from previous step. More infor can be found from this [mailing list post](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/zcHuWtmw-LE/M4kdIInQ7TQJ)

	```bash
	$ chainMergeSort chain/*.chain > all.chain
	$ chainPreNet all.chain ce10.chrom.sizes  cb3.chrom.sizes all.pre.chain
	```

12.  Load chains into UCSC genome browser. **Skip the UCSC genome visualization,** since it is installed on server.
	> ``` bash
	$ hgLoadChain ci2 chainCioSav2 all.pre.chain
	```
	>Adding a section similar to the following to trackDb.ra, running make alpha DBS=<yourdbname> ZOO= there, the track should then appear on the browser:
`track chainCioSav2
 shortLabel C. savignyi chain
 longLabel C. savignyi chain
 group compGeno
 priority 125
 visibility hide
 color 100,50,0
 altColor 255,240,200
 spectrum on
 type chain cioSav2
 otherDb cioSav2`

13.  Netting. Then the netting itself, combine the chains into nets and add synteny information
	```bash
	$ chainNet all.pre.chain -minSpace=1  ce10.chrom.sizes cb3.chrom.sizes stdout /dev/null | netSyntenic stdin noClass.net
	
	# This step is not necessary for getting same alignment result. It just provide more information from database for UCSC genome browser display, so alternatively, `cp noClass.net cb3.net` works.
	$ netClass -noAr noClass.net ce10 cb3 cb3.net
	```
`noClass.net` is the in.net, and `cioSav2.net` is out.net, `ci2` is the target databse tDb - database to fetch target repeat masker table information, and `cioSav2` is query database  qDb - database to fetch query repeat masker table information. Details do `netClass` get man page.

14. Maffing. Generate MAF file for pairwise alignment.  Note that in the prefix, the '.' sign is necessary and important.

	```bash
	$ netToAxt cb3.net all.pre.chain ce10/ cb3/ stdout | axtSort stdin ce10.cb3.axt
	$ axtToMaf ce10.cb3.axt ce10.chrom.sizes cb3.chrom.sizes ce10.cb3.maf -tPrefix=ce10. -qPrefix=cb3.
	```
<br>



**`work_steps_part1.sh`**
```bash
#!/bin/bash

WORKDIR=$PWD
tspec=$1
qspec=$2
#---------------------------------
# Prepare FASTA files
#----------------------------------
for spec in ${tspec} ${qspec}; do
    mkdir -p $spec
    for f in $WORKDIR/genome/$spec/chr*.fa.masked; do 
        f=`basename $f`
        cp $WORKDIR/genome/$spec/$f $WORKDIR/$spec/${spec}_${f%.masked}
    done
done

#---------------------------------
# Create chrom.sizes file
#---------------------------------
for spec in $tspec $qspec; do 
	cat $spec/*.fa > ../${spec}.fa; faSize -detailed ../${spec}.fa > ${spec}.chrom.sizes; 
done 

#---------------------------------
# convert FASTA to .2bit and .nib
#----------------------------------
for i in ${tspec}/*.fa ${qspec}/*.fa; do
    faToTwoBit $i `echo $i | sed -e s/.fa/.2bit/`;
    faToNib $i `echo $i | sed -e s/.fa/.nib/ `;
done

#-----------------
# lastz to align
#-----------------
mkdir -p $WORKDIR/lav
batch=batch_lastz.bat
echo ''> $batch
for i in $( ls -v ${tspec}/${tspec}_chr*2bit ); do
    id=`basename $i .2bit`
cat <<EOF > lastz_$id.sh
for j in ${qspec}/${qspec}_chr*.2bit; do
    echo "lastzing $i onto \$j ..."
    lastz $i \$j  H=2000 Y=3400 L=6000 K=2200 Q=HoxD55.q > lav/\`basename $i .2bit\`-\`basename \$j .2bit\`.lav;
done
EOF
chmod +x lastz_$id.sh
echo lastz_$id.sh '&' >> $batch
done
chmod +x $batch

if [ ! -f HoxD55.q ]; then
    echo "matrix scoring file HoxD55.q is missing"
    exit 0
fi

$batch
exit
```

**`work_steps_part2.sh`**
```
#!/bin/bash

WORKDIR=$PWD
tspec=$1
qspec=$2

#start this part 2 of script when part 1 finished.
#--------------------
# convert lav to psl
#---------------------
cd $WORKDIR/lav
for i in *.lav; do echo $i; lavToPsl $i `basename $i .lav`.psl; done; 
rm -f *.lav
#remove the line starts with '#', 
cat *.psl > $WORKDIR/all.psl.1
grep -v '^#' $WORKDIR/all.psl.1 > $WORKDIR/all.psl

mkdir -p $WORKDIR/psl
pslSplitOnTarget $WORKDIR/all.psl $WORKDIR/psl/ 

rm $WORKDIR/all.psl.1
#--------------------
# Chaining
#---------------------
mkdir -p $WORKDIR/chain

cd $WORKDIR
# The following step seems it complains .nib file prefix with species, to remove the prefix.

# for f in ${tspec}/*.nib ${qspec}/*.nib; do 
#     cp $f `dirname $f`/${f#*_};
# done
for i in psl/*.psl; do 
    echo $i; 
    axtChain $i ${tspec} ${qspec} chain/`basename $i .psl`.chain -linearGap=loose -psl; 
done

chainMergeSort chain/*.chain > all.chain

chainPreNet all.chain  ${tspec}.chrom.sizes ${qspec}.chrom.sizes all.pre.chain

#--------------------
# Netting
#---------------------
chainNet all.pre.chain -minSpace=1  ${tspec}.chrom.sizes ${qspec}.chrom.sizes stdout /dev/null | netSyntenic stdin noClass.net

#this step is taking ucsc genome browser database for extra information to display on browser, not mandatory. noClass.net can be used as ${qspec}.net.

#netClass -noAr noClass.net ${tspec} ${qspec} ${qspec}.net
cp noClass.net ${qspec}.net

#---------------------
# Mafing
#---------------------
netToAxt ${qspec}.net all.pre.chain ${tspec}/ ${qspec}/ stdout | axtSort stdin ${tspec}.${qspec}.axt

axtToMaf ${tspec}.${qspec}.axt ${tspec}.chrom.sizes ${qspec}.chrom.sizes ${tspec}.${qspec}.maf -tPrefix=${tspec}. -qPrefix=${qspec}.

```



As long as all the adaptations are changed, then simply run 

```bash
#usage: work_step_part1.sh tspec qspec
$ work_step_part1 ce10 cb3
#wait for part 1 to finish
$ work_step_part2 ce10 cb3
```


<br>
#### **Build Pairwise Whole Genome Alignment Approach 2** 
--- adapt the scripted method on a local machine and use  for _C. elegans-ce10_ and  _C. briggsae- cb3_

> **Date: 7/17/2015** 
> **work log: sandbox/align_howto5.2, 6.2, analysis/work_3a/align_ce10.caePb1, align_ce10.cb3** @usav1svBIOd1
align_howto5.2: test only chrI.
align_howto6.2: test whole genome.
align_ce10.caePb1: align whole genome for ce10 and cb3
align_ce10.caePb1: align whole genome for ce10 and caePb1

**OBJECTIVE:**
Run the whole genome alignment method and tools, following the one shell script: [Whole genome alignment howto](http://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt). Do a pairwise genome alignment of _C. elegans-ce10_ and  _C. briggsae- cb3_ 

**MATERIAL AND METHODS:**
Install mandatory tools: kenUtils, LASTZ, UCSC genome browser as instructed in approach 1. 

**Packages for the alignment:**
Perl scripts: (optional since KentUtils also includes these tools)
`constructLiftFile.pl`
`partitionSequence.pl`

Score matrix:
`human_chimp.v2.q`
`HoxD55.q`

**`modified core RunLastzChain.sh`**
```csh
#!/bin/sh

echo "example script to run lastz and chaining on two genomes in 2bit files"
echo "adjust this script for your local example, you have a couple of choices"
echo "for parameter sets.  You will need a parasol cluster computer system"
echo "to run the large number of lastz instances."
echo "requires companion script constructLiftFile.pl and"
echo "partitionSequence.pl"
echo 
echo "The point is to illustrate the steps of:"
echo "1. partitioning the two genomes into:"
echo "   a. 10,000,000 overlapping 10,000 chunks for the target sequence"
echo "   b. 20,000,000 no overlap chunks for the query sequence"
echo "2. setup cluster run target.list query.list lastz run script"
echo "3. chaining the psl results from the lastz procedure"

#exit 255

# typical axtChain and lastz parameter sets:
#export chainNear="-minScore=5000 -linearGap=medium"
#export chainMedium="-minScore=3000 -linearGap=medium"
#export chainFar="-minScore=5000 -linearGap=loose"
## To adapt.
#export lastzNear="B=0 C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 Q=human_chimp.v2.q T=2 Y=15000"
#export lastzMedium="B=0 C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400"
#export lastzMedium="B=2 C=0 E=30 H=2000 K=3000 L=3000 M=0 O=400 T=1 Y=9400"
## To adapt.
#export lastzFar="B=2 C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 Q=HoxD55.q T=2 Y=3400"

# select one of three different parameter sets
# Near == genomes close to each other
# Medium == genomes at middle distance from each other
# Far == genomes distant from each other

#export chainParams="$chainFar"
#to adapt
#export lastzParams="$lastzFar"
#export lastzParams="$lastzMedium"
#take an input file for align parameters
export lastzParams=`cat $3`
export chainParams=`cat $4`

#  WRKDIR is where your 2bit files are and where you want this to work
export WRKDIR=$PWD
cd ${WRKDIR}
## To adapt
export TNAME=$1
## To adapt
export QNAME=$2
export TARGET=${WRKDIR}/${TNAME}.2bit
export QUERY=${WRKDIR}/${QNAME}.2bit

ls -ld $TARGET $QUERY

if [ ! -s ${TNAME}.chrom.sizes ]; then
twoBitInfo ${TARGET} stdout | sort -k2nr > ${TNAME}.chrom.sizes
rm -fr ${TNAME}PartList ${TNAME}.part.list
mkdir ${TNAME}PartList
fi
if [ ! -s ${QNAME}.chrom.sizes ]; then
twoBitInfo ${QUERY} stdout | sort -k2nr > ${QNAME}.chrom.sizes
rm -fr ${QNAME}PartList ${QNAME}.part.list
mkdir ${QNAME}PartList
fi

if [ ! -s ${TNAME}.part.list ]; then
partitionSequence.pl 10000000 10000 ${TARGET} ${TNAME}.chrom.sizes 1 \
    -lstDir ${TNAME}PartList > ${TNAME}.part.list
fi
if [ ! -s ${QNAME}.part.list ]; then
partitionSequence.pl 20000000 0 ${QUERY} ${QNAME}.chrom.sizes 1 \
    -lstDir ${QNAME}PartList > ${QNAME}.part.list
fi

grep -v PartList ${TNAME}.part.list > target.list
for F in ${TNAME}PartList/*.lst
do
    cat ${F}
done >> target.list

grep -v PartList ${QNAME}.part.list > query.list
for F in ${QNAME}PartList/*.lst
do
    cat ${F}
done >> query.list

constructLiftFile.pl ${TNAME}.chrom.sizes target.list > target.lift
constructLiftFile.pl ${QNAME}.chrom.sizes query.list > query.lift

echo "#LOOP" > template
echo 'runLastz $(path1) $(path2) $(file1) $(file2) {check out exists+ psl/$(file1).$(file2).psl.gz}' >> template
echo "#ENDLOOP" >> template

cat <<_EOF_ > runLastz
#!/bin/csh -fe
set T = \$1
set Q = \$2
set FT = \$3
set FQ = \$4
## To adapt
set tmpDir = ./scratch/tmp/\${FT}
mkdir -p raw psl \${tmpDir}
twoBitToFa \${T} \${tmpDir}/\${FT}.fa
twoBitToFa \${Q} \${tmpDir}/\${FQ}.fa
## To adapt
lastz \${tmpDir}/\${FT}.fa \
    \${tmpDir}/\${FQ}.fa \
    ${lastzParams} \
    > raw/\${FT}.\${FQ}.lav
lavToPsl raw/\${FT}.\${FQ}.lav stdout \
    | liftUp -type=.psl stdout target.lift error stdin \
    | liftUp -nohead -pslQ -type=.psl stdout query.lift error stdin \
    | gzip -c > psl/\${FT}.\${FQ}.psl.gz
#rm -f \${tmpDir}/\${FT}.fa \${tmpDir}/\${FQ}.fa
#rmdir --ignore-fail-on-non-empty \${tmpDir}
_EOF_

echo "ready to run lastz kluster job:"
echo "gensub2 target.list query.list template jobList"
echo "para make jobList"

echo "when finished, run the commands in chainJobs.csh to perform the chaining"

mkdir -p chain
echo "#!/bin/csh -fe" > chainJobs.csh
for T in `cat target.list | sed -e "s#${WRKDIR}/##"`
do
echo "zcat psl/${T}.*.psl.gz \\"
echo "    | axtChain -psl -verbose=0 ${chainParams} \\"
echo -e "\tstdin ${TARGET} ${QUERY} stdout \\"
echo "   | chainAntiRepeat ${TARGET} ${QUERY} stdin chain/${T}.chain"
done >> chainJobs.csh

echo "find ./chain -name \"*.chain\" | chainMergeSort -inputList=stdin | gzip -c > ${TNAME}.${QNAME}.all.chain.gz" >> chainJobs.csh


```
**Notes of the adaptions:**
1. Comment out exit 255, we want the following code to be excuted.
2. Change the matrix score locations, since they are saved in the workdir, change Q option to  `Q=human_chimp.v2.q` and `Q=HoxD55.q` , remove the paths in the original script.
3. Change the target species to $1, and query speces to be  $2. `export TNAME=$1` and `export QNAME=$2`
4. Change the path of `lastz`, since it is installed and the path has been added into `$PATH`, there is no need for the path, or put absolute path in front of it.
5. change the tmpDir to current workDir, `set tmpDir =./scratch/tmp/\${FT}`
6. For lastz parameter, change B=2 for far species. This will allow both +/- strands of query sequence to be aligned.


**`align_part1.sh`**
```bash
#!/bin/bash
#Authored by Xiaoju (Ju) Zhang
#

WORKDIR=$PWD
GENOME=/home/xzhang/Projects/avalign/analysis/work_3a/genome
PACKAGE=/home/xzhang/Projects/avalign/analysis/work_3a/alignPackage
tspec=$1
qspec=$2
lastzParams=$3
chainParams=$4

#check inputs
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 tSpec qSpec lastzParams chainParams"
    exit 128 
fi

#----------------------------------------
# Prepare FASTA, 2bit, nib, and sizes file
#----------------------------------------


for spec in $tspec $qspec; do
    mkdir -p $spec
    echo "Preparing files for $spec ..."
    find $GENOME/$spec -name '*.fa*' | xargs cat > $WORKDIR/$spec.fa
    faToTwoBit $WORKDIR/$spec.fa $WORKDIR/$spec.2bit
    # split the fasta file into small pieces, 
    #+ also this step will name all fasta with .fa extension
    for f in $GENOME/$spec/*.fa*; do
        faSplit byname $f $WORKDIR/$spec/
    done


    # generate chromsome sizes file if it does not exist
    if [ ! -s $GENOME/$spec/$spec.chrom.sizes ]; then
        faSize -detailed $WORKDIR/$spec.fa > $WORKDIR/$spec.chrom.sizes
    else
        cp $GENOME/$spec/$spec.chrom.sizes $WORKDIR/$spec.chrom.sizes
    fi  

    # convert fast to .2bit and .nib for some following steps.
    for f in $WORKDIR/$spec/*.fa; do
        faToTwoBit $f ${f%fa}2bit
        faToNib $f ${f%fa}nib
        rm $f
    done
    rm $spec.fa
done


cp $PACKAGE/RunLastzChain.sh $WORKDIR
cp $PACKAGE/*.q $WORKDIR
#---------------------------------
# Run the job script
#----------------------------------
#define the two species and run the job.
RunLastzChain.sh $tspec $qspec $lastzParams $chainParams

# chmod +x is used for runLastZ and jobList, since it will be executed locally.
chmod +x runLastz 

#---------------------------------
# Genearte the jobList
#----------------------------------
gensub2 target.list query.list template jobList

# split the jobList to multiple jobs, based on the total counts of the jobs. One can type `cat jobList | wc -l` to get the total jobs. 
#get total job numbers
num_jobs=`cat jobList | wc -l`
#define the number of parallel jobs
num_cores=8
#compute number of jobs for each sub job script.
subjobs_size=$((num_jobs / num_cores + 1))

#use sed to get the fragment of jobs commands.
for f in `seq 1 $subjobs_size $num_jobs`; do 
    if [ "$((f + subjobs_size))" -le  "$num_jobs" ]; then 
        sed -n $f,$((f + subjobs_size - 1))p jobList > jobList${f}-$((f + subjobs_size - 1)); 
        chmod +x jobList${f}-$((f + subjobs_size -1));  
        jobList${f}-$((f + subjobs_size -1)) &
    else    
        sed -n $f,${num_jobs}p jobList > jobList${f}-${num_jobs};
        chmod +x jobList${f}-${num_jobs};
        jobList${f}-${num_jobs} &
    fi
done
exit 0
```


**`align_part2.sh`**
```bash
#!/bin/bash
#Authored by Xiaoju (Ju) Zhang

WORKDIR=$PWD
GENOME=/home/xzhang/Projects/avalign/analysis/work_3a/genome
tspec=$1
qspec=$2

####!!This section can only run after the previous 
#+ steps are done.
chmod +x chainJobs.csh && chainJobs.csh


################################################
# Netting, same steps with step by step method
################################################

chainPreNet ${tspec}.${qspec}.all.chain.gz  $WORKDIR/${tspec}.chrom.sizes $WORKDIR/${qspec}.chrom.sizes all.pre.chain
#--------------------
# Netting
#---------------------
chainNet all.pre.chain -minSpace=1  $WORKDIR/${tspec}.chrom.sizes $WORKDIR/${qspec}.chrom.sizes stdout /dev/null | netSyntenic stdin noClass.net


#this step is taking ucsc genome browser database for extra information to display on browser, not mandatory. noClass.net can be used as ${qspec}.net.
#netClass -noAr noClass.net ${tspec} ${qspec} ${qspec}.net

cp noClass.net ${qspec}.net

#---------------------
# Mafing
#---------------------
netToAxt ${qspec}.net all.pre.chain $WORKDIR/${tspec}/  $WORKDIR/${qspec}/ stdout | axtSort stdin ${tspec}.${qspec}.axt

axtToMaf ${tspec}.${qspec}.axt $WORKDIR/${tspec}.chrom.sizes $WORKDIR/${qspec}.chrom.sizes ${tspec}.${qspec}.maf -tPrefix=${tspec}. -qPrefix=${qspec}.

rm -rf $WORKDIR/scratch
```

Part 3 is still under testing. (This part seems not to affect the final result of MAF)
```bash


#---------------------
# Phastcons
#---------------------
#split the maf for better phastCons performance, if the ${tspec}-${qspec}.maf size is reasonable, 
mkdir -p maf
mafSplit -byTarget dummy.bed maf/ ${tspec}.${qspec}.maf  

#generate background model file, phyloFit.mod
#in the tutorial it uses the largest peace to get the background model, we can use the whole maf
phyloFit --tree "($tspec, (${tspec}, ${qspec}))" -i MAF ${tspec}.${qspec}.maf 

#
mkdir -p wig
mkdir -p mostCons
for i in maf/*.maf; do \
    x=`basename $i .maf`;
    phastCons --target-coverage 0.25 --expected-length 12 \
    --rho 0.4 --msa-format MAF $i phyloFit.mod \
    --seqname `cat $i | head -n 3  | tail -n 1 | tr -s ' ' | cut -f 2 -d ' ' | cut -d. -f2` \
    --most-conserved mostCons/$x.bed > wig/$x.wig; \
done

# create one piece .wig
```

As long as all the adaptations are changed, then simply run 

```bash
#usage: work_flow.sh tspec qspec
$ align_part1.sh ce10 cb3
#when part1 is done, then
$ align_part2.sh ce10 cb3
```




<br>
### SECTION 3: **Build Multiple Whole Genome Alignment** 

> **Date: 7/21/2015** 
> **work log: sandbox/align_howto7.1** @usav1svBIOd1

**OBJECTIVE:**
Test and run the multiple whole genome alignment, following a collective of resources. 

**MATERIAL AND METHODS:**
Mandatory tools: multiz
Following this discussion post: https://www.biostars.org/p/2882/
>I have a set of 10-12 very closely related chromosome sequences (from different strains) aligned to a "single" reference chromosome. Now I need to generate multiple sequence alignment of these without afftecting individual alignments to the reference. All that I need is to add relative inserts at respective sequence positions, so that I get a global alignment with respect to reference.
Pairwise alignments:


**RESULTS:**
1. Install multiz

	```bash
	$ wget http://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz
	$ tar -xzvf multiz-tba.012109.tar.gz
	$ mv multiz-tba.012109 ../Tools/multiz-tba
	$ cd  ../Tools/multiz-tba
	$ make
	```
2. Rename each species's genome fasta file to its species name first. The name needs to be specifically and restrictively in the following pattern, only replacing the species names.

	```bash
	$ cp ce10.fa ../tba_howto1/ce10
	$ cp cb3.fa ../tba_howto1/cb3
	$ cp caePb1.fa ../tba_howto1/caePb1
	$ cp ce10-cb3.maf ../tba_howto1/ce10.cb3.sing.maf
	$ cp ce10-caePb1.maf ../tba_howto1/ce10.caePb1.sing.maf
	```
3. Run tba, aligning multiple pairwise alignments.

	```bash
	$ tba "((ce10 cb3) caePb1)" *.*.maf tba.maf
	```
	In this case, it actually generates a two-way multiple alignment of cb3 and caePb1 taking ce10 as reference. The `tba.maf` is the final multiple alignment.

4. MAF project
Shift the reference in the alignment maf file. (ref:http://genomeview.org/manual/Preparing_whole_genome_alignments)

	```bash
	$ maf_project tba.maf cb3 > tba_project_cb3.maf
	```

<br>
### SECTION 4: **Understand Parameters**
There are two stages that need to set up parameters for optimal alignment. The first one is LASTZ aligning stage, and the second one is the chaining stage. Results alignment is sensitive to the parameters. To reference what parameters UCSC used for their multiple alignment dataset, go to these webpages. (http://genomewiki.ucsc.edu/index.php/Mm9_multiple_alignment)
or (http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters)



<br>
### SECTION 5: **Test run for hg19-turTru2**
pairwise align human and dolphin

Referencing to [Hg19_100way_Genome_size_statistics](http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics), hg19 chr counts is 93, and turTru2 is 240,901. The genomes are downloaded from UCSC database. Raw turTru2 fasta file is one single file, by using `faSplit` by name, the fasta file can be split into 240,901 pieces.

1. Prepare the genome.
	Download genomes from UCSC goldenPath
	
	```bash
	$ mkdir -p genome/hg19 genome/turTru2
	$ cd genome/turTru2 && wget http://hgdownload.soe.ucsc.edu/goldenPath/turTru2/bigZips/turTru2.fa.masked.gz
	$ cd genome/hg19 wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
	```
	Decompress genome,	Split the fasta file if it is in a single file.

	```bash
	$ faSplit byname turTru2.fa.masked
	$ rm turTru2.fa.masked
	```
2. Run `work_flow_part1.sh` .

	```bash
	$ work_flow_part1.sh hg19 turTru2
	```
