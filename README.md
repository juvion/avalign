---
title: Variants in Alignment
author: Xiaoju (Ju) Zhang
---


AVAlign - Protocol for Building whole genome alignment
----
Edited / Authored by Xiaoju (Ju) Zhang
Updated: 8/12/2015

### SECTION 1: **Install Tools**

#### **MySQL Server:**
Introduction: MySQL Server is required by UCSC genome Browser mirror, and it may complain MySQL version incompatibility.  
http://askubuntu.com/questions/489815/cannot-install-mysql-server-5-5-the-following-packages-have-unmet-dependicies
http://www.rackspace.com/knowledge_center/article/installing-mysql-server-on-ubuntu  

1. Downgrade MySQL from 5.6 to 5.5, installing mysql-server  

	```bash
    $ sudo apt-get purge mysql-client-core-5.6
    $ sudo apt-get autoremove
    $ sudo apt-get autoclean
    $ sudo rm  -r /var/lib/mysql # this step is added based on an error message
    $ sudo apt-get install mysql-client-core-5.5
    $ sudo apt-get install mysql-server   
	```
	for ubuntu, replace `yum` with `apt-get`


#### **kentUtils**: 
Introduction: kentUtils is the key tool box which supports the majority of the steps to perform pairwise alignment.
Download: https://github.com/ENCODE-DCC/kentUtils 
1. Download the source code, decompress the package.
2. Install required packages
	
	```bash
	$ sudo yum install git libssl-dev openssl
	```
3. cd into the folder, run following commands to compile the program

	```bash
	$ git clone git://github.com/ENCODE-DCC/kentUtils.git 
	$ cd kentUtils
	$ make # install general tools
	$ cd src/hg #make hg utils
	$ make
	$ cd src/utils #make other utils
	$ make
	```
4. If during the aligning steps, some command is missing, firstly make sure the tool is in the executable path, and then check kenUtils has it or not. It might not be compiled by default, and you need to manually run the makefile for it.


#### **LASTZ**: 
Introduction: LASTZ is the tool run pairwise alignment
Download: http://www.bx.psu.edu/miller_lab/dist/lastz-1.02.00.tar.gz 
1. Download the source code, decompress the tarball package.
	
	```bash
	$ wget http://www.bx.psu.edu/miller_lab/dist/lastz-1.02.00.tar.gz
	$ tar -zxvf lastz-1.02.00.tar.gz
	```
	
2. Change the folder from `lastz-xxx.xx` to `lastz`; edit file `make‑include.mak` to define `installDir` which is your installation directory. Also, make sure to put the directory `installDir` that you choose to your $PATH in `~/.bash_profile`.

3. Run make to install	
	```bash
	$ cd lastz
	$ make
	$ make install
	```


#### **Tandem Repeat Finder:**
Introduction: Tandem Repeat Finder is tool to find [tandem repeat](https://en.wikipedia.org/wiki/Tandem_repeat) of the sequence, which is needed to be masked before running alignment.
Download:http://tandem.bu.edu/trf/trf.download.html
1. Download Tandem Repeat Finder and put it into executable path (for example /usr/local/bin), and name it as "trf", 
2. make it executable `chmod +x trf`.


#### **UCSC Genome Browser**: 
Introduction: The UCSC source tree includes tools for the command line and the CGIs for a local UCSC genome browser mirror. Details for installation: http://genomewiki.ucsc.edu/index.php/Source_tree_compilation_on_Debian/Ubuntu
Download installation script: https://github.com/maximilianh/browserInstall
This shell script makes installation much easier, and it also offers command line to download and import database you are referring. The details instruction can be found in the repo's README.
>The fastest way ever to get a genome browser up and running on Ubuntu, Fedora, Centos, OSX
1. Download the installing script
2. Run the script to install UCSC genome browser mirror or import database

	``` bash
	$ sudo -i
	$ wget https://raw.githubusercontent.com/maximilianh/browserInstall/master/browserInstall.sh
	$ sudo bash browserInstall.sh
	# use the same script to download and import assemblies database.
	$ sudo bash  browserInstall.sh hg19 ce10 cb3
	#this may complains about the size of disk if the partition of the UCSC genome browser and mysql folder are on root mount. Do following.
	#create those two folders in my working folder.
	$ cd /home/xzhang/Tools/browser
	$ mkdir gbdb
	$ mkdir mysql
	# create symlink 
	$ ln -s /gbdb    /home/xzhang/Tools/browser/gbdb
	$ ln -s /var/lib/mysql      /home/xzhang/Tools/browser/mysql
	# To download and import database of any species
	$ bash browserInstall.sh ce10 cb3
	```


#### **Phast**:
Introduction: This package contains phastCons to compute Phastcon scores.
Download: 
Clapack: http://www.netlib.org/clapack/clapack.tgz
Phast: http://compgen.bscb.cornell.edu/phast/


**Part 1 - Installing Clapack** - (If you already have Clapack installed, skip to Part 2)

1. Download Clapack from the following URL http://www.netlib.org/clapack/clapack.tgz
2. Decompress the package
3. Go to the newly created Clapack directory install the package
    
    ```bash
    $ wget http://www.netlib.org/clapack/clapack.tgz
	$ tar -xvzf clapack.tgz
	$ mv CLAPACK-3.2.1 CLAPACK && cd CLAPACK
	$ cp make.inc.example make.inc && make f2clib && make blaslib && make lib
    ```
    
**Part 2 - Installing Phast**  

4. Download a copy of Phast from http://compgen.bscb.cornell.edu/phast/
5. Decompress the package 
5. Change directory to `phast/src/` and run `make CLAPACKPATH=/usr/local/software/clapack` replacing `/usr/local/software/clapack` with the path of you installed Clapack (e.g., `CLAPACKPATH=/home/username/CLAPACK-3.2.1`)
6. The Phast binaries should be created in the `../bin/` directory
7. add `$HOME/phast/bin` to `$PATH` in `~/.bash_profile`
   
    ```bash
	$ wget http://compgen.bscb.cornell.edu/phast/downloads/phast.v1_3.tgz
	$ tar -xvzf phast.v1_3.tgz
	$ mv  phast-1.3  phast && cd phast/src
	$ make CLAPACKPATH=/home/xzhang/Tools/CLAPACK
    ```  
    
####  **MULTIZ:**   
Introduction: MULTIZ is the program perform multiple alignment based on multiple pairwise alignments referencing to the same target genome.
Download: http://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz  
1. Download the package 
2. Decompress the package
3. Run make to install

	```bash
	$ wget http://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz
	$ tar -xzvf multiz-tba.012109.tar.gz
	$ mv multiz-tba.012109 ../Tools/multiz-tba
	$ cd  ../Tools/multiz-tba
	$ make
	```


<br>

###  SECTION 2: **Pairwise alignment and multiple alignment**  

Using `ce10` and `cb3` two genomes aligning as example.
`C. elegans-ce10`: http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/
Mandatory:
`cd10.chromFaMasked.tar.gz`: http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/chromFaMasked.tar.gz
Optional:
`ce10.chrom.sizes`: http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.chrom.sizes
`C. briggsae- cb3`:  http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/
Mandatory:
`cb3.chromFaMasked.tar.gz`: http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/chromFaMasked.tar.gz
Optional:
`cb3.chrom.sizes`: http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/cb3.chrom.sizes

1. Prepare the settings of working directory. The directory names need to follow the rules restrictively. 

	>Folder `genome` contains the genome FASTA files, in each sub folder named as assemblie id. eg. hg19, ce10, etc  
	 Genome data need to be downloaded and decompressed saved in each species subfolder 
	
	>Folder `tools` contains the align pipeline scripts: `align_part1.sh` and `align_part2.sh`. 
	 `align_part1.sh` : aligning 
	 `align_part2.sh` : after aligning is done, run the chaining, netting and maffing
	 
	>Folder `alignPackage` contains all the other needed scripts, scoring matrix templates, parameter input templates.
	`*.pl`: perl scripts used to manipulate the sequence file and liftovers, these scripts are also available in kentUtils.
	`RunLastzChain.sh`: Core script to run LASTZ
	`lastz*.in`: templates for lastz parameters, considering near medium and far distance between the pairs of species.
	`chain*.in`: templates for chain parameters, considering near medium and far distance between the pairs of species.
	
	>
	
```text    

|-- genome
|   |- ce10
|   |   |-  *.fa
|   |- cb3
|   |   |-  *.fa
|   |- caePb1
|   |   |-  *.fa
|
|-- alignPckage
|   |-  RunLastzChain.sh
|   |-  *.pl
|   |-  lastz*.in
|   |-  chain*.in
|
|-- tools
|   |-  align_part1.sh
|   |-  align_part2.sh
|
|-- align_ce10.cb3
|   |-  lastzParams.in
|   |-  chainParams.in
|
|-- align_ce10.caePb1
|   |-  lastzParams.in
|   |-  chainParams.in
|
|-- multiAlign_ce10.cb3.caePb1

```


2. Go to aligning working directory  `align_ce10.cb3`, edit `lastzParams.in` and `chainParams.in`  setting up parameters. 
Template `lastzParams.in` and `chainParams.in`:
	`lastzParams.in`
	```
	B=2 C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400
	```
	`chainParams.in`
	```
	-minScore=1000 -linearGap=loose
	```
3. Reference to this website to choose the right parameters.
Choose parameters for `lastzParams.in`: Refer to [UCSC hg19 100way alignment lastz parameters](http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters)
Choose parameters for `chainParams.in`: Refer to [UCSC hg19 100way alignment size statistics](http://genomewiki.ucsc.edu/index.php/Hg19_100way_Genome_size_statistics)

4. Run `align_part1.sh `
	
	```bash
	$ cd align_ce10.cb3
	# Usage: align_part1.sh tspec qspec lastzParams.in chainParams.in
	$ ../tools/align_part1.sh ce10 cb3 lastzParams.in chainParams.in
	```
5.  After all the LASTZ jobs are done, run `align_part2.sh `, and this step's final output will be `ce10.cb3.maf`
	
	```bash
	$ cd align_ce10.cb3
	# Usage: align_part2.sh tspec qspec
	$ ../tools/align_part2.sh ce10 cb3
	```
6.  Run MULTIZ program, assuming another pairwise alignment is done, and `ce10.caePb1.maf` is ready. Rename each species’s genome fasta file (a whole genome sequence single file) to its species name first. The name needs to be specifically and restrictively in the following pattern, only replacing the species names. 

	```bash
	$ cp ./align_ce10.cb3/ce10.fa ./align_ce10.caePb1/ce10
	$ cp ./align_ce10.cb3/cb3.fa ./align_ce10.caePb1/cb3
	$ cp ./align_ce10.caePb1/caePb1.fa ./align_ce10.caePb1/caePb1
	$ cp ./align_ce10.cb3/ce10.cb3.maf ./align_ce10.caePb1/ce10.cb3.sing.maf
	$ cp ./align_ce10.caePb1/ce10.caePb1.maf ./align_ce10.caePb1/ce10.caePb1.sing.maf
	```

7. Run tba, aligning multiple pairwise alignments. This multiple alignment is using default and preliminary input. Referencing to UCSC [in-house script](http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob;f=src/hg/makeDb/doc/anoCar2.txt) for advanced parameters and inputs. 

	```bash
	$ cd multiAlign_ce10.cb3.caePb1
	$ tba "((ce10 cb3) caePb1)" *.*.maf tba.maf
	```
	In this case, it actually generates a two-way multiple alignment of cb3 and caePb1 taking ce10 as reference. The `tba.maf` is the final multiple alignment.

8. MAF project
Shift the reference in the alignment maf file. (ref:http://genomeview.org/manual/Preparing_whole_genome_alignments)

	```bash
	$ maf_project tba.maf cb3 > tba_project_cb3.maf
	```


<br>

### SECTION 3: **Scripts and Input Files**
**Packages for the alignment:**
**Perl scripts**
Optional since KentUtils also includes these tools
`constructLiftFile.pl`
`partitionSequence.pl`

**Score matrix**
`human_chimp.v2.q`
`HoxD55.q`
`*.q` files are scoring matrices and they can be retrieved from [here](http://genomewiki.ucsc.edu/index.php/Ce10_conservation_lastz_parameters). 


**Core unLastzChain.sh** script initially is retrieved from [here](http://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt). The original script has some local settings and [mistakes](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/BZjUyNgRKnM/gDQoqe9PEgAJ), and it is modified for the purpose of this local aligning work.
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
**Notes of the modifications:**
1. Comment out exit 255, we want the following code to be excuted.
2. Change the matrix score locations, since they are saved in the workdir, change Q option to  `Q=human_chimp.v2.q` and `Q=HoxD55.q` , remove the paths in the original script.
3. Change the target species to `$1`, and query speces to be  `$2`, `export TNAME=$1` and `export QNAME=$2`
4. Change the path of `lastz`, since it is installed and the path has been added into `$PATH`, there is no need for the path, or put absolute path in front of it.
5. change the tmpDir to current workDir, `set tmpDir =./scratch/tmp/\${FT}`
6. For lastz parameter, change B=2 for far species. This will allow both +/- strands of query sequence to be aligned.

<br>

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

<br>

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
