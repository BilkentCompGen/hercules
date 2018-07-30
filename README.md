# Hercules: a profile HMM-based hybrid error correction algorithm for long reads

## Installing Hercules

* Make sure you have a compiler that has support for C++14.
* Download the code from its GitHub repository.

```bash
git clone https://github.com/BilkentCompGen/hercules.git
```
*  Change directory to `hercules/src/` and run the Makefile. If everything goes well, you will have a binary called `hercules` inside the `bin` folder.

```bash
cd hercules/src/
make
cd ../bin/
```
Now you can copy this binary wherever you want (preferably under a directory that is included in your `$PATH`). Assuming that you are in the directory that the binary exists, you may run the command below to display the help message.

```bash
./hercules -h
```

## Running Preprocessing step

To display the help message for the preprocessing step, you may run:

```bash
./hercules -1 -h
```
Assume that you have paired-end short reads and a long read (`short_1.fastq`, `short_2.fastq`, `long.fasta`). Then you may simply run:

```bash
mkdir preprocessing
./hercules -1 -li long.fasta -si short_1.fastq -si short_2.fastq -o preprocessing/
```
Note that the output folder `preprocessing` must exists prior to run the following command. The output of command above will give you the necessary information to proceed for the next steps until the correction step. You should just simply align compressed short reads to the compressed long reads (i.e. both are located in the output folder `preprocessing`). You must also sort them and preferably remove the duplicates. If you have `bowtie2` installed in one of your `$PATH` directories, then you may simply run:

```bash
../utils/runBowtieRmDup.sh preprocessing/compressed_long.fasta preprocessing/compressed_short.fasta bowtie 30
```
This will run `bowtie2` in `30` threads to align compressed short reads to the compressed long reads. The resulting alignment file will be stored in the directory `bowtie` with a file name `alignment.bam`. Note that this file will already be sorted and its duplicates will be removed. You do not need to run `afteralignment.sh` after `runBowtieRmDup.sh`. However, if you want to use another aligner without sorting its output file, then you must call `afteralignment.sh` to sort and remove its duplicates unless you want to do it by yourself:

```bash
../utils/afteralignment.sh alignment.bam output_alignment.bam 30 8G
```
Resulting alignment file will be `output_alignment.bam`. Note that the command above will use `30` threads and `8G` of your memory while sorting.

## Correction step

To get information about the parameters for the preprocessing step, you may run:

```bash
./hercules -2 -h
```
Assume that you have your alignment file `alignment.bam`, original long reads `long.fasta`, short reads (uncompressed, generated during preprocessing step) `preprocessing/short.fasta` and you would like to store corrected reads inside `corrected_long.fasta`. The command below will use `30` threads while correcting the original long reads:

```bash
./hercules -2 -li long.fasta -ai alignment.bam -si preprocessing/short.fasta -t 30 -o corrected_long.fasta
```
Resulting fasta file `corrected_long.fasta` will be the final output of Hercules.
