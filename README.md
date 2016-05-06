# Gerbil: A fast and memory-efficient k-mer counter with GPU-support

A basic task in bioinformatics is the counting of k-mers in
genome strings. The k-mer counting problem is to build a histogram of
all substrings of length k in a given genome sequence. 
Gerbil is a k-mer counter that is specialized for high effiency when counting large k-mers for k.

For a description of its functionality, see ...

## Install

Gerbil is developed and tested at Linux operating systems. Migrating it to other OS like Windows is a current issue. It follows a description of the installation process at Ubuntu 16.04.

1. Install 3rd-party libraries and neccessary software:

        sudo apt-get install git g++ libboost-all-dev

2. Download the Source Files. 

        git clone https://github.com/uni-halle/gerbil.git
        cd gerbil
        
3. Compile the Sources. Gerbil comes with a Makefile, that has been pre-configured to work in most cases. 
  
  a) To compile gerbil with CUDA support, run
        
        make
  
  b) At systems without CUDA-capable GPU's, Gerbil can be used without CUDA support:

        make GPU=false


## Usage

Gerbil can be controlled by several command line options and flags.

| Option                  | Description   | Default |
|:------------------------|:--------------| -------:|
| `-k <length>`        | Set the value of k, i.e. the length of k-mers to be counted. Supported k currently ranges from 8 to 200. Larger k can easily be activated. | 28 |
| `-m <length>`       | Set the length m of minimizers.      |   7 |
| `-e <size>(MB|GB)` | Restrict the maximal size of main memory in `MB` or `GB` that Gerbil is allowed to use.      |    auto |
| `-f <number>` | Set the number of temporary files.      |    512 |
| `-t <number>` | Set the maximal number of parallel threads to use.      |    auto |
| `-l <number>` | Set the minimal occurrence of a k-mer to be outputted.      |    3 |
| `-i` | Enable additional output.      |    |
| `-g` | Enable GPU mode. Gerbil will automatically detect CUDA-capable devices and will use them for counting in the second phase.      |     |
| `-v` | Show version number.      |     |
| `-d` | Disable normalization of k-mers. If normalization is disabled, a k-mer and its reverse complement are considered as different k-mers. If normalization is enabled, we map both k-mer and its reverse complement to the same k-mer.       |     |
| `-s` | Perform a system check and display information about your system.     |     |
