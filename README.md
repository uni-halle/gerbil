# Gerbil: A fast and memory-efficient k-mer counter with GPU-support

A basic task in bioinformatics is the counting of k-mers in
genome strings. The k-mer counting problem is to build a histogram of
all substrings of length k in a given genome sequence. 
Gerbil is a k-mer counter that is specialized for high effiency when counting k-mers for large k.

Developed by: Marius Erbert, Steffen Rechner, Martin Luther University Halle-Wittenberg

## Install

Gerbil is developed and tested at Linux operating systems. Migrating it to other OS like Windows is a current issue. It follows a description of the installation process at Ubuntu 16.04.

1. Install 3rd-party libraries and neccessary software:

        sudo apt-get install git g++ libboost-all-dev libz3-dev

2. Download the Source Files. 

        git clone https://github.com/uni-halle/gerbil.git
        cd gerbil
        
3. Compile the Sources. Gerbil comes with a Makefile, that has been pre-configured to work in most cases. 
  
  a) To compile gerbil with CUDA support, run
        
        make
  
  b) At systems without CUDA-capable GPU's, Gerbil can be used without CUDA support:

        make GPU=false

Problems while compiling the code are possible due to different locations of libraries or the `nvcc` Compiler. In these cases, you have to modify the first few lines of the Makefile.

## Usage

```gerbil [option|flag]* input temp output```

Gerbil can be controlled by several command line options and flags.

| Option                  | Description   | Default |
|:------------------------|:--------------| -------:|
| `-k int`        | Set the value of k, i.e. the length of k-mers to be counted. Supported k currently ranges from 8 to 300. | 28 |
| `-m int`       | Set the length m of minimizers.      |   auto |
| `-e int(MB|GB)` | Restrict the maximal size of main memory in `MB` or `GB` that Gerbil is allowed to use.      |    auto |
| `-f int` | Set the number of temporary files.      |    512 |
| `-t int` | Set the maximal number of parallel threads to use.      |    auto |
| `-l int` | Set the minimal occurrence of a k-mer to be outputted.      |    3 |
| `-i` | Enable additional output.      |    |
| `-g` | Enable GPU mode. Gerbil will automatically detect CUDA-capable devices and will use them for counting in the second phase.      |     |
| `-v` | Show version number.      |     |
| `-d` | Disable normalization of k-mers. If normalization is disabled, a k-mer and its reverse complement are considered as different k-mers. If normalization is enabled, we map both k-mer and its reverse complement to the same k-mer.       |     |
| `-s` | Perform a system check and display information about your system.     |     |
| `-x 1` | Stop execution after Phase One. Do not remove temporary files and `binStatFile` (with statistical information). Watch out: no `output` allowed. |     |
| `-x 2` | Only execute Phase Two. Requires temporary files and `binStatFile`. No `input` allowed. |     |
| `-x b` | Do not remove `binStatFile`. |     |
| `-x h` | Create a histogram of k-mers in a human readable format in output directory. |     |

## Input Formats

Gerbil supports the following input formats of genome read data in raw and compressed format: 
 * fastq
 * fasta
 * staden
 * Textfile containing paths to one or more input files

## Output Format

Gerbil uses an output format that is easy to parse and requires little space. The counter of each occuring k-mer is stored in binary form, followed by the corresponding byte-encoded k-mer. Each four bases of a k-mer are encoded in one single byte. We encode `A` with `00`, `C` with `01`, `G` with `10` and `T` with `11`. Most counters of k-meres are slightly smaller than the coverage of the genome data. We exploit this property by using only one byte for counters less than 255. A counter greater than or equal to 255 is encoded in five bytes. In the latter case, all bits of the first byte are set to 1. The remaining four bytes contain the counter in a conventional 32-bit unsigned integer.

Examples (`X` means undefined):

| Counter | k-mer   | Encoding                      |
|:--------|:--------|:------------------------------|
| 67      | AACGTG  | `01000011 00000110 1110XXXX` |
| 345     | TGGATC  | `11111111 00000000 00000000 00000001 01011001 11101000 1101XXXX` |

