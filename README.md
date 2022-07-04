# Binmarkers version 2

The is a new version for the script to create bin markers. Previous version: [BinMarkers](https://github.com/lileiting/BinMarkers)

* Author: [Leiting Li](https://github.com/lileiting)
* Email: lileiting@gmail.com
* LICENCE: [BSD](http://opensource.org/licenses/bsd-license.php)

## USAGE     

    USAGE
        binmarkers-v2.3.pl [Options] <input.markers.txt>
    
    Description
        Asssume the input data is a tab-seperated matrix,
        with first line as title [optional] and first column as marker name.
        Markers name format: scaffold_start:end:markers.
        For example, scaffold0001_10:10000:10 means a bin markers
            spanned scaffold0001, from 10 bp to 10000 bp, included
            10 original markers.
    
        Assume genotypes were encoded as a, A, b, B, h, H, `-`, u, `.`, `..`, `--`.
        A is equivalent to a; B is equivalent to b;
        H is equivalent to h; `u`, `--`, `.`, `..` are equivalent to `-`.
    
        This protocol required five main steps:
        1. Bin markers by 10 kb window using majority rules [-m 1 -w 10_000]
        2. Fill missing genotypes using majority rules:
            round 1: [-m 2 -w 3]
            round 2: [-m 2 -w 5]
            round 3: [-m 2 -w 7]
        3. Fill missing genotypes by in the breakpoints [-m 3 -w 3]
        4. Correct misscored genotypes with strict criteria [-m 4 -w 5]
        5. Merge 100% identical markers [-m 5]
        * Step 2 and step 4 will process the first and last [w] markers
    
        Majority rules ([w] markers above, [w] markers below):
            if two genotypes were equal weights, treat as missing data
    
        Strict criteria for correcting misscored genotypes:
            A genotype is different with surroundings genotypes;
            No missing data in surroundings genotypes;
            Surroundings genotypes are the same.
    
        Example commands:
            perl binmarkers-v2.3.pl example.input.markers.txt -m bin -w 10_000 -o out.step1.txt
            perl binmarkers-v2.3.pl out.step1.txt   -o out.step2.1.txt -m fill -w 3
            perl binmarkers-v2.3.pl out.step2.1.txt -o out.step2.2.txt -m fill -w 5
            perl binmarkers-v2.3.pl out.step2.2.txt -o out.step2.3.txt -m fill -w 7
            perl binmarkers-v2.3.pl out.step2.3.txt -o out.step3.txt   -m fill2 -w 3
            perl binmarkers-v2.3.pl out.step3.txt   -o out.step4.txt   -m correct -w 5
            perl binmarkers-v2.3.pl out.step4.txt   -o out.step5.txt   -m merge
    
        Or use pipe to run them in one command, such as
            perl binmarkers-v2.3.pl example.input.markers.txt -m bin -w 10_000 |
            perl binmarkers-v2.3.pl -m fill -w 3 |
            perl binmarkers-v2.3.pl -m fill -w 5 |
            perl binmarkers-v2.3.pl -m fill -w 7 |
            perl binmarkers-v2.3.pl -m fill2 -w 3 |
            perl binmarkers-v2.3.pl -m correct -w 5 |
            perl binmarkers-v2.3.pl -m merge -o out.step5.txt
    
        For manual checking, just copy data in the `example.input.markers.txt` file
            and the `out.step4.txt` file in an Excel file and then sort by
            marker names
    
    OPTIONS
        -t, --title       Treat first line as title [default: no title line]
        -w, --window NUM  Window size [defaults are based on differnt modes]
        -m, --mode NUM    Select modes:
                                bin: Bin markers for a larger block [default: -w 10_000]
                               fill: Fill missing data (majority rules) [default: -w 3]
                              fill2: Fill missing data in the breakpoints [default: -w 3]
                            correct: Correct misscored genotypes [default: -w 5]
                            correct2: Correct misscored genotypes [default: -w 5]
                            correct3: Correct misscored genotypes [default: -w 5]
                              merge: Merge adjacent 100% identical markers
        -o, --out FILE    Output file [default: stdout]
        --no_edge         Do not process the first and last [w] markers,
                          valid to -m 2 and -m 4
        --minimum NUM     Minimum block size for -m 2 and -m 4
                          [default: (w) * 2 + 1]
        -h, --help        Print help
    


## Citation

When using this program please cite the following paper:

Qin M-F, Li L-T, Singh J, Sun M-Y, Bai B, Li S-W, Ni J-P, Zhang J-Y, Zhang X, Wei W-L, Zhang M-Y, Li J-M, Qi K-J, Zhang S-L, Khan A, Wu J. [Construction of a high-density Bin-map and identification of the fruit quality related quantitative trait loci and functional genes in pear](https://doi.org/10.1093/hr/uhac141). _Horticulture Research_, 2022, uhac141. 

