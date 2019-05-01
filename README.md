## About this software
This software solves the [Seriation](http://www.jstatsoft.org/v25/i03) problem finding a suitable linear order for a set of proteins. The result is a list of proteins ordered in one dimension such that functionally associated proteins are closer.

## Authors
The software was developed by [Felipe Kuentzer](http://lattes.cnpq.br/1979213773480902), in collaboration with 
Douglas G. Ávila, Alexandre Pereira, Gabriel Perrone, Samoel da Silva, [Alexandre Amory](http://lattes.cnpq.br/2609000874577720), and [Rita de Almeida](http://lattes.cnpq.br/4672766298301524).

The version provided here was modified by [Clovis Ferreira dos Reis](http://lattes.cnpq.br/5487049518249525) to improve the textual feedback and to avoid bugs like:
* Duplication of identifiers on the ordering output.
* Segmentation fault while reading an input file containing many nodes.

## Download and Compilation
Compilation requires GCC. To compile this software invoke the following commands on the shell:
<pre>
wget https://github.com/arthurvinx/seriation/archive/master.zip
unzip master.zip
cd seriation-master/
gcc ordering1D.c -o ordering1D -lm
</pre>

## How to Use
type 'cfm-seriation' to show the options:

<pre>
cfm-seriation

Usage: cfm-seriation [OPTION...]

 Seriation Parameters:
   f=[NETWORK FILE]       Network file path name
   o=[ORDER FILE]         Apply initial order
   i=[INTERVAL]               Number of isothermal steps
   m=[STEPS]                  Number of steps
   c=[FACTOR]                 Cooling factor
   a=[ALPHA]                  Alpha value
   p=[PERCENTUAL]             Percentual energy for initial temperature
   s=[SEDD]                   Random seed
   P                          Plot graphs
   v                          Generate video
</pre>

type to execute the seriation. This process can take about 12 minutes, depending on the CPU.

> ./cfm-seriation f=data/Homo_sapiens.dat m=3000 P

In case you want a video of the process, type to execute the seriation. 

> ./cfm-seriation f=/usr/share/cfm-seriation/data/Homo_sapiens.dat m=3000 P v

This will consume some extra time.

<pre>
Usage: cfm-seriation [OPTION...]

 Seriation Parameters:
   f=[NETWORK FILE]       Network file path name
   o=[ORDER FILE]         Apply initial order
   i=[INTERVAL]               Number of isothermal steps
   m=[STEPS]                  Number of steps
   c=[FACTOR]                 Cooling factor
   a=[ALPHA]                  Alpha value
   p=[PERCENTUAL]             Percentual energy for initial temperature
   s=[SEDD]                   Random seed
   P                          Plot graphs
   v                          Generate video

Reading file...
	Proteins: 9684
	Interactions: 163509
Applying random order...
Saving and plotting initial order...
INITIAL Energy: 4123514310
Ordering...
100% [====================================================================================================]
</pre>

For a quicker test you can execute smaller dataset, like the Escherichia Coli.

> ./cfm-seriation f=Escherichia_coli.dat

<pre>
Reading file...
	Proteins: 3598
	Interactions: 13687
Applying random order...
Saving initial order...
INITIAL Energy: 129449102
Ordering...
100% [====================================================================================================]
FINAL Energy: 129025784
Saving final order...
Done!
</pre>

The results are save in the input file directory. 
## Inputs

The input file is a textual file describing an undirected network of nodes (in our examples the nodes are labeled by ENSEMBL Peptide IDs). Example:
<pre>
362663.ECP_0002	362663.ECP_4228
362663.ECP_0002	362663.ECP_0939
362663.ECP_4228 362663.ECP_0002
362663.ECP_0939 362663.ECP_0002
</pre>

This repository contains an example input file from Escherichia coli named [362663.protein.links.900.v11.0.txt](362663.protein.links.900.v11.0.txt).

## Outputs

The output is a text file with the order of the network nodes. Example:

<pre>
Protein	dim1
Z5822	0
Z5823	1
Z2911	2
Z2910	3
Z2909	4
Z4123	5
Z4124	6
Z3105	7
Z3106	8
...
</pre>

The following image represents the Homo sapiens network with a *random ordering*.

![initial](initial.png)

The next image represents the Homo sapiens network **'seriated'**.

![final](final.png)

## License
The source code is distributed under the terms of the GNU General Public License v3 [GPL](http://www.gnu.org/copyleft/gpl.html).

## How to Cite this software
If you are using this package on your research, please cite:

* [Kuentzer, F. A. et al. (2014). Optimization and analysis of seriation algorithm for ordering protein networks.
IEEE International Conference on Bioinformatics and Bioengineering, 231-237.](https://doi.org/10.1109/BIBE.2014.43)

## Similar softwares
* [Seriation R Package](http://www.jstatsoft.org/v25/i03), available at [CRAN](http://cran.r-project.org/web/packages/seriation/index.html).