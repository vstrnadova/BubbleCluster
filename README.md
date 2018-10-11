TO RUN BUBBLECLUSTER, VERSION 4:

Usage: bc datadirectory/inputfile <lodthreshold> <nonmissingthreshold> <lowqualc> <sigma> <clusterfilename>

Where the parameters are defined as follows:

-> inputfile: the file containing input data. The
input file should have line per marker, and be in 
the following format:
<MARKERNAME> <TAB> <GENOTYPES>
where GENOTYPES are either -1 or 1, or 0 for missing
data. The genotypes should also be tab-separated. 

-> lodthreshold: ("tau" in our paper) describes the 
LOD score that a marker must achieve with a 
representative point in order to join the cluster 
which the representative point belongs to

-> nonmissingthreshold: ("eta" in our paper)
the minimal number of non-missing entries that a marker must 
have in order to be included in Phase I, the sketch-
building phase of our algorithm. In other words, this 
limits the number of markers included in the high-quality
marker set based on the amount of missing data in each
marker

-> lowqualc: ("c" in our paper) the minimal difference
between the LOD score between a low-quality marker and 
the highest-scoring representative point and the LOD 
score between the marker and the second-highest scoring
representative point. If the difference is greater than
or equal to lowqualc, the marker gets assigned to the 
cluster containing the highest-scoring representative
point.

-> sigma: ("sigma" in our paper) the minimum size of a
"real-life" cluster. The program will attempt to merge
all cluster of size less than sigma to a cluster with
size greater than sigma

-> clusterfilename: this is just a label for the 
output files containing clusters and representative
points. The user may wish to include the species name 
or other identifying features in this filename  

----------------------------------------------------------------
The program generates one file per cluster, with one 
marker name per line, which for our simulated data is 
just "M[markernumber]", with file name:
        [clusterfilname]v4cluster[clusternumber]
The program also generates one file per cluster for 
representative points, with one representative point 
marker name per line, with file name:
        [clusterfilename]v4reppoints[clusternumber]
Finally, the program will also output a file with all
the markers that could not be placed in any cluster 
(based on the lowqualc threshold), with file name:
         singletons_[clusterfilename]v4

To get a log of useful output data, please redirect 
output to a file. For example:

~/mapdata/mapping/v4 ../markerstabseparated.txt 7.5 200 2 2000 lod7.5nm200 > out.txt

--------------------------------------------------------------------

