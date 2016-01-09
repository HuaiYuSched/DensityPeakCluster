# DensityPeakCluster
An implementation for 'Clustering by fast search and find of density peaks' in science 2014 with gnuplot to visualize the result.
##Prerequisite
This program needs gnuplot to generate the picture of result. You may need to install it in you machine if you need this function.Also need g++ and make.
##How to use
The program use the data from this website http://cs.joensuu.fi/sipu/datasets/, which is stored in data directory. The data consist of three column, x, y and the classification.
So run "make run" in the program directory, and firstly the program request a path to the data file. You can type any of file in the data directory like data/flame.txt

  Then it will run the program , you need to assign the minrho and min delta to control the number of cluster. After all calculation, it store the result(The classification of the points) in Clusterfile. Then invoke gnuplot to draw the result and store the result in a jpeg file which named by the time you invoke it.

If you do not want to generate jpeg file by gnuplot, you can run the command "make srun"
