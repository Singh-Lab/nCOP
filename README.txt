# --------------------------------------------------------------------------- #
# nCOP - Network Coverage of Patients
# --------------------------------------------------------------------------- #
#
# Author: Borislav Hristov (borislav@cs.princeton.edu)
#
# I. Input 
#
# There are two required inputs:
# 1) a network file and 2) a mutational file containing genes with list of
# individuals that have variants in these genes.
#
# Additionally, the user may provide:
# 3) a weight file specifying a weight for each node in the network
# 4) a value for alpha (in which case the program skips the step selecting 
# alpha and uses the user specified value)
# 5) output prefix which is used in the beginning of the name of 
# the output file
#
# II. Output
#
# output_prefix_results.txt is written in the Outputs directory. The file
# contain a list of candidate genes ranked by how frequently they appear
# in the randomized runs.
#
# III. How to run
#
# 1. To run with basic inputs:
# ./run_nCOP network_file.txt mutational_file.txt 
#
# 2. If you want to specify any additional parameter add "param=value"
# ./run_nCOP network_file.txt mutational_file.txt weights=weights_file.txt alpha=0.5 output_prefix=my_output
#
# Note that nCOP is implemented in Ruby. You need an installed Ruby which can
# be simply done via: sudo apt-get install ruby-full
#
# IV. Input File Formats
#
# 1. Network file: each line specifies an edge, white space delimited:
# GENE_ID GENE_ID
#
# 2. Mutational file: each line is white space delimited
# GENE_ID INDIVIDUAL_1 INDIVIDUAL_3
# GENE_ID INDIVIDUAL_5
#
# 3. Weights file: 
# GENE_ID WEIGHT
# i.e TP53 0.59
#