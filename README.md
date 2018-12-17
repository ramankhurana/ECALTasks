# ECALTasks

## To print the masked maps use: 

- python PrintMaskedTower.py

before running the macro do change the input file parameter, inputrootfilename

The rootfile can be obtained by runing the script run.sh which eventually run the L1Prefiring.C and fills the 2D occupancy maps of TT. 

once the text file with maps are ready, you can plug in this new text file in the L1Prefiring.C and run it again for the case 3 results when only an emulated TP exists. 

