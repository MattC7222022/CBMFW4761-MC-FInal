# CBMFW4761-MC-Final
Running the test set:

Download the entire folder.  Open the terminal and navigate to the downloaded folder. Run this command to compile the Cython file into C:
python hashtable_setup.py build_ext –inplace
This should create a .so file in the folder. Open MC_Final_Test_Script.py. Change the string of the PATH variable to the filepath of your downloaded folder. Run the script. You should end up with the species list in SRR32477976_test_hits.xlsx.

General Usage:

Download the entire folder and follow the instructions above to compile the Cyton. Open MC_Final_Test_Script.py. Change the string of the PATH variable to the filepath of your downloaded folder. Then, change ‘SRR32477976_test.fasta’ to whatever your fasta file is. It should be decompressed and in the folder.
The main script for preprocessing raw 16s rRNA sequences, filtering, generating graphs, and generating novel sequences is MC_CompGen_Final.py.
