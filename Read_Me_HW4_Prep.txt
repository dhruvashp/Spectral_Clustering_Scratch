The final files on which all the relevant steps will be performed have been 
stored as X_2, X_5, s_2, s_5 extracted from aca2 and aca5.


The code, in HW4_Prep, actually stores these files with different names when 
the code is run.

The names of the files have been changed to ensure that if the code in HW4_Prep
is run again, files that are relevant, aren't overwritten.

The reason for this is that we've used the same code for both aca2 and aca5, 
importing them one after the other and running the same code to prepare the 
final data.



NOTE :
Even though running the HW4_Prep code again won't really affect anything (since
the output files have been renamed, and those renamed ones will be subsequently 
used), I advise not running them again, since they will need -
1. Importing aca2 or aca5 (X and s)
2. Changing the name of the writematrix file output based on which file has been 
   imported

These steps have already been performed with the file outputs stored in X_2, X_5, s_2, s_5 (renaming the output after relevant imports)

The writematrix command in HW4_Prep has files corresponding to aca5, since that
was done after aca2. When aca2 was imported, aca2 filenames were used.