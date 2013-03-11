i) About LASAGNA
This is a small note describing the installation of the code LASAGNA.
LASAGNA was developed by Thomas Tram and Rasmus Sloth Hansen for solving
the Quantum Kinetic Equations for sterile neutrinos in the early Universe.
However, one can also see LASAGNA as a generic differential equation 
solver suited for large, stiff systems. If you use this code in a
scientific publication, please cite the paper

arXiv:1302.7279
"Can active-sterile neutrino oscillations lead to chaotic behavior of 
the cosmological lepton asymmetry?"

The appendix of this paper also describes some of the algorithms involved,
although not in great detail. If you depend on SuperLU for obtaining your 
results, you must also cite the two SuperLU_MT_2.0 papers. 
(See arXiv:1302.7279 or the SuperLU homepage for the references)

The overall structure of the code, some of the macros and the parser for
reading input files are adopted (with permission) from Julien Lesgourgues' 
Boltzmann code CLASS (http://class-code.net/).

ii) Installation
If you do not want to use SuperLU:
Search for "use_superlu" in the Makefile and set it to "no".
Then, choose the C-compiler and corresponding compiler flags of your 
choice.
write make and you are done
You can run the code by writing
./lasagna parameters.ini

If you plan to use SuperLU:
Download SuperLU_MT_2.0 from 
http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt
and extract it to some directory of your choice.
Copy the file "SuperLUpatch.tar.gz" from the LASAGNA root directory
to the SuperLU root directory and extract it. It overwrites a few files,
most notably "make.inc" in the root of the SuperLU directory. If you use
Intel compiler suite and MKL, you should not need to change this file,
otherwise you might need to.
write make to compile the SuperLU_MT library. Beware that on some systems
the SuperLU tests generate segmentation faults when they are executed by
the Makefile. These faults can be ignored, the library works well. (If it
makes you nervous, execute the test scripts manually after make is complete.)

Now edit the Makefile in the LASAGNA directory:
Search for "use_superlu" in the Makefile and set it to "yes",
and write the correct path for the SuperLU directory.
write make and you are done
You can run the code by writing
./lasagna parameters.ini

iii) Trouble shooting
SuperLU contains a single file which is a symbolic link,
/SuperLU_MT_2.0/TESTING/MATGEN/slu_mt_Cnames.h
If you are on a filesystem that does not support symbolic links
(For instance VirtualBox shared folder), you have to make a hard copy
of this file from /SuperLU_MT_2.0/SRC/slu_mt_Cnames.h.

The maximal amount of memory used by SuperLU must be set manually, and 
a too low or a too high value might result in a segmentation fault.
The amount of allocated memory can be increased by changing case 6, 7 
and 8 in SRC/sp_ienv.c on line 114 to 116. A value of (-100) was found
appropriate for vres=3200 while (-200) was needed for vres=6400.
Remember to recompile SuperLU after applying the changes.

iv) Bug reports
Please write any bug reports and comments to
rshansen@phys.au.dk or
tram@phys.au.dk
also feel free to ask about how to use the code, since a user guide 
is obviously missing. (If we get lots of emails we will probably be 
motivated for writing one.)
