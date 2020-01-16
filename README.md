# idX2
Second version of the idX algorithm. Implemented in C++ 17, using the idX Python project as a starting point.

A proteomics PSM assignment engine replacing most of the math with operations involving indexed collections. Identification confidence scoring based on Fisher's noncentral hypergeometric distribution, in conjunction with the number of ions matched with a sequence fragmentation pattern and the fraction of the ion intensity corresponding to those matched ions.

Current version 2020.01 (std)

##Kernels:

idX was designed to use a particular type of peptide sequence index file, which is referred to as a "kernel" file to avoid confusion with other projects that use indexing strategies. Kernel files are prepared prior to run time, so they can be edited and validated prior to use. All masses in kernels are in rounded integer millidaltons, so 1000.1231 Da would be represented as 1000123 mDa. All masses correspond to the neutral mass, not the ion mass. All reference to specific residues are in protein coordinates. 

These kernel files contain 3 types of entries and they are written in JSON Lines format.

1. The first JSON line must be like:

{"format":"jsms 1.0","source":"some description","created":"2019-12-17 21:55:36.492350"}

2. The lines containing peptide sequence information must be like:

{"lv":0,<br />
   "pm":parent mass,<br />
   "lb":"protein accession",<br />
   "pre":"residue - 1",<br />
   "post":"residue N+1",<br />
   "beg":start of peptide,<br />
   "end":end of peptide,<br />
   "seq":"peptide sequence",<br />
   "ns":[list of historical observations],<br />
   "bs":[list of b-ion masses],<br />
   "ys":[list of y-ion masses],<br />
   "mods":[lit of modifications [residue,pos,mass]],<br />
   "u":uid for line,<br />
   "h":uid for peptide}<br />

3. The last line must be like:

  {"validation":"sha256","value":"a21cd37e5ea05e36ddd264d0ede9abe108274ee7e58ad35dc08a468ff23d60bd"}

where the "value" is the SHA256 hash of all of the preceeding lines, with leading and trailing white splace removed.
  
##Design:

The software runs using the `main()` method in `idx.cpp` to control the workflow of PSM identification. This method is also responsible for most of the logging text output that is directed to stdout. The order of operations in this method are as follows:

1. Command line parameters are checked and stored. There are 4 required parameters:
- input spectrum file path (SPECTRUM_PATH)
- input kernel file path  (KERNEL_PATH)
- output file path (OUTPUT_PATH)
- fragment ion accuracy ("high" = 20 mDa, "medium" = 50 mDa and "low" = 400 mDa)

2. A `load_kernel` object is created and passed to a `load_spectra` object, that creates a list of spectra in the `load_spectra` object and an index of those spectra in the `load_kernel` object.

3. The `load_kernel` object is then loaded with kernel information from the kernel file specified on the command line.

4. A `create_results` object is created and it uses the `load_kernel` and `load_spectra` objects to generate a preliminary set of PSM assignments.

5. A `create_output` object is created which takes the preliminary list of PSM assignments and applies physical and statistical models to the list, generating a final list of PSM assignments that is recorded in a tab-separated value output file as specified on the command line. An additional metadata file in JSON format is also generated and stored at "OUTPUT_PATH.meta".

NOTE: As much as possible, neutral masses in integer millidaltons (type `int32_t`) are used throught the idX code, with as few variable type conversions as possible. C++ STL object indexes and sizes use type `size_t` as often as possible.
