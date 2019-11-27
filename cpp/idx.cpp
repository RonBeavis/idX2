/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/
/*
This C++ version of idX is the second iteration of the algorithm.

The first version, written in Python, explored the possibilities of
using a JSON Lines formatted collection of peptide sequence information
as a practical method of identifying peptides from peptide MS/MS data.
In this approach, each line of the JSON-formatted peptide sequence information was
called a "kernel". An example of a kernel is

{"lv":0,"pm":1665672,"lb":"sp|A5A3E0|","u":23973,"h":23973,
"pre":"K","post":"T","beg":916,"end":928,"seq":"LCYVALDFEQEMA",
"ns":[35,11339,2578,0],
"bs":[113084,...,1576625],
"ys":[89048,...,1552588],
"mods":[["C",917,119004],["M",927,15995]]
}

where all of the masses are expressed in millidaltons.

The name "kernel" is used both to signify that these are the basic building blocks of the
identification process and to distinguish them from other indexing approaches that
use different variables and file formats.

The arrays "bs" and "ys" correspond to the b-and y-fragment neutral masses of the sequence ("seq")
with the modifications described in "mods" (residue, position & mass).
The "ns" array is the number of times this sequence has been observed in GPMDB 
in the +1,+2,+3 and +4 charge states "u" is a unique serial number for the kernel 
within a file and "h" is the "u" value for the first time a kernel is listed in the file 
(allows for sequence redundancy). The "lv" value "0" indicates that all of the masses
represent neutral (uncharged) molecular species.

The first version indexed all values according to the parent neutral masses of the 
querying set of spectra and recorded the kernels associated with those masses. 
Those kernels where then compared with the fragments in each spectrum and 
an identification was generated.

The current version uses a similar approach but uses a much more aggressive indexing method. 
The fragment species in the spectra are indexed as pairs with the parent mass 

(pm,fm1):(pm,fm2):...:(pm,fmN). 

Then, as the kernels are being read, each of their y- and b- fragments are paired with 
the corresponding kernel "pm" and they are checked against all of the spectrum pairs.
Only kernel pairs that match spectrum pairs are taken forward into the identification 
phase, significantly reducing the amount of memory & calculation required to generate 
a spectrum-to-peptide match.
*/

/*
C++ STL elements are used throughout the code.
The 8-byte variable types "int64_t" and "double" are used as much as possible.
"size_t" is used when necessary for STL compatibility.
*/

/*
As mentioned above, all masses used are in millidaltons, represented as "int64_t" 
variables. When floating point masses are converted to integers, the following method is used:

	double dm = 1234.546789; //floating point mass in daltons
	int64_t im = (int64_t)(0.5+dm*1000.0); //integer mass in millidaltons

*/

#include "pch.h" //included for compatibility with Visual C++ precompiled header methods
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <exception>
#include <chrono> //used to keep track of method timing
#include "parallel_hashmap/phmap.h" //fast maps and sets used for indexing
using namespace std; //namespace used throughout
using namespace std::chrono; //namespace only used in this file
#include "load_spectra.hpp" //defines the load_spectra object
#include "load_kernel.hpp" //defines the load_kernels object
#include "create_results.hpp" //defines the create_results object
#include "create_output.hpp" //defines the create_output object

//
//	Simple method to serve as a cross-platform check for the existence of a file
//
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}
//
//	main procedure for performing tasks
//
int main(int argc, char* argv[])	{
	// checks for command line arguments
	if(argc < 4)	{
		cout << "usage:\t>idX SPECTRA_FILE KERNEL_FILE OUTPUT_FILE (high|medium|low*) (max_spectra*)" << endl;
		return 0;
	}
	map<string,string> params; //used to store command line and other constant values
	string version = "idX, 2020.1";
	params["version"] = version;

	int64_t fragment_tolerance = 400; // default fragment mass tolerance
	try	{
		if(argc > 4 and strcmp(argv[4],"high") == 0)	{
			fragment_tolerance = 20;
		}
		else if (argc > 4 and strcmp(argv[4],"medium") == 0)	{
			fragment_tolerance = 50;
		}
	}
	catch (...)	{
		cout << "Error (idx:0020): exception thrown trying to assign fragment tolerance" << endl;
		return 1;
	}
	ostringstream strStream;
	try	{
		strStream.clear();
		strStream.str("");
		strStream << (long)fragment_tolerance;
		params["fragment tolerance"] = strStream.str();
	}
	catch (...)	{
		cout << "Error (idx:0021): exception thrown trying to assign fragment tolerance" << endl;
		return 1;
	}
	string spectrum_file = argv[1]; //file containing spectra in MGF format
    	if(!exists(spectrum_file))	{
		cout << "Error (idx:0001): spectrum file \"" << spectrum_file << "\" does not exist" << endl;
		return 2;
	}
	params["spectrum file"] = spectrum_file;

	string kernel_file = argv[2]; //file containing kernels in JSON Lines format
    	if(!exists(kernel_file))	{
		cout << "Error (idx:0002): kernel file \"" << kernel_file << "\" does not exist" << endl;
		return 2;
	}
	params["kernel file"] = kernel_file;

	string output_file = argv[3]; //file that will contain the identifications in TSV format
	params["output file"] = output_file;

	int64_t maximum_spectra = -1; //if not -1, determines the number of spectra to consider
	try	{
		if(argc > 5)	{
			maximum_spectra = atoi(argv[5]);
		}
		strStream.clear();
		strStream.str("");
		strStream << (long)maximum_spectra;
		params["maximum spectra"] = strStream.str();
	}
	catch (...)	{
		cout << "Error (idx:0022): exception thrown trying to assign maximum spectra" << endl;
		return 1;
	}
	int64_t parent_tolerance = 20; //parent ion mass tolerance is fixed at 20 mDa
	try	{
		strStream.clear();
		strStream.str("");
		strStream << (long)parent_tolerance;
		params["parent tolerance"] = strStream.str();
	}
	catch (...)	{
		cout << "Error (idx:0023): exception thrown trying to assign parent tolerance" << endl;
		return 1;
	}
	cout << "\nstart ...\nidX parameters" << "\n"; //output the interpreted command line values for logging
	if(maximum_spectra != -1)	{
		cout << "\t   max spectra: " << maximum_spectra << "\n";
	}
	else	{
		cout << "\t   max spectra: unlimited" << "\n";
	}
	cout << "\t  fragment tol: " << fragment_tolerance << " mDa" << endl;
	cout << "\t    parent tol: " << params["parent tolerance"] << " ppm" << endl;
	cout << "\t spectrum file: " << spectrum_file << endl;
	cout << "\t   kernel file: " << kernel_file << endl;
	cout << "\t   output file: " << output_file << endl;
	cout << "\t       version: " << version << endl;
	cout << "load & index spectra" << endl;
	cout.flush();
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //begin timing spectrum loading
	load_spectra ls; 
	try	{
		if(!ls.load(params))	{ // load spectra into the ls object
			cout << "Error (idx:0003): failed to load spectrum file \"" << spectrum_file << "\"" << endl;
			return 1;
		}
		if(maximum_spectra != -1)	{ // use the maximum_spectra value, if specified
			ls.set_max(maximum_spectra);
		} 
	}
	catch (...)	{
		cout << "Error (idx:0024): failed to load spectrum file \"" << spectrum_file << "\"" << endl;
		return 1;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now(); //end timing spectrum loading and report
	cout << "	   spectra = " << ls.spectra.size() << "\n";
	cout << "	spectra &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s" << endl;
	strStream.str("");
	strStream.clear();
	strStream << (long)ls.spectra.size();
	params["spectra"] = strStream.str();
	t1 = high_resolution_clock::now(); //begin timing kernel loading
	cout << "load & index kernel"  << "\n";
	cout.flush();
	kernels kindex; //object that will contain kernel information
	map<int64_t,int64_t> mindex; //object that will contain (pm,fmN) index
	load_kernel lk; //object that will load kernel information from the specified file
	try	{
		if(!lk.load(params,ls,kindex,mindex))	{ //load kernel information 
			cout << "Error (idx:0005): failed to load kernel file \"" << kernel_file << "\"" << endl;
			return 1;
		}
	}
	catch (...)	{
		cout << "Error (idx:0025): failed to load kernel file \"" << kernel_file << "\"\n";
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing kernel loading and report
	cout << "	   kernels = " << kindex.size() << "\n";
	cout << "	kernels &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s" << endl;
	t1 = high_resolution_clock::now(); //begin timing peptide-to-spectrum matching
	cout << "perform ids"  << "\n";
	cout.flush();
	create_results cr; //object that will contain match information
	try	{
		if(!cr.create(params,ls,kindex,mindex))	{ //create peptide-to-spectrum matches
			cout << "Error (idx:0006): failed to create results " << "\n";
			return 1;
		}
	}
	catch (...)	{
		cout << "Error (idx:0026): failed to load create results" << endl;
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing peptide-to-spectrum matching and report
	cout << "	   results = " << cr.size() << "\n";
	cout << "	results &Delta;T = " << duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s" << endl;
	cout << "create report"  << "\n";
	cout.flush();
	t1 = high_resolution_clock::now(); //begin timing output file creation
	try	{
		create_output co;
		if(!co.create(params,cr))	{ //create output file, based on the matches in the cr object
			cout << "Error (idx:0007): failed to create output " << endl;
			return 1;
		}
	}
	catch (...)	{
		cout << "Error (idx:0027): failed to load create output" << endl;
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing output file creation
	cout << "... done\n";
	return 0;
}

