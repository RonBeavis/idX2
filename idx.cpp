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
The 8-byte variable types "int32_t" and "double" are used as much as possible.
"size_t" is used when necessary for STL compatibility.
*/

/*
As mentioned above, all masses used are in millidaltons, represented as "int32_t" 
variables. When floating point masses are converted to integers, the following method is used:

	double dm = 1234.546789; //floating point mass in daltons
	int32_t im = (int32_t)(0.5+dm*1000.0); //integer mass in millidaltons

*/
#include "pch.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <exception>
#include <chrono> //used to keep track of method timing
#include "parallel_hashmap/phmap.h" //fast maps and sets used for indexing
using namespace std; //namespace used throughout
using namespace std::chrono; //namespace only used in this file
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp" //defines the load_kernels object
#include "load_spectra.hpp" //defines the load_spectra object
#include "create_results.hpp" //defines the create_results object
#include "create_output.hpp" //defines the create_output object
#ifdef MSVC
	#define _TIMESPEC_DEFINED
#endif

//
//	Simple method to serve as a cross-platform check for the existence of a file
//
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

int load_params(map<string,string>& params,int argc,char* argv[])	{
	params["version"] = "idX, 2020.1 (std)";
	params["fragmentation"] = "";
	int32_t fragment_tolerance = 400; // default fragment mass tolerance
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
	ostringstream strStream; //using ostringstream to avoid potentially unsafe sprintf
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

	int32_t maximum_spectra = -1; //if not -1, determines the number of spectra to consider
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
	int32_t parent_tolerance = 20; //parent ion mass tolerance is fixed at 20 mDa
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
	return 0;

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
	int ret = load_params(params,argc,argv);
	if(ret != 0)	{
		return ret;
	}
	int32_t maximum_spectra = atol(params["maximum spectra"].c_str()); //if not -1, determines the number of spectra to consider
	cout << "\nstart ...\nidX parameters" << "\n"; //output the interpreted command line values for logging
	if(maximum_spectra != -1)	{
		cout << "\t   max spectra: " << maximum_spectra << endl;
	}
	else	{
		cout << "\t   max spectra: unlimited" << "\n";
	}
	cout << "\t  fragment tol: " << params["fragment tolerance"] << " mDa" << endl;
	cout << "\t    parent tol: " << params["parent tolerance"] << " ppm" << endl;
	cout << "\t spectrum file: " << params["spectrum file"] << endl;
	cout << "\t   kernel file: " << params["kernel file"] << endl;
	cout << "\t   output file: " << params["output file"] << endl;
	cout << "\t       version: " << params["version"] << endl;
	cout << endl << "load & index spectra" << endl;
	cout.flush();
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //begin timing spectrum loading
	high_resolution_clock::time_point t_origin = t1;
	load_spectra ls; 
	load_kernel lk_main;
	string strK = params["kernel file"];
	bool isB = false;
	if (strK.find(".b") == strK.size()-2) {
		cout.flush();
		isB = true;
	}
	lk_main.kfile = strK; //initialize lk variable
	lk_main.fragment_tolerance = (double)atof(params["fragment tolerance"].c_str()); //initialize lk variable
	try	{
		if(!ls.load(params,lk_main))	{ // load spectra into the ls object
			cout << "Error (idx:0003): failed to load spectrum file \"" << params["spectrum file"] << "\"" << endl;
			return 1;
		}
		if(maximum_spectra != -1)	{ // use the maximum_spectra value, if specified
			ls.set_max(maximum_spectra);
		} 
	}
	catch (...)	{
		cout << "Error (idx:0024): failed to load spectrum file \"" << params["spectrum file"] << "\"" << endl;
		return 1;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now(); //end timing spectrum loading and report
	cout << endl << "  spectra = " << ls.spectra.size() << endl;
	cout << "  spectra &Delta;T = " 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << endl;
	ostringstream strStream; //using ostringstream to avoid potentially unsafe sprintf
	strStream.str("");
	strStream.clear();
	strStream << (long)ls.spectra.size();
	params["spectra"] = strStream.str();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, spectra load (s)"] = strStream.str();
	t1 = high_resolution_clock::now(); //begin timing kernel loading
	cout << endl << "load & index kernel"  << endl;
	cout.flush();

//	lk_main.spectrum_pairs(ls);
	try	{
		if (isB) {
			lk_main.load_binary();
		}
		else {
			lk_main.load();
		}
	}
	catch (...)	{
		cout << "Error (idx:0028): failed to load kernels" << endl;
		return 1;
	}

	t2 = high_resolution_clock::now(); //end timing kernel loading and report
	cout << endl << "  kernel pairs = " << lk_main.kerns.size() << endl;
	cout << "  kernels &Delta;T = "  << fixed << setprecision(3) 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << endl;
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, kernel load (s)"] = strStream.str();
	t1 = high_resolution_clock::now(); //begin timing peptide-to-spectrum matching
	cout << endl << "perform ids"  << endl;
	cout.flush();
	create_results cr; //object that will contain match information
	try	{
		if(!cr.create(params,ls,lk_main))	{ //create peptide-to-spectrum matches
			cout << "Error (idx:0006): failed to create results " << endl;
			return 1;
		}
	}
	catch (...)	{
		cout << "Error (idx:0026): failed to load create results" << endl;
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing peptide-to-spectrum matching and report
	cout << endl << "  results = " << cr.size() << endl;
	cout << "  results &Delta;T = "  << fixed << setprecision(3)
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << endl;
	cout << endl << "create models & report"  << endl;
	cout.flush();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, id process (s)"] = strStream.str();
	lk_main.clean_up();
	ls.clean_up();
	t1 = high_resolution_clock::now(); //begin timing output file creation
	create_output co;
	try	{
		if(isB)	{
			if(!co.create_binary(params,cr,lk_main.hu_set))	{ //create output file, based on the matches in the cr object
				cout << "Error (idx:0007): failed to create output " << endl;
				return 1;
			}
		}
		else	{
			if(!co.create(params,cr,lk_main.hu_set))	{ //create output file, based on the matches in the cr object
				cout << "Error (idx:0007): failed to create output " << endl;
				return 1;
			}
		}
	}
	catch (...)	{
		cout << "Error (idx:0027): failed to load create output" << endl;
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing output file creation
	cout << "  reporting &Delta;T = " << fixed << setprecision(3) 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s" << endl;
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, output file creation (s)"] = strStream.str();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t_origin).count()/1000.0;
	params["time, total (s)"] = strStream.str();
	co.dump_meta(params);
	cout << endl << "... done" << endl;
	return 0;
}

