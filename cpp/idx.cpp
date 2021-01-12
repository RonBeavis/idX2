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
/*
The text output to cout is meant to serve as logging information for the process. The
suggested use for this text is to redirect it to files on the command line, something like:

>idx ... 1>session.log 2>error.log

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
#include <ctime>
#include "parallel_hashmap/phmap.h" //fast maps and sets used for indexing
using namespace std; //namespace used throughout
using namespace std::chrono; //namespace only used in this file
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp" //defines the load_kernels object
#include "load_spectra.hpp" //defines the load_spectra object
#include "create_results.hpp" //defines the create_results object
#include "create_output.hpp" //defines the create_output object
#ifdef MSVC //define MSVC if you are compiling with Visual Studio
	#define _TIMESPEC_DEFINED
#endif

// Simple method to serve as a cross-platform check for the existence of a file

inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void print_help(void)	{
	cout << "\nusage:\t>idX -s FILE -k FILE -o FILE -f (low) -m (-1) -p (20)" << "\n";
	cout << "\n1. Required parameters.\n";
	cout << "-s FILE           spectrum peak list file path\n";
	cout << "                  FILE:  MGF or CMN format peak list\n";
	cout << "-k FILE           kernel file path\n";
	cout << "                  FILE:  JSON lines format kernel list\n";
	cout << "-o FILE           output file path\n";
	cout << "                  FILE:  TSV format result file\n";
	cout << "\n2. Optional parameters.\n";
	cout << "-f VALUE          fragment mass tolerance\n";
	cout << "                  VALUE: low, medium or high\n";
	cout << "                  default: low\n";
	cout << "-m VALUE          maximum number of spectra to use\n";
	cout << "                  VALUE:  positive integer or -1 if use all\n";
	cout << "                  default: -1\n";
	cout << "-p VALUE          parent mass tolerance (ppm)\n";
	cout << "                  VALUE: positive integer\n";
	cout << "                  default: 20\n";
	cout << "\n3. Other parameters.\n";
	cout << "-v, --version     displays current version of idX\n";
	cout << "-h, --help        displays this page\n";
	return;
}

// takes command line values and processes them into a list of process
// parameters, stored in a map<string,string> object

int load_params(map<string,string>& params,int argc,char* argv[])	{
	params["version"] = "idX, 2021.1 (std)";
	params["fragmentation"] = "";
	params["parent tolerance"] = "20";
	params["fragment tolerance"] = "300";
	params["maximum spectra"] = "-1";
	params["spectrum file"] = "";
	params["kernel file"] = "";
	params["output file"] = "";
	char flag[16] = "";
	if(argc > 1)	{
		if(strcmp(argv[1],"--version") == 0 || strcmp(argv[1],"-v") == 0)	{
			cout << params["version"] << "\n";
			exit(0);
		} 
		if(strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-h") == 0)	{
			print_help();
			exit(0);
		} 
	}
	if(argc == 0)	{
		print_help();
		exit(0);
	}
	char *pvalue = NULL;
	ostringstream strStream; //using ostringstream to avoid potentially unsafe sprintf
	for(int i = 1; i < argc; i++)	{
		if(argv[i][0] == '-')	{
			flag[0] = argv[i][0];
			flag[1] = argv[i][1];
			flag[2] = '\0';
			continue;
		}
		pvalue = argv[i];
		if(strcmp(flag,"-f") == 0)	{
			int32_t fragment_tolerance = 300; // default fragment mass tolerance
			try	{
				if(strcmp(pvalue,"high") == 0)	{
					fragment_tolerance = 20;
				}
				else if (strcmp(pvalue,"medium") == 0)	{
					fragment_tolerance = 50;
				}
			}
			catch (...)	{
				cout << "Error (idx:0020): exception thrown trying to assign fragment tolerance" << '\n';
				return 1;
			}
			try	{
				strStream.clear();
				strStream.str("");
				strStream << (long)fragment_tolerance;
				params["fragment tolerance"] = strStream.str();
			}
			catch (...)	{
				cout << "Error (idx:0021): exception thrown trying to assign fragment tolerance" << '\n';
				return 1;
			}
		}
		if(strcmp(flag,"-s") == 0)	{
			string spectrum_file = pvalue; //file containing spectra in MGF format
		    	if(!exists(spectrum_file))	{
				cout << "Error (idx:0001): spectrum file \"" << spectrum_file << "\" does not exist" << '\n';
				return 2;
			}
			params["spectrum file"] = spectrum_file;
		}
		if(strcmp(flag,"-k") == 0)	{
			string kernel_file = pvalue; //file containing kernels in JSON Lines format
		    	if(!exists(kernel_file))	{
				cout << "Error (idx:0002): kernel file \"" << kernel_file << "\" does not exist" << '\n';
				return 2;
			}
			params["kernel file"] = kernel_file;
		}
		if(strcmp(flag,"-o") == 0)	{
			string output_file = pvalue; //file that will contain the identifications in TSV format
			params["output file"] = output_file;
		}
		if(strcmp(flag,"-m") == 0)	{
			try	{
				strStream.clear();
				strStream.str("");
				strStream << atoi(pvalue);
				params["maximum spectra"] = strStream.str();
			}
			catch (...)	{
				cout << "Error (idx:0022): exception thrown trying to assign maximum spectra" << '\n';
				return 1;
			}
		}
		if(strcmp(flag,"-p") == 0)	{
			try	{
				strStream.clear();
				strStream.str("");
				strStream << atoi(pvalue);
				params["parent tolerance"] = strStream.str();
			}
			catch (...)	{
				cout << "Error (idx:0023): exception thrown trying to assign parent tolerance" << '\n';
				return 1;
			}
		}
	}
	return 0;

}

// main procedure for performing tasks

int main(int argc, char* argv[])	{
	// checks for command line arguments
	if(argc < 2)	{
		print_help();
		return 0;
	}
	map<string,string> params; //used to store command line and other constant values
	int ret = load_params(params,argc,argv);
	if(ret != 0)	{
		return ret;
	}
	int32_t maximum_spectra = atol(params["maximum spectra"].c_str()); //if not -1, determines the number of spectra to consider
	cout << "started ..." << '\n';
	ostringstream strStream; //using ostringstream to avoid potentially unsafe sprintf
	strStream.str("");
	strStream.clear();
	// get current time information for logging into the .meta file
	auto now = std::chrono::system_clock::now();
    	std::time_t end_time = std::chrono::system_clock::to_time_t(now);
	strStream << std::ctime(&end_time);
	string strTemp = strStream.str();
	// deal with annoying trailing whitespace
	strTemp.erase(std::remove(strTemp.begin(), strTemp.end(), '\n'), strTemp.end());
	strTemp.erase(std::remove(strTemp.begin(), strTemp.end(), '\r'), strTemp.end());
	cout << "start time: " << strTemp << '\n';
	params["time, started"] = strTemp;
	// record the epoch timestamp too
	strStream.str("");
	strStream.clear();
	strStream << fixed << setprecision(3);
	strStream << std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count()/1000.0;
	params["time, started timestamp"] = strStream.str();
	//output the interpreted command line values for logging
	cout << "idX parameters" << "\n"; 
	if(maximum_spectra != -1)	{
		cout << "\t   max spectra: " << maximum_spectra << '\n';
	}
	else	{
		cout << "\t   max spectra: unlimited" << "\n";
	}
	cout << "\t  fragment tol: " << params["fragment tolerance"] << " mDa" << '\n';
	cout << "\t    parent tol: " << params["parent tolerance"] << " ppm" << '\n';
	cout << "\t spectrum file: " << params["spectrum file"] << '\n';
	cout << "\t   kernel file: " << params["kernel file"] << '\n';
	cout << "\t   output file: " << params["output file"] << '\n';
	cout << "\t       version: " << params["version"] << '\n';
	cout << '\n' << "load & index spectra" << '\n';
	cout.flush();
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //begin timing spectrum loading
	high_resolution_clock::time_point t_origin = t1;
	load_spectra ls; 
	load_kernel lk;
	// set up the load_kernel object
	string strK = params["kernel file"];
	bool isB = false;
	if (strK.find(".b") == strK.size()-2) {
		cout.flush();
		isB = true;
	}
	lk.kfile = strK; //initialize lk variable
	lk.fragment_tolerance = (double)atof(params["fragment tolerance"].c_str()); //initialize lk variable
	try	{
		if(!ls.load(params,lk))	{ // load spectra into the ls object
			cout << "Error (idx:0003): failed to load spectrum file \"" << params["spectrum file"] << "\"" << '\n';
			return 1;
		}
		if(maximum_spectra != -1)	{ // use the maximum_spectra value, if specified
			ls.set_max(maximum_spectra);
		} 
	}
	catch (...)	{
		cout << "Error (idx:0024): failed to load spectrum file \"" << params["spectrum file"] << "\"" << '\n';
		return 1;
	}
	// record times and report process information
	high_resolution_clock::time_point t2 = high_resolution_clock::now(); //end timing spectrum loading and report
	cout << '\n' << "  spectra = " << ls.spectra.size() << '\n';
	cout << "  skipped = " << ls.skipped << '\n';
	cout << "  spectra &Delta;T = " 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << '\n';
	params["spectrum file validation"] = ls.validation;
	strStream.str("");
	strStream.clear();
	strStream << (long)ls.spectra.size();
	params["spectra"] = strStream.str();
	strStream.str("");
	strStream.clear();
	strStream << (long)ls.spectra.size()+(long)ls.skipped;
	params["spectra, total"] = strStream.str();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, spectra load (s)"] = strStream.str();
	t1 = high_resolution_clock::now(); //begin timing kernel loading
	cout << '\n' << "load & index kernel"  << '\n';
	cout.flush();
	try	{
		if (isB) {
			lk.load_binary(); // load the binary JSON kernel
		}
		else {
			lk.load(); // load the JSON kerenl
		}
	}
	catch (...)	{
		cout << "Error (idx:0028): failed to load kernels" << '\n';
		return 1;
	}
	params["input kernel validation"] = lk.validation;
	t2 = high_resolution_clock::now(); //end timing kernel loading and report
	cout << '\n' << "  kernel pairs = " << lk.kerns.size() << '\n';
	cout << "  kernels &Delta;T = "  << fixed << setprecision(3) 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << '\n';
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, kernel load (s)"] = strStream.str();
	t1 = high_resolution_clock::now(); //begin timing peptide-to-spectrum matching
	cout << '\n' << "perform ids"  << '\n';
	cout.flush();
	create_results cr; //object that will contain match information
	try	{
		if(!cr.create(params,ls,lk))	{ //create peptide-to-spectrum matches
			cout << "Error (idx:0006): failed to create results " << '\n';
			return 1;
		}
	}
	catch (...)	{
		cout << "Error (idx:0026): failed to load create results" << '\n';
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing peptide-to-spectrum matching and report
	cout << '\n' << "  results = " << cr.size() << '\n';
	cout << "  results &Delta;T = "  << fixed << setprecision(3)
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 
				<< " s" << '\n';
	cout << '\n' << "create models & report"  << '\n';
	cout.flush();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, id process (s)"] = strStream.str();
	lk.clean_up();
	ls.clean_up();
	t1 = high_resolution_clock::now(); //begin timing output file creation
	create_output co;
	try	{
		if(isB)	{
			if(!co.create_binary(params,cr,lk.hu_set))	{ //create output file, based on the matches in the cr object
				cout << "Error (idx:0007): failed to create output " << '\n';
				return 1;
			}
		}
		else	{
			if(!co.create(params,cr,lk.hu_set))	{ //create output file, based on the matches in the cr object
				cout << "Error (idx:0007): failed to create output " << '\n';
				return 1;
			}
		}
	}
	catch (...)	{
		cout << "Error (idx:0027): failed to load create output" << '\n';
		return 1;
	}
	t2 = high_resolution_clock::now(); //end timing output file creation
	// record additional parameters for the .meta file
	params["output kernel validation"] = co.validation;
	cout << "  reporting &Delta;T = " << fixed << setprecision(3) 
				<< duration_cast<milliseconds>(t2 - t1).count()/1000.0 << " s" << '\n';
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t1).count()/1000.0;
	params["time, output file creation (s)"] = strStream.str();
	strStream.str("");
	strStream.clear();
	strStream << duration_cast<milliseconds>(t2 - t_origin).count()/1000.0;
	params["time, total (s)"] = strStream.str();
	strStream.str("");
	strStream.clear();
	now = std::chrono::system_clock::now();
    	end_time = std::chrono::system_clock::to_time_t(now);
	strStream << std::ctime(&end_time);
	strTemp = strStream.str();
	// deal with annoying trailing whitespace
	strTemp.erase(std::remove(strTemp.begin(), strTemp.end(), '\n'), strTemp.end());
	strTemp.erase(std::remove(strTemp.begin(), strTemp.end(), '\r'), strTemp.end());
	params["time, completed"] = strTemp;
	// record the epoch time too
	strStream.str("");
	strStream.clear();
	strStream << fixed << setprecision(3);
	strStream << std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count()/1000.0;
	params["time, completed timestamp"] = strStream.str();
	// create the .meta file
	co.dump_meta(params);
	cout << "\ncompleted: " << strTemp << '\n';
	cout << "... done" << '\n';
	return 0;
}


