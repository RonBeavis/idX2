/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include <fstream>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h" //fast maps and sets used for indexing
#include "picosha2/picosha2.h"
using namespace std; //namespace used throughout
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp"
#include "load_spectra.hpp"

// loads spectra using the MGF file specified in _params

bool load_spectra::load(map<string,string>& _params,load_kernel& _lk)	{
	// check for MGF format file via file extension
	size_t f = _params["spectrum file"].find(".mgf");
	if(f != _params["spectrum file"].npos && f == _params["spectrum file"].length() - 4)	{
		cout << "mgf format file detected\n";
		cout.flush();
		return load_mgf(_params,_lk);
	}
	// check for CMN format file via file extension
	f = _params["spectrum file"].find(".cmn");
	if(f != _params["spectrum file"].npos && f == _params["spectrum file"].length() - 4)	{
		return load_cmn(_params,_lk);
	}
	return false;
}

bool load_spectra::load_mgf(map<string,string>& _params,load_kernel& _lk)	{
	// open the spectrum file
	ifstream istr;
	istr.open(_params["spectrum file"]); // opens input file stream
	if(istr.fail())	{ // bail on fail
		return false;
	}
	// set up parameters and temporary values
	int32_t res = (int32_t)atoi(_params["fragment tolerance"].c_str());
	double pt = atoi(_params["parent tolerance"].c_str());
	size_t len = 1024*4; // maximum line length
	size_t size = 0;  // size of current line
	char *line = new char[len]; //char buffer used to read file
	string temp = ""; // temporary string object
	string desc = ""; //description information (if available)
	char *pos = NULL; // char pointer for parsing
	size_t equals = 0;
	double parent = 0.0; //parent mass
	double charge = 1.0; //parent charge
	string run_time = ""; //chromatographic retention time information (if available)
	vector<double> masses; // vector of fragment masses
	vector<double> intensities; // vector of fragment intensities
	spectrum sp; // spectrum object
	const double proton = 1.007276; //constant used to recalculate neutral masses
	int32_t scan = 0; //scan number of the spectrum
	int32_t s = 1; // spectrum cound
	double rt = 0.0; // run time (seconds)
	// process the file, one line at a time
	while(istr.good() && !istr.eof())	{
		istr.getline(line,len-1,'\n');
		temp = line;
		size = temp.size();
		equals = temp.find("=");
		if(atof(line) > 0.0)	{ //if no tag present, interpret as mass intensity pair
			masses.push_back(atof(line));
			pos = line;
			while(*pos != '\0' && isspace(*pos) != 0)	{
				pos++;
			}
			while(*pos != '\0' && isspace(*pos) == 0)	{
				pos++;
			}
			intensities.push_back(atof(pos));
		}
		else if(temp.find("PEPMASS=") != temp.npos)	{ //tag corresponding to the parent ion mass
			temp = temp.substr(equals+1,size-equals+1);
			parent = atof(temp.c_str());
		}
		else if(temp.find("BEGIN IONS") != temp.npos)	{ //tag starting a single spectrum
		}
		else if(temp.find("TITLE=") != temp.npos)	{ //tag for the description information
			desc = temp.substr(equals+1,size-equals+1);
		}
		else if(temp.find("RTINSECONDS=") != temp.npos)	{ //tag for retention time information
			run_time = temp.substr(equals+1,size-equals+1);
			rt = (double)atof(run_time.c_str());
		}
		else if(temp.find("CHARGE=") != temp.npos)	{ //tag for parent ion charge
			temp = temp.substr(equals+1,size-equals+1);
			charge = atof(temp.c_str());
		}
		else if(temp.find("SCANS=") != temp.npos)	{ //tag for spectrum scan number (if available)
			temp = temp.substr(equals+1,size-equals+1);
			scan = atoi(temp.c_str());
		}
		else if(temp.find("END IONS") != temp.npos)	{ //tag for the end of a spectrum
			//carry out calculations and type conversions
			sp.clear();
			// skip spectrum if the parent ion mass <= 600 Da
			if(parent*charge > 600.0)	{
				sp.pm = (int32_t)(0.5 + 1000*(parent*charge - proton*charge)); // convert to mDa
				sp.pz = (int32_t)charge; // record charge
				sp.pi = 100; // parent intensity
				sp.pt = pt; // record parent tolerance
				sp.desc = desc; // record spectrum description
				sp.rt = rt; // record retention time
				//substitute ordinal value for scan number, if no scan available
				if(scan == 0)	{
					sp.sc = s;
					equals = desc.find("scan=");
					if(equals != desc.npos)	{
						equals += 5;
						sp.sc = atol(desc.substr(equals,size-equals).c_str());
					}
				}
				else	{
					sp.sc = scan;
				}
				size_t i = 0;
				pair<int32_t,int32_t> p;
				while(i < masses.size())	{
					p.first = (int32_t)(0.5+1000.0*(masses[i]-proton));
					p.second = (int32_t)(0.5+intensities[i]);
					sp.mis.push_back(p);
					i++;
				}
				// clean up spectrum pairs for later use
				sp.condition(res,50);
				// add information to _lk
				_lk.sp_set.insert(sp.pm);
				_lk.spairs.insert(sp.spairs.begin(),sp.spairs.end());
				// remove pair information
				sp.spairs.clear();
				// record spectrum info spectra vector
				spectra.push_back(sp);
			}
			else	{
				// track number of spectra that didn't pass mass test
				skipped++;
			}
			//clean up to be ready for the next spectrum
			desc = "";
			masses.clear();
			intensities.clear();
			parent = 0.0;
			charge = 0.0;
			s++;
			//output progress text for logging
			if(s % 2500 == 0)	{
				cout << '.';
				cout.flush();
			}
			if(s != 0 and s % 50000 == 0)	{
				cout << ' ' << s << '\n';
				cout.flush();
			}
			scan = 0;
		}
	}
	cout << "\n";
	cout.flush();
	delete line;
	istr.close();
	//calculate a SHA256 hash value for the spectrum file
	//for validation use
	std::ifstream ifile(_params["spectrum file"], std::ios::binary);
	std::vector<unsigned char> sfile(picosha2::k_digest_size);
	picosha2::hash256(ifile, sfile.begin(), sfile.end());
	ifile.close();
	std::string hash_hex_str;
	picosha2::hash256_hex_string(sfile, hash_hex_str);
	validation = hash_hex_str;
	return true;
}

bool load_spectra::load_cmn(map<string,string>& _params,load_kernel& _lk)	{
	// open the spectrum file
	std::fstream inFile(_params["spectrum file"].c_str(),std::ios_base::binary|std::ios_base::in);
	if(inFile.fail())	{
		return false;
	}
	char *pLine = new char[1024];
	inFile.read((char*)pLine,256);
	if(inFile.fail())	{
		cout << "File not in CMN format (1)\n";
		cout.flush();
		return false;
	}
	unsigned int tLength = 255;
	pLine[tLength] = '\0';
	size_t version = 1;
	if(strstr(pLine,"CMN ") != pLine)	{
		inFile.close();
		cout << "File not in CMN format (2)\n";
		cout.flush();
		delete pLine;
		return false;
	}
	if(pLine[64] != 0)	{
		version = 2;
	}
	else	{
		version = 1;
	}
	// set up parameters and temporary values
	int32_t res = (int32_t)atoi(_params["fragment tolerance"].c_str());
	double pt = atoi(_params["parent tolerance"].c_str());
	size_t len = 1024*4; // maximum line length
	size_t size = 0;  // size of current line
	char *line = new char[len]; //char buffer used to read file
	string temp = ""; // temporary string object
	string desc = ""; //description information (if available)
	size_t equals = 0;
	double parent = 0.0; //parent mass
	double charge = 1.0; //parent charge
	string run_time = ""; //chromatographic retention time information (if available)
	vector<double> masses; // vector of fragment masses
	vector<double> intensities; // vector of fragment intensities
	spectrum sp; // spectrum object
	const double proton = 1.007276; //constant used to recalculate neutral masses
	uint32_t scan = 0; //scan number of the spectrum
	uint32_t s = 1; // spectrum count
	double rt = 0.0; // run time (seconds)
	unsigned short sValue = 0;
	unsigned char cValue = 0;
	uint32_t iValue = 0;
	float fValue = 0.0;
	double dValue = 0.0;
	float pi = 100.0;
	// process the file, one line at a time
	while(inFile.good())	{
		iValue = 0;
		inFile.read((char *)&iValue,4);
		scan = iValue;
		inFile.read((char *)&dValue,sizeof(double));
		parent = dValue - proton;
		inFile.read((char *)&cValue,1);
		charge = (double)cValue;
		if(version == 2)	{
			uint32_t tValue = 0;
			inFile.read((char *)&tValue,4);
			if(iValue > tLength)	{
				tLength = tValue + 255;
				delete pLine;
				pLine = new char[tLength];
			}
			inFile.read((char *)pLine,tValue);
			pLine[tValue] = '\0';
		}
		else	{	
			inFile.read((char *)&cValue,1);
			inFile.read((char *)pLine,(int)cValue);

			pLine[cValue] = '\0';
		}
		desc = pLine;
		sp.clear();
		sp.pm = (int32_t)(0.5 + 1000*parent); // convert to mDa
		sp.pz = (int32_t)charge; // record charge
		sp.pi = (int32_t)pi; // parent intensity
		sp.pt = pt; // record parent tolerance
		sp.desc = desc; // record spectrum description
		sp.rt = rt; // record retention time
		//substitute ordinal value for scan number, if no scan available
		sp.sc = scan;
		equals = desc.find("scan=");
		if(equals != desc.npos)	{
			equals += 5;
			sp.sc = atol(desc.substr(equals,size-equals).c_str());
		}
		cout << sp.sc << "\n";
		fValue = 0.0;
		inFile.read((char *)&fValue,sizeof(float));
		cValue = 0;
		inFile.read((char *)&cValue,1);
		size_t tSize = (size_t)cValue;
		fValue = 0.0;
		inFile.read((char *)&fValue,sizeof(float));
		inFile.read((char *)&cValue,1);
		float fScale = fValue;
		size_t a = 0;
		inFile.read((char *)&sValue,2);
		iValue = (unsigned int)sValue;
		masses.push_back((float)iValue/fScale);
		a++;
		while(a < tSize)	{
			inFile.read((char *)&sValue,2);
			iValue += (unsigned int)sValue;
			masses.push_back((float)iValue/fScale);
			a++;
		}
		a = 0;
		double dSum = 0.0;
		char cMax = 0;
		while(a < tSize)	{
			inFile.read((char *)&cValue,1);
			intensities.push_back((float)cValue);
			if(cMax < cValue)	{
				cMax = cValue;
			}
			dSum += (double)cValue;
			a++;
		}
		if(parent > 600.0)	{
			size_t i = 0;
			pair<int32_t,int32_t> p;
			while(i < masses.size())	{
				p.first = (int32_t)(0.5+1000.0*(masses[i]-proton));
				p.second = (int32_t)(0.5+intensities[i]);
				sp.mis.push_back(p);
				i++;
			}
			// clean up spectrum pairs for later use
			sp.condition(res,50);
			// add information to _lk
			_lk.sp_set.insert(sp.pm);
			_lk.spairs.insert(sp.spairs.begin(),sp.spairs.end());
			// remove pair information
			sp.spairs.clear();
			// record spectrum info spectra vector
			spectra.push_back(sp);
		}
		else	{
			// track number of spectra that didn't pass mass test
			skipped++;
		}
		//clean up to be ready for the next spectrum
		desc = "";
		masses.clear();
		intensities.clear();
		parent = 0.0;
		charge = 0.0;
		s++;
		//output progress text for logging
		if(s % 2500 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(s != 0 and s % 50000 == 0)	{
			cout << ' ' << s << '\n';
			cout.flush();
		}
		scan = 0;
	}
	cout << "\n";
	cout.flush();
	delete line;
	inFile.close();
	//calculate a SHA256 hash value for the spectrum file
	//for validation use
	std::ifstream ifile(_params["spectrum file"], std::ios::binary);
	std::vector<unsigned char> sfile(picosha2::k_digest_size);
	picosha2::hash256(ifile, sfile.begin(), sfile.end());
	ifile.close();
	std::string hash_hex_str;
	picosha2::hash256_hex_string(sfile, hash_hex_str);
	validation = hash_hex_str;
	return true;
}

