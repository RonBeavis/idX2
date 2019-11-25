/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include "rapidjson/document.h"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <unordered_map>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
#include "load_spectra.hpp"
#include "load_kernel.hpp"

load_kernel::load_kernel(void)	{
}

load_kernel::~load_kernel(void)	{
}
//
// Loads information from a kernel file specified in _params based on the spectra in the _load_spectra object.
// The information is returned in the _kernesl object and _mindex map
//
bool load_kernel::load(map<string,string>& _params,load_spectra& _l,kernels& _kernels,map<int64_t,int64_t>& _mindex)	{
	ifstream istr;
	istr.open(_params["kernel file"]); //open kernel file stream
	if(istr.fail())	{
		return false;
	}
	string line;
	using namespace rapidjson; //namespace for the rapidjson methods
	const double ft = 1.0/atof(_params["fragment tolerance"].c_str()); //fragment tolerance
	const double pt = 1.0/70.0; //maximum parent tolerance in millidaltons
	const double ppm = 2.0E-5; //parent tolerance in ppm
	set<int64_t> sp_set; //set of spectrum parent masses
	phmap::parallel_flat_hash_set<sPair> spairs; //map of spectrum (parent:fragment) mass pairs
	for(size_t a = 0; a < _l.spectra.size();a++)	{
		sp_set.insert(_l.spectra[a].pm);
		spairs.insert(_l.spectra[a].spairs.begin(),_l.spectra[a].spairs.end());
	}
	auto itsp = sp_set.end();
	auto itppm = sp_set.end();
	int64_t skipped = 0;
	int64_t hmatched = 0;
	int64_t pm = 0;
	int64_t mv = 0;
	int64_t u = 0;
	const int64_t c13 = 1003; //difference between the A0 and A1 peaks
	int64_t val = 0;
	int64_t lines = 0;
	bool skip = true;
	int64_t lower = 0;
	int64_t delta = 0;
	kPair pr;
	//loop through kernel lines
	while(getline(istr,line))	{
		//print keep-alive text for logging
		if(lines != 0 and lines % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(lines != 0 and lines % 200000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		lines++;
		Document js; //rapidjson main object
   		js.Parse(line.c_str(),line.length()); //add information for a single JSON Lines object
		if(!js.HasMember("pm"))	{ //bail out if the JSON object does not have a parent mass
			continue;
		}
		if(js["u"].GetInt() != js["h"].GetInt())	{ //bail out if the JSON object is not the first instance of the peptide & modifications
			hmatched++;
			continue;
		}
		pm = (int64_t)js["pm"].GetInt(); //parent mass
		mv = (int64_t)(0.5+(double)pm*pt); //reduced parent mass
		delta = (int64_t)(0.5+(double)pm*ppm); //parent mass tolerance based on ppm
		//check parent mass for ppm tolerance
		lower = pm-delta;
		itppm = sp_set.lower_bound(lower);
		skip = true;
		if(itppm != itsp and (*itppm-lower) <= 2*delta)	{
			skip = false;
		}
		//check A1 mass for ppm tolerance
		lower = pm+c13-delta;
		itppm = sp_set.lower_bound(lower);
		if(itppm != itsp and (*itppm-lower) <= 2*delta)	{
			skip = false;
		}
		if(skip)	{ //bail out if the parent mass doesn't match
			skipped++;
			continue;
		}
		u = (int64_t)js["u"].GetInt();  //record the unique kernel id
		_mindex[u] = pm;
		const Value& jbs = js["bs"]; //retrieve reference to the b-type fragment fragments                                                           
//		size_t vpos = 0;
		pr.first = (int64_t)mv; //initialize the parent mass element of the (parent:fragment) pair
		for(SizeType a = 0; a < jbs.Size();a++)	{
			val = (int64_t)(0.5+jbs[a].GetDouble()*ft); //reduced fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(_kernels.kindex.find(pr) == _kernels.kindex.end())	{ 
				_kernels.add_pair(pr); //create a new vector for the object
			}
			_kernels.mvindex.insert((int64_t)mv); //add parent mass to set
			_kernels.kindex[pr].push_back(u); //add kernel id to vector
		}
		const Value& jys = js["ys"]; //retrieve reference to the y-type fragments
		for(SizeType a = 0; a < jys.Size();a++)	{
			val = (int64_t)(0.5+jys[a].GetDouble()*ft); //reducted fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(_kernels.kindex.find(pr) == _kernels.kindex.end())	{
				_kernels.add_pair(pr); //create a new vector for the object
			}
			_kernels.mvindex.insert((int64_t)mv); //add parent mass to set
			_kernels.kindex[pr].push_back(u);//add kernel id to vector
		}
	}
	cout << "\n";
	cout.flush();
	istr.close();
	return true;
}

