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

load_kernel::load_kernel()	{
	kfile = "";
	fragment_tolerance = 0;
	thread = 0;
	threads = 1;
}

load_kernel::~load_kernel(void)	{
}
//
// Loads information from a kernel file specified in _params based on the spectra in the _load_spectra object.
// The information is returned in the kerns and pmindex member objects
//
bool load_kernel::load(void)	{
	FILE *pFile = ::fopen(kfile.c_str(),"r");
	if(pFile == NULL)	{
		return false;
	}
	string line;
	const double ft = 1.0/fragment_tolerance; //fragment tolerance
	const double pt = 1.0/70.0; //maximum parent tolerance in millidaltons
	const double ppm = 2.0E-5; //parent tolerance in ppm
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
	const int max_buffer = 1024*16-1;
	char *buffer = new char[max_buffer+1];
	char *pPm = NULL;
	//loop through kernel lines
	while(::fgets(buffer,max_buffer,pFile) != NULL)	{
		//print keep-alive text for logging
		if(thread == 0 and lines != 0 and lines % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(thread == 0 and lines != 0 and lines % 200000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		if(lines % threads != thread)	{
			lines++;
			continue;
		}
		lines++;
		// quickly get the pm value without loading the rapidjson::Document
		pPm = std::strstr(buffer,"\"pm\":");
		if(pPm != NULL)	{
			pm = atoi(pPm + 5);
		}
		else	{
			continue;
		}

//		pm = (int64_t)js["pm"].GetInt(); //parent mass
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
		Document js; //rapidjson main object
		js.ParseInsitu(buffer);
		if(!js.HasMember("pm"))	{ //bail out if the JSON object does not have a parent mass
			continue;
		}
		add_hu(js["h"].GetInt(),js["u"].GetInt());
		if(js["u"].GetInt() != js["h"].GetInt())	{ //bail out if JSON is not 1st instance of the peptide + mods
			hmatched++;
			continue;
		}
		u = (int64_t)js["u"].GetInt();  //record the unique kernel id
		pmindex[u] = pm;
		const Value& jbs = js["bs"]; //retrieve reference to the b-type fragments                                                           
		pr.first = (int64_t)mv; //initialize the parent mass element of the (parent:fragment) pair
		for(SizeType a = 0; a < jbs.Size();a++)	{
			val = (int64_t)(0.5+jbs[a].GetDouble()*ft); //reduced fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(kerns.kindex.find(pr) == kerns.kindex.end())	{ 
				kerns.add_pair(pr); //create a new vector for the object
			}
			kerns.mvindex.insert((int64_t)mv); //add parent mass to set
			kerns.kindex[pr].push_back(u); //add kernel id to vector
		}
		const Value& jys = js["ys"]; //retrieve reference to the y-type fragments
		for(SizeType a = 0; a < jys.Size();a++)	{
			val = (int64_t)(0.5+jys[a].GetDouble()*ft); //reducted fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(kerns.kindex.find(pr) == kerns.kindex.end())	{
				kerns.add_pair(pr); //create a new vector for the object
			}
			kerns.mvindex.insert((int64_t)mv); //add parent mass to set
			kerns.kindex[pr].push_back(u);//add kernel id to vector
		}
	}
	if(thread == 0)	{
		cout << "\n";
		cout.flush();
	}
	::fclose(pFile);
	delete buffer;
	return true;
}

bool load_kernel::get_next(FILE *_pFile,jsObject& _js)
{
	_js.reset();
	int jsl = 0;
	size_t ret = fread(&jsl,4,1,_pFile);
	if(ret != 1)	{
		cout << "failed to get json size" << endl;
		return false;
	}
	int count = 0;
	int klen = 0;
	char element = '\0';
	int tlen = 0;
	int itemp = 0;
	size_t i = 0;
	int *pI = 0;
	while(count < jsl)	{
		ret = fread(&klen,4,1,_pFile);
		ret = fread(_js.pKey,klen,1,_pFile);
		_js.pKey[klen] = '\0';
		_js.key = _js.pKey;
		if(_js.key == "validation")	{
			return false;
		}
		ret = fread(&element,1,1,_pFile);
		switch(element)	{
			case 'm':
				ret = fread(&tlen,4,1,_pFile);
				for(i = 0; i < (size_t)tlen;i++)	{
					ret = fread(_js.pBuffer,9,1,_pFile);
				}
				break;
			case 'l':
				ret = fread(&tlen,4,1,_pFile);
				ret = fread(_js.pBuffer,4*tlen,1,_pFile);
				pI = (int *)_js.pBuffer;
				if(_js.key == "bs")	{
					_js.bs.insert(_js.bs.end(),pI,pI+tlen);
				}
				else if(_js.key == "ys")	{
					_js.ys.insert(_js.ys.end(),pI,pI+tlen);
				}
				break;
			case 's':
				ret = fread(&tlen,4,1,_pFile);
				ret = fread(_js.pBuffer,tlen,1,_pFile);
				_js.pBuffer[tlen] = '\0';
				break;
			case 'i':
				ret = fread(&itemp,4,1,_pFile);
				if(_js.key == "pm")	{
					_js.pm = itemp;
				}
				else if(_js.key == "h")	{
					_js.h = itemp;
				}
				else if(_js.key == "u")	{
					_js.u = itemp;
				}
				break;
			default:
				cout << "bad element value" << endl;
		}
		count++;
	}
	return true;
}

bool load_kernel::load_binary(void)	{
	FILE *pFile = ::fopen(kfile.c_str(),"rb");
	if(pFile == NULL)	{
		return false;
	}
	string line;
	const double ft = 1.0/fragment_tolerance; //fragment tolerance
	const double pt = 1.0/70.0; //maximum parent tolerance in millidaltons
	const double ppm = 2.0E-5; //parent tolerance in ppm
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
	jsObject js;
	//loop through kernel lines
	while(get_next(pFile,js))	{
		//print keep-alive text for logging
		if(thread == 0 and lines != 0 and lines % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(thread == 0 and lines != 0 and lines % 200000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		if(lines % threads != thread)	{
			lines++;
			continue;
		}
		lines++;
		// quickly get the pm value without loading the rapidjson::Document
		pm = js.pm;
		if(pm == 0)	{
			continue;
		}
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
		add_hu(js.h,js.u);
		if(js.u != js.h)	{ //bail out if JSON is not 1st instance of the peptide + mods
			hmatched++;
			continue;
		}
		u = js.u;  //record the unique kernel id
		pmindex[u] = pm;
		pr.first = (int64_t)mv; //initialize the parent mass element of the (parent:fragment) pair
		for(size_t a = 0; a < js.bs.size();a++)	{
			val = (int64_t)(0.5+js.bs[a]*ft); //reduced fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(kerns.kindex.find(pr) == kerns.kindex.end())	{ 
				kerns.add_pair(pr); //create a new vector for the object
			}
			kerns.mvindex.insert((int64_t)mv); //add parent mass to set
			kerns.kindex[pr].push_back(u); //add kernel id to vector
		}
		for(size_t a = 0; a < js.ys.size();a++)	{
			val = (int64_t)(0.5+js.ys[a]*ft); //reducted fragment mass
			pr.second = val;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			if(kerns.kindex.find(pr) == kerns.kindex.end())	{
				kerns.add_pair(pr); //create a new vector for the object
			}
			kerns.mvindex.insert((int64_t)mv); //add parent mass to set
			kerns.kindex[pr].push_back(u);//add kernel id to vector
		}
	}
	if(thread == 0)	{
		cout << "\n";
		cout.flush();
	}
	::fclose(pFile);
	return true;
}

