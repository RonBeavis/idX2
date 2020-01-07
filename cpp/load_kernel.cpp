/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
#include "pch.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unordered_map>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp"
#include "load_spectra.hpp"

load_kernel::load_kernel()	{
	const int32_t c13 = 1003; //difference between the A0 and A1 peaks
	kfile = "";
	fragment_tolerance = 0;
	channel ch;
	ch.delta = 0;
	ch.lower = 0;
	channels.push_back(ch);
	ch.delta = c13;
	ch.lower = 1000 * 1000;
	channels.push_back(ch);
	ch.delta = 2 * 	c13;
	ch.lower = 1500 * 1000;
	channels.push_back(ch);
	clength = channels.size();
	validation = "";
}

load_kernel::~load_kernel(void)	{
}
//
// Loads information from a kernel file specified in _params based on the spectra in the _load_spectra object.
// The information is returned in the kerns and pmindex member objects
//
bool load_kernel::load(void)	{
	ifstream ifs;
	ifs.open(kfile,std::ifstream::in);
	if(!ifs.good())	{
		cout << "Kernel file \"" << kfile << "\" would not open" << endl;
		return false;
	}
	string line;
	const double ft = 1.0/fragment_tolerance; //fragment tolerance
	const double pt = 1.0/70.0; //maximum parent tolerance in millidaltons
	const double ppm = 2.0E-5; //parent tolerance in ppm
	auto itsp = sp_set.end();
	auto itppm = sp_set.end();
	int32_t skipped = 0;
	int32_t hmatched = 0;
	int32_t pm = 0;
	int32_t mv = 0;
	int32_t u = 0;
	int32_t lines = 0;
	bool skip = true;
	int32_t lower = 0;
	int32_t ppm_delta = 0;
	kPair pr;
	const int max_buffer = 1024*16-1;
	char *buffer = new char[max_buffer+1];
	char *pPm = NULL;
	SizeType a = 0;
	size_t b = 0;
	//loop through kernel lines
	const size_t clength = channels.size();
	while(ifs.getline(buffer,max_buffer))	{
		//print keep-alive text for logging
		if(lines != 0 and lines % 5000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(lines != 0 and lines % 100000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		lines++;
		// quickly get the pm value without loading the rapidjson::Document
		pPm = std::strstr(buffer,"\"pm\":");
		if(pPm != NULL)	{
			pm = atoi(pPm + 5);
		}
		else	{
			Document js; //rapidjson main object
			js.ParseInsitu(buffer);
			if(js.HasMember("validation"))	{
				validation = js["value"].GetString();
			}
			continue;
		}
		ppm_delta = (int32_t)(0.5+(double)pm*ppm); //parent mass tolerance based on ppm
		//check parent mass for ppm tolerance
		skip = true;
		for(b = 0; b < clength; b++)	{
			channels[b].mv = (int32_t)(0.5+(double)(pm + channels[b].delta)*pt); //reduced parent mass
			lower = pm - ppm_delta + channels[b].delta;
			itppm = sp_set.lower_bound(lower);
			if(pm > channels[b].lower)	{
				if(itppm != itsp and (*itppm-lower) <= 2*ppm_delta)	{
					skip = false;
					channels[b].dead = false;
				}
				else	{
					channels[b].dead = true;
				}
			}
		}
		//check A1 mass for ppm tolerance
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
		u = (int32_t)js["u"].GetInt();  //record the unique kernel id
		const Value& jbs = js["bs"]; //retrieve reference to the b-type fragments                                                           
		for(a = 0; a < jbs.Size();a++)	{
			pr.second = (int32_t)(0.5+jbs[a].GetDouble()*ft); //reduced fragment mass
			pr.first = mv;
			check_and_update(pr,u);
		}
		const Value& jys = js["ys"]; //retrieve reference to the y-type fragments
		for(a = 0; a < jys.Size();a++)	{
			pr.second = (int32_t)(0.5+jys[a].GetDouble()*ft); //reducted fragment mass
			pr.first = mv;
			check_and_update(pr,u);
		}
	}
	cout.flush();
	ifs.close();
	delete buffer;
	sp_set.clear();
	return true;
}

bool load_kernel::get_next(ifstream &_ifs,jsObject& _js)
{
	_js.reset();
	int32_t jsl = 0;
	if(_ifs.fail())	{
		return false;
	}
	_ifs.read((char *)(&jsl),4);
	if(_ifs.bad())	{
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
	while(count < jsl && !(_ifs.bad() or _ifs.eof()))	{
		_ifs.read((char *)(&klen),4);
		_ifs.read((char *)(_js.pKey),klen);
		_js.pKey[klen] = '\0';
		_js.key = _js.pKey;
		_ifs.read((char *)(&element),1);
		switch(element)	{
			case 'm':
				_ifs.read((char *)(&tlen),4);
				for(i = 0; i < (size_t)tlen;i++)	{
					_ifs.read((char *)(_js.pBuffer),9);
				}
				break;
			case 'l':
				_ifs.read((char *)(&tlen),4);
				_ifs.read((char *)(_js.pBuffer),4*tlen);
				pI = (int *)_js.pBuffer;
				if(_js.key == "bs")	{
					_js.bs.insert(_js.bs.end(),pI,pI+tlen);
				}
				else if(_js.key == "ys")	{
					_js.ys.insert(_js.ys.end(),pI,pI+tlen);
				}
				break;
			case 's':
				_ifs.read((char *)(&tlen),4);
				_ifs.read((char *)(_js.pBuffer),tlen);
				_js.pBuffer[tlen] = '\0';
				break;
			case 'i':
				_ifs.read((char *)(&itemp),4);
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
		if(_js.key == "value")	{
			validation = (char *)_js.pBuffer;
			return false;
		}
	}
	return true;
}

bool load_kernel::load_binary(void)	{
	ifstream ifs(kfile, ios::in | ios::binary);
	if(ifs.fail())	{
		return false;
	}
	string line;
	const double ft = 1.0/fragment_tolerance; //fragment tolerance
	const double pt = 1.0/70.0; //maximum parent tolerance in millidaltons
	const double ppm = 2.0E-5; //parent tolerance in ppm
	auto itsp = sp_set.end();
	auto itppm = sp_set.end();
	int32_t skipped = 0;
	int32_t hmatched = 0;
	int32_t pm = 0;
	int32_t mv = 0;
	int32_t u = 0;
	int32_t lines = 0;
	bool skip = true;
	int32_t lower = 0;
	int32_t ppm_delta = 0;
	kPair pr;
	jsObject js;
	//loop through kernel lines
	size_t a = 0;
	size_t b = 0;
	//loop through kernel lines
	const size_t clength = channels.size();
	while(get_next(ifs,js))	{
		//print keep-alive text for logging
		if(lines != 0 and lines % 5000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(lines != 0 and lines % 100000 == 0)	{
			cout << " " << lines << endl;
			cout.flush();
		}
		lines++;
		// quickly get the pm value without loading the rapidjson::Document
		pm = js.pm;
		if(pm == 0)	{
			continue;
		}
		ppm_delta = (int32_t)(0.5+(double)pm*ppm); //parent mass tolerance based on ppm
		//check parent mass for ppm tolerance
		skip = true;
		for(b = 0; b < clength; b++)	{
			channels[b].mv = (int32_t)(0.5+(double)(pm + channels[b].delta)*pt); //reduced parent mass
			lower = pm - ppm_delta + channels[b].delta;
			itppm = sp_set.lower_bound(lower);
			if(pm > channels[b].lower)	{
				if(itppm != itsp and (*itppm-lower) <= 2*ppm_delta)	{
					skip = false;
					channels[b].dead = false;
				}
				else	{
					channels[b].dead = true;
				}
			}
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
		for(a = 0; a < js.bs.size();a++)	{
			pr.second = (int32_t)(0.5+js.bs[a]*ft); //reduced fragment mass
			pr.first = mv;
			check_and_update(pr,u);
		}
		for(a = 0; a < js.ys.size();a++)	{
			pr.second  = (int32_t)(0.5+js.ys[a]*ft); //reduced fragment mass
			pr.first = mv;
			check_and_update(pr,u);
		}
	}
	cout << "\n";
	cout.flush();
	ifs.close();
	sp_set.clear();
	return true;
}

