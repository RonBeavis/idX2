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
#include <map>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp"
#include "load_spectra.hpp"
#include "create_results.hpp"

create_results::create_results(void) {
}

create_results::~create_results(void) {
}
//
//create_results takes the information recorded from the spectrum file and kernel file
//and generates spectrum-to-peptide matches
//
bool create_results::create(map<string, string>& _p, //parameters
		load_spectra& _l, //spectrum informaiton
		load_kernel& _k) { //kernel information
	int32_t z = 1; //serial number for PSM: used if spectrum::sc is not available
	const double pt = 1.0/70.0;
	vector<int32_t> dvals{ -1,0,1 }; //values to compensate for integer rounding effects
//	vector<int32_t> dvals{ 0 }; //values to compensate for integer rounding effects
	//initialize variables for the identification process
	int32_t d = 0;
	int32_t m = 0;
	int32_t idi = 0;
	int32_t cpm = 0;
	vector<int32_t> idx;
	vector<vector<int32_t> > ident;
	vector<int32_t> ims;
	auto itk = _k.kerns.kindex_a[0].end(); //store value for subsequent map finds
	kPair pv; //a kernel (parent,fragment) pair
	const size_t clength = _k.channels.size();
	size_t a = 0;
	size_t n = 0;
	size_t b = 0;
	int32_t mi = 0;
	//loop through spectra
	for (size_t s = 0; s < _l.spectra.size(); s++) {
		//output keep-alive text for logging
		if (s != 0 and s % 500 == 0) {
			cout << '.';
			cout.flush();
		}
		if (s != 0 and s % 5000 == 0) {
			cout << " " << s << endl;
			cout.flush();
		}
		ident.clear(); //initialize the temporary identification vector
		ims.clear(); //initialize the temporary mass vector
		mi = _l.spectra[s].pm;
		idx.clear(); //initialize temporary list of kernel identifiers
		idi = 0;
		for(a = 0; a < clength; a++)	{
			if(mi > _k.channels[a].lower)	{
				_k.channels[a].mv = mi - _k.channels[a].delta; //add the parent mass
				_k.channels[a].dead = false;
			}
			else	{
				_k.channels[a].dead = true;
			}
		}
		//loop though masses
		for (a = 0; a < clength; a++) {
			if(_k.channels[a].dead)	{
				continue;
			}
			cpm = (int32_t)(0.5 + (double)_k.channels[a].mv * pt); //current reduced parent mass
			//loop through possible parent mass values to compensate for rounding errors
			for (n = 0; n < dvals.size(); n++) {
				d = dvals[n];
				pv.first = cpm + d; //initialize (parent,fragment) pair
				if(_k.kerns.mvindex_a[a].find(pv.first) == _k.kerns.mvindex_a[a].end())	{ //bail if parent had no fragments in any kernel
					continue;
				}
				//loop through spectrum fragment (mass,intensity) pairs
				for(b = 0; b < _l.spectra[s].mis.size(); b++)	{
					m = _l.spectra[s].mis[b].first; //fragment mass
					pv.second = m; //set fragment mass in (mass,intensity) pair
					itk = _k.kerns.kindex_a[a].find(pv); //
					idx.clear();
					//check if pair was in kernels and record result
					if(itk !=  _k.kerns.kindex_a[a].end())	{ 
						idx.insert(idx.end(),_k.kerns.kindex_a[a][pv].begin(),_k.kerns.kindex_a[a][pv].end());
						idi = _l.spectra[s].mis[b].second;
					}
					else	{
						continue;
					}
					//record kernel matching information
					ims.push_back(idi);
					ident.push_back(idx);
				}
			}
		}
		if(ident.empty())	{ //bail out if no identifications were found
			continue;
		}
		//check all fragment identifications and find best matches
		map<int32_t,int32_t> ans;
		map<int32_t,int32_t> aint;
		//iterate through all information recorded in the ident vector
		//and find the kernel identifiers corresponding to the best identification
		for(b = 0; b < ident.size(); b++)	{
			set<int32_t> sv(ident[b].begin(),ident[b].end());
			set<int32_t>::iterator it = sv.begin();
			while(it != sv.end())	{
				if(ans.find(*it) != ans.end()) {
					ans[*it] += 1;
					aint[*it] += ims[b];
				}
				else	{
					ans[*it] = 1;
					aint[*it] = ims[b];
				}
				it++;
			}
		}
		int32_t mn = 0;
		vector<int32_t> mv;
		vector<int32_t> iv;
		map<int32_t,int32_t>::iterator itans = ans.begin();
		//arrange identification information
		while(itans != ans.end())	{
			if(itans->second > mn)	{
				mn = itans->second;
				mv.clear();
				mv.push_back(itans->first);
				iv.clear();
				iv.push_back(aint[itans->first]);
			}
			else if(itans->second == mn)	{
				mv.push_back(itans->first);
				iv.push_back(aint[itans->first]);
			}
			itans++;
		}
		//record identification information if at least 5 fragments were matched
		if(mn > 4)	{
			id r;
			int32_t max_i = *max_element(iv.begin(),iv.end());
			int32_t tot_i = 0;
			for(b = 0; b < mv.size(); b++)	{
				tot_i += iv[b];
				if(iv[b] < max_i)	{
					continue;
				}
				r.ks.push_back(mv[b]);
			}
			r.sn = z;
			r.peaks = mn;
			r.ri = (double)tot_i/(double)_l.spectra[s].isum;
			r.pm = _l.spectra[s].pm;
			r.pz = _l.spectra[s].pz;
			r.sc = _l.spectra[s].sc;
			r.rt = _l.spectra[s].rt;
			r.ions = _l.spectra[s].pks;
			ids.push_back(r);
		}
		z += 1;
	}
	cout << endl;
	cout.flush();
	return true;
}
