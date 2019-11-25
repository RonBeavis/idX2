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
#include "load_spectra.hpp"
#include "load_kernel.hpp"
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
		kernels& _k, //kernel information
		map<int64_t, int64_t>& _m) { //allowed (parent,fragment) mass pairs
	int64_t z = 1; //serial number for PSM: used if spectrum::sc is not available
	double pt = 1.0 / 70.0;
	double ppm = atof(_p["parent tolerance"].c_str()) / 1.0e06; //fractional parent mass tolerance
	vector<int64_t> dvals{ -1,0,1 }; //values to compensate for integer rounding effects
	//initialize variables for the identification process
	int64_t d = 0;
	int64_t m = 0;
	int64_t m2 = 0;
	double c13 = 1003.0; //mass difference between A0 and A1 peaks, mass(13C)-mass(12C)
	double rm = 0.0;
	int64_t pm = 0;
	int64_t pz = 0;
	int64_t idi = 0;
	int64_t cpm = 0;
	vector<int64_t> idx;
	vector<vector<int64_t> > ident;
	vector<int64_t> ims;
	bool use2 = false; //use z=2 fragments if fragment mass large enough
	bool use3 = false; //use z=2 fragments if parent z=3
	auto itk = _k.kindex.end(); //store value for subsequent map finds
	kPair pv; //a kernel (parent,fragment) pair
	vector<double> mvs;
	vector<double> ms;
	//loop through spectra
	for (size_t s = 0; s < _l.spectra.size(); s++) {
		//output keep-alive text for logging
		if (s != 0 and s % 500 == 0) {
			cout << '.';
			cout.flush();
		}
		if (s != 0 and s % 10000 == 0) {
			cout << " " << s << endl;
			cout.flush();
		}
		ident.clear(); //initialize the temporary identification vector
		ims.clear(); //initialize the temporary mass vector
		rm = (double)_l.spectra[s].pm; //parent mass
		pm = (int64_t)(0.5 + rm * pt); //reduced parent mass
		pz = (int64_t)_l.spectra[s].pz; //parent charge
		idx.clear(); //initialize temporary list of kernel identifiers
		idi = 0;
		use2 = (pm > 1500000 and pz == 2); //check if parent mass large enough to consider z=2 fragments
		use3 = (pz > 2); //check if charge large enough to consider z=2 fragments
		ms.clear();
		mvs.clear(); //clear list of masses to check
		mvs.push_back(rm); //add the parent mass
		if (rm > 1500000.0) {
			mvs.push_back(rm - c13); //add the A1 mass
		}
		//loop though masses
		for (size_t a = 0; a < mvs.size(); a++) {
			rm = mvs[a];
			cpm = (int64_t)(0.5 + rm * pt); //current reduced parent mass
			//loop through possible parent mass values to compensate for round errors
			for (size_t n = 0; n < dvals.size(); n++) {
				d = dvals[n];
				pv.first = (int64_t)(cpm + d); //initialize (parent,fragment) pair
				if(_k.mvindex.find(pv.first) == _k.mvindex.end())	{ //bail if parent had no fragments in any kernel
					continue;
				}
				//loop through spectrum fragment (mass,intensity) pairs
				for(size_t b=0; b < _l.spectra[s].mis.size(); b++)	{
					m = _l.spectra[s].mis[b].first; //fragment mass
					m2 = m * 2; //fragment mass if fragment charge was 2
					pv.second = m; //set fragment mass in (mass,intensity) pair
					itk = _k.kindex.find(pv); //
					idx.clear();
					//check if pair was in kernels and record result
					if(itk !=  _k.kindex.end())	{ 
						idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
						idi = _l.spectra[s].mis[b].second;
					}
					else if(use2 and m2 > 100000)	{
						pv.second = m2;
						itk = _k.kindex.find(pv);
						if(itk !=  _k.kindex.end())	{
							idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
							idi = _l.spectra[s].mis[b].second;
						}
					}
					else if(use3)	{
						pv.second = m2;
						itk = _k.kindex.find(pv);
						if(itk !=  _k.kindex.end())	{
							idx.insert(idx.end(),_k.kindex[pv].begin(),_k.kindex[pv].end());
							idi = _l.spectra[s].mis[b].second;
						}
					}
					else	{
						continue;
					}
					//record kernel matching information
					ims.push_back(idi);
					ident.push_back(idx);
					ms.push_back(rm);
				}
			}
		}
		if(ident.empty())	{ //bail out if no identifications were found
			continue;
		}
		//check all fragment identifications and find best matches
		int64_t isum = 0;
		//find the sum of spectrum fragment intensities
		for(size_t a = 0; a < _l.spectra[s].mis.size(); a++)	{
			isum += _l.spectra[s].mis[a].second;
		}
		double total = (double)isum/3; //compensate for the method used to compensate for rounding errors
		set<int64_t> aok;
		//test to be sure that parent masses are within the parent mass ppm tolerance
		//mark kernel match as not ok if parent mass outside of allowed tolerance
		for(size_t b = 0; b < ident.size(); b++)	{
			for(size_t c = 0; c < ident[b].size(); c++)	{
				int64_t a = ident[b][c];
				if(fabs(ms[b]-_m[a]) < ppm*ms[b])	{
					aok.insert(a);
				}
			}
		}
		map<int64_t,int64_t> ans;
		map<int64_t,int64_t> aint;
		//iterate through all information recorded in the ident vector
		//and find the kernel identifiers corresponding to the best identification
		for(size_t b = 0; b < ident.size(); b++)	{
			set<int64_t> sv(ident[b].begin(),ident[b].end());
			set<int64_t>::iterator it = sv.begin();
			while(it != sv.end())	{
				if(aok.find(*it) == aok.end())	{
					ans[*it] = 0;
					aint[*it] = 0;
				}
				else if(ans.find(*it) != ans.end()) {
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
		int64_t mn = 0;
		vector<int64_t> mv;
		vector<int64_t> iv;
		map<int64_t,int64_t>::iterator itans = ans.begin();
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
			int64_t max_i = *max_element(iv.begin(),iv.end());
			for(size_t b = 0; b < mv.size(); b++)	{
				if(iv[b] < max_i)	{
					continue;
				}
				r.ks.push_back(mv[b]);
			}
			r.sn = z;
			r.peaks = mn;
			r.ri = max_i/total;
			r.pm = _l.spectra[s].pm;
			r.pz = _l.spectra[s].pz;
			r.sc = _l.spectra[s].sc;
			r.ions = (int64_t)_l.spectra[s].mis.size()/3;
			ids.push_back(r);
		}
		z += 1;
	}
	cout << endl;
	cout.flush();
	return true;
}

