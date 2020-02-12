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
#include <set>
#include <vector>
#include <deque>
#include "parallel_hashmap/phmap.h"
using namespace std;
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs
#include "load_kernel.hpp"
#include "load_spectra.hpp"
#include "create_results.hpp"

//
// create_results takes the information recorded from the spectrum file and kernel file
// and generates spectrum-to-peptide matches
//

bool create_results::create(map<string, string>& _p, // parameter strings
		const load_spectra& _l, //spectrum informaiton
		load_kernel& _k) { // kernel information
	int32_t z = 1; // serial number for PSM: used if spectrum::sc is not available
	const double pt = 1.0/70.0;
	vector<int32_t> dvals{ -1,0,1 }; // values to compensate for integer rounding effects
	// initialize variables for the identification process
	int32_t d = 0;
	int32_t m = 0;
	int32_t idi = 0;
	int32_t cpm = 0;
	vector<int32_t> idx;
	vector<vector<int32_t> > ident;
	vector<int32_t> ims;
	auto itk = _k.kerns.kindex_a[0].end(); // store value for subsequent map finds
	kPair pv; // a kernel (parent,fragment) pair
	const size_t clength = _k.channels.size();
	size_t a = 0;
	size_t n = 0;
	size_t b = 0;
	int32_t mi = 0;
	// loop through spectra
	for (size_t s = 0; s < _l.spectra.size(); s++) {
		// output keep-alive text for logging
		if (s != 0 and s % 500 == 0) {
			cout << '.';
			cout.flush();
		}
		if (s != 0 and s % 5000 == 0) {
			cout << " " << s << '\n';
			cout.flush();
		}
		ident.clear(); // initialize the temporary identification vector
		ims.clear(); // initialize the temporary mass vector
		mi = _l.spectra[s].pm; // record the parent mass in a local integer
		idx.clear(); // initialize temporary list of kernel identifiers
		idi = 0; // initialize the fragment ion intensity record
		for(a = 0; a < clength; a++)	{
			if(mi > _k.channels[a].lower)	{
				_k.channels[a].mv = mi - _k.channels[a].delta; //add the parent mass
				_k.channels[a].dead = false;
			}
			else	{
				_k.channels[a].dead = true;
			}
		}
		// loop though masses
		for (a = 0; a < clength; a++) {
			if(_k.channels[a].dead)	{ //skip channel if it is not active
				continue;
			}
			cpm = (int32_t)(0.5 + (double)_k.channels[a].mv * pt); //current reduced parent mass
			// loop through possible parent mass values to compensate for rounding errors
			for (n = 0; n < dvals.size(); n++) {
				d = dvals[n];
				pv.first = cpm + d; //initialize (parent,fragment) pair
				// bail if parent had no fragments in any kernel
				if(_k.kerns.mvindex_a[a].find(pv.first) == _k.kerns.mvindex_a[a].end())	{
					continue;
				}
				//loop through spectrum fragment (mass,intensity) pairs
				for(b = 0; b < _l.spectra[s].mis.size(); b++)	{
					m = _l.spectra[s].mis[b].first; //fragment mass
					pv.second = m; //set fragment mass in (mass,intensity) pair
					itk = _k.kerns.kindex_a[a].find(pv); //
					idx.clear();
					// check if pair was in kernels and record result
					if(itk !=  _k.kerns.kindex_a[a].end())	{ 
						idx.insert(idx.end(),_k.kerns.kindex_a[a][pv].begin(),_k.kerns.kindex_a[a][pv].end());
						idi = _l.spectra[s].mis[b].second;
					}
					else	{
						continue;
					}
					// record kernel matching information
					ims.push_back(idi); // matched fragment ion intensity
					ident.push_back(idx); // vector of matched kernels
				}
			}
		}
		if(ident.empty())	{ //bail out if no identifications were found, do not update the PSM serial number "z"
			continue;
		}

		// ATTENTION: This next section is where the peptide identifications are assigned.
		//		The point is to create a set of the ids with the maximum ion count
		//		and from that set record the ids that have the maximum summed intensity.

		// map of kernel numbers and the number of times the kernel occurs for this spectrum (ion count)
		map<int32_t,int32_t> ans;
		// map of kernel numbers and the summed intensity of fragment ions observed for that kernel (total intensity)	
		map<int32_t,int32_t> aint;	

		// loop through all information recorded in the ident (kernel list) and ims (intensity list) vectors
		// and find the kernel identifiers corresponding to the best identification
		for(b = 0; b < ident.size(); b++)	{
			// create a set of the kernels present in the ident vector for each fragment ion
			// the set prevents double counting if a kernel has been assigned to both
			// +1 and +2 fragments
			set<int32_t> sv(ident[b].begin(),ident[b].end());
			// iterate through the kernel identifiers (*it) and map the kernel numbers with
			// the number of times they occur for the spectrum (ans) and their associated
			// fragment ion intensity (aint)
			set<int32_t>::iterator it = sv.begin(); // get the beginning of the set of kernels
			while(it != sv.end())	{
				if(ans.find(*it) != ans.end()) {
					ans[*it] += 1; // update the ion count
					aint[*it] += ims[b]; // update the intensity sum
				}
				else	{
					ans[*it] = 1; // initialize the ion count
					aint[*it] = ims[b]; // initialize the intensity sum
				}
				it++;
			}
		}
		int32_t mn = 0; // best ion count
		vector<int32_t> mv; // record of kernels with the best ion count
		vector<int32_t> iv; // record of summed intensities with the best ion count

		// NOTE: the records in mv and iv with the same index correspond to each other and could be
		//       represented in a single vector with a pair of values 

		// iterate through the map of kernels and ion counts to find the maximum ion count
		map<int32_t,int32_t>::iterator itans = ans.begin(); // get an iterator to the first kernel record
		//arrange identification information
		while(itans != ans.end())	{
			// if the ion count (itans->second) beats the current record, reset the vectors
			if(itans->second > mn)	{
				mn = itans->second; // update the best ion count
				mv.clear(); // reset the kernel vector
				mv.push_back(itans->first); // add the current kernel to the vector
				iv.clear(); // clear the intensity vector
				iv.push_back(aint[itans->first]); // add the current summed intensity to the vector
			}
			// if the ion count (itans->second) equals the current record, record it
			else if(itans->second == mn)	{
				mv.push_back(itans->first); // add the current kernel to the vector
				iv.push_back(aint[itans->first]); // add the current summed intensity to the vector
			}
			itans++;
		}
		// record identification information if at least 5 fragments were matched
		// mn = 5 is below what would normally be considered appropriate for a match
		// low quality ids will be filtered out by the physical model used in create_output
		if(mn > 4)	{
			id r;
			// find the maximum intensity in the set of solutions recorded in mv and iv
			int32_t tot_i = *max_element(iv.begin(),iv.end());
			// filter out solutions with less than tot_i intensity
			for(b = 0; b < mv.size(); b++)	{
				if(iv[b] < tot_i)	{
					continue;
				}
				r.ks.push_back(mv[b]); // record kernel numbers that pass this test
			}
			r.sn = z; // record serial number of the PSM
			r.peaks = mn; // record the best ion count
			r.ri = (double)tot_i/(double)_l.spectra[s].isum; // record the best summed fragment ion intensity
			r.pm = _l.spectra[s].pm; // record the spectrum parent mass
			r.pz = _l.spectra[s].pz; // record the spectrum parent ion charge
			r.sc = _l.spectra[s].sc; // record the spectrum scan number
			r.rt = _l.spectra[s].rt; // record the spectrum retention time
			r.ions = _l.spectra[s].pks; // record the number of peaks to be used in the physical model
			ids.push_back(r); // save the id record for use in create_output
		}
		z += 1; // increment PSM serial number
	}
	// guarantee that the log text output is terminated with a new line and flushed
	cout << '\n';
	cout.flush();
	return true;
}

