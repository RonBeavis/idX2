/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/

#include <algorithm>
#include <cmath>

class spectrum
{
public:
	spectrum(void)	{
		pm = 0;
		pi=0;
		pz=0;
		sc=0;
		pt=0.0;
		rt=0.0;
	}
	virtual ~spectrum(void)	{
		clear(); 
	}
	int32_t pm; //parent mass
	int32_t pi; //parent intensity
	int32_t pz; //parent charge
	int32_t sc; //scan number
	double rt; //run time
	int32_t isum; //sum of fragment intensities
	int32_t pks; //number of peaks
	double pt; //parent mass tolerance
	vector<pair<int32_t,int32_t>> mis; //fragment mass/intensity pairs
	phmap::flat_hash_set<sPair> spairs; //index of (parent,fragment) masses 
	string desc;//description of spectrum
	// reset object to its inital state
	bool clear()	{
		sc = 0;
		pm = 0; 
		pi = 0; 
		pz = 0; 
		sc=0; 
		desc = ""; 
		mis.clear();
		spairs.clear();
		return true;
	}
	//copy operator
	spectrum& operator=(const spectrum &rhs)	{ 
		mis.clear();
		pm = rhs.pm;
		pi = rhs.pi;
		pz = rhs.pz;
		pt = rhs.pt;
		rt = rhs.rt;
		pks = rhs.pks;
		desc = rhs.desc;
		sc = rhs.sc;
		isum = rhs.isum;
		pair<int32_t,int32_t> p;
		for(size_t i = 0; i < rhs.mis.size(); i++)	{
			p.first = rhs.mis[i].first;
			p.second = rhs.mis[i].second;
			mis.push_back(p);
		}
		spairs.clear();
		spairs.insert(rhs.spairs.begin(),rhs.spairs.end());
		return *this;
	}

	// condition converts the information in the MGF file into the data types and indexes
	// used by idX. _ires is the fragment ion mass tolerance and _l is the maximum number
	// of fragment ions that may be considered.

	bool condition(const int32_t _ires, const int32_t _l)	{
		double i_max = 0.0;
		const double res = 1.0/(double)_ires;
		int32_t m = 0;
		int32_t i = 0;
		//make sure the fragment ions are sorted by mass
		sort(mis.begin(), mis.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );
		//find maximum intensity of relevant ions
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			i = mis[a].second;		
			if(m < 150000 or (int32_t)fabs(pm-m) < 45000 or (int32_t)fabs(pm/pz- m) < 2000)	{
				continue;
			}
			if(i > i_max)	{
				i_max = (double)i;
			}
		}
		vector<int32_t> sMs;
		vector<int32_t> sIs;
		i_max = i_max/100.0;
		//clean out fragments that are too small, too close to the parent mass 
		//or have intensities < 1% of the maximum
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			if(m < 150000 or fabs(pm-m) < 45000.0 or fabs(pm/pz-m) < 2000.0)	{
				continue;
			}
			if((double)mis[a].second/i_max >= 1.00)	{
				i = mis[a].second;
				if(!sMs.empty())	{
					if(fabs(sMs.back() - m) < 500.0)	{
						if(sIs.back() < i)	{
							sMs.pop_back();
							sIs.pop_back();
							sMs.push_back(m);
							sIs.push_back(i);
						}
					}
					else	{
						sMs.push_back(m);
						sIs.push_back(i);
					}
				}
				else	{
					sMs.push_back(m);
					sIs.push_back(i);
				}
			}
		}
		//make a temporary vector of mass,intensity pairs
		vector<pair<int32_t,int32_t> > pMs;
		pair<int32_t,int32_t> p(0,0);
		for(size_t a = 0; a < sMs.size();a++)	{
			p.first = sMs[a];
			p.second = sIs[a];
			pMs.push_back(p);
		}
		//sort the temporary vector by fragment ion intensity
		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) {
							if(x.second > y.second) return true; 
							if(x.second == y.second) return x.first > y.first;
							return false;} );
		//estimate the maximum number of relevant ions
		int32_t max_l = 2 * (int32_t)((0.5 + (double)pm)/100000.0);
		if(max_l > _l)	{
			max_l = _l;
		}
		//remove ions below the relevance cutoff
		while(pMs.size() > (size_t)max_l)	{
			pMs.pop_back();
		}
		//re-sort the temporary vector by mass
		sort(pMs.begin(), pMs.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );
		//generate a normalized set of spectrum masses
		isum = 0;
		mis.clear();
		sPair spr;
		const double ptd = 1.0/70.0; 
		int32_t spr1 = (int32_t)(0.5+(double)pm*ptd); //calculate the reduced value used for indexing the parent ion mass
		spr.first = spr1;
		//create a new vector of fragment ion mass,intensity pairs using
		//the reduced value corresponding to the fragment ion mass tolerance.
		//3 values are recorded, to compensate for errors introduced
		//by the fragment mass reduction rounding process.
		pks = 0;
		int32_t pks2 = 0;
		int32_t pm_limit = (pm - 100000)/2;
		int32_t j = 0;
		int32_t p1 = 0;
		for(size_t a=0; a < pMs.size();a++)	{
			m = pMs[a].first;
			isum += pMs[a].second;
			pks++;
			p1 = (int32_t)(0.5+(double)m*res);
			p.first = p1;
			p.second = pMs[a].second;
			for(i = -1; i < 2; i++)	{
				p.first = p1 + i;
				mis.push_back(p);
				for(j = -1; j < 2; j++)	{
					spr.first = spr1 + j;
					spr.second = p.first;
					spairs.insert(spr);
				}
			}

			// add z=2 fragment masses as a possibility if appropriate
			if(pz > 1 and m > 400000 and m < pm_limit)	{
				pks2++;
				p1 = (int32_t)(0.5+(double)(m*2-1003)*res);
				p.first = p1;
				p.second = pMs[a].second;
				for(i = -1; i < 2; i++)	{
					p.first = p1 + i;
					mis.push_back(p);
					for(j = -1; j < 2; j++)	{
						spr.first = spr1 + j;
						spr.second = p.first;
						spairs.insert(spr);
					}
				}
			}
		}
		// add 3 to pks if there are +2 peaks recorded and z = 2
		// add pks2/3 if z > 2
		if(pks2 > 0 and pz == 2)	{
			pks += 3;
		}
		else if(pks2 > 0 and pz > 2)	{
			pks += (int32_t)(0.5 + (float)pks2/3.0);
		}
		// re-sort pairs by accending fragment mass
		sort(mis.begin(), mis.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );

		auto itmis = mis.begin();
		int32_t erased = 0;
		// remove any repeated values that may have gotten recorded
		while(itmis != mis.end())	{
			if((itmis + 1) != mis.end())	{
				if(itmis->first == (itmis+1)->first)	{
					if(itmis->second >= (itmis+1)->second)	{
						itmis = mis.erase(itmis+1);
						itmis--;
						erased++;
					}
					else	{
						itmis = mis.erase(itmis);
						erased++;
					}
				}
				else	{
					itmis++;
				}
			}
			else	{
					itmis++;
			}
		}
		return true;
	}
};

// load_spectra takes a spectrum file path, reads the spectra and stores the relevant information
// in a vector of spectrum objects

class load_spectra
{
public:
	load_spectra(void)	{
		validation = "";
		skipped = 0;
	}
	virtual ~load_spectra(void) {}

	// load spectra using the information in _p
	bool load(map<string,string>& _p,load_kernel& _lk);
	bool load_mgf(map<string,string>& _p,load_kernel& _lk);
	bool load_cmn(map<string,string>& _p,load_kernel& _lk);
	// enforces the maximum number of spectra to use by truncating
	// the spectra vector
	void set_max(const int32_t _max)	{
		if(spectra.size() <= (size_t)_max)	{
			return;
		}
		spectra.erase(spectra.begin()+_max,spectra.end());
		return;
	}
	// recover memory when spectra no longer required
	void clean_up(void)	{
		spectra.clear();
		spectra.shrink_to_fit();
	}
	vector<spectrum> spectra; // spectrum information for identifications
	string validation; // SHA256 hash of input file
	int32_t skipped;
};



