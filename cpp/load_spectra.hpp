/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/

/*
load_spectra is a specialty class for taking a file in MGF format and
creating a vector of spectrum objects containing the information
relevant for the purposes of the idX algorithm. All masses recorded
correspond to neutral molecules and are recorded in millidaltons.
*/
#include <algorithm>
#include <cmath>

typedef std::pair <int64_t,int64_t> sPair; //type used to record (parent,fragment) pairs

class spectrum
{
public:
	spectrum(void)	{pm = 0;pi=0;pz=0;sc=0;pt=0.0;rt=0.0;}
	virtual ~spectrum(void)	{ clear(); }
	int64_t pm; //parent mass
	int64_t pi; //parent intensity
	int64_t pz; //parent charge
	int64_t sc; //scan number
	double rt; //run time
	int64_t isum; //sum of fragment intensities
	double pt; //parent mass tolerance
	vector<pair<int64_t,int64_t>> mis; //fragment mass/intensity pairs
	phmap::flat_hash_set<sPair> spairs; //index of (parent,fragment) masses 
	string desc;//description of spectrum
	bool clear()	{sc = 0;pm = 0; pi = 0; pz = 0; sc=0; desc = ""; mis.clear();spairs.clear();return true;}
	spectrum& operator=(const spectrum &rhs)	{ //copy operator
		mis.clear();
		pm = rhs.pm;
		pi = rhs.pi;
		pz = rhs.pz;
		pt = rhs.pt;
		rt = rhs.rt;
		desc = rhs.desc;
		sc = rhs.sc;
		pair<int64_t,int64_t> p;
		for(size_t i = 0; i < rhs.mis.size(); i++)	{
			p.first = rhs.mis[i].first;
			p.second = rhs.mis[i].second;
			mis.push_back(p);
		}
		spairs.clear();
		spairs.insert(rhs.spairs.begin(),rhs.spairs.end());
		return *this;
	}
	//condition() converts the information in the MGF file into the data types and indexes
	//used by idX. _ires is the fragment ion mass tolerance and _l is the maximum number
	//of fragment ions that may be considered.
	bool condition(int64_t _ires, int64_t _l)	{
		double i_max = 0.0;
		const double res = 1.0/(double)_ires;
		int64_t m = 0;
		int64_t i = 0;
		//make sure the fragment ions are sorted by mass
		sort(mis.begin(), mis.end(), 
               		[](const auto& x, const auto& y) { return x.first < y.first; } );
		//find maximum intensity of relevant ions
		for(size_t a = 0; a < mis.size(); a++)	{
			m = mis[a].first;
			i = mis[a].second;		
			if(m < 150000 or (int64_t)fabs(pm-m) < 45000 or (int64_t)fabs(pm/pz- m) < 2000)	{
				continue;
			}
			if(i > i_max)	{
				i_max = (double)i;
			}
		}
		vector<int64_t> sMs;
		vector<int64_t> sIs;
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
		//make a temporary vector of mass,intensity paris
		vector<pair<int64_t,int64_t> > pMs;
		pair<int64_t,int64_t> p(0,0);
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
		int64_t max_l = 2 * (int64_t)((0.5 + (double)pm)/100000.0);
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

//
//		generate a normalized set of spectrum masses
//
		int64_t is = 0;
		mis.clear();
		sPair pr;
		const double ptd = 1.0/70.0; 
		pr.first = (int64_t)(0.5+(double)pm*ptd); //calculate the reduced value used for indexing the parent ion mass
		//create a new vector of fragment ion mass,intensity pairs using
		//the reduced value corresponding to the fragment ion mass tolerance.
		//3 values are recorded, to compensate for errors introduced
		//by the fragment mass reduction rounding process.
		for(size_t a=0; a < pMs.size();a++)	{
			m = pMs[a].first;
			p.first = (int64_t)(0.5+(double)m*res);
			pr.second = p.first;
			spairs.insert(pr);
			pr.first += 1;
			spairs.insert(pr);
			pr.first -= 2;
			spairs.insert(pr);
			pr.first += 1;
			p.second = pMs[a].second;
			mis.push_back(p);
			p.first -= 1;
			pr.second = p.first;
			spairs.insert(pr);
			pr.first += 1;
			spairs.insert(pr);
			pr.first -= 2;
			spairs.insert(pr);
			pr.first += 1;
			mis.push_back(p);
			p.first += 2;
			pr.second = p.first;
			spairs.insert(pr);
			pr.first += 1;
			spairs.insert(pr);
			pr.first -= 2;
			spairs.insert(pr);
			pr.first += 1;
			mis.push_back(p);
			is += 3*p.second;
		}
/*		if(sc == 3485)	{
			auto itp = spairs.begin();
			while(itp != spairs.end())	{
				cout << itp->first << "\t" << itp->second << "\t" << pm << endl;
				itp++;
			}
		} */
		isum = is; //note: this value is 3x the actual sum of intensities
		return true;
	}
};

class load_spectra
{
public:
	load_spectra(void);
	virtual ~load_spectra(void);
//
//	load spectra using the information in _p
//
	bool load(map<string,string>& _p);
//
//	enforces the maximum number of spectra to use by truncating
//	the spectra vector
//
	void set_max(const int64_t _max)	{
		if(spectra.size() <= (size_t)_max)	{
			return;
		}
		spectra.erase(spectra.begin()+_max,spectra.end());
		return;
	}
	vector<spectrum> spectra; // spectrum information for identifications

};



