/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/


#include <algorithm>
typedef std::pair <int64_t,int64_t> kPair; // type used to record kernel (parent,fragment) mass pairs

/*
kernels is a class that stores the indexes derived from
reading the kernel file with load_kernel. kindex
is the primary record of information. mvindex records the
parent ion masses that have at least one entry in kindex.
*/
class kernels
{
public:
	kernels(void)	{}
	virtual ~kernels(void)	{}
	phmap::flat_hash_map<kPair,vector<int64_t> > kindex; //records the kernels containing the specified pair 
	phmap::flat_hash_set<int64_t> mvindex; //set of parent masses with at least one entry in kindex 
	void add_pair(kPair _v) {kindex[_v] = vector<int64_t>();} //addes a new vector to kindex
	int64_t size(void)	{ return (int64_t)kindex.size();} // returns the size of kindex
	void clear(void) { kindex.clear(); mvindex.clear();} //resets kindex and mvindex
};

/*
load_kernel is a specialty class that loads kernels into memory from
a file that was specified on the command line. 
*/

class load_kernel
{
public:
	load_kernel(void);
	virtual ~load_kernel(void);
	//load uses spectrum information and the kernel file specified in _p
	//to retrieve information about candidate kernels that will be used in the peptide-sequence matching process
	bool load(void);
	string kfile; //path to the kernel file
	double fragment_tolerance; //fragment mass tolerance in mDa
	kernels kerns; //object that will contain kernel information
	map<int64_t,int64_t> pmindex; //object that will contain (pm,fmN) index
	int thread;
	int threads;
	load_kernel& operator+=(const load_kernel &rhs)	{ //copy operator
		kfile = rhs.kfile;
		pmindex.insert(rhs.pmindex.begin(),rhs.pmindex.end());
		kerns.mvindex.insert(rhs.kerns.mvindex.begin(),rhs.kerns.mvindex.end());
		auto itk = rhs.kerns.kindex.begin();
		kPair kp;
		while(itk != rhs.kerns.kindex.end())	{
			kp = itk->first;
			if(kerns.kindex.find(kp) == kerns.kindex.end())	{
				kerns.add_pair(kp);
			}
			kerns.kindex[kp].insert(kerns.kindex[kp].end(),itk->second.begin(),itk->second.end());
			itk++;
		}
		auto ithu = rhs.hu_set.begin();
		while(ithu != rhs.hu_set.end())	{
			if(hu_set.find(ithu->first) == hu_set.end())	{
				hu_set[ithu->first] = set<int64_t>();
			}
			hu_set[ithu->first].insert(ithu->second.begin(),ithu->second.end());
			ithu++;
		}
		return *this;
	}
	set<int64_t> sp_set; //set of spectrum parent masses
	map<int64_t, set<int64_t> > hu_set;
	phmap::flat_hash_set<sPair> spairs; //map of spectrum (parent:fragment) mass pairs 
	void spectrum_pairs(load_spectra& _l)	{ //creates independent spairs and sp_set used in each thread
		for(size_t a = 0; a < _l.spectra.size();a++)	{
			sp_set.insert(_l.spectra[a].pm);
			spairs.insert(_l.spectra[a].spairs.begin(),_l.spectra[a].spairs.end());
		}
	}
	void add_hu(const int64_t _h,const int64_t _u)	{
		if(hu_set.find(_h) == hu_set.end())	{
			hu_set[_h] = set<int64_t>();
		}
		hu_set[_h].insert(_u);
	}
};



