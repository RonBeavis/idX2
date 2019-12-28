/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/


#include <algorithm>
#include "rapidjson/document.h"
using namespace rapidjson; //namespace for the rapidjson methods

typedef std::pair <int32_t,int32_t> kPair; // type used to record kernel (parent,fragment) mass pairs

/*
kernels is a class that stores the indexes derived from
reading the kernel file with load_kernel. kindex
is the primary record of information. mvindex records the
parent ion masses that have at least one entry in kindex.
*/

class jsObject
{
public:
	jsObject(void)	{pm = 0; u = 0; h = 0; 
			pBuffer = new unsigned char[1024*16-1];
			pKey = new char[256];
	}
	virtual ~jsObject(void)	{delete pKey;delete pBuffer;}
	int32_t pm;
	int32_t u;
	int32_t h;
	unsigned int size;
	vector<int32_t> bs;
	vector<int32_t> ys;
	char *pKey;
	unsigned char *pBuffer;
	string key;
	void reset(void)	{
		pm = 0;
		u = 0;
		h = 0;
		bs.clear();
		ys.clear();
	}
};

class channel
{
public:
	channel(void) {lower = 0; delta = 0; dead = true; mv = 0;}
	virtual ~channel(void) {}
	int32_t lower;
	int32_t delta;
	bool dead;
	int32_t mv;
};

class kernels
{
public:
	kernels(void)	{ 
		clength = 3;
		for(size_t a = 0; a < clength; a++)	{
			kindex_a.push_back(phmap::flat_hash_map<kPair,vector<int32_t> >());
			mvindex_a.push_back(phmap::flat_hash_set<int32_t>());
		}

	}
	virtual ~kernels(void)	{}
	vector<phmap::flat_hash_map<kPair,vector<int32_t> > > kindex_a; //records the kernels containing the specified pair 
	vector<phmap::flat_hash_set<int32_t> > mvindex_a; //set of parent masses with at least one entry in kindex 
	void add_pair(kPair _v,size_t _a) {kindex_a[_a][_v] = vector<int32_t>();} //addes a new vector to kindex
	size_t clength;
	int32_t size(void)	{
		int32_t sz = 0;
		for(size_t a = 0; a < clength; a++)	{
			sz += kindex_a[a].size();
		}
		return sz;
	} // returns the size of kindex
	void clear(void) {
		for(size_t a = 0; a < clength; a++)	{
			kindex_a[a].clear();
			mvindex_a[a].clear();
		}
	} //resets kindex and mvindex
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
	bool load_binary(void);
	bool get_next(FILE *_pFile,jsObject& _js);
	string kfile; //path to the kernel file
	double fragment_tolerance; //fragment mass tolerance in mDa
	kernels kerns; //object that will contain kernel information
	map<int32_t,int32_t> pmindex; //object that will contain (pm,fmN) index
	set<int32_t> sp_set; //set of spectrum parent masses
	map<int32_t, set<int32_t> > hu_set;
	phmap::flat_hash_set<sPair> spairs; //map of spectrum (parent:fragment) mass pairs 
	void spectrum_pairs(load_spectra& _l)	{ //creates independent spairs and sp_set used in each thread
		for(size_t a = 0; a < _l.spectra.size();a++)	{
			sp_set.insert(_l.spectra[a].pm);
			spairs.insert(_l.spectra[a].spairs.begin(),_l.spectra[a].spairs.end());
		}
	}
	void add_hu(const int32_t _h,const int32_t _u)	{
		if(hu_set.find(_h) == hu_set.end())	{
			hu_set[_h] = set<int32_t>();
		}
		hu_set[_h].insert(_u);
	}
	vector<channel> channels;
};



