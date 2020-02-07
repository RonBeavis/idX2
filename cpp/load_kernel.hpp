/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/


#include <algorithm>
#include "rapidjson/document.h"
#include <deque>
using namespace rapidjson; //namespace for the rapidjson methods

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
	size_t clength; // number of mass channels (a constant)
	int32_t size(void)	{ // retrieve the number of values in kindex_a
		size_t sz = 0;
		for(size_t a = 0; a < clength; a++)	{
			sz += kindex_a[a].size();
		}
		return (int32_t)sz;
	}
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
	bool load(void); // load from a text JSON file
	bool load_binary(void); // load from a binary JSON file
	bool get_next(ifstream& _ifs,jsObject& _js); //retrieve the next JSON object from a binary stream
	string kfile; //path to the kernel file
	double fragment_tolerance; //fragment mass tolerance in mDa
	kernels kerns; //object that will contain kernel information
	set<int32_t> sp_set; //set of spectrum parent masses
	map<int32_t, set<int32_t> > hu_set; //mapping unique kernel ids to homologous kernel ids
	string validation; // hex-encoded SHA256 validation string for kernel file
	void clean_up(void)	{ // free up memory and reset containers to zero length
		sp_set.clear();
		kerns.kindex_a.clear();
		kerns.mvindex_a.clear();
		spairs.clear();
	}
	phmap::flat_hash_set<sPair> spairs; //map of spectrum (parent:fragment) mass pairs 
	void add_hu(const int32_t _h,const int32_t _u)	{ // add a homologous kernel id/unique kernel id pair
		if(hu_set.find(_h) == hu_set.end())	{
			hu_set[_h] = set<int32_t>();
		}
		hu_set[_h].insert(_u);
	}
	void check_and_update(kPair& _pr,int32_t _u)	{ //check a new pair/kernel unique id combination
		for(size_t b = 0; b < clength;b++)	{
			if(channels[b].dead)	{
				continue;
			}
			_pr.first = channels[b].mv;
			if(spairs.find(_pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			_pr.first = channels[0].mv;
			if(kerns.kindex_a[b].find(_pr) == kerns.kindex_a[b].end())	{ 
				kerns.add_pair(_pr,b); //create a new vector for the object
			}
			kerns.mvindex_a[b].insert(_pr.first); //add parent mass to set
			kerns.kindex_a[b][_pr].push_back(_u); //add kernel id to vector
		}
		return;
	}
	vector<channel> channels; // parent ion mass derived set of masses, e.g. [A0,A1,A2]
	size_t clength;
};



