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

// object to store kernel information from a binary JSON file object
class jsObject
{
public:
	jsObject(void)	{
		pm = 0;
		u = 0;
		h = 0; 
		pBuffer = new unsigned char[1024*16-1];
		pKey = new char[256];
	}
	virtual ~jsObject(void)	{delete pKey;delete pBuffer;}
	int32_t pm; // parent mass
	int32_t u; // uid for kernel
	int32_t h; // uid for first instance of a specifically modified peptide sequence
	unsigned int size; // temporary value of an element's size
	vector<int32_t> bs; // b fragment information
	vector<int32_t> ys; // y fragment information
	char *pKey; // buffer to store the key value of a binary JSON element
	unsigned char *pBuffer; // buffer to read a binary JSON object
	string key; // key value for a binary JSON element
	// reset object to initial state, leaving pBuffer and pKey unchanged
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
	channel(void) {
		lower = 0;
		delta = 0;
		dead = true;
		mv = 0;
	}
	virtual ~channel(void) {}
	int32_t lower; // lower parent mass limit for a channel to be operative
	int32_t delta; // mass offset of a channel
	bool dead; // if dead = true, channel not available for use
	int32_t mv; // reduced parent mass value of a channel
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
	size_t clength; // number of mass channels (a constant)
	// retrieve the number of values in kindex_a
	int32_t size(void)	{
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
	load_kernel(void)		{
		const int32_t c13 = 1003; //difference between the A0 and A1 peaks
		kfile = "";
		fragment_tolerance = 0;
		// set up the A0, A1 and A2 parent ion mass channels
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
	virtual ~load_kernel(void) {}
	//load uses spectrum information and the kernel file specified in _p
	//to retrieve information about candidate kernels that will be used in the peptide-sequence matching process

	// load from a text JSON file
	bool load(void);
	// load from a binary JSON file
	bool load_binary(void);
	//retrieve the next JSON object from a binary stream
	bool get_next(ifstream& _ifs,jsObject& _js); 
	string kfile; //path to the kernel file
	double fragment_tolerance; //fragment mass tolerance in mDa
	kernels kerns; //object that will contain kernel information
	set<int32_t> sp_set; //set of spectrum parent masses
	map<int32_t, set<int32_t> > hu_set; //mapping unique kernel ids to homologous kernel ids
	string validation; // hex-encoded SHA256 validation string for kernel file
	phmap::flat_hash_set<sPair> spairs; //map of spectrum (parent:fragment) mass pairs 
	vector<channel> channels; //parent ion mass derived set of masses, e.g. [A0,A1,A2]
	size_t clength; //number of elements in channels
	// free up memory and reset containers to zero length
	void clean_up(void)	{ 
		sp_set.clear();
		kerns.kindex_a.clear();
		kerns.kindex_a.shrink_to_fit();
		kerns.mvindex_a.clear();
		kerns.mvindex_a.shrink_to_fit();
		spairs.clear();
	}
	// add a homologous kernel id/unique kernel id pair
	inline void add_hu(const int32_t _h,const int32_t _u)	{
		if(hu_set.find(_h) == hu_set.end())	{
			hu_set[_h] = set<int32_t>();
		}
		hu_set[_h].insert(_u);
	}
	//check a new pair/kernel unique id combination and record the information in the kernels object
	inline void check_and_update(const kPair& _pr,const int32_t _u)	{
		kPair pr = _pr;
		for(size_t b = 0; b < clength;b++)	{
			if(channels[b].dead)	{
				continue;
			}
			pr.first = channels[b].mv;
			if(spairs.find(pr) == spairs.end())	{ //bail if pair not in spectrum pairs
				continue;
			}
			pr.first = channels[0].mv;
			if(kerns.kindex_a[b].find(pr) == kerns.kindex_a[b].end())	{ 
				kerns.kindex_a[b][pr] = vector<int32_t>();
			}
			kerns.mvindex_a[b].insert(pr.first); //add parent mass to set
			kerns.kindex_a[b][pr].push_back(_u); //add kernel id to vector
		}
		return;
	}
};



