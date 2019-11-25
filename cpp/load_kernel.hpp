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
	phmap::parallel_flat_hash_map<kPair,vector<int64_t> > kindex; //records the kernels containing the specified pair
	phmap::parallel_flat_hash_set<int64_t> mvindex; //set of parent masses with at least one entry in kindex
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
	bool load(map<string,string>& _p,load_spectra& _l,kernels& _k,map<int64_t,int64_t>& _m);
};



