/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/
/*
id records the information about a spectrum-to-peptide match in enough detail
to generate the output for that match
*/
class id
{
public:
	id(void)	{
		sn = 0; //scan number
		peaks = 0; //number of fragment ions identified
		ri = 0.0; //fraction of ion intensity identified
		pm = 0; //parent ion mass
		pz = 0; //parent ion charge
		sc = 0; //scan number
		ions = 0; //number of fragment ions in the spectrum
		rt = 0.0; //retention time
	}
	virtual ~id(void)	{}
	int32_t sn; //scan number
	int32_t peaks; //number of fragment ions identified
	vector<int32_t> ks; //vector of kernel unique identifiers
	double ri; //fraction of ion intensity identified
	double rt; //chromatographic retention time
	int32_t pm; //parent ion mass
	int32_t pz; //parent ion charge
	int32_t sc; //scan number
	int32_t ions; //number of fragment ions in the spectrum
	//copy operator
	id& operator=(const id &rhs)	{
		sn = rhs.sn;
		peaks = rhs.peaks;
		ri = rhs.ri;
		pm = rhs.pm;
		pz = rhs.pz;
		sc = rhs.sc;
		rt = rhs.rt;
		ions = rhs.ions;
		ks.clear();
		for(size_t k = 0; k < rhs.ks.size(); k++)	{
			ks.push_back(rhs.ks[k]);
		}
		return *this;
	}
};

class create_results
{
public:
	create_results(void) {}
	virtual ~create_results(void) {}
	//the create method checks the information from the spectrum file with
	//the information from the kernel file
	bool create(map<string,string>& _p,
			const load_spectra& _l,
			load_kernel& _k);
	// get size of ids vector
	inline int32_t size(void) { return (int32_t)ids.size(); }
	
	vector<id> ids; //identification information
};



