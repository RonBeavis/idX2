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
	id(void)	{}
	virtual ~id(void)	{}
	int64_t sn; //scan number
	int64_t peaks; //number of fragment ions identified
	vector<int64_t> ks; //vector of kernel unique identifiers
	double ri; //fraction of ion intensity identified
	int64_t pm; //parent ion mass
	int64_t pz; //parent ion charge
	int64_t sc; //scan number
	int64_t ions; //number of fragment ions in the spectrum
	id& operator=(const id &rhs)	{
		sn = rhs.sn;
		peaks = rhs.peaks;
		ri = rhs.ri;
		pm = rhs.pm;
		pz = rhs.pz;
		sc = rhs.sc;
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
	create_results(void);
	virtual ~create_results(void);
	//the create method checks the information from the spectrum file with
	//the information from the kernel file
	bool create(map<string,string>& _p,
			load_spectra& _l,
			kernels& _k,
			map<int64_t,int64_t>& _m);
	int64_t size(void) { return (int64_t)ids.size(); }
	
	vector<id> ids; //identification information
};



