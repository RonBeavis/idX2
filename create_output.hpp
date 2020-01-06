/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Identifies kernels corresponding to spectra
#
*/
/*
hypergeom facilitates the calculation of the PDF of a specified hypergeometric distribution
*/
#include "rapidjson/document.h"
using namespace rapidjson;

class mod
{
public:
	mod(void)	{
		pos = 0; //residue position in protein coordinates
		res = ""; //the single character abbreviation for the genome-encoded residue at pos
		mass = 0; //the mass change associated with the PTM
	}
	virtual ~mod(void)	{}
	int32_t pos; //residue position in protein coordinates
	string res; //the single character abbreviation for the genome-encoded residue at pos
	int32_t mass; //the mass change associated with the PTM
	mod& operator=(const mod &rhs)	{ //copy operator
		res = rhs.res;
		pos = rhs.pos;
		mass = rhs.mass;
		return *this;
	}
	bool operator<( const mod& rhs ) const { //less than operator, used for formating output
		return pos < rhs.pos; 
	}
};

class osObject
{
public:
	osObject(void)	{pm = 0; u = 0; h = 0; pz = 0;
			pBuffer = new unsigned char[1024*16-1];
			pKey = new char[256];
	}
	virtual ~osObject(void)	{delete pKey;delete pBuffer;}
	int32_t pm;
	int32_t u;
	int32_t h;
	int32_t pz;
	int32_t beg;
	int32_t end;
	vector<int32_t> ns;
	mod mod_temp;
	vector<mod> mods;
	vector<mod> savs;
	unsigned int size;
	char *pKey;
	unsigned char *pBuffer;
	string key;
	string seq;
	string pre;
	string post;
	string lb;
	void reset(void)	{
		pm = 0;
		u = 0;
		h = 0;
		pz = 0;
		ns.clear();
		mods.clear();
		savs.clear();
		key.clear();
		seq.clear();
		pre.clear();
		post.clear();
		lb.clear();
	}
};

class hypergeom	{
public:
	// specify hypergeometric distribution parameters
	hypergeom(int32_t _n,int32_t _r,int32_t _N)	{n = _n; r = _r; N = _N;}
	virtual ~hypergeom() {}
	int32_t n;
	int32_t N;
	int32_t r;
	double pdf(int32_t k)	{ //calculate the PDF
		double top = st(n)+st(r)+st(N-n)+st(N-r);
		double bottom = st(N)+st(k)+st(n-k)+st(r-k)+st(N-n-r+k);
		double lp = top - bottom;
		return exp(lp);
	}
	double st(int32_t _n)	{ //either directly calculate the log of a factorial or use Sterling's approximation
		if(_n < 50)	{
			double f = 1.0;
			for(int32_t i = 1; i <= _n; i++)    {
        			f *= (double)i;
   			 }
			return log(f);
		}
		double n = (double)_n;
		double pi = 3.1416;
		return (n*log(n) - n + log(sqrt(2*pi*n)));
	}
};
/*
mod records a PTM as a residue, its protein coordinates and the mass of the modification
*/
/*
create_output takes peptide-to-spectrum matches and generates a tab-separated value output
file describing those identifications.
*/
class create_output
{
public:
	create_output(void);
	virtual ~create_output(void);
	// generates a TSV formated output file from results using a JSON kernel
	bool create(map<string,string>& _p,create_results& _cr, map<int32_t, set<int32_t> >& _hu);
	// generates a TSV formated output file from results using a JSON binary kernel
	bool create_binary(map<string,string>& _p,create_results& _cr, map<int32_t, set<int32_t> >& _hu);
	bool dump_meta(map<string,string>& _p);
	int32_t roundf(double _x)	{ 
		if(_x < 0.0) return (int32_t)(_x - 0.5);
		return (int32_t)(_x + 0.5);
	}
private:
	bool load_mods(void); //loads a map with strings associated with particular modification masses
	bool find_window(void); //determines a valid window for results, in parent mass ppm
	int32_t get_cells(double _pm,int32_t _res); //retrieves one of the parameters necessary for a hypergeometric model
	//calculates a probability model for a particular identification
	bool apply_model(int32_t _r,string& _s,double _pm,int32_t _ions,int32_t _lspectrum,double& pscore,double& p);
	bool create_line(id& _s, double _pm, double _d, double _ppm, double _score, Document& _js, string& _line);
	bool create_line_binary(id& _s, double _pm, double _d, double _ppm, double _score, osObject& _js, string& _line);
	bool create_header_line(string& _h); // generates a TSV formated header line for output
	bool get_next(FILE *_pFile,osObject& _os); // gets the next osObject from a JSON binary file
	bool dump_lines(string& _ofile,double _tp); // serializes odict lines into a file
	int32_t low; // lower value for the ppm window calculated in find_window
	int32_t high; // upper value for the ppm window calculated in find_window
	map<int32_t,id> sv;
	map<int32_t,set<int32_t> > sdict;
	map<int32_t,string> mt;
	map<int32_t,vector<string> > odict;
	vector<int32_t> ppms;
	map<string,string> info;
	map<int32_t,int32_t> ppm_map;
	string fragmentation;
	const int32_t c13 = 1003; //mass difference between the A1 and A0 peaks
	map<int32_t,vector<double> > distribution;
	//retrieves the ppm column from a formatted output string
	int32_t get_scan(string& t)	{
		size_t s = t.find("\t");
		s = t.find("\t",s+1);
		return (int32_t)atoi(t.c_str());
	}
	double get_ppm(string& t)	{
		size_t s = t.find("\t");
		s = t.find("\t",s+1);
		s = t.find("\t",s+1);
		s = t.find("\t",s+1);
		return atof((t.substr(s+1,t.size()-1)).c_str());
	}
	//loads a structure with the number of cells available for a hypergeometic model, based one
	//parent mass and fragment mass tolerance, either 20, 50 or 400 mDa.
	bool load_distribution(void)	{
		distribution.clear();
		vector<double> v = {1.3,1.0,0.6};
		distribution[800] = v;
		v = {3.2,2.0,0.8};
		distribution[900] = v;
		v = {4.9,2.7,1.0};
		distribution[1000] = v;
		v = {5.6,3.1,1.0};
		distribution[1100] = v;
		v = {6.2,3.3,1.1};
		distribution[1200] = v;
		v = {6.6,3.5,1.1};
		distribution[1300] = v;
		v = {7.1,3.8,1.2};
		distribution[1400] = v;
		v = {7.5,4.0,1.2};
		distribution[1500] = v;
		v = {7.8,4.1,1.3};
		distribution[1600] = v;
		v = {8.1,4.3,1.3};
		distribution[1700] = v;
		v = {8.2,4.4,1.3};
		distribution[1800] = v;
		v = {8.4,4.5,1.3};
		distribution[1900] = v;
		v = {8.5,4.6,1.3};
		distribution[2000] = v;
		v = {8.4,4.7,1.4};
		distribution[2100] = v;
		v = {8.6,4.7,1.4};
		distribution[2200] = v;
		v = {8.3,4.7,1.4};
		distribution[2300] = v;
		v = {8.2,4.6,1.4};
		distribution[2400] = v;
		v = {8.2,4.7,1.4};
		distribution[2500] = v;
		v = {8.1,4.8,1.4};
		distribution[2600] = v;
		v = {7.8,4.7,1.4};
		distribution[2700] = v;
		v = {7.5,4.7,1.4};
		distribution[2800] = v;
		v = {7.3,4.6,1.4};
		distribution[2900] = v;
		v = {7.3,4.7,1.5};
		distribution[3000] = v;
		v = {6.9,4.5,1.4};
		distribution[3100] = v;
		v = {6.5,4.4,1.4};
		distribution[3200] = v;
		v = {6.4,4.4,1.4};
		distribution[3300] = v;
		v = {5.8,4.1,1.4};
		distribution[3400] = v;
		v = {5.0,3.8,1.4};
		distribution[3500] = v;
		v = {5.5,4.1,1.4};
		distribution[3600] = v;
		v = {4.7,3.6,1.4};
		distribution[3700] = v;
		v = {4.2,3.3,1.4};
		distribution[3800] = v;
		v = {3.8,3.1,1.3};
		distribution[3900] = v;
		v = {3.7,3.0,1.3};
		distribution[4000] = v;
		return true;
	}
};



