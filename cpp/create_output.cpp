/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Creates output files based on a PSM identification session
#

create_output is an object to perform all of the tasks necessary to generate
a useful output file specified on the command line

create_output takes a list of identified kernels and reads through the
original kernel file to find them. It then constructs a line of text
describing the results in tab-separated value format and stores that
line in a text file.

create_output has methods that use either JSON formatted kernels or
binary JSON formatted kernels. The binary format was made necessary
because none of the JSON reading modules tested ran quickly enough
on Microsoft Windows. Because of this, using binary JSON kernels is
recommended when running the software on Windows.
*/

#include "pch.h"

#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include "parallel_hashmap/phmap.h"
#include "picosha2/picosha2.h"
using namespace std;
typedef std::pair <int32_t, int32_t> sPair; //type used to record (parent,fragment) pairs
typedef std::pair <int32_t, int32_t> kPair; //type used to record (parent,fragment) pairs

#include "load_kernel.hpp"
#include "load_spectra.hpp"
#include "create_results.hpp"
#include "create_output.hpp"

//  calculates the number of cells to use in the hypergeometric distribution
//  used to model a result, based on the parent mass and fragment mass tolerance
//  using the information in the class member map distribution.

int32_t create_output::get_cells(double _pm,int32_t _res)	{
	int32_t pm = 100*(int32_t)(_pm/100000);
	int32_t r = 2;
	if(_res == 50)	{
		r = 1;
	}
	else if(_res == 20)	{
		r = 0;
	}
	if(pm < 800)	{
		return (int32_t)(pm*distribution[800][r]);
	}
	if(pm > 4000)	{
		return (int32_t)(pm*distribution[4000][r]);
	}
	return (int32_t)(pm*distribution[pm][r]);
}

// generates a probability of random assignment for a particular peptide-to-spectrum match
// using a version of the hypergeometric distribution customized for the experiment situation

bool create_output::apply_model(int32_t _res, // fragent ion mass resolution in mDa
				string& _seq, // best match peptide sequence
				double _pm, // peptide mass
				int32_t _ions, // # of ions identified
				int32_t _lspectrum, // # of ions in spectrum being identified
				int32_t  _z, // parent ion charge
				double& pscore, // score to be returned
				double& p)	{ // probability to be returned
	p = 0.0001; //initialize output value
	pscore = 400; //initialize output value
	int32_t sfactor = 20;
	int32_t sadjust = 3;
	// make adjustments based on mass resolution
	if(_res > 100)	{
		sfactor = 40;
	}
	if(fragmentation == "hcd")	{
		sfactor = 20;
	}
	else if(fragmentation == "cid")	{
		sfactor = 40;
	}
	// calculate parameters for hypergeometric model
	int32_t cells = get_cells(_pm,_res);
	int32_t total_ions = 2*(int32_t)(_seq.size() - 1);
	if(_z > 2)	{
		total_ions *= 2;
	}
	if(total_ions > sfactor)	{
		total_ions = sfactor;
	}
	if(total_ions < _ions)	{
		total_ions = _ions + 1;
	}
	int32_t sc = _lspectrum * sadjust;
	if(_ions >= sc)	{
		sc = _ions + 2;
	}
	// initialize the hypergeometric model calculation object
	hypergeom hp(sc,total_ions,cells);
	// retrieve the probability, based on the hypergeometric model
	p = hp.pdf(_ions);
	// convert probability into a score
	pscore = -100.0*log(p)/2.3;
	return true;
}

// load_mods uses the information stored in "report_mods.txt" to create
// user-friendly strings in place of modification masses

bool create_output::load_mods(void)	{
	ifstream istr;
	// note: report_mods.txt should be in the cwd that idx is run from
	// if this is a problem for your installation, you may want to
	// either add it to the parameter strings object or
	// change the path used here
	istr.open("report_mods.txt");
	if(istr.fail())	{
		return false;
	}
	mt.clear();
	string line;
	while(getline(istr,line))	{
		size_t tab = line.find('\t');
		mt[atoi(line.substr(0,tab).c_str())] = line.substr(tab+1,line.size()-1);
	}
	istr.close();
	return true;
}

//  finds a window (in parent mass ppm) to exclude identifications that
//  are outside of the main distribution of results. This window is
//  found by creating a histogram of parent mass differences from the assigned
//  peptide sequence. The largest bin in the histogram is then located and
//  the lower and upper ppm values at which the histogram has dropped below 1% of the largest bin
//  is determined.

bool create_output::find_window(void)	{
	map<int32_t,int32_t> vs;
	int32_t i = 0;
	int32_t ilow = -1*parent_tolerance;
	int32_t ihigh = parent_tolerance;
	// initialize a histogram between ilow ppm to ihigh ppm
	for(i = ilow-1; i <= ihigh+1; i++)	{
		vs[i] = 0;
	}
	auto it = ppms.begin();
	// iterate through ids to populate the histogram
	while(it != ppms.end())	{
		i = roundf(*it);
		if(i >= ilow-1 and i <= ihigh+1)	{
			vs[i] += 1;
		}
		it++;
	}
	ppm_map.clear();
	int32_t max = 0;
	int32_t center = ilow-1;
	// put the histogram into map form
	for(i = ilow; i <= high;i++)	{
		if(vs[i] > max)	{
			max = vs[i];
			center = i;
		}
		ppm_map[i] = vs[i];
	}
	// bail out if there aren't enough ids
	// and use a +/- 10 ppm window
	if(max < 100)	{
		low = center - 10;
		high = center + 10;
		if(low < ilow-1)	low = ilow-1;
		if(high > ihigh+1)	high = high+1;
		return true;
	}
	// bail out if there aren't enough ids
	// and set no window
	else if(max < ihigh)	{
		low = ilow-1;
		high = ihigh+1;
		return true;
	}
	double ic = (double)max;
	int32_t l = ilow;
	// find the lower edge of the histogram (1%)
	for(i = center;i >= ilow ; i--)	{
		if(l == ilow and vs[i]/ic < 0.01)	{
			l = i;
		}
		else if(l > ilow and vs[i]/ic > 0.05)	{
			l = ilow;
		}
	}
	int32_t h = ihigh;
	// find the upper edge of the histogram (1%)
	for(i = center; i <= ihigh ; i++)	{
		if(h == ihigh and vs[i]/ic < 0.01)	{
			h = i;
		}
		else if(l < ihigh and vs[i]/ic > 0.05)	{
			h = ihigh;
		}
	}
	low = l;
	high = h;
	return true;
}

// creates a single line of TSV formatted output
// change this method if you want to alter the output format

bool create_output::create_line(id& _s, // an id object
				double _pm, // parent mass
				double _d, // assigned mass error (mDa)
				double _ppm, // assigned mass error (ppm)
				double _score, // assigned id score
				rapidjson::Document& _js, // kernel object
				int32_t _u, // unique id of kernel
				string& _line)	{ // output line
	ostringstream oline;
	oline << _s.sc << "\t"; // spectrum scan number
	//write TSV string
	if(_s.rt != 0)	{ // spectrum retention time
		oline << fixed << setprecision(3)  << _s.rt << "\t";
	}
	else	{
		oline << "\t";
	}
	oline << fixed << setprecision(3);
	oline << _pm/1000.0 << "\t"; // spectrum parent mass
	oline << _d << "\t"; // kernel-spectrum parent mass
	oline << setprecision(1);
	oline << _ppm << "\t"; // kernel-spectrum parent mass in ppm
	oline << setprecision(3);
	oline << _s.pz << "\t"; // spectrum parent charge
	oline  << _js["lb"].GetString() << "\t"; // kernel label
	oline << _js["beg"].GetInt() << "\t"; // kernel start residue#
	oline << _js["end"].GetInt() << "\t"; // kernel end residue #
	oline << _js["pre"].GetString() << "\t"; // residue N-terminal to kernel
	oline << _js["seq"].GetString() << "\t"; // kernel sequence
	oline << _js["post"].GetString() << "\t"; // residue C-terminal to kernel
	oline << _s.peaks << "\t"; // # of identified spectrum signals
	const Value& jbs = _js["ns"];
	int32_t lns = 0;
	for(SizeType a = 0; a < jbs.Size();a++)	{ // sum potentials z for the kernel
		lns += jbs[a].GetInt();
	}
	oline << fixed << setprecision(2); 
	if(lns > 0)	{ 
		oline << (float)round(100.0*_s.ri)/100.0 << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << (float)round(10.0*log(lns)/2.3)/10.0 << "\t"; // # of times kernel sequence observed (GPMDB)
		oline << (float)round(-0.01*_score*10.0)/10.0 << "\t"; // score
	}
	else	{
		oline << (float)round(100.0*_s.ri)/100.0 << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << "-" << "\t";// kernel sequence not observed (GPMDB)
		oline << (float)round(10.0*-0.01*_score)/10.0 << "\t"; // score
	}
	//deal with recoding modifications
	if(_js.HasMember("mods"))	{
		const Value& jmods = _js["mods"];
		vector<mod> mods;
		mod tmod;
		for(SizeType a = 0; a < jmods.Size();a++)	{
			const Value& lmods = jmods[a];
			tmod.pos = lmods[1].GetInt();
			tmod.res = lmods[0].GetString();
			tmod.mass = lmods[2].GetInt();
			mods.push_back(tmod);
		}
		sort(mods.begin(),mods.end());
		for(size_t a = 0; a < mods.size(); a++)	{
			if(mt.find(mods[a].mass) != mt.end())	{
				oline << mods[a].res << mods[a].pos << "~" 
					<< mt.find(mods[a].mass)->second << ";";
			}
			else	{
				oline << mods[a].res << mods[a].pos << "#" 
					<< mods[a].mass/1000.0 << ";";
			}
		}
	}
	oline << "\t" << _u;
	_line = oline.str();
	_line = oline.str();
	return true;
}

// creates a single line of TSV formatted output
// change this method if you want to alter the output format

bool create_output::create_line_binary(id& _s, // an id object
					double _pm, // parent mass
					double _d, // assigned mass error (mDa)
					double _ppm, // assigned mass error (ppm)
					double _score, // assigned id score
					osObject& _js, // kernel binary JSON object
					int32_t _u, // unique id of kernel
					string& _line)	{ // output line
	ostringstream oline;
	oline << _s.sc << "\t"; // spectrum scan number
	//write TSV string
	if(_s.rt != 0)	{ // spectrum retention time
		oline << fixed << setprecision(3)  << _s.rt << "\t";
	}
	else	{
		oline << "\t";
	}
	oline << fixed << setprecision(3);
	oline << _pm/1000.0 << "\t"; // spectrum parent mass
	oline << _d << "\t"; // kernel-spectrum parent mass
	oline << setprecision(1);
	oline << _ppm << "\t"; // kernel-spectrum parent mass in ppm
	oline << setprecision(3);
	oline << _s.pz << "\t"; // spectrum parent charge
	oline  << _js.lb << "\t"; // kernel label
	oline << _js.beg << "\t"; // kernel start residue#
	oline << _js.end << "\t"; // kernel end residue #
	oline << _js.pre << "\t"; // residue N-terminal to kernel
	oline << _js.seq << "\t"; // kernel sequence
	oline << _js.post << "\t"; // residue C-terminal to kernel
	oline << _s.peaks << "\t"; // # of identified spectrum signals
	int32_t lns = 0;
	for(size_t i = 0; i < _js.ns.size();i++)	{
		lns += _js.ns[i];
	}
	oline << fixed << setprecision(2); 
	if(lns > 0)	{ 
		oline << (float)round(100.0*_s.ri)/100.0 << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << (float)round(10.0*log(lns)/2.3)/10.0 << "\t"; // # of times kernel sequence observed (GPMDB)
		oline << (float)round(-0.01*_score*10.0)/10.0 << "\t"; // score
	}
	else	{
		oline << (float)round(100.0*_s.ri)/100.0 << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << "-" << "\t";// kernel sequence not observed (GPMDB)
		oline << (float)round(-0.01*_score*10.0)/10.0 << "\t"; // score
	}
	//deal with recoding modifications
	if(!_js.mods.empty())	{
		sort(_js.mods.begin(),_js.mods.end());
		for(size_t a = 0; a < _js.mods.size(); a++)	{
			if(mt.find(_js.mods[a].mass) != mt.end())	{
				oline << _js.mods[a].res << _js.mods[a].pos << "~" 
					<< mt.find(_js.mods[a].mass)->second << ";";
			}
			else	{
				oline << _js.mods[a].res << _js.mods[a].pos << "#" 
					<< _js.mods[a].mass/1000.0 << ";";
			}
		}
	}
	oline << "\t" << _u;
	_line = oline.str();
	return true;
}

// create the header line for the output TSV
// be sure to align this with the results of create_line

bool create_output::create_header_line(string& _h)	{
	_h = "Id\tSub\tScan\tRT(s)\tPeptide mass\tDelta\tppm\tz\tProtein acc\t";
	_h += "Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tlog(f)\tlog(p)\tModifications\tKernel";
	return true;
}

// gets the next kernel from a binary formatted JSON file with an open ifstream

bool create_output::get_next(ifstream& _ifs,osObject& _js)
{
	_js.reset();
	int32_t jsl = 0;
	_ifs.read((char *)&jsl,4); // get the number of entries corresponding to a binary JSON object
	if(_ifs.fail())	{
		cout << "failed to get json size" << '\n';
		return false;
	}
	int32_t count = 0;
	int32_t klen = 0;
	char element = '\0';
	int32_t tlen = 0;
	int32_t itemp = 0;
	size_t i = 0;
	int32_t *pI = 0;
	while(count < jsl and !(_ifs.bad() or _ifs.eof()))	{
		_ifs.read((char *)&klen,4);
		_ifs.read((char *)_js.pKey,klen);
		_js.pKey[klen] = '\0';
		_js.key = _js.pKey;
		_ifs.read(&element,1);
		switch(element)	{
			case 'm':	// special for modification-like entries
				_ifs.read((char *)&tlen,4);
				if(_js.key == "mods")	{	//deal with modification specifications
					_js.mods.clear();
					char c = '\0';
					for(i = 0; i < (size_t)tlen;i++)	{
						_ifs.read((char *)&c,1);
						_js.mod_temp.res = c;
						_ifs.read((char *)&itemp,4);
						_js.mod_temp.pos = itemp;
						_ifs.read((char *)&itemp,4);
						_js.mod_temp.mass = itemp;
						_js.mods.push_back(_js.mod_temp);
					}
				}
				if(_js.key == "savs")	{	//deal with amino acid variant specifications
					_js.savs.clear();
					char c = '\0';
					for(i = 0; i < (size_t)tlen;i++)	{
						_ifs.read((char *)&c,1);
						_js.mod_temp.res = c;
						_ifs.read((char *)&itemp,4);
						_js.mod_temp.pos = itemp;
						_ifs.read((char *)&itemp,4);
						_js.mod_temp.mass = itemp;
						_js.savs.push_back(_js.mod_temp);
					}
				}
				break;
			case 'l':	// arrays of integers
				_ifs.read((char *)&tlen,4);
				_ifs.read((char *)_js.pBuffer,4*tlen);
				pI = (int32_t *)_js.pBuffer;
				if(_js.key == "ns")	{
					_js.ns.insert(_js.ns.end(),pI,pI+tlen);
				}
				break;
			case 's':	// strings
				_ifs.read((char *)&tlen,4);
				_ifs.read((char *)_js.pBuffer,tlen);
				_js.pBuffer[tlen] = '\0';
				if(_js.key == "seq")	{
					_js.seq = (char *)_js.pBuffer;
				}
				else if(_js.key == "pre")	{
					_js.pre = (char *)_js.pBuffer;
				}
				else if(_js.key == "post")	{
					_js.post = (char *)_js.pBuffer;
				}
				else if(_js.key == "lb")	{
					_js.lb = (char *)_js.pBuffer;
				}
				break;
			case 'i':	//	integers
				_ifs.read((char *)&itemp,4);
				if(_js.key == "pm")	{
					_js.pm = itemp;
				}
				else if(_js.key == "h")	{
					_js.h = itemp;
				}
				else if(_js.key == "u")	{
					_js.u = itemp;
				}
				else if(_js.key == "pz")	{
					_js.pz = itemp;
				}
				else if(_js.key == "beg")	{
					_js.beg = itemp;
				}
				else if(_js.key == "end")	{
					_js.end = itemp;
				}
				break;
			default: // should never happen
				cout << "ERROR: bad element value, binary JSON file corrupt" << '\n';
				exit(1);
		}
		count++;
		if(_js.key == "value")	{ // retrieve the file validation string
			validation = (char *)_js.pBuffer;
			return false;
		}
	}
	return true;
}

//creates an output file, as specified in _params for a binary JSON kernel file

bool create_output::create_binary(map<string,string>& _params,
				const create_results& _cr,
				map<int32_t, set<int32_t> >& _hu)	{
	spectrum_count = atoi(_params.find("spectra")->second.c_str());
	parent_tolerance = atoi(_params.find("parent tolerance")->second.c_str());
	low = parent_tolerance - 1;
	high = parent_tolerance + 1;
	ifstream ifs(_params.find("kernel file")->second,ios::in | ios::binary);
	if(ifs.fail())	{
		return false;
	}
	if(!load_mods())	{ //warns if "reports_mods.txt" is not present
		cout << "Warning (idx:1001): annotation file \"report_mods.txt\" was not present" << '\n';
	}
//	if(_params.find("kernel file")->second.find("_tmt") != _params.find("kernel file")->second.npos)	{
//	}
	int32_t k = 0;
	//loops through ids and creates some structures to facilitate generating report lines
	for(size_t j = 0; j < _cr.ids.size(); j++)	{
		sv[_cr.ids[j].sn] = _cr.ids[j];
		for(size_t i = 0; i < _cr.ids[j].ks.size(); i++)	{
			k = _cr.ids[j].ks[i];
			if(sdict.find(k) != sdict.end())	{
				sdict[k].insert(_cr.ids[j].sn);
			}
			else	{
				sdict.insert(pair<int32_t,set<int32_t> >(k,set<int32_t>()));
				sdict[k].insert(_cr.ids[j].sn);
			}
		}
	}
	auto itsd = sdict.begin();
	set<int32_t> hu_set;
	while(itsd != sdict.end())	{
		if(_hu.find(itsd->first) != _hu.end())	{
			hu_set.insert(_hu[itsd->first].begin(),_hu[itsd->first].end());
		}
		itsd++;
	}
	//initialize variables
	int32_t res = atoi(_params.find("fragment tolerance")->second.c_str());
	int32_t inferred = 0;
	double score_min = 200.0;
	double total_prob = 0.0; //sum of all assigned probabilities
	int32_t min_c = 8; //minimum number of assignments necessary for a spectrum-to-kernel match
	//updated min_c value based on instrument resolution
	if(res == 50)	{
		min_c = 7;
	}
	else if(res == 20)	{
		min_c = 6;
	}
	int32_t min_l = min_c;
	string line; //will contain a JSON Lines JSON object
	int32_t c = 0; //lines read counter
	int32_t h = 0; //for checking peptide homology
	double ri_limit = 0.2;
	double dparent = (double)parent_tolerance;
	//loop through the JSON Lines entries in the kernel file to find
	//the information about individual ids
	osObject js;
	while(get_next(ifs,js))	{
		//output keep-alive text for logging
		if(c != 0 and c % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(c != 0 and c % 200000 == 0)	{
			cout << " " << c << '\n';
			cout.flush();
		}
		c++; //increment line count
		if(hu_set.find(c-1) == hu_set.end())	{ //bail out if the line was not in solution
			continue;
		}
		if(js.pm == 0)	{ //bail out if the JSON object does not have a parent mass
			continue;
		}
		h = (int32_t)js.h;
		if(sdict.find(h) == sdict.end())	{ //bail out if the kernel was not identified
			continue;
		}
		if(h != (int32_t)js.u)	{ //track # of peptide alternate solutions
			inferred += 1;
		}
		//initialize variables
		double max_prob = 0.0; //maximum probability for a set of spectrum-to-kernel assignments
		double prob = 0.0; //probability of a spectrum-to-kernel assignment
		double delta = 0.0; //parent mass difference spectrum parent - kernel parent
		double pm = 0.0; //kernel parent mass
		double ppm = -10000000.0; //delta in ppm
		double score = 0.0; //score for a spectrum-to-kernel assignment
		auto it = sdict[h].begin(); //iterator for sdict vectors
		int32_t s = 0; //spectrum index
		string seq = js.seq; //kernel peptide sequence
		min_l = min_c;
		int32_t slen = (int32_t)(0.5 + (float)seq.length()*0.70);
		if(min_c > slen)	{
			min_l = slen;
		}
		//iterate through the spectra corresponding the index value h
		while(it != sdict[h].end())	{
			s = *it;
			it++;
			pm = (double)js.pm;
			//calculated mass difference, model - observed
			delta = (sv[s].pm-pm)/1000.0;
			// convert to parts-per-million
			ppm = 1.0e6*(sv[s].pm-pm)/pm;
			// deal with A1 peaks
			if(delta > 0.5 and delta < 1.5)	{
				ppm = 1.0e6*(sv[s].pm-c13-pm)/pm;
			}
			else if(delta > 1.5)	{
				ppm = 1.0e6*(sv[s].pm-2*c13-pm)/pm;
			}
			if(fabs(ppm) > dparent)	{
				continue;
			}
			seq = js.seq;
			// apply math model to result to obtain a probability of stochastic assignment
			// get the total modification mass for the peptide and remove it from the peptide mass
			// in the apply_model method to get the right density from the model
			double mod_mass = 0.0;
			for(size_t a = 0; a < js.mods.size();a++)	{
				mod_mass += (double)js.mods[a].mass;
			}
			apply_model(res,seq,pm-mod_mass,sv[s].peaks,sv[s].ions,sv[s].pz,score,prob);
			ri_limit = get_ri_limit(seq.size());
			if(score < score_min or sv[s].ri < ri_limit or sv[s].peaks < min_l)	{
				continue; //bail out if match does not pass conditions
			}
			//keep track of the probabilities
			if(prob > max_prob)	{
				max_prob = prob;
			}
			// construct a line of output for this spectrum object
			string new_line;
			create_line_binary(sv[s],pm,delta,ppm,score,js,c-1,new_line);
			//add output file line to mapped vector for spectrum index s
			if(odict.find(s) != odict.end())	{
				odict[s].push_back(new_line);			
			}
			else	{
				odict.insert(pair<int32_t,vector<string> >(s,vector<string>()));
				odict[s].push_back(new_line);
				//add ppm value to ppms histogram
				if(delta < 0.5)	{
					ppms.push_back(roundf(ppm));
				}
			}
		}
		total_prob += max_prob;
	}
	ifs.close();
	//open a file stream to output information in odict
	dump_lines(_params.find("output file")->second,total_prob);
	cout << '\n';
	cout.flush();
	return true;
}

// Serializes the formulated lines to a specified output file

bool create_output::dump_lines(string& _ofile,double _tp)	{
	//determine parent mass delta ppm acceptance window
	string temp;
	find_window();
	int32_t err = 0;
	int32_t sub = 0;
	int32_t tot = 0;
	string t;
	//loop through result lines and record the information
	int32_t low_t = low;
	int32_t high_t = high;
	int32_t ps_t = 0;
	int32_t scan = 0;
	multimap<int32_t,string> ostrings;
	multimap<int32_t,string> estrings;
	pair<int32_t,string> opair;
	char *pString =  new char[1024*8];
	set<int32_t> scans;
	for(int32_t a = 0; a < (int32_t)odict.size(); a++)	{
		sub = 1;
		for(size_t b = 0; b < odict[a].size(); b++)	{
			t = odict[a][b];
			ps_t = roundf(get_ppm(t));
			sprintf(pString,"%li\t%li\t%s",(long)a,(long)sub,t.c_str());
			temp = pString;
			scan = get_scan(t);
			scans.insert(scan);
			opair.first = scan;
			opair.second = temp;
			if(ps_t < high_t and ps_t > low_t)	{ //apply the parent mass window
				ostrings.insert(opair);
				sub++;
				tot++; //record number of parent mass window hits
			}
			else	{
				estrings.insert(opair);
				err++; // record number of parent mass window misses
			}
		}
	}
	delete pString;
	//define headers for the output TSV file
	temp = "";
	create_header_line(temp);

	ofstream ofs;
	// write output file without line ending characters
	// so that the hash value can be calculated in 
	// a platform independent manner
	ofs.open(_ofile); //open output stream
	ofs << temp;
	auto itS = ostrings.begin();
	while(itS != ostrings.end())	{
		ofs << itS->second;
		itS++;
	}
	ofs.close();
	if(errors_tsv)	{
		string efile = _ofile;
		efile.insert(_ofile.length()-4,".err");
		ofs.open(efile); //open output stream
		create_header_line(temp);
		ofs << temp << '\n';
		itS = estrings.begin();
		while(itS != estrings.end())	{
			ofs << itS->second << '\n';
			itS++;
		}
		ofs.close();
	}
	// create validation hash value
	std::ifstream ifs(_ofile, std::ios::binary);
	std::vector<unsigned char> s(picosha2::k_digest_size);
	picosha2::hash256(ifs, s.begin(), s.end());
	ifs.close();
	// re-write output file with platform-specific line ending characters
	ofs.open(_ofile);
	ofs << temp << '\n';
	itS = ostrings.begin();
	while(itS != ostrings.end())	{
		ofs << itS->second << '\n';
		itS++;
	}
	ofs.close();

	std::string hash_hex_str;
	picosha2::hash256_hex_string(s, hash_hex_str);
	info["output file validation"] = hash_hex_str;
	//output some additional information for logging
	char *str = new char[1024];
	sprintf(str,"%i",tot);
	info["lines"] = str;
	double dtot = 0.0;
	for(int32_t i = low_t+1; i < high_t; i++)	{
		dtot += (double)ppm_map[i];
	}
	double evalue = 0.0;
	if(dtot > 0.0 and low_t > 5-parent_tolerance)	{
		evalue = 0.0;
		for(int32_t i = -1*parent_tolerance; i <= 5-parent_tolerance; i++)	{
			evalue += (double)ppm_map[i];
		}
		evalue /= 6;
		double ble = 100.0*(evalue*(high_t - low_t)/dtot);
		sprintf(str,"%.1f",ble);
		info["baseline error (%)"] = str;
		sprintf(str,"%.1f",100.0*get_ble_error(evalue*(high_t - low_t),dtot));
		info["baseline error (+/- %)"] = str;
		sprintf(str,"%.1f",evalue);
		info["baseline error (per ppm)"] = str;
		if(spectrum_count - scans.size() > 0)	{
			sprintf(str,"%.1e",evalue*6.0/(spectrum_count - scans.size()));
			info["fpr"] = str;
		}
		else	{
			info["fpr"] = "";
		}
	}
	else if(dtot > 0.0 and high_t < parent_tolerance-5)	{
		evalue = 0.0;
		for(int32_t i = parent_tolerance-5; i <= parent_tolerance; i++)	{
			evalue += (double)ppm_map[i];
		}
		evalue /= 6.0;
		double ble = 100.0*(evalue*(double)(high_t - low_t)/dtot);
		sprintf(str,"%.1f",ble);
		info["baseline error (%)"] = str;
		sprintf(str,"%.1f",100.0*get_ble_error(evalue*(high_t - low_t),dtot));
		info["baseline error (+/- %)"] = str;
		sprintf(str,"%.1f",evalue);
		info["baseline error (ppm)"] = str;
		if(spectrum_count - scans.size() > 0)	{
			sprintf(str,"%.1e",evalue*6.0/(spectrum_count - scans.size()));
			info["fpr"] = str;
		}
		else	{
			info["fpr"] = "";
		}
	}
	else	{
		info["baseline error (%)"] = "";
		info["baseline error (per ppm)"] = "";
		info["fpr"] = "";
	}
	sprintf(str,"%i,%i",low_t,high_t);
	info["parent ion tolerance"] = str;
	return true;
}

// creates a file containing information about the PSM identification process
// as well as file validation records

bool create_output::dump_meta(const map<string,string>& _p)	{
	auto itp = _p.find("output file");
	string mpath = itp->second;
	mpath += ".meta";
	ofstream mfs;
	mfs.open(mpath);
	mfs << "{\"input\" : {" << '\n';
	string first;
	string second;
	itp = _p.begin();
	while(itp != _p.end())	{
		first = itp->first;
		second = itp->second;
		while(first.find("\\") != first.npos)	{
			first.replace(first.find("\\"),1,"/");
		}
		while(second.find("\\") != second.npos)	{
			second.replace(second.find("\\"),1,"/");
		}
		mfs << "\"" << first << "\" : \"" << second << "\"";
		itp++;
		if(itp != _p.end())	{
			mfs << " , " << '\n';
		}
		else	{
			mfs << '\n';
		}
	}
	mfs << "},\n\"data\" : {" << '\n'; 
	itp = info.begin();
	while(itp != info.end())	{
		first = itp->first;
		second = itp->second;
		while(first.find("\\") != first.npos)	{
			first.replace(first.find("\\"),1,"/");
		}
		while(second.find("\\") != second.npos)	{
			second.replace(second.find("\\"),1,"/");
		}
		mfs << "\"" << first << "\" : \"" << second << "\"";
		cout << "  " << first << " : " << second << '\n';
		itp++;
		mfs << " , " << '\n';
	}
	mfs << "\"ppms\" : {";
	auto itppm = ppm_map.begin();
	while(itppm != ppm_map.end())	{
		mfs << "\"" << itppm->first << "\" : " << itppm->second;
		itppm++;
		if(itppm != ppm_map.end())	{
			mfs << ",";
		}
	}
	mfs << "} " << '\n';
	mfs << "}}" << '\n';
	mfs.close(); 
	return true;
}

// creates an output file, as specified in _params for a JSON kernel

bool create_output::create(map<string,string>& _params, 
				const create_results& _cr, 
				map<int32_t, set<int32_t> >& _hu)	{
	spectrum_count = atoi(_params["spectra"].c_str());
	parent_tolerance = atoi(_params["parent tolerance"].c_str());
	low = parent_tolerance - 1;
	high = parent_tolerance + 1;
	ifstream ifs;
	ifs.open(_params.find("kernel file")->second,std::ifstream::in);
	if(!ifs.good())	{
		cout << "Kernel file \"" << _params.find("kernel file")->second << "\" would not open" << '\n';
		return false;
	}
	if(!load_mods())	{ //warns if "reports_mods.txt" is not present
		cout << "Warning (idx:1001): annotation file \"report_mods.txt\" was not present" << '\n';
	}
	int32_t k = 0;
//	if(_params.find("kernel file")->second.find("_tmt") != _params.find("kernel file")->second.npos)	{
//	}
	//loops through ids and creates some structures to facilitate generating report lines
	for(size_t j = 0; j < _cr.ids.size(); j++)	{
		sv[_cr.ids[j].sn] = _cr.ids[j];
		for(size_t i = 0; i < _cr.ids[j].ks.size(); i++)	{
			k = _cr.ids[j].ks[i];
			if(sdict.find(k) != sdict.end())	{
				sdict[k].insert(_cr.ids[j].sn);
			}
			else	{
				sdict.insert(pair<int32_t,set<int32_t> >(k,set<int32_t>()));
				sdict[k].insert(_cr.ids[j].sn);
			}
		}
	}
	auto itsd = sdict.begin();
	set<int32_t> hu_set;
	while(itsd != sdict.end())	{
		if(_hu.find(itsd->first) != _hu.end())	{
			hu_set.insert(_hu[itsd->first].begin(),_hu[itsd->first].end());
		}
		itsd++;
	}
	//initialize variables
	int32_t res = atoi(_params.find("fragment tolerance")->second.c_str());
	int32_t inferred = 0;
	double score_min = 200.0;
	double total_prob = 0.0; //sum of all assigned probabilities
	double limit = 0.2;
	int32_t min_c = 8; //minimum number of assignments necessary for a spectrum-to-kernel match
	//updated min_c value based on instrument resolution
	if(res == 50)	{
		min_c = 7;
	}
	else if(res == 20)	{
		min_c = 6;
	}
	int32_t min_l = min_c;
	string line; //will contain a JSON Lines JSON object
	int32_t c = 0; //lines read counter
	int32_t h = 0; //for checking peptide homology
	//loop through the JSON Lines entries in the kernel file to find
	//the information about individual ids
	const int max_buffer = 1024*8-1;
	char *buffer = new char[max_buffer+1];
	char *last = new char[max_buffer+1];
	double dparent = (double)parent_tolerance;
	while(ifs.getline(buffer,max_buffer))	{
		//output keep-alive text for logging
		if(c != 0 and c % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(c != 0 and c % 200000 == 0)	{
			cout << " " << c << '\n';
			cout.flush();
		}
		memcpy(last,buffer,max_buffer);
		c++; //increment line count
		if(hu_set.find(c-1) == hu_set.end())	{ //bail out if the line was not in solution
			continue;
		}
		Document js; //rapidjson main object
   		js.ParseInsitu(buffer); //add information for a single JSON Lines object
		if(!js.HasMember("pm"))	{ //bail out if the JSON object does not have a parent mass
			if(js.HasMember("validation"))	{
				validation = js["value"].GetString();
			}
			continue;
		}
		h = (int32_t)js["h"].GetInt();
		if(sdict.find(h) == sdict.end())	{ //bail out if the kernel was not identified
			continue;
		}
		if(h != (int32_t)js["u"].GetInt())	{ //track # of peptide alternate solutions
			inferred += 1;
		}
		//initialize variables
		double max_prob = 0.0; //maximum probability for a set of spectrum-to-kernel assignments
		double prob = 0.0; //probability of a spectrum-to-kernel assignment
		double delta = 0.0; //parent mass difference spectrum parent - kernel parent
		double pm = 0.0; //kernel parent mass
		double ppm = -10000000.0; //delta in ppm
		double score = 0.0; //score for a spectrum-to-kernel assignment
		auto it = sdict[h].begin(); //iterator for sdict vectors
		int32_t s = 0; //spectrum index
		string seq= js["seq"].GetString(); //kernel peptide sequence
		//iterate through the spectra corresponding the index value h
		min_l = min_c;
		int32_t slen = (int32_t)(0.5 + (float)seq.length()*0.70);
		if(min_c > slen)	{
			min_l = slen;
		}
		while(it != sdict[h].end())	{
			s = *it;
			it++;
			pm = js["pm"].GetFloat();
			//calculated mass difference, model - observed
			delta = (sv[s].pm-pm)/1000.0;
			// convert to parts-per-million
			ppm = 1.0e6*(sv[s].pm-pm)/pm;
			// deal with A1 peaks
			if(delta > 0.5 and delta < 1.5)	{
				ppm = 1.0e6*(sv[s].pm-c13-pm)/pm;
			}
			else if(delta > 1.5)	{
				ppm = 1.0e6*(sv[s].pm-2*c13-pm)/pm;
			}
			if(fabs(ppm) > dparent)	{
				continue;
			}
			seq = js["seq"].GetString();
			// get the total modification mass for the peptide and remove it from the peptide mass
			// in the apply_model method to get the right density from the model
			double mod_mass = 0.0;
			if(js.HasMember("mods"))	{
				const Value& jmods = js["mods"];
				for(SizeType a = 0; a < jmods.Size();a++)	{
					const Value& lmods = jmods[a];
					mod_mass += (double)(lmods[2].GetInt());
				}
			}	
			// apply math model to result to obtain a probability of stochastic assignment
			apply_model(res,seq,pm-mod_mass,sv[s].peaks,sv[s].ions,sv[s].pz,score,prob);
			limit = get_ri_limit(seq.size());
			if(score < score_min or sv[s].ri < limit or sv[s].peaks < min_l)	{
				continue; //bail out if match does not pass conditions
			}
			//keep track of the probabilities
			if(prob > max_prob)	{
				max_prob = prob;
			}
			// construct a line of output for this spectrum object
			string new_line;
			create_line(sv[s],pm,delta,ppm,score,js,c-1,new_line);
			//add output file line to mapped vector for spectrum index s
			if(odict.find(s) != odict.end())	{
				odict[s].push_back(new_line);			
			}
			else	{
				odict.insert(pair<int32_t,vector<string> >(s,vector<string>()));
				odict[s].push_back(new_line);
				//add ppm value to ppms histogram
				if(delta < 0.5)	{
					ppms.push_back(roundf(ppm));
				}
			}
		}
		total_prob += max_prob;
	}
	if(strstr(last,"\"validation\"") != NULL)	{
		Document js; //rapidjson main object
		js.ParseInsitu(last); //add information for a single JSON Lines object
		if(js.HasMember("validation"))	{
			validation = js["value"].GetString();
		}
	}
	ifs.close();
	delete buffer;
	delete last;
	//open a file stream to output information in odict
	dump_lines(_params.find("output file")->second,total_prob);
	cout << '\n';
	cout.flush();
	return true;
}

