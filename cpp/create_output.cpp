/*
#
# Copyright © 2019 Ronald C. Beavis
#
# Loads a spectrum file into a vector of spectrum objects
#
*/
/*
create_output is an object to perform all of the tasks necessary to generate
a useful output file specified on the command line
*/
#include "pch.h"

#include "rapidjson/document.h"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <unordered_map>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "parallel_hashmap/phmap.h"
using namespace std;
using namespace rapidjson; 
//simplify calling rapidjson methods
#include "load_spectra.hpp"
#include "load_kernel.hpp"
#include "create_results.hpp"
#include "create_output.hpp"

//
// initialize values
//
create_output::create_output(void)	{
	low = -20.0;
	high = 20.0;
	load_distribution();
}

create_output::~create_output(void)	{
}
//
//calculates the number of cells to use in the hypergeometric distribution
//used to model a result, based on the parent mass and fragment mass tolerance
//using the information in the class member map distribution.
//
int64_t create_output::get_cells(double _pm,int64_t _res)	{
	int64_t pm = 100*(int64_t)(_pm/100000);
	int64_t r = 2;
	if(_res == 50)	{
		r = 1;
	}
	else if(_res == 20)	{
		r = 0;
	}
	if(pm < 800)	{
		return (int64_t)(pm*distribution[800][r]);
	}
	if(pm > 4000)	{
		return (int64_t)(pm*distribution[4000][r]);
	}
	return (int64_t)(pm*distribution[pm][r]);
}
//
//generates a probability of random assignment for a particular peptide-to-spectrum match
//using a version of the hypergeometric distribution customized for the experiment situation
//
bool create_output::apply_model(int64_t _res,string& _seq,double _pm,int64_t _ions,int64_t _lspectrum,double& pscore,double& p){
	p = 0.0001; //initialize output value
	pscore = 400; //initialize output value
	int64_t sfactor = 20;
	int64_t sadjust = 3;
	if(_res > 100)	{
		sfactor = 40;
	}
	if(fragmentation == "hcd")	{
		sfactor = 20;
	}
	else if(fragmentation == "cid")	{
		sfactor = 40;
	}
	int64_t cells = get_cells(_pm,_res);
	int64_t total_ions = 2*(int64_t)(_seq.size() - 1);
	if(total_ions > sfactor)	{
		total_ions = sfactor;
	}
	if(total_ions < _ions)	{
		total_ions = _ions + 1;
	}
	int64_t sc = _lspectrum * sadjust;
	if(_ions >= sc)	{
		sc = _ions + 2;
	}
	hypergeom hp(sc,total_ions,cells);
	p = hp.pdf(_ions);
	pscore = -100.0*log(p)/2.3;
	return true;
}
//
//load_mods uses the information stored in "report_mods.txt" to create
//user-friendly strings in place of modification masses
//
bool create_output::load_mods(void)	{
	ifstream istr;
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
//
//finds a window (in parent mass ppm) to exclude identifications that
//are outside of the main distribution of results. This window is
//found by creating a histogram of parent mass differences from the assigned
//peptide sequence. The largest bin in the histogram is then located and
//the lower and upper ppm values at which the histogram has dropped below 1% of the largest bin
//is determined.
//
bool create_output::find_window(void)	{
	map<int64_t,int64_t> vs;
	int64_t i = 0;
	for(i = -21; i < 22; i++)	{
		vs[i] = 0;
	}
	auto it = ppms.begin();
	int64_t max = -22;
	int64_t center = -1;
	while(it != ppms.end())	{
		i = (int64_t)(0.5 + *it);
		if(i >= -21 and i <= 21)	{
			vs[i] += 1;
		}
		if(vs[i] > max)	{
			max = vs[i];
			center = i;
		}
		it++;
	}
//	cout << endl << max << endl;
//	for(int64_t j = -20; j <= 20; j++)	{
//		cout << j << "\t" << vs[j] << endl;
//	}
	double ic = (double)max;
	if(ic < 100.0)	{
		low = -20.0;
		high = 20.0;
		return true;
	}
	int64_t l = -20;
	for(i = center;i >= -20 ; i--)	{
		if(vs[i]/ic < 0.01 and vs[i-1]/ic < 0.01)	{
			l = i;
			break;
		}
	}
	int64_t h = 20;
	for(i = center; i <= 20 ; i++)	{
		if(vs[i]/ic < 0.01 and vs[i+1]/ic < 0.01)	{
			h = i;
			break;
		}
	}
	low = (double)l;
	high = (double)h;
	return true;
}
//
// creates a single line of TSV formatted output
// change this method if you want to alter the output format
//
bool create_output::create_line(id& _s,double _pm,double _d,double _ppm,double _score,rapidjson::Document& _js,string& _line)	{
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
	int64_t lns = 0;
	for(SizeType a = 0; a < jbs.Size();a++)	{ // sum potentials z for the kernel
		lns += jbs[a].GetInt();
	}
	oline << fixed << setprecision(2); 
	if(lns > 0)	{ 
		oline << _s.ri << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << log(lns)/2.3 << "\t"; // # of times kernel sequence observed (GPMDB)
		oline << -0.01*_score << "\t"; // score
	}
	else	{
		oline << _s.ri << "\t"; // fraction of spectrum intensity identified
		oline << setprecision(1);
		oline << "-" << "\t";// kernel sequence not observed (GPMDB)
		oline << -0.01*_score << "\t"; // score
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
	_line = oline.str();
	return true;
}
//
// create the header line for the output TSV
// be sure to align this with the results of create_line
//
bool create_output::create_header_line(string& _h)	{
	_h = "Id\tSub\tScan\tRT(s)\tPeptide mass\tDelta\tppm\tz\tProtein acc\t";
	_h += "Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tlog(f)\tlog(p)\tModifications";
	return true;
}
//
//create coordinates the creation of an output file, as specified in _params
//
bool create_output::create(map<string,string>& _params,create_results& _cr, map<int64_t, set<int64_t> >& _hu)	{
	FILE *pFile = ::fopen(_params["kernel file"].c_str(),"r");
	if(pFile == NULL)	{
		return false;
	}
	if(!load_mods())	{ //warns if "reports_mods.txt" is not present
		cout << "Warning (idx:1001): annotation file \"report_mods.txt\" was not present" << endl;
	}
	int64_t k = 0;
	//loops through ids and creates some structures to facilitate generating report lines
	for(size_t j = 0; j < _cr.ids.size(); j++)	{
		sv[_cr.ids[j].sn] = _cr.ids[j];
		for(size_t i = 0; i < _cr.ids[j].ks.size(); i++)	{
			k = _cr.ids[j].ks[i];
			if(sdict.find(k) != sdict.end())	{
				sdict[k].insert(_cr.ids[j].sn);
			}
			else	{
				sdict.insert(pair<int64_t,set<int64_t> >(k,set<int64_t>()));
				sdict[k].insert(_cr.ids[j].sn);
			}
		}
	}
	auto itsd = sdict.begin();
	set<int64_t> hu_set;
	while(itsd != sdict.end())	{
		if(_hu.find(itsd->first) != _hu.end())	{
			hu_set.insert(_hu[itsd->first].begin(),_hu[itsd->first].end());
		}
		itsd++;
	}
	//initialize variables
	int64_t res = atoi(_params["fragment tolerance"].c_str());
	int64_t inferred = 0;
	int64_t specs = atoi(_params["spectra"].c_str());
	double score_min = 200.0;
	if(specs > 0)	{
		score_min += 100.0 * log(specs)/2.3;
	}
	double total_prob = 0.0; //sum of all assigned probabilities
	int64_t min_c = 8; //minimum number of assignments necessary for a spectrum-to-kernel match
	//updated min_c value based on instrument resolution
	if(res == 50)	{
		min_c = 7;
	}
	else if(res == 20)	{
		min_c = 6;
	}
	string line; //will contain a JSON Lines JSON object
	int64_t c = 0; //lines read counter
	int64_t h = 0; //for checking peptide homology
	//loop through the JSON Lines entries in the kernel file to find
	//the information about individual ids
	const int max_buffer = 1024*16-1;
	char *buffer = new char[max_buffer+1];
	while(::fgets(buffer,max_buffer,pFile) != NULL)	{
		//output keep-alive text for logging
		if(c != 0 and c % 10000 == 0)	{
			cout << '.';
			cout.flush();
		}
		if(c != 0 and c % 200000 == 0)	{
			cout << " " << c << endl;
			cout.flush();
		}
		c++; //increment line count
		if(hu_set.find(c-1) == hu_set.end())	{
			continue;
		}
		Document js; //rapidjson main object
   		js.ParseInsitu(buffer); //add information for a single JSON Lines object
		if(!js.HasMember("pm"))	{ //bail out if the JSON object does not have a parent mass
			continue;
		}
		h = (int64_t)js["h"].GetInt();
		if(sdict.find(h) == sdict.end())	{ //bail out if the kernel was not identified
			continue;
		}
		if(h != (int64_t)js["u"].GetInt())	{ //track # of peptide alternate solutions
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
		int64_t s = 0; //spectrum index
		string seq; //kernel peptide sequence
		//iterate through the spectra corresponding the index value h
		while(it != sdict[h].end())	{
			s = *it;
			it++;
			pm = js["pm"].GetFloat();
			//calculated mass difference, model - observed
			delta = (sv[s].pm-pm)/1000.0;
			// convert to parts-per-million
			ppm = 1.0e6*(sv[s].pm-pm)/pm;
			// deal with A1 peaks
			if(delta > 0.5)	{
				ppm = 1.0e6*(sv[s].pm-c13-pm)/pm;
			}
			seq = js["seq"].GetString();
			// apply math model to result to obtain a probability of stochastic assignment
			apply_model(res,seq,pm,sv[s].peaks,sv[s].ions,score,prob);
			if(score < score_min or sv[s].ri < 0.20 or sv[s].peaks < min_c)	{
				continue; //bail out if match does not pass conditions
			}
			//keep track of the probabilities
			if(prob > max_prob)	{
				max_prob = prob;
			}
			// construct a line of output for this spectrum object
			string new_line;
			create_line(sv[s],pm,delta,ppm,score,js,new_line);
			//add output file line to mapped vector for spectrum index s
			if(odict.find(s) != odict.end())	{
				odict[s].push_back(new_line);			
			}
			else	{
				odict.insert(pair<int64_t,vector<string> >(s,vector<string>()));
				odict[s].push_back(new_line);
			}
			//add ppm value to ppms histogram
			ppms.push_back((int64_t)(0.5+ppm));
		}
		total_prob += max_prob;
	}
	::fclose(pFile);
	delete buffer;
	//open a file stream to output information in odict
	ofstream ofs;
	ofs.open(_params["output file"]); //open output stream
	//define headers for the output TSV file
	string header;
	create_header_line(header);
	ofs << header << endl;
	//determine parent mass delta ppm acceptance window
	find_window();
	int64_t err = 0;
	int64_t sub = 0;
	int64_t tot = 0;
	string t;
	//loop through result lines and record the information
	int64_t low_t = (int64_t)(0.5 + low);
	int64_t high_t = (int64_t)(0.5 + high);
	int64_t ps_t = 0;
	for(int64_t a = 0; a < (int64_t)odict.size(); a++)	{
		sub = 1;
		for(size_t b = 0; b < odict[a].size(); b++)	{
			t = odict[a][b];
			ps_t = (int64_t)(0.5+get_ppm(t));
			if(ps_t <= high_t and ps_t >= low_t)	{ //apply the parent mass window
				ofs << a + 1 << "\t";
				ofs << sub << "\t";
				ofs << t << endl;
				sub++;
				tot++; //record number of parent mass window hits
			}
			else	{
				err++; // record number of parent mass window misses
			}
		}
	}
	//output some additional information for logging
	cout << "\n     lines = " << tot << endl;
	if(total_prob > 0)	{
		cout << "     fpr = " << scientific << setprecision(1) << total_prob << endl;
	}
	if(low != -20.0 and high != 20.0)	{
		double ble = (20.0 - high) + (low + 20.0) + 2;
		ble = err/ble;
		ble = 100.0*(ble*(high - low - 1)/tot);
		cout << "     baseline error = " << fixed << setprecision(1) << ble << "% (" << err << ")" << endl;
	}
	else	{
		cout << "     baseline error = n/a" << endl;
	}
	cout << "     parent ion tolerance = " << fixed << setprecision(0) << low+1.0 << "," << high-1.0 << endl;
	

	ofs.close();
	cout << endl;
	cout.flush();
	return true;
}

