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
	for(i = -20; i < 21; i++)	{
		vs[i] = 0;
	}
	auto it = ppms.begin();
	int64_t max = -21;
	int64_t center = -1;
	while(it != ppms.end())	{
		i = (int64_t)(0.5 + it->second);
		if(i >= -20 and i <= 20)	{
			vs[i] += 1;
		}
		if(vs[i] > max)	{
			max = vs[i];
			center = i;
		}
		it++;
	}
	double ic = (double)max;
	if(ic < 100.0)	{
		low = -20.0;
		high = 20.0;
		return true;
	}
	int64_t l = -20;
	for(i = -20;i < center; i++)	{
		if(vs[i]/ic >= 0.01)	{
			l = i;
			break;
		}
	}
	int64_t h = 20;
	for(i = 20; i > center; i--)	{
		if(vs[i]/ic >= 0.01)	{
			h = i;
			break;
		}
	}
	low = (double)l;
	high = (double)h;
	return true;
}
//
//create coordinates the creation of an output file, as specified in _params
//
bool create_output::create(map<string,string>& _params,create_results& _cr)	{
	ifstream istr;
	istr.open(_params["kernel file"]); //opens a stream associated with the kernel file
	if(istr.fail())	{
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
	fragmentation = _params["fragmentation"];
	if(res == 50)	{
		min_c = 7;
	}
	else if(res == 20)	{
		min_c = 6;
	}
	if(fragmentation == "hcd")	{
		min_c = 6;
	}		
	else if(fragmentation == "cid")	{
		min_c = 8;
	}		
	string line; //will contain a JSON Lines JSON object
	using namespace rapidjson; //simplify calling rapidjson methods
	int64_t c = 0; //lines read counter
	int64_t h = 0; //for checking peptide homology
	//loop through the JSON Lines entries in the kernel file to find
	//the information about individual ids
	while(getline(istr,line))	{
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
		Document js; //rapidjson main object
   		js.Parse(line.c_str(),line.length()); //add information for a single JSON Lines object
		if(!js.HasMember("pm"))	{ //bail out if the JSON object does not have a parent mass
			continue;
		}
		h = (int64_t)js["h"].GetInt();
		if(sdict.find(h) == sdict.end())	{ //bail out if the kernel was not identified
			continue;
		}
		if(h != (int64_t)js["u"].GetInt())	{ //keep track of the number of peptide alternate solutions
			inferred += 1;
		}
		//initialize variables
		double max_prob = 0.0; //maximum probability for a set of spectrum-to-kernel assignments
		double prob = 0.0; //probability of a spectrum-to-kernel assignment
		double delta = 0.0; //parent mass difference spectrum parent - kernel parent
		double pm = 0.0; //kernel parent mass
		double ppm = 0.0; //delta in ppm
		double score = 0.0; //score for a spectrum-to-kernel assignment
		auto it = sdict[h].begin(); //iterator for sdict vectors
		int64_t s = 0; //spectrum index
		string seq; //kernel peptide sequence
		//iterate through the spectra corresponding the index value h
		while(it != sdict[h].end())	{
			s = *it;
			it++;
			pm = js["pm"].GetFloat();
			delta = (sv[s].pm-pm)/1000.0;
			ppm = 1.0e6*(sv[s].pm-pm)/pm;
			if(delta > 0.5)	{
				ppm = 1.0e6*(sv[s].pm-c13-pm)/pm;
			}
			seq = js["seq"].GetString();
			apply_model(res,seq,pm,sv[s].peaks,sv[s].ions,score,prob);
			if(score < score_min or sv[s].ri < 0.20 or sv[s].peaks < min_c)	{
				continue; //bail out if match does not pass conditions
			}
			//keep track of the probabilities
			if(prob > max_prob)	{
				max_prob = prob;
			}
			//initialize a string-based stream to hold one line of the output file
			ostringstream oline;
			//write TSV string
			oline << sv[s].sc << "\t" << fixed << setprecision(3) << pm/1000.0 << "\t";
			oline << delta << "\t" << setprecision(1) << ppm << setprecision(3) << "\t";
			oline << sv[s].pz << "\t" << js["lb"].GetString() << "\t";
			oline << js["beg"].GetInt() << "\t" << js["end"].GetInt() << "\t";
			oline << js["pre"].GetString() << "\t" << js["seq"].GetString() << "\t";
			oline << js["post"].GetString() << "\t" << sv[s].peaks << "\t";
			const Value& jbs = js["ns"];
			int64_t lns = 0;
			for(SizeType a = 0; a < jbs.Size();a++)	{
				lns += jbs[a].GetInt();
			}
			oline << fixed << setprecision(2);
			if(lns > 0)	{
				oline << sv[s].ri << "\t" << setprecision(1) << log(lns)/2.3 << "\t" << -0.01*score << "\t";
			}
			else	{
				oline << sv[s].ri << "\t" << setprecision(1) << "-" << "\t" << -0.01*score << "\t";
			}
			//deal with modifications
			if(js.HasMember("mods"))	{
				const Value& jmods = js["mods"];
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
						oline << mods[a].res << mods[a].pos << "~" << mt.find(mods[a].mass)->second << ";";
					}
					else	{
						oline << mods[a].res << mods[a].pos << "#" << mods[a].mass/1000.0 << ";";
					}
				}
			}
			//add output file line to mapped vector for spectrum index s
			if(odict.find(s) != odict.end())	{
				odict[s].push_back(oline.str());			
			}
			else	{
				odict.insert(pair<int64_t,vector<string> >(s,vector<string>()));
				odict[s].push_back(oline.str());
			}
			//add ppm value to ppms histogram
			ppms[s] = (int64_t)(0.5+ppm);
		}
		total_prob += max_prob;
	}
	//open a file stream to output information in odict
	ofstream ofs;
	ofs.open(_params["output file"]); //open output stream
	//define headers for the output TSV file
	string header = "Id\tSub\tScan\tPeptide mass\tDelta\tppm\tz\tProtein acc\t";
	header += "Start\tEnd\tPre\tSequence\tPost\tIC\tRI\tlog(f)\tlog(p)\tModifications";
	ofs << header << endl;
	//determine parent mass delta ppm acceptance window
	find_window();
	int64_t err = 0;
	int64_t sub = 0;
	int64_t tot = 0;
	string t;
	//loop through result lines and record the information
	for(int64_t a = 0; a < (int64_t)odict.size(); a++)	{
		sub = 1;
		for(size_t b = 0; b < odict[a].size(); b++)	{
			t = odict[a][b];
			double ps = get_ppm(t);
			if(ps <= high and ps >= low)	{
				ofs << a << "\t" << sub << "\t" << t << endl;
				sub++;
				tot++;
			}
			else	{
				err++;
			}
		}
	}
	//output some additional information for logging
	cout << "\n     lines = " << tot << endl;
	if(total_prob > 0)	{
		cout << "     fpr = " << scientific << setprecision(1) << total_prob << endl;
	}
	if(low != -20 and high != 20)	{
		double ble = 100.0*((high-low)/41.0)*(double)err/(double)tot;
		cout << "     baseline error = " << fixed << setprecision(1) << ble << "%" << endl;
	}
	else	{
		cout << "     baseline error = n/a" << endl;
	}
	cout << "     parent ion tolerance = " << fixed << setprecision(0) << low << "," << high << endl;
	

	ofs.close();
	cout << endl;
	cout.flush();
	istr.close();
	return true;
}

