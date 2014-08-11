//
//  StopLight.cpp
//  
//
//  Created by Kim Ngo Aug 11, 2014.

//  StopLight class creates & stores parallel vector of values to %Min-Max sequence that denote regions that correspond with fast, slow, and slowest translational speed
//	1:Fast--region with high codon usage frequency
//	0:GAPS
//	-1:Slower--region with relatively lesser codon usage frequency
//	-2:Slow--region with rare codon clusters

#include "StopLight.h"
#include "TTest.h"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <math.h>
#include <limits> // NAN
#include <sstream>

using namespace std;

StopLight :: StopLight (const string& homologID, const vector<float>& pos, const vector<float>& pv, const vector<pair<string,vector<string> > >& seqs) {

	set_groupID(homologID);
	set_position(pos);
	set_pvalue(pv);
	set_MM_map(seqs);
}

// Creates output file with all seqID, gapped-aligned %Min-Max sequence, and stoplight sequence for homolog group
// In columns
void StopLight :: create_output_file(const string& filename) {

}



void StopLight :: set_groupID(const string& ID) {
	groupID = ID;
}
string StopLight :: get_groupID() {
	return groupID;
}


// Sets position vector by parsing string of positions
void StopLight :: set_position(const vector<float>& pos_v) {
	position = pos_v;
}
// Returns vector of sequentional positions
vector<float> StopLight :: get_position() {
	return position;
}


// Sets pvalue vector by parsing string of values
void StopLight :: set_pvalue(const vector<float>& pv_v) {
	pvalue = pv_v;
}
// Returns vector of sequential pvalues
vector<float> StopLight :: get_pvalue() {
	return pvalue;
}


// creates map of seqID & its %Min-Max sequence
void StopLight :: set_MM_map(const vector<pair<string, vector<string> > >& seqs) {
	map<string, vector<float> > seq_map;
	
	for (int n = 0; n < seqs.size(); n++) {
		string ID = seqs[n].first;
		vector<float> mm_v (seqs[n].second.size(), 0);
		
		vector<string>::const_iterator it = seqs[n].second.begin();
		vector<float>::iterator mm_it = mm_v.begin();
		for (it; it != seqs[n].second.end(); ++it, ++mm_it) {
		
			if (*it != "-") {
				stringstream ss(*it);
				float i;
				while (ss >> i) {
					*mm_it = i;
				}
			} else *mm_it = numeric_limits<float>::quiet_NaN();
		}
		
		seq_map.insert(pair<string, vector<float> >(ID,mm_v));
	}
	MM = seq_map;
}
// Returns map of seqID & its %Min-Max sequence
map<string, vector<float> > StopLight :: get_MM() {
	return MM;
}


// Returns map of seqID & its stop light sequence
map<string, vector<int> > StopLight :: get_SL() {
	return SL;
}
// Creates parallel mapping of %MM & stoplight values
void StopLight :: create_StopLight_map(const map<string, vector<float> >& seqs) {
	map<string, vector<int> > stoplight;
	
	map<string, vector<float> >::const_iterator it = seqs.begin();
	for (it; it != seqs.end(); ++it) {
	
		string ID = it->first;
		vector<int> sl_v (it->second.size(), 0);
		vector<int>::iterator sl_it = sl_v.begin();
		
		bool first_region = true;
		vector<float>::const_iterator mm_it = it->second.begin();
		for (mm_it, sl_it; mm_it != it->second.end(); ++mm_it) {
			if (isnan(*mm_it)) continue;	// Moves on to next mm, sl @ that pos remains 0 for gap
			
			// IF FIRST REGION
			if (first_region) {
				if (*mm_it < 0) {
					while(*mm_it < 0) {
						++mm_it;
					}
				
				} else {
				
				}
				first_region = false;
			}
			
			// NOT FIRST REGION
			else {
			
			
			
			
			
			}
		}
	}
	SL = stoplight;
}

