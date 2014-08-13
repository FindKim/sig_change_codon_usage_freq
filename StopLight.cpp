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

// Can't declare const for seqs b/c prevents from random access with key map[]
StopLight :: StopLight (const string& homologID, const vector<float>& pos, const vector<float>& pv, vector<pair<string,vector<string> > > seqs) {

	set_groupID(homologID);
	set_position(pos);
	set_pvalue(pv);
	set_MM_map(seqs);
	create_StopLight_regions(MM);
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
void StopLight :: create_StopLight_map(map<string, vector<float> >& seqs, const map<string, vector<pair<int,int> > >& stoplight_regions_homolog) {

	map<string, vector<int> > stoplight_map;

	// Iterates through each orf in homolog group
	map<string, vector<pair<int, int> > >::const_iterator orf = stoplight_regions_homolog.begin();
	for (orf; orf != stoplight_regions_homolog.end(); ++orf) {
	
		string orf_id = orf->first;
		vector<float> current_mm_seq = seqs[orf_id];
		// Iterates through regions within orf & creates a parallel sequence of stoplight values
		vector<int> stoplight_seq (current_mm_seq.size(), 0);	// Instantiate with same size
		vector<pair<int,int> >::const_iterator region = orf->second.begin();
		for (region; region+1 != orf->second.end(); ++region) {
			
			TTest check_sig(current_mm_seq, *region, *(region+1));

			// Rare codon cluster
			if (current_mm_seq[(region+1)->first] < 0) {
				for (int i = (region+1)->first; i < (region+1)->second; i++)
					stoplight_seq[i] = -2;

			// Not rare codon cluster, but significant
			} else if (check_sig.significant()) {
				
				// First region is greater
				if (check_sig.greater() == 1) {
						
					for (int i = (region+1)->first; i < (region+1)->second; i++) {
						stoplight_seq[i] = -1;	
					}
				
				// First region is less
				} else {
					for (int i = (region+1)->first; i < (region+1)->second; i++) {
						stoplight_seq[i] = 1;
					}
				}
			}
		}
		
		stoplight_map.insert(pair<string, vector<int> > (orf_id, stoplight_seq));
	}
	SL = stoplight_map;
}




// Creates parallel mapping of %MM & stoplight values
void StopLight :: create_StopLight_regions(map<string, vector<float> >& seqs) {
	map<string, vector<pair<int,int> > > stoplight_regions_homolog;
	
	// Iterate through all sequences in homolog group
	map<string, vector<float> >::const_iterator it = seqs.begin();
	for (it; it != seqs.end(); ++it) {
	
		vector<pair<int,int> > stoplight_regions;
		string ID = it->first;
		cout << ID << endl;
		
		// Iterate through %MM sequence for each orfeome
		// Stores start & end of each region
		bool first_region = true;
		int start = 0;
		int end = 0;
		pair<int, int> regionA, regionB;
		vector<float> mm = it->second;
		for (int i = 0; i < mm.size(); i++) {
//			cout << mm[i] << endl;
			// Moves on to next mm, sl @ that pos remains 0 for gap
			if (isnan(mm[i])) continue;
			
			int highest_local_mm = 0;
			// IF FIRST REGION
			if (first_region) {
				start = i;
				cout << "FIRST REGION" << endl;
				
				// First region begins with rare codon cluster
				if (mm[i] < 0) {
					while(mm[i] < 0) {
						if (i < mm.size() || isnan(mm[i])) i++;
						else break;
					}
					
					int j = 1;
					bool rcclust = true;
					while (rcclust) {
						bool still_rcclust = false;
						int last_rcclust = 0;
						for (int k = 0; j+k <= 10; j++) {
							if (j+i > mm.size()) break;
							else if (isnan(mm[j+i])) k--;
							else if (mm[j+i] < 0) {
								last_rcclust = i+j;
								still_rcclust = true;
							}
						}
						if (still_rcclust) i = last_rcclust;
						else rcclust = false;
						i++;
					}
					
					end = i;
					cout << mm[start] << ", " << mm[end] << endl;
					cout << "First rare codon cluster: " << start << ", " << end << endl;
				
				// First region not a rcc -- find highest local peak
				} else {
					highest_local_mm = i;
					int j = 0;
					for (int k = 0; j+k < 10; j++) {	// 10 b/c of rcclust length
						if (j+i > mm.size()) break;			// Out of bounds
						else if (isnan(mm[j+i])) k--;		// Counter balances j so only iterate over 10 windows that aren't gaps
						else if (mm[highest_local_mm] < mm[j+i]) highest_local_mm = j+i;
					}
					end = highest_local_mm;
					i = end;	// Moves iterator to end of region
					cout << "First region: " << start << ", " << end << endl;
				}
				
				regionB = make_pair (start,end);	// Will swap to A for convention
				stoplight_regions.push_back(regionB);
				first_region = false;
				// Finished dealing with first region
			
			
			
			// NOT FIRST REGION
			} else {
				if (isnan(mm[i])) continue;
				cout << "ANOTHER REGION" << endl;
				
				regionA = regionB;	// B -> A, find new B
				start = i;

				// Region begins with rare codon cluster
				if (mm[i] < 0) {
					while(mm[i] < 0) {
						if (i < mm.size()) i++;
						else break;
					}
					end = i;
					cout << start << ", " << end << endl;

					pair<int,int> temp_region = make_pair (start, end);
					TTest check_sig (mm, regionA, temp_region);
					bool sig = check_sig.significant();
					cout << "Is it sig? " << sig << endl;
					while (sig) {
						temp_region.second++;
						TTest check_next_sig (mm, regionA, temp_region);
						sig = check_next_sig.significant();
						cout << "Is next one sig? " << sig << endl;
					}
					cout << "Consecutive rcclust: " << temp_region.first << ", " << temp_region.second-1 << endl;
					end = temp_region.second -1;	// Moves interator to end of temp_region (-1 because last one made it not signficant)
					i = end;

				// Region doesn't begin with rcc
				} else {
					int j = 0;
					bool cluster = true;
					highest_local_mm = i;
					while (cluster) {
						bool same_cluster = false;
						for (int k = 0; j+k < 10; j++) {	// 10 b/c of rcclust length
							cout << i+j << ") " << mm[i+j] << endl;
							if (i+j > mm.size()) break;			// Out of bounds
							else if (isnan(mm[i+j])) k--;		// Counter balances j so only iterate over 10 windows that aren't gaps
							else if (mm[highest_local_mm] < mm[i+j]) {
								highest_local_mm = i+j;
								same_cluster = true;
								cout << "highest local mm " << mm[i+j] << endl;
							}
						}
						if (same_cluster) i = highest_local_mm;
						else cluster = false;
						i++;
					}
					end = highest_local_mm;
					i = end;	// Moves iterator to end of region
					
				}

				regionB = make_pair (start,end);
				cout << "Region A: " << regionA.first << "," << regionA.second << endl << "\t";
				for (int n = regionA.first; n <= regionA.second; n++) cout << mm[n] << " ";
				cout << endl;
				cout << "Region B: " << regionB.first << "," << regionB.second << endl << "\t";
				for (int n = regionB.first; n <= regionB.second; n++) cout << mm[n] << " ";
				cout << endl;
				TTest check_sig (mm, regionA, regionB);
				bool sig = check_sig.significant();
				if (sig) stoplight_regions.push_back (regionB);
				else {
					cout << "MERGING PREVIOUS REGION" << endl;
					cout << regionA.first << "," << regionA.second << "\t" << regionB.first << "," << regionB.second << endl;
					regionB.first = regionA.first;	// regionB takes start of prev region
					vector<pair<int, int> >::iterator prev_region = stoplight_regions.end()-1;
					prev_region->second = regionB.second;	// Replace previous region with extended region
					cout << "Merged to: " << prev_region->first << "," << prev_region->second << endl;
				}
			}
		}
		stoplight_regions_homolog.insert (pair<string, vector<pair<int,int> > > (ID, stoplight_regions));
	}
	create_StopLight_map (seqs, stoplight_regions_homolog);
}

