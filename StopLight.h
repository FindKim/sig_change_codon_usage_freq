//
//  StopLight.h
//  
//
//  Created by Kim Ngo Aug 11, 2014.

//  StopLight class creates & stores parallel vector of values to %Min-Max sequence that denote regions that correspond with fast, slow, and slowest translational speed
//	1:Fast--region with high codon usage frequency
//	0:GAPS
//	-1:Slower--region with relatively lesser codon usage frequency
//	-2:Slow--region with rare codon clusters

#ifndef STOPLIGHT_H
#define STOPLIGHT_H

#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace std;

class StopLight {

	public:
		StopLight(const string&, const vector<float>&, const vector<float>&, vector<pair<string,vector<string> > >);
			// string homolog group ID
			// position #
			// p-values
			// seqID, gapped-aligned %Min-Max sequence
			
		void create_output_file(const string&);
			// Creates output file with all seqID, gapped-aligned %Min-Max sequence, and stoplight sequence for homolog group
			// In columns
		string get_groupID();
			// Returns homolog group ID
		vector<float> get_position();
			// Returns vector of sequentional positions
		vector<float> get_pvalue();	
			// Returns vector of sequential pvalues
		map<string, vector<float> > get_MM();
			// Returns map of seqID & its %Min-Max sequence
		map<string, vector<int> > get_SL();
			// Returns map of seqID & its stop light sequence
	
	private:
		string groupID;						// Homolog group ID
		vector<float> position;		// Sequential Postions for homolog group
		vector<float> pvalue;			// Sequential p-values for homolog group
		map<string, vector<float> > MM;
		map<string, vector<int> > SL;
		void set_groupID(const string&);
			// Sets groupID by stripping filename
		void set_position(const vector<float>&);
			// Sets position vector by parsing string of positions
		void set_pvalue(const vector<float>&);
			// Sets pvalue vector by parsing string of values
		void set_MM_map(const vector<pair<string, vector<string> > >&);
			// creates map of seqID & its %Min-Max sequence
		void create_StopLight_map(map<string, vector<float> >&, const map<string, vector<pair<int,int> > >&);
			// Creates parallel mapping of %MM & stoplight values
		void create_StopLight_regions(map<string, vector<float> >&);
			// breaks up sequences into regions for T-Test
};

#endif

