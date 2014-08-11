//
//  Parse_pvalue_ortholog_map.h
//  
//
//  Created by Kim Ngo July 23, 2014.

//  Parses pvalue and ortholog mapping files (ultimate/ult_mm_splice_pvalue_mapped) into an array
//	Array for position, pvalue, and listed %MinMax values

#ifndef PARSE_PVALUE_ORTHOLOG_MAP_H
#define PARSE_PVALUE_ORTHOLOG_MAP_H

//#include "Parse_pvalue_ortholog_map.cpp" // for templating

#include <string>
#include <vector>
#include <utility>

using namespace std;

class Parse_pvalue_ortholog_map {

	public:
		Parse_pvalue_ortholog_map(const string&);
		static bool valid_file(const string&);
			// Returns T if ultimateORFS[0-9]+.fasta.mm_spliced_pvalue_mapped
		string get_groupID();
			// Returns group id number
		vector<float> get_position();
			// Returns vector of sequentional positions
		vector<float> get_pvalue();	
			// Returns vector of sequential pvalues
		vector<pair<string,vector<string> > > get_seqs();
			// Pairs of Seq ID and %MinMax values for the ortholog group
	
	
	private:
		string groupID;
		vector<float> position;
		vector<float> pvalue;
		vector< pair<string,vector<string> > > seqs;
		void set_groupID(const string&);
			// Sets groupID by stripping filename
		void set_position(const string&);
			// Sets position vector by parsing string of positions
		void set_pvalue(const string&);
			// Sets pvalue vector by parsing string of values
		pair<string,vector<string> > create_seq_pair(const string&, const string&);
			// Parses string of %MinMax values into vector and creates pair with seqID
		void add_pair(const pair<string,vector<string> >&);
			// Adds pair to seqs vector
		vector<float> parse_string_vector(const string&, float);
			// Parses string with delim=, into vector of floats
		vector<string> parse_string_vector(const string&, string);
			// Parses string with delim=, into vector of strings (overloaded)
		void parse_file(const string&);
			// Parses file line by line
			// Sets groupID, position v, pvalue v, and seqs v of pairs

};

#endif

