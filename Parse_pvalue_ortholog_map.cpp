//
//  Parse_pvalue_ortholog_map.cpp
//  
//
//  Created by Kim Ngo July 23, 2014.

//  Parses pvalue and ortholog mapping files (ultimate/ult_mm_splice_pvalue_mapped) into an array
//	Array for position, pvalue, and listed %MinMax values

#include "Parse_pvalue_ortholog_map.h"

#include <iostream>
#include <string>
#include <fstream>	// infile
#include <sstream>	// stringstream
#include <vector>
#include <utility>	// pair
#include <algorithm>// replace
#include <iterator>

const string FILE_BEGIN = "ultimateORFS";
const string FILE_END = ".fasta.mm_spliced_pvalue_mapped";

using namespace std;

Parse_pvalue_ortholog_map :: Parse_pvalue_ortholog_map (const string& filename) {

	set_groupID(filename);
	parse_file(filename);
}


// Sets groupID by stripping filename
void Parse_pvalue_ortholog_map :: set_groupID(const string& filename) {
	int id_len = filename.size() - FILE_BEGIN.size() - FILE_END.size();
	string ID = filename.substr(FILE_BEGIN.size(), id_len);
	groupID = ID;
}
// Returns group ID number
string Parse_pvalue_ortholog_map :: get_groupID() {
	return groupID;
}

void Parse_pvalue_ortholog_map :: set_position(const string& pos_str) {
	float type;
	vector<float> pos_v = parse_string_vector(pos_str, type);
	position = pos_v;
}
// Returns vector of sequential positions
vector<float> Parse_pvalue_ortholog_map :: get_position() {
	return position;
}

void Parse_pvalue_ortholog_map :: set_pvalue(const string& pvalue_str) {
	float type;
	vector<float> pvalue_v = parse_string_vector(pvalue_str, type);
	pvalue = pvalue_v;
}
// Returns vector of sequential pvalues
vector<float> Parse_pvalue_ortholog_map :: get_pvalue() {
	return pvalue;
}
// Returns pairs of seq ID and %MinMax values for that ortholog group
vector<pair<string,vector<string> > > Parse_pvalue_ortholog_map :: get_seqs() {
	return seqs;
}


// Returns T if filename is ultimate[0-9]+.fasta.mm_splice_pvalue_mapped
bool Parse_pvalue_ortholog_map :: valid_file(const string& filename) {
	size_t found_beginning = filename.find(FILE_BEGIN);
	size_t found_ending = filename.find(FILE_END);
	if (found_beginning != string::npos && found_ending != string::npos)
		return true;
	else
		return false;
}

// Parses minmax sequence string into vector of floats=
//template <class T>
vector<float> Parse_pvalue_ortholog_map :: parse_string_vector(const string& str, float type) {

	float i;
	vector<float> v;
	
	stringstream ss(str);
	while (ss >> i) {
		v.push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return v;
}
// overloaded
// Parses minmax sequence string into vector of floats=
vector<string> Parse_pvalue_ortholog_map :: parse_string_vector(const string& str, string type) {

	string new_str = str;
	replace( new_str.begin(), new_str.end(), ',', ' ');	// replaces all , to space
	stringstream ss(new_str);
	
	istream_iterator<string> begin(ss);
	istream_iterator<string> end;
	vector<string> v(begin, end);
	return v;
}

// Parses string of %MinMax values into vector and creates pair with seqID
pair<string,vector<string> > Parse_pvalue_ortholog_map :: create_seq_pair(const string& id, const string& minmax_str) {
	string s;
	vector<string> minmax = parse_string_vector(minmax_str, s);
	pair<string,vector<string> > seq_pair = make_pair(id, minmax);
	vector<string>::iterator it = seq_pair.second.begin();
	return seq_pair;
}

void Parse_pvalue_ortholog_map :: add_pair(const pair<string,vector<string> >& seq_pair) {
	seqs.push_back(seq_pair);
}

void Parse_pvalue_ortholog_map :: parse_file(const string& filename) {

	ifstream file(filename.c_str());
	string header, pos, pv, id, mm_str;
	getline(file, header);	// "position"
	getline(file, pos);			// string of positions
	getline(file, header);	// "pvalue"
	getline(file, pv);			// string of pvalues	
//	cout << pos << endl << endl << pv << endl << endl;
	
	set_position(pos);
	set_pvalue(pv);
	
	while (!file.eof()) {
		getline(file, id);
		getline(file, mm_str);
		if (id.size() > 0 && mm_str.size() > 0) {
			add_pair( create_seq_pair(id, mm_str) );
		}
	}
}
