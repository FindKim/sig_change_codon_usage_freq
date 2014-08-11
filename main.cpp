//
//  Main.cpp
//  
//	Driver program to create parallel data for %Min-Max sequences
//	Reads file; parses file; prints name, description, and sequence
//
//  Created by Kim Ngo on Aug 11, 2014.
//
//

#include "Parse_pvalue_ortholog_map.h"
#include "StopLight.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <fstream>	// Output file
#include <dirent.h>	// Traverse directory
#include <algorithm> // Sort
#include <ctime>

using namespace std;

//const string DIRECTORY = "/afs/crc.nd.edu/user/k/kngo/co-occ_rcclust/subset/";
const string DIRECTORY = "/afs/crc.nd.edu/group/NDBL/data/rcclust/ultimate/ult_mm_splice_pvalue_mapped/";
const string DIRECTORY_SIGORFS = "/afs/crc.nd.edu/group/NDBL/data/rcclust/ultimate/prunedMaskedSig/";
//const string SIGORF_FILE = "/afs/crc.nd.edu/group/NDBL/data/rcclust/ultimate/prunedMaskedSig/sigOrfs_p.0001.txt";
//const string OUTPUT_FILE = "coords_co-occuring_p.0001.txt";
//const string OUTPUT_FILE = "coords_co-occuring_all.txt";

// Checks if provided directory exists
bool directory_exists(const char* pzPath) {
	if(pzPath == NULL) return false;
	
	DIR *pDir;
	bool bExists = false;
	
	pDir = opendir(pzPath);
	
	if(pDir != NULL) {
		bExists = true;
		(void) closedir (pDir);
	}
	return bExists;
}


// Traverses directory and creates a vector of filenames
vector<string> traverse_directory(string& directory) {

	vector<string> co_occ_files;
	DIR *dpdf;
	struct dirent *epdf;
	dpdf = opendir(directory.c_str());
	if (dpdf) {
		while(true) {
			epdf = readdir(dpdf);
			if (epdf == NULL) break;
			string file = directory;
			file.append(string(epdf->d_name));	// Absolute path
			co_occ_files.push_back(file);
		}
		closedir(dpdf);
		sort(co_occ_files.begin(), co_occ_files.end());
	}
	return co_occ_files;
}


// Initializes outputfile with two column format
void initialize_outputfile(const string& outputfile) {

//	cout << "Creating " << outputfile << "..." << endl;
	ofstream ofile;
	ofile.open (outputfile.c_str());
	
	if (ofile.is_open()) {
	
		string x_coord = "co-occ rcclust length";
		string y_coord = "co-occ rcclust %min sum";
		ofile << x_coord << "," << y_coord << endl;
		ofile.close();

	} else cout << "Unable to open " << outputfile << endl;
}

void timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    printf("%s",asctime( localtime(&ltime) ) );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -------------------------------------CALCULATES LENGTH,SUM COORD

int main() {

	string temp = "/afs/crc.nd.edu/group/NDBL/data/rcclust/ultimate/prunedMaskedSig/sigOrfs_p";
	string directory = DIRECTORY;
	string directory_sigOrfs = DIRECTORY_SIGORFS;
//	string output_file = OUTPUT_FILE;
//	string sigOrf_file = SIGORF_FILE;

	vector<string> mmfiles;
	vector<string> sigOrf_files;
	vector<string> co_occ_files;
/*	
	string output_file = "results/coords_co-occuring_p";
	string output_file_nsig = "results/coords_not_co-occuring_p";
	output_file.append(sigOrf_file_it->substr(temp.length(), sigOrf_file_it->length()));
	output_file_nsig.append(sigOrf_file_it->substr(temp.length(), sigOrf_file_it->length()));
	initialize_outputfile(output_file);
	initialize_outputfile(output_file_nsig);
*/
	if (directory_exists(directory.c_str())) {

		// Traverses directory and returns vector of filenames
		co_occ_files = traverse_directory(directory);
		cout << "Computing Stop Light values at ";
		timestamp();
		cout << endl;

		// Iterates through each file
		vector<string>::iterator file_it = co_occ_files.begin();
		for (file_it; file_it != co_occ_files.end(); ++file_it) {

		if (*file_it == "/afs/crc.nd.edu/group/NDBL/data/rcclust/ultimate/ult_mm_splice_pvalue_mapped/ultimateORFS13888.fasta.mm_spliced_pvalue_mapped") {
		cout << *file_it << endl;
			if (Parse_pvalue_ortholog_map::valid_file(*file_it)) {
				Parse_pvalue_ortholog_map map(*file_it);
				StopLight stoplight(map.get_groupID(), map.get_position(), map.get_pvalue(), map.get_seqs());
			}
		}
		}
	}
	cout << "Finished computing Stop Light values at ";
	timestamp();
	cout << endl;
}


