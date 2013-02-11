/*
 * driver.cpp
 *
 *	This is the driver file for creating a position weight matrix
 *  for coding locations from a genbank file.
 *
 *	Typical use:
 *		pwm genbankFile
 *
 *  Created on: 2-7-13
 *      Author: tomkolar
 */
#include "GenbankFile.h"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main( int argc, char *argv[] ) {

	// Check that file name was  entered as argument
	if (argc < 2) {
		cout << "Invalid # of arguments\n";
		cout << "usage: pwm genbankFile\n";
		return -1;
	}

	cout << "Starting\n";

	// Get Genbank File name
	string genbankFileName = argv[1];
/*
	// Test file name
	string genbankFileName = "c:/Users/kolart/Documents/Genome540/Assignment4/NC_004317.gbk";

	// Set Genbank File name
//	string genbankFileName = "c:/Users/kolart/Documents/Genome540/Assignment4/NC_004353.gbk";
*/
	// Create the Genbank File object
	GenbankFile* genbankFile = new GenbankFile(genbankFileName);
	cout << genbankFile->firstLineResultString();
	cout << genbankFile->backgroundCountsResultString();
	cout << genbankFile->backgroundFrequenciesResultString();
	cout << genbankFile->spliceSiteCountsResultString();
	cout << genbankFile->spliceSiteFrequenciesResultString();
	cout << genbankFile->spliceSiteWeightsResultString();
	cout << genbankFile->spliceSiteScoreCountsResultString();
	cout << genbankFile->siteScoreCountsResultString();
	cout << genbankFile->highNonCodingSiteScoresResultString();

	cout << "Genbank's done\n";

	genbankFile = NULL;
}
