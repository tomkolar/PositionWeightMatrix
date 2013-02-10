/*
 * GenbankCodingSequence.cpp 
 *
 *	This is the cpp file for the GenbankCodingSequence object. This 
 *  object represents the coding sequence information found in a Genbank
 *  file.  At this point in time it only contains information pertaining
 *  to the location for the coding sequence.  Future implementations will
 *  include other data elements associated with the Genbank CDS element.
 *
 *  Created on: 2-6-13
 *      Author: tomkolar
 */

#include "GenbankCodingSequence.h"
#include "StringUtilities.h"
#include <string>
#include <vector>
#include <iostream>
using namespace std;


// Constuctors
// ==============================================
GenbankCodingSequence::GenbankCodingSequence(){
}

GenbankCodingSequence::GenbankCodingSequence(string locationInfo) {
	complement = false;
	populate(locationInfo);
}

// Destructor
// =============================================
GenbankCodingSequence::~GenbankCodingSequence() {
}
// Public Methods
// =============================================

// spliceSiteLocations()
//  Purpose:
//		Returns the locations for the start of a coding sequence
//		+/- 10
//		The map has keys -10,-9,...0...9,10 and values of the
//		location of the residue in the sequence or complement
map<int, int> GenbankCodingSequence::spliceSiteLocations() {
	map<int, int> spliceSiteLocations;
	if (!isComplement()) {
		int start = locationMap.at(0);
		spliceSiteLocations[0] = start;
		for (int i = 1; i <= 10; i++) {
			spliceSiteLocations[-i] = start - i;
			spliceSiteLocations[i] = locationMap.at(i);
		}
	}
	else {
		int mapEnd = locationMap.size() -1;
		int start = locationMap.at(mapEnd);
		spliceSiteLocations[0] = start;
		for (int i = 1; i <= 10; i++) {
			spliceSiteLocations[-i] = start + i;
			spliceSiteLocations[i] = locationMap.at(mapEnd - i);
		}
	}

	return spliceSiteLocations;
}

// start()
//  Purpose:
//		Returns the start location of this coding sequence
const int GenbankCodingSequence::start() {
	if (!isComplement()) 
		return locationMap.at(0);
	
	int mapEnd = locationMap.size() -1;
	return locationMap.at(mapEnd);
}

// Public Accessors
// =============================================
const bool GenbankCodingSequence::isComplement() {
	return complement;
} 

const map<int, int>& GenbankCodingSequence::locations() {
	return locationMap;
}

// Private Methods
// =============================================

// populate()
//  Purpose:
//		Parses the locationInfo string to populate the object with locations
//  Postconditions:
//		isComplement - set to true if complement is first word in
//					   locationInfo
//		locaitons - populated with locations defined in locationInfor
void GenbankCodingSequence::populate(string& locationInfo) {

	int locationsStartPos = 0;
	int locationsEndPos = locationInfo.length(); 

	// Check for complement
	if (locationInfo.substr(0,11) == "complement(") {
		complement = true;
		locationsStartPos +=11;
		locationsEndPos --;
	}

	if (locationInfo.substr(locationsStartPos, 5) == "join(") {
		locationsStartPos += 5;
		locationsEndPos --;
	}

	// Split the remaining string into location tokens
	vector<string> infoTokens;
	StringUtilities::split(locationInfo.substr(locationsStartPos, locationsEndPos-locationsStartPos), ',', infoTokens);

	// Populate the locations map
	int mapLocation = 0;
	for (string infoToken : infoTokens) {
		// Split token string into start and stop locations
		vector<string> locationTokens;
		StringUtilities::split(infoToken, '.', locationTokens);
		
		int locationStart = atoi(locationTokens.front().c_str());
		int locationEnd = atoi(locationTokens.back().c_str());

		for (int location = locationStart; location <= locationEnd; location++) {
			locationMap[mapLocation] = location;
			mapLocation++;
		}
	}
}
