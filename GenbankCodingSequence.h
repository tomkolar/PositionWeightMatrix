/*
 * GenbankCodingSequence.h
 *
 *	This is the header file for the GenbankCodingSequence object. This 
 *  object represents the coding sequence information found in a Genbank
 *  file.  At this point in time it only contains information pertaining
 *  to the location for the coding sequence.  Future implementations will
 *  include other data elements associated with the Genbank CDS element.
 *
 *  Created on: 2-6-13
 *      Author: tomkolar
 */

#ifndef GENBANKCODINGSEQUENCE_H
#define GENBANKCODINGSEQUENCE_H

#include <string>
#include <map>
using namespace std;

class GenbankCodingSequence
{

public:

	// Constuctors
	// ==============================================
	GenbankCodingSequence();
	GenbankCodingSequence(string locationInfo);

	// Destructor
	// =============================================
	virtual ~GenbankCodingSequence();

	// Public Methods
	// =============================================

	// spliceSiteLocations()
	//  Purpose:
	//		Returns the locations for the start of a coding sequence
	//		+/- 10
	//		The map has keys -10,-9,...0...9,10 and values of the
	//		location of the residue in the sequence or complement
	map<int, int> spliceSiteLocations();

	// start()
	//  Purpose:
	//		Returns the start location of this coding sequence
	const int start();

	// Public Accessors
	// =============================================
	const bool isComplement(); 
	const map<int, int>& locations();




private:

	// Attributes
	// =============================================
	bool complement;
	map<int, int> locationMap;

	// Private Methods
	// =============================================

	// populate()
	//  Purpose:
	//		Parses the locationInfo string to populate the object with locations
	//  Postconditions:
	//		isComplement - set to true if complement is first word in
	//					   locationInfo
	//		locaitons - populated with locations defined in locationInfor
	void populate(string& locationInfo);
};

#endif //GENBANKCODINGSEQUENCE_H
