/*
 * GenbankFile.cpp
 *
 *	This is the cpp file for the GenbankFile object. The GenbankFile object
 *  is a utility object designed to read in a Genbank File and keep the 
 *  information for the file in memory.
 *
 *  ************* WARNING ***********************************
 *  *  There is no error handling in place for this object! *
 *  *  This means that if the file does not exist or is     *
 *  *  formatted incorrectly, you will get an error and     *
 *  *  will not be able to use this object.                 *
 *  *                                                       *
 *  *  If this code gets moved to a production setting      *
 *  *  appropriate error handling should be implemented!    *
 *  *********************************************************
 *
 * Typical use for the file would be to use the GenbankFile(aFileName)
 * constructor to create the object.  This will automatically open the
 * Genbank File specified by aFileName, and read its contents storing them
 * in the firstLine, and sequence, reverseComplement and codingSequencdes
 *  attributes.
 *
 *  Created on: 2-6-13
 *      Author: tomkolar
 */
#include "GenbankFile.h"
#include "StringUtilities.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <set>
using namespace std;

// Constuctors
// ==============================================
GenbankFile::GenbankFile() {
}


GenbankFile::GenbankFile(string aFileName) {
	
	parseFileName(aFileName);
	populate();
}

// Destructor
// ==============================================
GenbankFile::~GenbankFile() {
}

// Public Methods
// =============================================

// string firstLineResultString()
//  Purpose:
//		Returns the string value of an XML element representing the first line of 
//		the Genbank file.
//
//		format:
//			<result type='first line' file='<<fileName>>' >
//				<<firstLine>>
//			</result>
//  Preconditions:
//		Genbank File has been read and firstLine has been populated
string GenbankFile::firstLineResultString() {
	stringstream ss;

	ss << "    <result type=\"first line\" file=\"" << fileName << "\">\n";
	ss << "      " << firstLine << "\n";
	ss << "    </result>\n";

	return ss.str();
}


// string backgroundCountsResultString()
//  Purpose:
//		Returns the string value of an XML element representing the background
//		counts of the sequence.
//
//		format:
//			<result type='nucleotide histogram' file='<<fileName>>' >
//				<<residue>>=<<backgroundCountForResidue>>,
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::backgroundCountsResultString() {
	stringstream ss;

	// Header
	ss <<  "    <result type=\"nucleotide histogram\" file=\"" << fileName << "\"'>\n";

	// Counts
	ss << "      ";
	ss << "A=" << backgroundCounts['a'] << ",";
	ss << "C=" << backgroundCounts['c'] << ",";
	ss << "G=" << backgroundCounts['g'] << ",";
	ss << "T=" << backgroundCounts['t'] << ",";
	ss << "N=" << backgroundCounts['n'];
	ss << "\n";

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string backgroundFrequenciesResultString()
//  Purpose:
//		Returns the string value of an XML element representing the background
//		frequencies of the sequence.
//
//		format:
//			<result type='nucleotide histogram' file='<<fileName>>' >
//				<<residue>>=<<backgroundFrequencyForResidue>>,
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::backgroundFrequenciesResultString() {
	stringstream ss;
	ss.precision(3);

	// Header
	ss <<  "    <result type=\"background frequency\" file=\"" << fileName << "\"'>\n";

	// Counts
	ss << "      ";
	ss << "A=" << backgroundFrequencies['a'] << ",";
	ss << "C=" << backgroundFrequencies['c'] << ",";
	ss << "G=" << backgroundFrequencies['g'] << ",";
	ss << "T=" << backgroundFrequencies['t'];
	ss << "\n";

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string spliceSiteCountsResultString()
//  Purpose:
//		Returns the string value of an XML element representing the counts 
//		of the splice site residues.
//
//		format:
//			<result type='count matrix' file='<<fileName>>' >
//				(<<siteIndex>>,<<residue>>)=<<residueCount>>,
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::spliceSiteCountsResultString() {
	stringstream ss;

	// Header
	ss <<  "    <result type=\"count matrix\" file=\"" << fileName << "\">\n";

	for (int i = -10; i <= 10; i++) {
		ss << "      ";
		ss << "(" << i <<",A)=" << spliceSiteCounts[i]['a'] << ", ";
		ss << "(" << i <<",C)=" << spliceSiteCounts[i]['c'] << ", ";
		ss << "(" << i <<",G)=" << spliceSiteCounts[i]['g'] << ", ";
		ss << "(" << i <<",T)=" << spliceSiteCounts[i]['t'] << ", ";
		ss << "\n";
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string spliceSiteFrequenciessResultString()
//  Purpose:
//		Returns the string value of an XML element representing the  
//		frequency of the splice site residues.
//
//		format:
//			<result type='count matrix' file='<<fileName>>' >
//				(<<siteIndex>>,<<residue>>)=<<residueFrequency>>,
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::spliceSiteFrequenciesResultString() {
	stringstream ss;
	ss.setf(ios::fixed,ios::floatfield);
	ss.precision(3);

	// Header
	ss <<  "    <result type=\"frequency matrix\" file=\"" << fileName << "\">\n";

	for (int i = -10; i <= 10; i++) {
		ss << "      ";
		ss << "(" << i <<",A)=" << spliceSiteFrequencies[i]['a'] << ", ";
		ss << "(" << i <<",C)=" << spliceSiteFrequencies[i]['c'] << ", ";
		ss << "(" << i <<",G)=" << spliceSiteFrequencies[i]['g'] << ", ";
		ss << "(" << i <<",T)=" << spliceSiteFrequencies[i]['t'] << ", ";
		ss << "\n";
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string spliceSiteWeightsResultString()
//  Purpose:
//		Returns the string value of an XML element representing the  
//		weights of the splice site residues.
//
//		format:
//			<result type='count matrix' file='<<fileName>>' >
//				(<<siteIndex>>,<<residue>>)=<<residueWeight>>,
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::spliceSiteWeightsResultString() {
	stringstream ss;
	ss.setf(ios::fixed,ios::floatfield);
	ss.precision(3);

	// Header
	ss <<  "    <result type=\"weight matrix\" file=\"" << fileName << "\">\n";

	for (int i = -10; i <= 10; i++) {
		ss << "      ";
		ss << "(" << i <<",A)=" << spliceSiteWeights[i]['a'] << ", ";
		ss << "(" << i <<",C)=" << spliceSiteWeights[i]['c'] << ", ";
		ss << "(" << i <<",G)=" << spliceSiteWeights[i]['g'] << ", ";
		ss << "(" << i <<",T)=" << spliceSiteWeights[i]['t'] << ", ";
		ss << "\n";
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string spliceSiteScoreCountsResultString()
//  Purpose:
//		Returns the string value of an XML element representing the  
//		score counts of the splice site residues.
//
//		format:
//			<result type='score histogram' positions='true sites' file='<<fileName>>' >
//				(<<score>>,<<score count>>),
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::spliceSiteScoreCountsResultString() {
	stringstream ss;

	// Header
	ss <<  "    <result type=\"score histogram\" positions=\"true sites\" file=\"" << fileName << "\">\n";

	// Counts
	for (int i = -20; i <= 20; i++) {
		int count = spliceSiteScoreCounts[i] ;
		if (count > 0) {
			ss << "      ";
			ss << "(" << i <<"," << count << ")";
			ss << "\n";
		}
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string siteScoreCountsResultString()
//  Purpose:
//		Returns the string value of an XML element representing the  
//		score counts of all site residues.
//
//		format:
//			<result type='score histogram' positions='all' file='<<fileName>>' >
//				(<<score>>,<<score count>>),
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::siteScoreCountsResultString() {
	stringstream ss;

	// Header
	ss <<  "    <result type=\"score histogram\" positions=\"all\" file=\"" << fileName << "\">\n";

	// Counts
	for (int i = -20; i <= 20; i++) {
		int count = siteScoreCounts[i] ;
		if (count > 0) {
			ss << "      ";
			ss << "(" << i <<"," << count << ")";
			ss << "\n";
		}
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string highNonCodingSiteScoresResultString()
//  Purpose:
//		Returns the string value of an XML element representing the  
//		sites that are non coding yet scored >= 10.0
//
//		format:
//			<result type='position list' file='<<fileName>>' >
//				(<<position>>,<<strand>>,<<score>>),
//			</result>
//  Preconditions:
//		Genbank File has been read and sequence has been populated
string GenbankFile::highNonCodingSiteScoresResultString() {
	stringstream ss;

	// Header
	ss <<  "    <result type=\"position list\" file=\"" << fileName << "\">\n";

	// Positions
	int count = 1;
	for (SiteScore siteScore : highNonCodingSiteScores) {
		ss << "      ";
		ss << "(" << siteScore.site << ",";
		ss << siteScore.strand << ",";
		ss << siteScore.score << "), ";
		if (count % 5 == 0)
			ss << "\n";
		count++;
	}

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// Public Accessors
// =============================================
const int GenbankFile::getSequenceLength() {
	return sequence.length();
}

string& GenbankFile::getFileName() {
	return fileName;
}

string& GenbankFile::getSequence() {
	return sequence;
}

vector<GenbankCodingSequence>& GenbankFile::getCodingSequences() {
	return codingSequences;
}

// Private Methods
// =============================================

// parseFileName(string& aFileName)
//  Purpose:
//		Parse aFileName into the filePath and fileName attributes
//  Postconditions:
//		fileName and filePath have been set
void GenbankFile::parseFileName(string& aFileName) {

	// Check if in local file
	size_t lastSlashPosition = aFileName.find_last_of('/');
	if (lastSlashPosition == string::npos) {
		filePath = ".";
		fileName = aFileName;
	}
	else { // Path is in file name so parse it out
		filePath = aFileName.substr(0, lastSlashPosition);
		fileName = aFileName.substr(lastSlashPosition + 1, string::npos);
	}

}

// populate()
//  Purpose:
//		Reads in the Genbank File specified by filePath and fileName and populates
//		the object with its contents
//	Preconditions:
//		fileName and filePath have been set
//  Postconditions:
//		firstLine - populated with first line from file
//		sequence - populated with sequence from file
//		reverseComplement - populated with reverse complement of sequence
//		codingSequences - populated from all CDS elements in genbank file
void GenbankFile::populate() {

	ifstream inputFile(filePath + "/" + fileName);
	stringstream dnaSeqStream;
	stringstream codingSeqStream;
	string line;

	getline(inputFile, firstLine);
		
	while(getline(inputFile, line)) {
		// Look for Coding Sequences and parse
		if (line.substr(5,3) == "CDS") {
			// If coding sequence start or end is uncertain then ignore
			if (line.find('<') != string::npos || line.find('>') != string::npos)
				continue;

			// Add line to codingSeq stream
			codingSeqStream.str("");
			codingSeqStream << line.substr(21, string::npos);
			
			// Read rest of lines that define coding sequence location
			bool foundEndLocation = false;
			while (!foundEndLocation) {
				getline(inputFile, line);

				// Check if we fond end of location
				if (line.find('/') != string::npos)
					foundEndLocation = true;
				else {
					// If coding sequence start or end is uncertain then ignore this CDS
					if (line.find('<') != string::npos || line.find('>') != string::npos)
						break;

					// Add line to codingSeq stream
					codingSeqStream << line.substr(21, string::npos);
				}
			}
			
			if (foundEndLocation) {
				GenbankCodingSequence codingSequence(codingSeqStream.str());
				codingSequences.push_back(codingSequence);
			}
			
		}

		// Break if beginning of DNA Sequence
		if (line.substr(0,6) == "ORIGIN")
			break;
	}

	// Exited above loop so must be at first line of DNA Sequence
	while(getline(inputFile, line)) {
		// Check for end of sequence
		if (line.substr(0,2) == "//")
			break;

		if (line.length() == 75) {
			// Add the dna sequence to the stream
			dnaSeqStream << line.substr(10,10);
			dnaSeqStream << line.substr(21,10);
			dnaSeqStream << line.substr(32,10);
			dnaSeqStream << line.substr(43,10);
			dnaSeqStream << line.substr(54,10);
			dnaSeqStream << line.substr(65,10);
		}
		else {// special handling for last line
			vector<string> tokens;
			StringUtilities::split(line.substr(10, string::npos), ' ', tokens);
			for (string token : tokens) {
				dnaSeqStream << token;
			}
		}

	}
	
	sequence = dnaSeqStream.str();

	inputFile.close();

	createComplement();
	calculateBackgroundCounts();
	calculateBackgroundFrequencies();
	calculateSpliceSiteCounts();
	calculateSpliceSiteFrequencies();
	calculateSpliceSiteWeights();
	calculateSpliceSiteScoreCounts();
	calculateSiteScoreCounts();
}

// createComplment()
//  Purpose:
//		populates the complement attribute with the comlpement
//		of the sequence
//	Preconditions:
//		dnaSequence has been set
//  Postconditions:
//		complement - populated with complement of sequence
void GenbankFile::createComplement() {
	stringstream ss;

	// Append the complements in order
	for (string::size_type i = 0; i < sequence.length(); i++) {
			ss << complementOf(sequence[i]);
	}

	// Reverse the string and set it to the reverseComplement attribute
	complement = ss.str();
}

// char complementOf(char aChar)
//  Purpose:  returns the dna complement of aChar
char GenbankFile::complementOf(char aChar) {
	if (aChar == 'a')
		return 't';

	if (aChar == 't')
		return 'a';

	if (aChar == 'g')
		return 'c';

	if (aChar == 'c')
		return 'g';

	return aChar;
}

// calculateBackgroundCounts()
//  Purpose:
//		populates the background counts map with the counts for residues
//		in the sequence (and its complement)
void GenbankFile::calculateBackgroundCounts() {
	
	// initialzie backgroundCounts
	backgroundCounts['a'] = 0;
	backgroundCounts['c'] = 0;
	backgroundCounts['g'] = 0;
	backgroundCounts['t'] = 0;
	backgroundCounts['n'] = 0;

	calculateBackgroundCounts(sequence);
	calculateBackgroundCounts(complement);
}

// calculateBackgroundCounts(string& aSequence)
//  Purpose:
//		populates the background counts map with the counts for residues
//		in aSequence
void GenbankFile::calculateBackgroundCounts(string& aSequence) {
	for (string::size_type i=0; i < aSequence.length(); i++) {
		char residue = aSequence[i];
		if (residue == 'a' || residue == 'c' || residue == 'g' || residue == 't')
			backgroundCounts[residue]++;
		else
			backgroundCounts['n']++;
	}
}

// calculateBackgroundFrequencies()
//  Purpose:
//		populates the background frequencies map with the frequencies
//		of residues in the sequence (and its complement)
void GenbankFile::calculateBackgroundFrequencies() {
	
	// initialzie backgroundFrequencies
	backgroundFrequencies['a'] = 0;
	backgroundFrequencies['c'] = 0;
	backgroundFrequencies['g'] = 0;
	backgroundFrequencies['t'] = 0;

	// Get the total count for the sequence
	int backgroundTotal = 0;
	for (pair<const char, int> residueCount : backgroundCounts) {
		backgroundTotal += residueCount.second;
	}

	// Calculate the frequenies 
	for (pair<const char, int> residueCount : backgroundCounts) {
		char residue = residueCount.first;
		if (residue != 'n')
			backgroundFrequencies[residue] = (double) residueCount.second / backgroundTotal;
	}

}

// calculateSpliceSiteCounts
//  Purpose:
//		Calculate the splice site counts and populate the
//		spliceSiteCounts map with the values.
//
//		Counts are from 10 sites before the start of the
//		coding sequence throught 10 sites after the start
//		of the coding sequence.
//  Postconditions:
//		spliceSiteCounts will be populated
void GenbankFile::calculateSpliceSiteCounts() {
	// initialzie spliceSiteCounts
	for (int i = -10; i <= 10; i++) {
		map<char, int> countsMap;
		countsMap['a'] = 0;
		countsMap['c'] = 0;
		countsMap['g'] = 0;
		countsMap['t'] = 0;
		spliceSiteCounts[i] = countsMap;
	}

	// Iterate throught the coding sequences
	for (GenbankCodingSequence codingSequence : codingSequences) {
		
		// Get the location on the sequence (or complement) of the splice sites 
		map<int, int>& spliceSiteLocations = codingSequence.spliceSiteLocations();

		// Iterate through the splice sites
		for (int i = -10; i <= 10; i++) {

			// Find the residue that is at the site
			char residue;
			if (!codingSequence.isComplement())
				residue = sequence.at(spliceSiteLocations.at(i) -1);
			else
				residue = complement.at(spliceSiteLocations.at(i) -1);

			// Only increment if known residue
			if (residue == 'a' || residue == 'c' || residue == 'g' || residue == 't')
				spliceSiteCounts[i][residue]++;
		}
	}
}

// calculateSpliceSiteFrequencies
//  Purpose:
//		Calculate the splice site frequencies and populate the
//		spliceSiteFrequencies map with the values.
//  Postconditions:
//		spliceSiteFrequencies will be populated
void GenbankFile::calculateSpliceSiteFrequencies() {
	// initialzie spliceSiteFrequencies
	for (int i = -10; i <= 10; i++) {
		map<char, double> frequenciesMap;
		frequenciesMap['a'] = 0;
		frequenciesMap['c'] = 0;
		frequenciesMap['g'] = 0;
		frequenciesMap['t'] = 0;
		spliceSiteFrequencies[i] = frequenciesMap;
	}

	// Iterate through the splice sites
	for (int i = -10; i <= 10; i++) {

		// Get counts for this site
		map<char, int>& countsMap = spliceSiteCounts.at(i);

		// Get the total count for this site
		int siteTotal = 0;
		for (pair<const char, int> residueCount : countsMap) {
			siteTotal += residueCount.second;
		}

		// Calculate the frequenies for this site
		map<char, double>& frequenciesMap = spliceSiteFrequencies.at(i);
		for (pair<const char, int> residueCount : countsMap) {
			frequenciesMap[residueCount.first] = (double) residueCount.second / siteTotal;
		}
	}
}

// calculateSpliceSiteWeights()
//  Purpose:
//		Calculate the splice site weights and populate the
//		calculateSpliceSiteWieghts map with the values.
//  Postconditions:
//		spliceSiteWeights will be populated
void GenbankFile::calculateSpliceSiteWeights() {
	// initialzie spliceSiteWeights
	for (int i = -10; i <= 10; i++) {
		map<char, double> weightsMap;
		weightsMap['a'] = 0;
		weightsMap['c'] = 0;
		weightsMap['g'] = 0;
		weightsMap['t'] = 0;
		spliceSiteWeights[i] = weightsMap;
	}

	// Iterate through the splice sites
	for (int i = -10; i <= 10; i++) {

		// Get frequencies for this site
		map<char, double>& frequenciesMap = spliceSiteFrequencies.at(i);

		// Calculate the weights for this site
		map<char, double>& weightsMap = spliceSiteWeights.at(i);
		for (pair<const char, double> freqPair : frequenciesMap) {
			char residue = freqPair.first;
			double residueFreq = freqPair.second;

			// Calculate score
			double score;
			if (residueFreq == 0)
				score = -99;
			else
				score = (log(residueFreq / backgroundFrequencies.at(residue))) / log(2);

			// Set score on weights map
			weightsMap[residue] = score;
		}
	}
}

// calculateSpliceSiteScoreCounts
//  Purpose:
//		Calculate the splice site scores and populate the
//		calculateSpliceSiteScoreCounts map with the values.
//  Postconditions:
//		spliceSiteScoreCounts will be populated
void GenbankFile::calculateSpliceSiteScoreCounts() {

	// Initialize counts
	for (int i = -20; i <= 20; i++) {
		spliceSiteScoreCounts[i] = 0;
	}
	
	// Iterate the true splice sites
	for (GenbankCodingSequence codingSequence : codingSequences) {

		// Get the location on the sequence (or complement) of the splice sites 
		map<int, int>& spliceSiteLocations = codingSequence.spliceSiteLocations();

		// Iterate through the splice sites and calculate score
		double score = 0;
		for (int i = -10; i <= 10; i++) {

			// Find the residue that is at the site
			char residue;
			if (!codingSequence.isComplement())
				residue = sequence.at(spliceSiteLocations.at(i) -1);
			else
				residue = complement.at(spliceSiteLocations.at(i) -1);

			// Calculate the score (only for known residues)
			if (residue == 'a' || residue == 'c' || residue == 'g' || residue == 't') 
				score += spliceSiteWeights[i][residue];
		}

		int bin = floor(score);

		spliceSiteScoreCounts[bin]++;

	}
}

// calculateSiteScoreCounts
//  Purpose:
//		Calculate the site scores and populate the
//		siteScoreCounts map with the values.
//  Postconditions:
//		spliceSiteScoreCounts - will be populated
//		highNonCodingSiteScores - will be populated with scores
//							      that are >= 10.0
void GenbankFile::calculateSiteScoreCounts() {

	// Initialize counts
	for (int i = -20; i <= 20; i++) {
		siteScoreCounts[i] = 0;
	}

	// Create lookup for coding sequence starts
	set<int> forwardStartLocs;
	set<int> complementStartLocs;
	for (GenbankCodingSequence codingSequence : codingSequences) {
		if (codingSequence.isComplement())
			complementStartLocs.insert(codingSequence.start());
		else
			forwardStartLocs.insert(codingSequence.start());
	}

	for (int seqLoc = 10; seqLoc < (int)sequence.length() - 10; seqLoc++) {

		// Forward Strand
		double score = 0;
		for (int offset = -10; offset <= 10; offset++) {

			// Find the residue that is at the site
			char residue = sequence.at(seqLoc + offset);

			// Calculate the score (only for known residues)
			if (residue == 'a' || residue == 'c' || residue == 'g' || residue == 't') 
				score += spliceSiteWeights[offset][residue];
		}

		// Add to counts
		int bin = floor(score);
		siteScoreCounts[bin]++;

		// Add to highNonCodingSiteScores if > 10.0
		if (score >= 10.0) {
			if (forwardStartLocs.find(seqLoc) == forwardStartLocs.end()) {
				SiteScore siteScore;
				siteScore.site = seqLoc + 1;
				siteScore.strand = 0;
				siteScore.score = score;
				highNonCodingSiteScores.push_back(siteScore);
			}
		}

		// Complement
		score = 0;
		for (int offset = -10; offset <= 10; offset++) {

			// Find the residue that is at the site
			char residue = complement.at(seqLoc - offset);

			// Calculate the score (only for known residues)
			if (residue == 'a' || residue == 'c' || residue == 'g' || residue == 't') 
				score += spliceSiteWeights[offset][residue];
		}

		// Add to counts
		bin = floor(score);
		siteScoreCounts[bin]++;

		// Add to highNonCodingSiteScores if > 10.0
		if (score >= 10.0) {
			if (complementStartLocs.find(seqLoc) == complementStartLocs.end()) {
				SiteScore siteScore;
				siteScore.site = seqLoc + 1;
				siteScore.strand = 1;
				siteScore.score = score;
				highNonCodingSiteScores.push_back(siteScore);
			}
		}
	}

}

// bool isStartOfCodingSequence(int aLocation)
//  Purpose:
//		Returns true if aLocation is the start on one of the known
//		coding sequences for this genbank file
bool GenbankFile::isStartOfCodingSequence(int aLocation, bool isComplement) {
	bool result = false;

	for (GenbankCodingSequence codingSequence : codingSequences) {
		if (isComplement && codingSequence.isComplement() && (aLocation == codingSequence.start())) {
			result = true;
			break;
		}

		if (!isComplement && !codingSequence.isComplement() && (aLocation == codingSequence.start())) {
			result = true;
			break;
		}
	}

	return result;
}

