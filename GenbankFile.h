/*
 * GenbankFile.h
 *
 *	This is the header file for the GenbankFile object. The GenbankFile object
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

#ifndef GenbankFile_H
#define GenbankFile_H

#include "GenbankCodingSequence.h"
#include <string>
#include <vector>
using namespace std;

class GenbankFile {

public:

	// Constuctors
	// ==============================================
	GenbankFile();
	GenbankFile(string fileName);  

	// Destructor
	// =============================================
	virtual ~GenbankFile();

	// Public Methods
	// =============================================

	// string firstLineResultString()
	//  Purpose:
	//		Returns the string value of an XML element representing the first line of 
	//		the Fasta file.
	//
	//		format:
	//			<result type='first line' file='<<fileName>>' >
	//				<<firstLine>>
	//			</result>
	//  Preconditions:
	//		Fasta File has been read and firstLine has been populated
	string firstLineResultString();

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
	string backgroundCountsResultString();

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
	string backgroundFrequenciesResultString();

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
	string spliceSiteCountsResultString();

	// string spliceSiteFrequenciesResultString()
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
	string spliceSiteFrequenciesResultString();

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
	string spliceSiteWeightsResultString();

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
	string spliceSiteScoreCountsResultString();

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
	string siteScoreCountsResultString();

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
	string highNonCodingSiteScoresResultString();

	// Public Accessors
	// =============================================
	const int getSequenceLength();  // length of dnaSequence
	string& getFileName();
	string& getSequence();
	vector<GenbankCodingSequence>& getCodingSequences();

private:
	// Attributes
	// =============================================
	string filePath;
	string fileName;
	string firstLine;
	string sequence;
	string complement;
	vector<GenbankCodingSequence> codingSequences;
	map<char, int> backgroundCounts;
	map<char, double> backgroundFrequencies;
	map<int, map<char, int>> spliceSiteCounts;
	map<int, map<char, double>> spliceSiteFrequencies;
	map<int, map<char, double>> spliceSiteWeights;
	map<int, int> spliceSiteScoreCounts;
	map<int, int> siteScoreCounts;

	struct SiteScore {
		int site;
		int strand;
		double score;
	};

	vector<SiteScore> highNonCodingSiteScores;


	// Private Methods
	// =============================================

	// parseFileName(string& aFileName)
	//  Purpose:
	//		Parse aFileName into the filePath and fileName attributes
	//  Postconditions:
	//		fileName and filePath have been set
	void parseFileName(string& aFileName);

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
	void populate();

	// createComplment()
	//  Purpose:
	//		populates the complement attribute with the comlpement
	//		of the sequence
	//	Preconditions:
	//		dnaSequence has been set
	//  Postconditions:
	//		complement - populated with complement of sequence
	void createComplement();

	// char complementOf(char aChar)
	//  Purpose:  returns the dna complement of aChar
	char complementOf(char aChar);

	// calculateBackgroundCounts()
	//  Purpose:
	//		populates the background counts map with the counts for residues
	//		in the sequence (and its complement)
	void calculateBackgroundCounts();

	// calculateBackgroundCounts(string& aSequence)
	//  Purpose:
	//		populates the background counts map with the counts for residues
	//		in aSequence
	void calculateBackgroundCounts(string& aSequence);

	// calculateBackgroundFrequencies()
	//  Purpose:
	//		populates the background frequencies map with the frequencies
	//		of residues in the sequence (and its complement)
	void calculateBackgroundFrequencies();

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
	void calculateSpliceSiteCounts();

	// calculateSpliceSiteFrequenciess
	//  Purpose:
	//		Calculate the splice site frequencies and populate the
	//		spliceSiteFrequencies map with the values.
	//  Postconditions:
	//		spliceSiteFrequenciess will be populated
	void calculateSpliceSiteFrequencies();

	// calculateSpliceSiteWeights
	//  Purpose:
	//		Calculate the splice site weights and populate the
	//		spliceSiteWieghts map with the values.
	//  Postconditions:
	//		spliceSiteWeights will be populated
	void calculateSpliceSiteWeights();

	// calculateSpliceSiteScoreCounts
	//  Purpose:
	//		Calculate the splice site scores and populate the
	//		spliceSiteScoreCounts map with the values.
	//  Postconditions:
	//		spliceSiteScoreCounts will be populated
	void calculateSpliceSiteScoreCounts();

	// calculateSiteScoreCounts
	//  Purpose:
	//		Calculate the site scores and populate the
	//		siteScoreCounts map with the values.
	//  Postconditions:
	//		spliceSiteScoreCounts - will be populated
	//		highNonCodingSiteScores - will be populated with scores
	//							      that are >= 10.0
	void calculateSiteScoreCounts();
	
};

#endif /* GenbankFile_H */
