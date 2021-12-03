/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#ifndef _KMER_FILE_HEADER_H
#define _KMER_FILE_HEADER_H
#include "defs.h"
#include <string>
#include <iostream>

//************************************************************************************************************
// CKmerFileHeader - represents header of k-mer database.
//************************************************************************************************************

enum class KmerFileType { KMC1, KMC2 };

struct CKmerFileHeader
{
public:
	uint32_t kmer_len = 0;
	uint32_t mode = 0;
	uint32_t counter_size = 0;
	uint32_t lut_prefix_len = 0;
	uint32_t signature_len = 0; //only for kmc2
	uint32_t min_count = 0;
	uint64_t max_count = 0;
	uint64_t total_kmers = 0;
	bool both_strands = true;
	uint32_t db_version = 0;
	uint32_t header_offset = 0;
	
	uint32_t no_of_bins = 0; //only for kmc2
	KmerFileType kmer_file_type;	
	KmerFileType GetType() const
	{
		return kmer_file_type;
	}
	
	CKmerFileHeader(std::string file_name);
	

private:
	template<typename T> void load_uint(FILE* file, T& res)
	{
		res = 0;
		for (uint32_t i = 0; i < sizeof(T); ++i)
			res += (T)getc(file) << (i << 3);
	}
};

#endif


// ***** EOF