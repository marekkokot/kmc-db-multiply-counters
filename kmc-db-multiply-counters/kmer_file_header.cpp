/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Marek Kokot
  
  Version: 3.1.1
  Date   : 2019-05-19
*/

#include "kmer_file_header.h"
#include <cstring>
#include <set>

/*****************************************************************************************************************************/
/******************************************************** CONSTRUCTOR ********************************************************/
/*****************************************************************************************************************************/

CKmerFileHeader::CKmerFileHeader(std::string file_name)
{
	file_name += ".kmc_pre";
	FILE* file = my_fopen(file_name.c_str(), "rb");
	if (!file)
	{
		std::cerr << "Error: Cannot open file " << file_name << "\n";
		exit(1);
	}
	char marker[4];
	if (fread(marker, 1, 4, file) != 4)
	{
		std::cerr << "Error while reading start marker in " << file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCP", 4) != 0)
	{
		std::cerr << "Error: wrong start marker in " << file_name << "\n";
		exit(1);
	}

	my_fseek(file, -4, SEEK_END);
	if (fread(marker, 1, 4, file) != 4)
	{
		std::cerr << "Error while reading end marker in " << file_name << "\n";
		exit(1);
	}

	if (strncmp(marker, "KMCP", 4) != 0)
	{
		std::cerr << "Error: wrong end marker in " << file_name << "\n";
		exit(1);
	}

	my_fseek(file, 0, SEEK_END);
	uint64_t file_size = my_ftell(file);

	my_fseek(file, -8, SEEK_END);
	load_uint(file, header_offset);

	my_fseek(file, -12, SEEK_END);
	load_uint(file, db_version);

	kmer_file_type = db_version == 0x200 ? KmerFileType::KMC2 : KmerFileType::KMC1;

	my_fseek(file, 0LL - (header_offset + 8), SEEK_END);
	load_uint(file, kmer_len);
	load_uint(file, mode);
	load_uint(file, counter_size);
	load_uint(file, lut_prefix_len);
	if (kmer_file_type == KmerFileType::KMC2)
		load_uint(file, signature_len);
	load_uint(file, min_count);
	uint32_t max_count_lo;
	load_uint(file, max_count_lo);
	load_uint(file, total_kmers);
	uint8_t both_s_tmp;
	load_uint(file, both_s_tmp);
	both_strands = both_s_tmp == 1;
	both_strands = !both_strands;

	fseek(file, 3, SEEK_CUR);
	uint32_t max_count_hi;
	load_uint(file, max_count_hi);
	max_count = (((uint64_t)max_count_hi) << 32) + max_count_lo;
	fclose(file);

	if (kmer_file_type == KmerFileType::KMC2)
	{
		uint32_t single_lut_size = (1ull << (2 * lut_prefix_len)) * sizeof(uint64_t);
		uint32_t map_size = ((1 << 2 * signature_len) + 1) * sizeof(uint32_t);
		no_of_bins = (uint32_t)((file_size - sizeof(uint64_t) - 12 - header_offset - map_size) / single_lut_size);
	}
}

// ***** EOF