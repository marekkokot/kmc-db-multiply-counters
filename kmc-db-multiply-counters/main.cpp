#include "kmer_file_header.h"
#include <string>
#include <iostream>
#include <limits>

uint64_t read_counter(FILE* file, uint64_t counter_size)
{
	uint64_t res{};
	char data[8];
	fread(data, 1, counter_size, file);
	for (uint64_t i = 0; i < counter_size; ++i)
	{
		res <<= 8;		
		res += data[counter_size - 1 - i];
	}
	return res;
}

void write_counter(FILE* file, uint64_t c, uint64_t counter_size)
{
	char data[8];
	for (uint64_t i = 0; i < counter_size; ++i)
	{
		data[i] = c >> (8 * i);		
	}
	fwrite(data, 1, counter_size, file);
}

uint64_t get_max_value_to_store_on_n_bytes(uint64_t n)
{
	uint64_t res = (1ull << 8 * n) - 1;
	return res;
}

int multiply_counters(const std::string& kmcDbPath, uint64_t total_kmers, uint64_t suffix_len, uint64_t counter_size, uint64_t min_count, uint64_t max_count, double multiplier)
{
	auto _kmcDbPath = kmcDbPath + ".kmc_suf";
	auto file = my_fopen(_kmcDbPath.c_str(), "r+b");
	if (!file)
	{
		std::cerr << "Error: cannot open file " << _kmcDbPath << "\n";
		return 1;
	}
	my_fseek(file, 0, SEEK_SET);
	std::string marker(4, 'X');
	fread((void*)marker.data(), 1, 4, file);
	if (marker != "KMCS")
	{		
		std::cerr << "Error: wrong start marker\n";
		return 1;
	}

	my_fseek(file, -4, SEEK_END);	
	fread((void*)marker.data(), 1, 4, file);
	if (marker != "KMCS")
	{
		std::cerr << "Error: wrong end marker\n";
		return 1;
	}

	my_fseek(file, 4, SEEK_SET); //skip start marker

	auto suffix_bytes = suffix_len / 4;
	
	uint64_t lowest_count = std::numeric_limits<uint64_t>::max();
	uint64_t highest_count = 0;

	uint64_t max_legal = get_max_value_to_store_on_n_bytes(counter_size);
	uint64_t n_corrected{};
	for (uint64_t i = 0; i < total_kmers; ++i)
	{
		my_fseek(file, suffix_bytes, SEEK_CUR); //shift to counter
		auto c = read_counter(file, counter_size);
		c *= multiplier;

		if (c > highest_count)
			highest_count = c;
		if (c < lowest_count)
			lowest_count = c;

		if (c > max_legal)
		{
			c = max_legal;
			++n_corrected;
		}

		my_fseek(file, -(int64_t)counter_size, SEEK_CUR); //shift back to set counter
		write_counter(file, c, counter_size);
	}

	//have we reached the end marker?
	fread((void*)marker.data(), 1, 4, file);
	if (marker != "KMCS")
	{
		std::cerr << "Error: something went wrong\n";
		return 1;
	}

	auto pos = my_ftell(file);
	my_fseek(file, 0, SEEK_END);
	//is this the end of file?
	if (pos != my_ftell(file))
	{
		std::cerr << "Error: something went wrong (EOF not reached)\n";
		return 1;
	}

	fclose(file);

	return 0;
}

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << "<kmcDbPath> <multiplier>\n";
		std::cerr << "Examlple: " << argv[0] << "A 5.5\n";
	}
	std::string kmcDbPath = argv[1];
	double multiplier = atof(argv[2]);

	CKmerFileHeader header(kmcDbPath);
	
	uint64_t suffix_len = header.kmer_len - header.lut_prefix_len;

	std::cerr << "kmer_len: " << header.kmer_len << "\n";
	std::cerr << "total_kmers: " << header.total_kmers << "\n";
	std::cerr << "suffix_len: " << suffix_len << "\n";
	std::cerr << "lut_prefix_len: " << header.lut_prefix_len << "\n";
	std::cerr << "counter_size: " << header.counter_size << "\n";
	std::cerr << "min_count: " << header.min_count << "\n";
	std::cerr << "max_count: " << header.max_count << "\n";
	std::cerr << "multiplier: " << multiplier << "\n";

	return multiply_counters(kmcDbPath, header.total_kmers, suffix_len, header.counter_size, header.min_count, header.max_count, multiplier);
}
