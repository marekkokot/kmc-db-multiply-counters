#pragma once

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#endif
