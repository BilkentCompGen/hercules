all:
	g++ -I /usr/local/include/ -std=c++14 -O3 -DNDEBUG -W -Wall -o ../bin/hercules main.cpp -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -pedantic -fopenmp -lpthread -lrt -lz -lbz2
