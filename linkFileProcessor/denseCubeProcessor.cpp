// Copyright 2024 NVIDIA Corporation.
// For Copyright and Licensing please refer to the LICENSE and
// THIRD_PARTY_LICENSES file in the top level directory of this package


#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
// Needed for easy determination of the 
// Filesize 
#include <sys/types.h>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include <limits>
#include <array>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

#include "raw_data_reader.hpp"
#include "gmy.h"
#include <rpc/rpc.h>
#include <zlib.h>

// This will be my map of blocks to file offsets. 

constexpr uint64_t blockDim = 8;
constexpr uint64_t blockSites = (blockDim*blockDim*blockDim);



// The block map is a map from block id to to a vector
// the vector stores pairs of (site_id, offset in local array)

uint64_t nBlocksX = 0;
uint64_t nBlocksY = 0;
uint64_t nBlocksZ = 0;
int this_rank=0;
int num_ranks=0;
std::vector<Site> rearranged_data;

struct ConvertedBlockInfo { 
	std::array<Site*,blockSites> sites;	
	NonEmptyHeaderRecord header;
	ConvertedBlockInfo() {
		for(int i=0; i < blockSites; i++) sites[i]=nullptr;
	} 
};

std::map<size_t, ConvertedBlockInfo> output_info;
std::vector<char> outbuf;		// Output buffer

void generateData(uint64_t cube_dimension, const std::string& output_filename)
{

	//Here we generate the sites for a N^3 cube of blocks for large domain scaling purposes
	//
	
	/*
	 * Step 1: Work out size of the Block Grid
	 */
	MPI_Comm_rank(MPI_COMM_WORLD,&this_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
	
	// This rank is responsible for the following blocks
	uint64_t sizeDiv = (cube_dimension * cube_dimension * cube_dimension) / num_ranks;
	uint64_t sizeRem = (cube_dimension * cube_dimension * cube_dimension) % num_ranks;

	uint64_t localBlockMin, localBlockMax, numlocalBlocks, numlocalSites;

	if (this_rank < sizeRem) {
		localBlockMin = this_rank * (sizeDiv + 1);
		localBlockMax = (this_rank + 1) * (sizeDiv + 1) - 1;
	} else {
		localBlockMin = sizeRem * (sizeDiv+1) + (this_rank - sizeRem)*sizeDiv;
		localBlockMax = sizeRem * (sizeDiv+1) + (this_rank + 1 - sizeRem)*sizeDiv - 1;
	}
	numlocalBlocks = localBlockMax - localBlockMin +1;

	//std::cout << "DEBUG: rank " << this_rank << " (sD " << sizeDiv << " sR " << sizeRem << ") has " << numlocalBlocks << " blocks with blockMin " << localBlockMin << ", blockMax " << localBlockMax << std::endl;

	nBlocksX = cube_dimension;
	nBlocksY = cube_dimension;
	nBlocksZ = cube_dimension;
	
	numlocalSites = numlocalBlocks * blockSites;
	
	// Now we want to set up buffers for output
	size_t outbuf_idx = 0;

	// I don't understand this magic number from mpivx2gmy but 
	// if I replace it with what I think it should be I seem to have weird problems.
	//	
	uint64_t max_site_size = (2 + 26 * 4 * sizeof(uint64_t));

	uint64_t max_buffer_size = blockSites * max_site_size;
	
	std::vector<unsigned char> outputBuffer( numlocalSites * max_site_size );
	

	double generate_starttime = MPI_Wtime();
	
	for (int b=localBlockMin; b < localBlockMax+1; b++) {

		uint64_t block_id = b;
		uint64_t blockX = block_id / (cube_dimension * cube_dimension);
		uint64_t blockY = (block_id / cube_dimension) % cube_dimension;
		//uint64_t blockY = (block_id % (cube_dimension * cube_dimension)) / cube_dimension;
		uint64_t blockZ = block_id % cube_dimension;	
		//uint64_t blockZ = (block_id % (cube_dimension * cube_dimension)) % cube_dimension;	
		
		ConvertedBlockInfo binfo;
		binfo.header.sites = blockSites;
		binfo.header.blockNumber = b;
		binfo.header.fileOffset = outbuf_idx; // we will need to add the length of the header block to this
		
		std::array<OutputSite*, blockSites> conv_sites;        // The converted sites

		// For now we set all the converted sites to bet the null pointer
		for(int i=0; i < blockSites; i++) conv_sites[i]=nullptr;

		// We will need to store the temporary converted sites
		std::vector<OutputSite> converted;
		
		std::vector<unsigned char> decompressedBuffer(max_buffer_size);
		std::vector<unsigned char> compressedBuffer(max_buffer_size);

		uint32_t blockUncompressedLen = 0;

		// First set up a map of converted blockinfo 
		for(int i=0; i < blockSites; i++) {	
			uint64_t site_id = i; // Within a block, sID = siteZ + blockDim * (siteY + blockDim * siteX)
			uint64_t siteX = site_id / (blockDim * blockDim); //local coords, modify these based on BJ slack 7/5
			//uint64_t siteY = (site_id % (blockDim * blockDim)) / blockDim;
			uint64_t siteY = (site_id / blockDim) % blockDim;
			//uint64_t siteZ = (site_id % (blockDim * blockDim)) % blockDim;
			uint64_t siteZ = site_id % blockDim;
			uint64_t x = 8 * blockX + siteX; //global coords
			uint64_t y = 8 * blockY + siteY;
			uint64_t z = 8 * blockZ + siteZ;
	
			OutputSite newsite; 
			newsite.hasWallNormal = false;
			newsite.normalX = 0.0;
			newsite.normalY = 0.0;
			newsite.normalZ = 0.0;
			uint32_t num_intersections = 0;

			newsite.x = x;
			newsite.y = y;
			newsite.z = z;
			uint64_t linkType, configID;
			std::vector<int> linkList;

			for(int linkdir = 0; linkdir < 26; linkdir++) {	
				
				float wallDistance=0.5, normalX=0.0, normalY=0.0, normalZ=0.0;

				linkType = 0;
				configID = 0;

				if (x==0 and linkdir<9) {
					linkType = 1;
					normalX = -1.0;
				}
				if (x==8*cube_dimension-1 and linkdir>16) {
					linkType = 1;
					normalX = 1.0;
				}
				linkList = {0,1,2,9,10,11,17,18,19};
				if (y==0 and (std::find(std::begin(linkList), std::end(linkList), linkdir) != std::end(linkList)) ) {
					linkType = 1;
					normalY = -1.0;
				}
			
				linkList = {6,7,8,14,15,16,23,24,25};
				if (y==8*cube_dimension-1 and (std::find(std::begin(linkList), std::end(linkList), linkdir) != std::end(linkList)) ) {
					linkType = 1;
					normalY = 1.0;
				}
			
				linkList = {0,3,6,9,12,14,17,20,23};
				if (z==0 and (std::find(std::begin(linkList), std::end(linkList), linkdir) != std::end(linkList)) ) {
					linkType = 2;
					normalZ = -1.0;
				}
			
				linkList = {2,5,8,11,13,16,19,22,25};
				if (z==8*cube_dimension-1 and (std::find(std::begin(linkList), std::end(linkList), linkdir) != std::end(linkList)) ) {
					linkType = 3;
					normalZ = 1.0;
				}

				newsite.links[linkdir].linkType = (uint32_t)linkType;
				newsite.links[linkdir].configID = (uint32_t) configID;
		
				newsite.links[linkdir].wallDistance = wallDistance;

				// if link is to a wall...
				if(newsite.links[linkdir].linkType == 1) {
					newsite.hasWallNormal = true;
					newsite.normalX += normalX;
					newsite.normalY += normalY;
					newsite.normalZ += normalZ;
					num_intersections++;
				}
			}

			if ( newsite.hasWallNormal ) {
				newsite.normalX /= (float) num_intersections;
				newsite.normalY /= (float) num_intersections;
				newsite.normalZ /= (float) num_intersections;
			}
               
			// do something with the newsite that has been created
			converted.push_back(newsite);
			conv_sites[ i ] = &converted[ converted.size()-1 ];	
		
			// At this point we can encode the conv_sites to XDR
			OutputSite* siteptr = conv_sites[i]; 	        	
		
			if( siteptr == nullptr ) { // if solid 
				blockUncompressedLen += sizeof(uint32_t); // SiteIsSimulated
			}
			else { // if fluid
				blockUncompressedLen += sizeof(uint32_t); // SiteIsSimulated
				for(uint32_t m = 0; m < 26; m++) {
					blockUncompressedLen += sizeof(uint32_t); // linkType
					uint32_t linkType = (*siteptr).links[m].linkType;
	
					switch ( linkType ) {
						case 0: // linkType = FLUID (no further data)
							break;
						case 1: // linkType = WALL (write distance to nearest obstacle)
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						case 2:
							blockUncompressedLen += sizeof(uint32_t); // configEntry
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						case 3: // linkType = INLET or OUTLET (write config ID and distance to nearest obstacle)
							blockUncompressedLen += sizeof(uint32_t); // configEntry
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						default:
							fprintf(stderr, "ERROR: Unrecognised linkType %u on line %d\n", linkType, __LINE__);
							MPI_Finalize();
							exit(1);
					}
				}
				blockUncompressedLen += sizeof(uint32_t); // hasWallNormal
				
				if ( siteptr->hasWallNormal ) {	
					blockUncompressedLen += sizeof(float); // normalX
					blockUncompressedLen += sizeof(float); // normalY
					blockUncompressedLen += sizeof(float); // normalZ
				}
			}

		}
	
		if( blockUncompressedLen > max_buffer_size ) { 
			std::cout << "Rank " << this_rank << " : blockUncompressedLen= " << blockUncompressedLen  << " < max_buffer_size = " 
			<< max_buffer_size <<"\n";
		}

		binfo.header.uncompressedBytes = blockUncompressedLen;

		XDR xdrbs;
		xdrmem_create(&xdrbs, (char *)decompressedBuffer.data(), blockUncompressedLen, XDR_ENCODE);
	        
		// Encode the block
		for(int i=0; i < blockSites; i++) {
			OutputSite* siteptr = conv_sites[i];
			uint32_t siteIsSimulated = ( siteptr != nullptr ) ? 1 : 0;
			xdr_u_int(&xdrbs, &siteIsSimulated);

			if( siteIsSimulated == 0 ) continue;

			for( uint32_t link=0; link < 26; link++) {

				// write type of link
				uint32_t linkType = converted[i].links[link].linkType;
				//uint32_t linkType = (*siteptr).links[link].linkType;
				xdr_u_int(&xdrbs, &linkType);
				switch (linkType) {
					case 0: // linkType = FLUID (no further data)
						break;
					case 1: // linkType = WALL (write distance to nearest obstacle)
						xdr_float(&xdrbs, &(converted[i].links[link].wallDistance));
						//xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					case 2: // linkType = INLET (write inletID and distance to nearest obstacle
						xdr_u_int(&xdrbs, &(converted[i].links[link].configID));
						xdr_float(&xdrbs, &(converted[i].links[link].wallDistance));
						//xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
						//xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					case 3: // linkType = OUTLET (write outletID and distance to nearest obstacle
						xdr_u_int(&xdrbs, &(converted[i].links[link].configID));
						xdr_float(&xdrbs, &(converted[i].links[link].wallDistance));
						//xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
						//xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					default:
						fprintf(stderr, "ERROR: Unrecognised linkType %u on line %d .\n", linkType, __LINE__);
						MPI_Finalize();
						exit(1);
				}
			}

				// state if there are wall normal coordinates to be read (1 for yes)
			uint32_t hasWallNormal = (converted[i].hasWallNormal == true)? 1: 0;
			//uint32_t hasWallNormal = (siteptr->hasWallNormal == true)? 1: 0;
			xdr_u_int(&xdrbs, &hasWallNormal);
			if (hasWallNormal == 1) {
				// write wall normal coordinates as separate floats
				xdr_float(&xdrbs, &(converted[i].normalX));
				xdr_float(&xdrbs, &(converted[i].normalY));
				xdr_float(&xdrbs, &(converted[i].normalZ));
				//xdr_float(&xdrbs, &(siteptr->normalX));
				//xdr_float(&xdrbs, &(siteptr->normalY));
				//xdr_float(&xdrbs, &(siteptr->normalZ));
			}

		}

		xdr_destroy(&xdrbs);
		// Now we compress
		z_stream strm;

		// setup zlib for compression
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;

		// input
		strm.avail_in = blockUncompressedLen;
		strm.next_in = decompressedBuffer.data();

		// output
		strm.avail_out = max_buffer_size;
		strm.next_out =compressedBuffer.data();

		uint32_t ret;
		ret = deflateInit(&strm, Z_BEST_COMPRESSION);
		if(ret != Z_OK) {
			fprintf(stderr, "ERROR: zlib deflation init.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflate(&strm, Z_FINISH);
		if (ret != Z_STREAM_END) {
			fprintf(stderr, "ERROR: Deflation error for block.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflateEnd(&strm);
		if (ret != Z_OK) {
			fprintf(stderr, "ERROR: Deflation end error for block.\n");
			MPI_Finalize();
			exit(1);
		}

		// get new compressed size
		uint32_t blockCompressedLen = (unsigned char*)strm.next_out - compressedBuffer.data();
		binfo.header.bytes = blockCompressedLen;
		memcpy(&outputBuffer[outbuf_idx], compressedBuffer.data(), blockCompressedLen);
		outbuf_idx += blockCompressedLen;
		
		output_info.insert( std::pair<uint64_t, ConvertedBlockInfo >( block_id, binfo) );

	} 
	
	double generate_endtime=MPI_Wtime(); 
	if (this_rank == 0 ) std::cout << "Generation and compression completed: " << generate_endtime-generate_starttime << " sec.\n";


	if (this_rank == 0 ) std::cout << "Gathering information for parallel I/O\n";

	uint64_t rank_compressed_bytes = outbuf_idx;  // This is the total number of compressed data on the local rank

	// Get the total compressed bytes in the file
	uint64_t total_compressed_bytes = 0;

	// Get the offsets of the compressed data in the file
	uint64_t compressed_offset = 0;

	// At this point we can get a total file size for the data portion.
	MPI_Allreduce( &rank_compressed_bytes, &total_compressed_bytes, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);


	// We can also get the offsets of the data (from the end of the headers) for each rank by performaing an exclusive prefix scan.
	MPI_Exscan(&rank_compressed_bytes, &compressed_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	// To fill out the header info we need the maximum compressed and uncompressed bytes as well as total number of non-empty blocks
	// We may have these already but they are easy to find by traversing out map.
	uint32_t global_max_uncompressed_bytes = 0;
	uint32_t global_max_compressed_bytes = 0;
	uint64_t global_uncompressed_bytes = 0;
	uint64_t global_compressed_bytes = 0;
	{
		// We can compute the maxes of both the uncompressed and copressed in a single loop
		// then we can do a length 2 MPI Allreduce 
		uint32_t local_maxes[2] = { 0, 0 };
		uint64_t local_bytes[2] = { 0, 0 };
		for( const auto& kv : output_info ) {
			const auto& header = kv.second.header;
			local_bytes[0] += (uint64_t)(header.bytes);
			local_bytes[1] += (uint64_t)(header.uncompressedBytes);
			if( header.bytes > local_maxes[0] ) local_maxes[0] = header.bytes;
			if( header.uncompressedBytes > local_maxes[1] ) local_maxes[1] = header.uncompressedBytes;
		} 
		uint32_t global_maxes[2] = {0,0}; 
		uint64_t global_bytes[2] = {(uint64_t)0,(uint64_t)0};
		MPI_Allreduce(local_maxes, global_maxes, 2, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
		global_max_compressed_bytes = global_maxes[0];
		global_max_uncompressed_bytes = global_maxes[1];
		MPI_Allreduce(local_bytes, global_bytes, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
		global_compressed_bytes = global_bytes[0];
		global_uncompressed_bytes = global_bytes[1];
	}

	// The number of non empty blocks locally is just the number of elements in output info
	uint64_t global_nonempty_blocks=0;
	{	
		uint64_t local_nonempty_blocks = output_info.size();
		MPI_Allreduce(&local_nonempty_blocks, &global_nonempty_blocks, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	}

	if ( this_rank == 0 ) { 
		std::cout << "Global Nonempty Blocks = " << global_nonempty_blocks << "\n";
		std::cout << "(INFO) Total Sites = " << global_nonempty_blocks *512 << "\n";
		std::cout << "Total Uncompressed bytes = " << global_uncompressed_bytes << " = "
				  << (double)global_uncompressed_bytes/(double)(1024*1024) << " MiB\n";

		std::cout << "Total Compressed bytes = " << total_compressed_bytes << " = "
				  << (double)total_compressed_bytes/(double)(1024*1024) << " MiB\n"; 
		std::cout << "Total Compressed bytes 2 = " << global_compressed_bytes << " = "
				  << (double)global_compressed_bytes/(double)(1024*1024) << " MiB\n";
		std::cout << "Max Uncomressed bytes (per block) = " << global_max_uncompressed_bytes << " = " 
				  << (double)global_max_uncompressed_bytes / (double)(1024) << " KiB\n";

		std::cout << "Max Compressed bytes (per block) = " << global_max_compressed_bytes << " = " 
				  << (double)global_max_compressed_bytes / (double)(1024) << " KiB\n";

		std::cout << "Preparing Preamble and headers\n";
	}

	// Now fill out the Preamble Info
	OutputPreambleInfo pinfo;
	pinfo.HemeLBMagic = HemeLbMagicNumber;         // From gmy.h
	pinfo.GmyNativeMagic = GmyNativeMagicNumber;   // From gmy.h  
	pinfo.Version = GmyNativeVersionNumber;        // From gmy.h
	pinfo.BlocksX = nBlocksX;					   // From earlier calculations
	pinfo.BlocksY = nBlocksY;
	pinfo.BlocksZ = nBlocksZ;
	pinfo.BlockSize = blockDim;				       // We set this as constexpr earlier in the file
	pinfo.MaxCompressedBytes = global_max_compressed_bytes;			// We computed this (these are per-block quantities so uint32_t is safe)
	pinfo.MaxUncompressedBytes = global_max_uncompressed_bytes;     // We computed this  (these are per-block quantities so uint32_t is safe)
	// Header offset is 56 -- default on construction
	pinfo.NonEmptyBlocks = global_nonempty_blocks;					// We computed this 

	// This is the offset to the data: HeaderOffset is fixed at 56. NonemptyHeaderRecordSize is 28. 
	pinfo.DataOffset = pinfo.HeaderOffset + pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize;
	if( this_rank == 0 ) { 
		std::cout << "Header Size = " << pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize << " bytes = " 
				   << (double)pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize/(double)(1024*1024) << " MiB\n";
	}	
	// Now the headers
	// We walk the map and unroll it into a vector for writing
	std::vector<NonEmptyHeaderRecord> local_header( output_info.size() );
	uint64_t header_idx = 0;
	for( const auto& kv : output_info )  {
		local_header[header_idx] = kv.second.header;
		local_header[header_idx].fileOffset += compressed_offset;

#if 0
		std::cout << "Rank " << this_rank << " : idx = "<< header_idx 
				     << " block = " << local_header[header_idx].blockNumber
		       	  	     << " sites = " << local_header[header_idx].sites
				     << " comp_bytes = " << local_header[header_idx].bytes 
				     << " uncomp_bytes = " << local_header[header_idx].uncompressedBytes 
				     << " offset = " << local_header[header_idx].fileOffset << "\n"; 
#endif
		header_idx++;
	}
	
	// Each rank needs to offset its write to not stomp on the other nodes headers. 
	// We can work out the displacements using an exclusive scan
	uint64_t header_offset = 0;
	MPI_Exscan( &header_idx, &header_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD );
	// Header offset now holds where my rank should wrirte in units of NonEmptyHeaderRecords. 
	// Let us convert this to bytes
	header_offset *= NonemptyHeaderRecordSize;  // NB: The size of the struct is 32 from sizeof because of alignment we don't want to pad

	// And let us add on the Preamble offset. 
	header_offset += pinfo.HeaderOffset;

	// This is just to convert to an MPI_Offset (from uint64_t)
	MPI_Offset my_header_offset = header_offset; // Offset to start writing our local array

	// Now we will make an MPI type for our header: 6 elements (blocks): 2 uint64_ts, 4 uint32_ts
	MPI_Datatype header_type;
	int blen[6] = {1,1,1,1,1,1};  // Number of elements in each 'block'
	MPI_Aint  bdisp[6] = {0,8,16,20,24,28}; // Displacements of the elements

	MPI_Datatype btypes[6] = { MPI_UINT64_T, MPI_UINT64_T, MPI_UINT32_T, MPI_UINT32_T, MPI_UINT32_T, MPI_UINT32_T };
	MPI_Type_create_struct( 6, blen, bdisp, btypes, &header_type );
	MPI_Type_commit(&header_type);

	size_t chunk_size = 2*1024*1024; // 2 MiB transaction size
	

	// OK. We need to open an output file.
	if( this_rank == 0 ) std::cout << "Opening and writing SGMY file\n";

	MPI_File out_fh;
	MPI_File_open(MPI_COMM_WORLD, output_filename.c_str(), MPI_MODE_CREATE |MPI_MODE_WRONLY, MPI_INFO_NULL, &out_fh);
	double startio = MPI_Wtime();
	
	// Now for the IO part
	// Rank 0 writes the preamble
	if( this_rank == 0 ) {
		MPI_Status status;
		MPI_File_write_at(out_fh, 0, &pinfo, pinfo.HeaderOffset, MPI_BYTE, &status);
	}

	// Everyone writes the headers using collective write_at_all()
	// FIXME:  Check return and status
	

	MPI_Status header_write_status;
	MPI_File_seek(out_fh, my_header_offset, MPI_SEEK_SET);

	// Use the number of records from the scan
	MPI_File_write_all(out_fh, local_header.data(), header_idx, header_type, &header_write_status); 
	// Now write the data. We are going to write a stream of bytes, but we may have that we have over 2GB locally to write.
	// So unless we make a large type that is hard. However making a large type is also hard because each block may have 
	// a different compressed length. So we will use a transaction size (2 MiB) in this case and loop with that until 
	// we reach the end of the data for each rank.
	// Compressed offset is from the exclusive scan earlier;
	// Writing loop: we will start writing locally at &outputBuffer[bytes_written] and write bytes_to_write bytes
	// bytes_to_write will be our chunk size, or a mop-up amount at the end of the data which is less.
	// Then we increase 'bytes_written' and the start offset by the amount we just wrote.
	// Finally in principle each process could have different amounts of data, so I will use the "MPI_File_write_at_all" collective MPI I/O routine
	size_t bytes_written = 0;
 	MPI_Offset current_offset = pinfo.DataOffset+compressed_offset;
	while( bytes_written < rank_compressed_bytes ) {

		// Determine bytes to write
		size_t bytes_to_write =  ( bytes_written + chunk_size > rank_compressed_bytes ) ? rank_compressed_bytes - bytes_written : chunk_size; 

		// Do the write:
		// FIXME: Check status
		MPI_Status mpi_data_write_status;
		MPI_File_write_at(out_fh, current_offset, &outputBuffer[bytes_written], bytes_to_write, MPI_BYTE, &mpi_data_write_status);

		// Increase bytes_written (start location in our local buffer) 
		bytes_written += bytes_to_write;

		// Increase offset (start_location in the file)
		current_offset += bytes_to_write;
	}
	// We are done	
	double endio = MPI_Wtime();

	MPI_File_close(&out_fh);
	if( this_rank == 0 ) {
		double io_time = endio-startio;
		uint64_t size_in_bytes = pinfo.DataOffset + global_compressed_bytes;
		double size_in_GiB = (double)(size_in_bytes)/(double)(1024*1024*1024);
		std::cout << "Output IO took: " << endio - startio << " sec. Average BW: " << size_in_GiB/io_time << " GiB/sec. \n";
		std::cout << "Expected file size = " << size_in_bytes  << " bytes = " << size_in_GiB <<" GiB\n";
	}	
}

int main(int argc, char *argv[]) 
{
	if( argc != 3 ) {
	  std::cout << "Usage: denseCubeProcessor <cube dimension in blocks> < SGMY filename>\n";
	  return -1;
	}
	uint64_t cube_dimension = std::stoi(argv[1]);
	std::string filename_out(argv[2]);

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&this_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
	if( this_rank == 0) std::cout << "MPI Initialized: I am rank " << this_rank << " out of " << num_ranks << "\n";
	
	generateData(cube_dimension,filename_out);

	MPI_Finalize();
}
