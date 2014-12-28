/*
 * RMQ.hpp
 *
 *  Created on: Oct 12, 2012
 *      Author: alexei
 */

#ifndef RMQ_HPP_
#define RMQ_HPP_

#include<iostream>
#include<cmath>

#include<omp.h>

template < class DT = int, class DTidx = uint >
class RMQ : public RangeMinimumQuery
{


private:

	ushort* typeTable;
	// a table storing the type for each individual micro-block, for fast lookup

	uchar** inBlock;
	// a 'predecessor' like structure, for solving in-block-queries fast

	DTidx inBlock_size;

	/**
	 *	M[i][j] -> stores the position of the minimum in the subarray
	 *			   data[ j * block_size, (j + 2^i) * block_size - 1 ]
	 */

	DTidx** M;
	DTidx M_sz;

	/**
	 *	sM[i][j] -> stores the position of the minimum in the subarray
	 *			    data[ j * block_size, (j + 2^i) * block_size - 1 ]
	 */

	DTidx** sM;
	DTidx sM_sz;
	// actual size of the two matrixes

	static DTidx** compute_catalan_table( DTidx size )
	{
		DTidx** Catalan = (DTidx**) calloc( size, sizeof(DTidx*) );
		// allocate memory for the full table

		Catalan[0] = (DTidx*) calloc( 1 , sizeof(DTidx) );
		Catalan[0][0] = 1;
		// init

		for( DTidx i = 1; i < size; ++i )
		{
			Catalan[i] = (DTidx*) calloc ( i+1, sizeof(DTidx) );

			Catalan[i][0] = Catalan[i-1][0];

			for( DTidx j = 1; j < i; ++j ){
				Catalan[i][j] = Catalan[i][j-1] + Catalan[i-1][j];
			}

			Catalan[i][i] = Catalan[i][i-1];
		}

		return Catalan;
	}

	DT* data;

	DTidx size;

	DTidx super_block_size;
	DTidx block_size;
	DTidx micro_block_size;

	DTidx super_block_count;
	DTidx block_count;
	DTidx micro_block_count;

	void serial_compute_typeTable()
	{
		DTidx** Catalan = compute_catalan_table( micro_block_size+1 );
		// precompute Catalan for fast lookup

		/**
		 * Compute type for every micro_block as seen in Fischer/Heun CPM'06
		 */

		typeTable  = (ushort*) calloc( micro_block_count,    sizeof( ushort ) );
		DT* right_path = (DT*) calloc( micro_block_size + 1, sizeof(DT) );

		inBlock = (uchar**) calloc ( Catalan[ micro_block_size ][ micro_block_size ], sizeof(uchar*) );
		// allocate space for Alstrup et al. in-block preprocessing

		inBlock_size = Catalan[ micro_block_size ][ micro_block_size ];
		// save max size ( so we can free Catalan matrix before inBlock )

		char* Guardian = (char*) calloc ( inBlock_size, sizeof(char) );

		/*
		 * Install a guarding array ( 0 by default, idx if type was found )
		 * We want to avoid allocating memory / computing types for types that don't
		 * exist in the initial array
		 */

		DTidx* stack = ( DTidx* ) calloc ( micro_block_size + 1, sizeof(DTidx) );
		uchar stack_size;

		// when going parallel, allocate as many rp ( and stacks ) as nr of threads
		// and don't forget to free the memory

		// TODO - parallel

		for( DTidx t = 0; t < micro_block_count; ++t )
		{
			DTidx start = t * micro_block_size;

			DTidx it = start;
			// position the iterator on the t'th micro_block

			DTidx mblock_end = it + micro_block_size;
			// compute micro-block end

			mblock_end = mblock_end > size ? size : mblock_end;
			// check boundary condition; if last micro-block is smaller than 'micro_block_size'
			// consider the array it's filled with INF values

			DTidx rp_it = 0;
			// initially, there are no elements in the tree ( implicitly, none on the right path )

			DTidx q = micro_block_size;
			DTidx p = q - 1;
			// start from the bottom of the table

			DT dcopy;
			// store value pointed by iterator to avoid table lookup

			ushort type = 0;
			// store block type

			for( ; it < mblock_end; ++it, --p ) // iterate through block / go 'westwards'
			{
				dcopy = data[it];
				// copy data at with offset i

				// as long as the current element is smaller than the one on top of right_path stack
				while( rp_it && right_path[ rp_it ] > dcopy )
				{
					type += Catalan[q][p]; // skip C(q,p) paths
					--rp_it; --q;	       // by going upwards one step
				}

				right_path[ ++rp_it ] = dcopy;
				// push the newest element on the right path
			}

			typeTable[t] = type;  // congrats, you've just computed the type for this new block

			/**
			 * Compute answers for in-block-queries as seen in Alstrup et al. SPAA '02
			 */

			bool flag = false;

			if( Guardian[type] == 0 ){
				Guardian[type] = 1; flag = true;
			}

			if( flag )
			{
				inBlock[type] = (uchar*) calloc ( micro_block_size, sizeof(uchar) );
				// smallest value yet

				stack_size = 0;
				// init stack

				for( it = start; it < mblock_end; ++it )
				{
					DT dcopy = data[it];

					while( stack_size && ( dcopy < data[ stack[ stack_size ] ] ) ){
						--stack_size;
					} // pop elements bigger than current element

					if( stack_size )
					{
						uchar smaller = stack[ stack_size ] - start;
						// offset of first element on your left that is smaller than you

						inBlock[type][ it-start ] = inBlock[type][ smaller ] | ( 1 << smaller );
					}
					else
						inBlock[type][ it-start ] = 0; // smallest element

					stack[ ++stack_size ] = it; // push new element
				}
			}
		}

		// === free memory

		free( right_path );

		for( DTidx i = 0; i <= micro_block_size; ++i )
			free( Catalan[i] );
		free(Catalan);

		// deallocate memory for the right_path stack and Catalan table


		// === free memory

		free(Guardian);
		free(stack);
	}

	void serial_init_DP()
	{
		M_sz = log2fast(block_count);

		M    = (DTidx**) calloc( M_sz + 1, sizeof(DTidx*) );
		M[0] = (DTidx*)  calloc( block_count, sizeof(DTidx) );
		// allocate memory for matrix M

		sM_sz = log2fast(super_block_count);

		sM   = (DTidx**) calloc( sM_sz + 1, sizeof(DTidx*) );
		sM[0]= (DTidx*)  calloc( super_block_count + 1, sizeof(DTidx) );
		// allocate memory for matrix sM

		DT smin_v, min_v; // current super-block-min and block-min

		DTidx smin_p, min_p; // .. and their respective positions

		DTidx alpha = super_block_size / block_size; // how many blocks are there in a super_block

		DTidx it; // iterator

		DTidx end; // for easy computation of block end

		DTidx blk_cnt = 0; // block index

		// == initialization step in the DP - compute min for each individual (super)block

		for( DTidx i = 0; i <= size; i+=super_block_size ) // iterate through every super-block
		{
			smin_p = i;       // save position of super-block min
			smin_v = data[i]; // and it's actual value

			it = i; // position iterator on super_block_start

			for( DTidx j = 0; j < alpha; ++j, ++blk_cnt )
			{
				end = it + block_size;
				if( end > size )
				{
					end = size; // compute block end
					j = alpha;  // limit reached, finish the computation
				}

				min_p = it;
				min_v = data[it]; // init block_min and skip to the next item

				while( ++it < end ) // iterate through block
				{
						if( data[it] < min_v )
							min_p = it, min_v = data[it]; // update min
				}

				M[0][blk_cnt] = min_p;

				if( min_v < smin_v )
					smin_v = min_v, smin_p = min_p; // update super block

			}

			sM[0][ i/super_block_size ] = smin_p; // init super-block

		} // end DP initialization


	}

	void compute_table( DTidx**& T, DTidx size, DTidx count )
	{

		for( DTidx i = 1; i <= size; ++i )
		{
			DTidx pw = 1 << (i-1);
			DTidx sz = count - (1<<i);

			T[i] = (DTidx*) calloc( sz + 1, sizeof(DTidx) );

			for( DTidx j = 0; j < sz; ++j )
			{
				if( data[ T[i-1][j] ] <= data[ T[i-1][j+pw] ] ) // strict order
				{
					T[i][j] = T[i-1][j];
				}
				else
				{
					T[i][j] = T[i-1][j+pw];
				}
			}
		}

	}

	void parallel_compute_table( DTidx** T, DTidx size, DTidx count )
	{
//		#pragma omp parallel
//		{
//			printf("%d\n",omp_get_num_threads());
//		}

		for( DTidx i = 1; i <= size; ++i )
		{
			DTidx pw = 1 << (i-1);
			DTidx sz = count - (1<<i);

			T[i] = (DTidx*) calloc( sz + 1, sizeof(DTidx) );

			DTidx j;
			for( j = 0; j <= sz; ++j )
			{
				if( data[ T[i-1][j] ] <= data[ T[i-1][j+pw] ] ) // strict order
				{
					T[i][j] = T[i-1][j];
				}
				else
				{
					T[i][j] = T[i-1][j+pw];
				}
			}

		}

	}


	/**
	 * Compute answers to queries that span over multiple (super)-blocks using the classic
	 * NlogN DP trick
	 */

	void serial_compute_Mtables()
	{
		serial_init_DP();

		compute_table(  M, M_sz, block_count );
		compute_table( sM, sM_sz, super_block_count );
	}

	void parallel_compute_Mtables()
	{
		serial_init_DP();

		parallel_compute_table( M, M_sz, block_count );
		parallel_compute_table( sM, sM_sz, super_block_count );
	}


	void init_serial()
	{
			serial_compute_typeTable();
			serial_compute_Mtables();
	}


	void init_parallel()
	{
			serial_compute_typeTable();
			parallel_compute_Mtables();
	}

public:

	RMQ( DT* data, DTidx size, bool serial )
	{
		this->size = size;
		this->data = data;

		/*
		 * Optimized fixed size ( as recommended here [1] )
		 * To make it dependent on input size, set
		 * block_size = log^(2+eps)(n)
		 * super_block_size = log n / (2 + delta)
		 * where eps and delta are arbitrary constants > 0
		 */

		micro_block_size = 1 << 3;	// one byte / word
		block_size 		 = 1 << 4;
		super_block_size = 1 << 8;  // multiple of block_size

		super_block_count = size / super_block_size + 1;
		block_count       = size / block_size + 1;
		micro_block_count = size / micro_block_size + 1;


		if( serial == false ){
			init_parallel();
		}
		else
			init_serial(); // debugging purpose for the moment
	}

	inline DTidx micro_query( DTidx block, DTidx blk_idx, DTidx a_offset, DTidx b_idx )
	{
		////std::cerr << "(Micro)Query between " << blk_idx + a_offset << " and " << b_idx << "\n";

		DTidx min_pos;

		if( a_offset != 0)
		{
			min_pos = clearbits( inBlock[ typeTable[ block ] ][ b_idx-blk_idx ], a_offset );
		}
		else
		{
			min_pos = inBlock[ typeTable[block] ][b_idx-blk_idx]; // no need to call clearbits(,)
		}

		if( min_pos == 0 )
			return b_idx;

		return blk_idx + lsb(min_pos);
	}

	inline DTidx block_query( DTidx blk_a, DTidx blk_b )
	{
		//std::cerr << "Query between block " << blk_a << " " << blk_b << "\n";

		if( blk_a == blk_b )
			return M[0][blk_a];

		DTidx log_ = log2fast( blk_b - blk_a );

		if( data[ M[log_][blk_a] ] <= data[ M[log_][blk_b - (1<<log_)] ] ) // ensure strict order
			return M[log_][blk_a];

		return M[log_][blk_b-(1<<log_)];
	}

	inline DTidx super_block_query( DTidx sblk_a, DTidx sblk_b )
	{
		//std::cerr << "Super query between sblocks " << sblk_a << "  " << sblk_b << "\n";

		if( sblk_a == sblk_b )
			return sM[0][sblk_a];

		DTidx log_ = log2fast( sblk_b - sblk_a );

		if( data[ sM[log_][sblk_a] ] <= data[ sM[log_][sblk_b - (1<<log_) ]] ) // ensure strict order
			return sM[ log_][ sblk_a ];

		return sM[log_][sblk_b-(1<<log_)];
	}

	/**
	 * Ensure strict order
	 *
	 * A < B iff value(A) < value(B) OR ( value(A) == value(B) AND pos(A) < pos(B) )
	 *
	 * @param: cmin   -> current min
	 *		   cmin_p -> min pos
	 *		   tmin   -> target min ( value we try to update min with )
	 *		   tmin_p -> target pos
	 */

	inline void update_min( DTidx& cmin_p, DT& cmin, DTidx& tmin_p, DT& tmin )
	{
		if( cmin > tmin || ( cmin == tmin && cmin_p > tmin_p ) )
 			cmin = tmin, cmin_p = tmin_p;
	}


	virtual DTidx query( DTidx a, DTidx b )
	{
		if( a > b ){
			a^=b; b^=a; a^=b; // swap positions
		}

		//std::cerr << "Query: " << a << "  " << b << "\n";

		DT dcopy;

		DTidx a_micro_block = a / micro_block_size;
		DTidx b_micro_block = b / micro_block_size;

		DTidx a_mb_start = a_micro_block * micro_block_size;
		DTidx a_offset   = a - a_mb_start;

		if( a_micro_block == b_micro_block ){
			return micro_query( a_micro_block, a_mb_start, a_offset, b );
		}

		//std::cerr << "Different micros\n";

		DTidx b_mb_start = b_micro_block * micro_block_size;

		DTidx min_a = micro_query( a_micro_block, a_mb_start, a_offset, a_mb_start + micro_block_size - 1 );
		DTidx min_b = micro_query( b_micro_block, b_mb_start, 0, b );

		//std::cerr << "Results: (values) " << data[min_a] << "  " << data[min_b] << "\n";

		dcopy = data[min_a];

		if( data[min_b] < dcopy )				// strict order is ensured
			dcopy = data[min_b], min_a = min_b; // min_a now stores the min of the two micro-queries

		//std::cerr << min_a << " <- current min\n";

		if( a_micro_block + 1 == b_micro_block ) // the micro_block are adjoin
			return min_a; // victory

		//std::cerr << "Not adjoin micros\n";

		DTidx a_block = a / block_size;
		DTidx b_block = b / block_size;

		DTidx a_blk_start = a_block * block_size;
		DTidx b_blk_start = b_block * block_size;

		DTidx min_tmp;

		if( a_blk_start + micro_block_size >= a )
		{
			DTidx next_mb_idx = a_blk_start + micro_block_size;
			min_tmp = micro_query( a_micro_block+1, next_mb_idx, 0,  next_mb_idx + micro_block_size - 1 );

			//std::cerr << "0. trying to improve with " << data[min_tmp] << "\n";

			update_min( min_a, dcopy, min_tmp, data[min_tmp] );

			//std::cerr << min_a << " <- current min " << data[min_a] << "\n";
		}

		if( b_blk_start + micro_block_size <= b )
		{
			DTidx next_mb_idx = b_blk_start;
			min_tmp = micro_query( b_micro_block-1, next_mb_idx, 0, next_mb_idx + micro_block_size - 1);

			//std::cerr << "1. trying to improve with " << data[min_tmp] << "\n";

			update_min( min_a, dcopy, min_tmp, data[min_tmp] );

			//std::cerr << min_a << " <- current min " << data[min_a] << "\n";
		}

		if( a_block + 1 == b_block ) // the two blocks are adjoin
			return min_a; // yet again, victory!

		//std::cerr << "Not adjoin blocks\n";

		if( b_blk_start - a_blk_start - block_size <= super_block_size )
		{
			DTidx blk_q = block_query( a_block + 1, b_block );

			//std::cerr << "2. trying to improve with " << data[blk_q] << "\n";

			update_min( min_a, dcopy, blk_q, data[blk_q] );

			//std::cerr << min_a << " <- current min " << data[min_a] << "\n";

			return min_a;
		}

		DTidx a_sblk = a / super_block_size;
		DTidx b_sblk = b / super_block_size;

		//=== block query
		DTidx last_block = (a_sblk + 1) * super_block_size / block_size;
		// last block in the super-block containing a

		//std::cerr << "a's block" << a_block << " starts on " << a_blk_start << "\n";

		DTidx blk_q = block_query( a_block + 1, last_block );

		//std::cerr << "value:" << data[blk_q] << " " << blk_q << "\n";

		update_min( min_a, dcopy, blk_q, data[blk_q] );

//		if( data[ blk_q ] < dcopy )
//			min_a = blk_q, dcopy = data[blk_q]; // update min

		//=== block query
		DTidx first_block = b_sblk * super_block_size / block_size;
		// first block in the super-block containing b

		blk_q = block_query( first_block, b_block );

		//std::cerr << "value:" << data[blk_q] << " " << blk_q << "\n";

		update_min( min_a, dcopy, blk_q, data[blk_q] );

		if( a_sblk + 1 == b_sblk ) // the two superblocks are adjoin
			return min_a;

		//== super block query

		//std::cerr << "super query (cmin ->)" << min_a << "\n";

		DTidx sblk_q = super_block_query( a_sblk + 1, b_sblk );

		update_min( min_a, dcopy, sblk_q, data[sblk_q] );

		return min_a;
	}

	virtual ~RMQ()
	{
		free(typeTable);

		for( DTidx i = 0; i < inBlock_size; ++i )
			free(inBlock[i]);
		free(inBlock);

		for( DTidx i = 0; i <= M_sz; ++i )
			free(M[i]);
		free(M);

		for( DTidx i = 0; i <= sM_sz; ++i )
			free(sM[i]);
		free(sM);

	}


	void init_deprecated()
	{
		DTidx super_block_count = size / super_block_size + 1;
					DTidx block_count 	    = size / block_size + 1;
					DTidx micro_block_count = size / micro_block_size + 1;

					DTidx** Catalan = compute_catalan_table( micro_block_size+1 );
					// precompute Catalan for fast lookup

					/**
					 * Compute type for every micro_block as seen in Fischer/Heun CPM'06
					 */

					typeTable  = (ushort*) calloc( micro_block_count,    sizeof( ushort ) );
					DT* right_path = (DT*) calloc( micro_block_size + 1, sizeof(DT) );

					inBlock = (uchar**) calloc ( Catalan[ micro_block_size ][ micro_block_size ], sizeof(uchar*) );
					// allocate space for Alstrup et al. in-block preprocessing

					inBlock_size = Catalan[ micro_block_size ][ micro_block_size ];
					// save max size ( so we can free Catalan matrix before inBlock )

					char* Guardian = (char*) calloc ( inBlock_size, sizeof(char) );

					/*
					 * Install a guarding array ( 0 by default, idx if type was found )
					 * We want to avoid allocating memory / computing types for types that don't
					 * exist in the initial array
					 */

					DTidx* stack = ( DTidx* ) calloc ( micro_block_size + 1, sizeof(DTidx) );
					uchar stack_size;

					// when going parallel, allocate as many rp ( and stacks ) as nr of threads
					// and don't forget to free the memory

					// TODO - parallel

					for( DTidx t = 0; t < micro_block_count; ++t )
					{
						DTidx start = t * micro_block_size;

						DTidx it = start;
						// position the iterator on the t'th micro_block

						DTidx mblock_end = it + micro_block_size;
						// compute micro-block end

						mblock_end = mblock_end > size ? size : mblock_end;
						// check boundary condition; if last micro-block is smaller than 'micro_block_size'
						// consider the array it's filled with INF values

						DTidx rp_it = 0;
						// initially, there are no elements in the tree ( implicitly, none on the right path )

						DTidx q = micro_block_size;
						DTidx p = q - 1;
						// start from the bottom of the table

						DT dcopy;
						// store value pointed by iterator to avoid table lookup

						ushort type = 0;
						// store block type

						for( ; it < mblock_end; ++it, --p ) // iterate through block / go 'westwards'
						{
							dcopy = data[it];
							// copy data at with offset i

							// as long as the current element is smaller than the one on top of right_path stack
							while( rp_it && right_path[ rp_it ] > dcopy )
							{
								type += Catalan[q][p]; // skip C(q,p) paths
								--rp_it; --q;	       // by going upwards one step
							}

							right_path[ ++rp_it ] = dcopy;
							// push the newest element on the right path
						}

						typeTable[t] = type;  // congrats, you've just computed the type for this new block

						/**
						 * Compute answers for in-block-queries as seen in Alstrup et al. SPAA '02
						 */

						bool flag = false;

						if( Guardian[type] == 0 ){
							Guardian[type] = 1; flag = true;
						}

						if( flag )
						{
							inBlock[type] = (uchar*) calloc ( micro_block_size, sizeof(uchar) );
							// smallest value yet

							stack_size = 0;
							// init stack

							for( it = start; it < mblock_end; ++it )
							{
								DT dcopy = data[it];

								while( stack_size && ( dcopy < data[ stack[ stack_size ] ] ) ){
									--stack_size;
								} // pop elements bigger than current element

								if( stack_size )
								{
									uchar smaller = stack[ stack_size ] - start;
									// position of first element on your left that is smaller than you
									// ..relative to the block start

									inBlock[type][ it-start ] = inBlock[type][ smaller ] | ( 1 << smaller );
								}
								else
									inBlock[type][ it-start ] = 0;

								stack[ ++stack_size ] = it;
							}
						}
					}

					// === free memory

					free( right_path );

					for( DTidx i = 0; i <= micro_block_size; ++i )
						free( Catalan[i] );
					free(Catalan);

					// deallocate memory for the right_path stack and Catalan table


					// === free memory

					free(Guardian);
					free(stack);

					/**
					 * Compute answers to queries that span over multiple (super)-blocks using the classic
					 * NlogN DP trick
					 */

					M_sz = log2fast(block_count);

					M    = (DTidx**) calloc( M_sz + 1, sizeof(DTidx*) );
					M[0] = (DTidx*)  calloc( block_count, sizeof(DTidx) );

					sM_sz = log2fast(super_block_count);

					sM   = (DTidx**) calloc( sM_sz + 1, sizeof(DTidx*) );
					sM[0]= (DTidx*)  calloc( super_block_count, sizeof(DTidx) );

					// allocate memory for the two matrixes

					DT smin_v, min_v;
					// current super-block-min and block-min

					DTidx smin_p, min_p;
					// .. and their respective positions

					DTidx alpha = super_block_size / block_size;
					// how many blocks are there in a super_block

					DTidx it;
					// iterator

					DTidx end;
					// for easy computation of block end

					DTidx blk_cnt = 0;

					// TODO - parallel

					// == initialization step in the DP - compute min for each individual (super)block

					for( DTidx i = 0; i < size; i+=super_block_size )
					{
						// //std::cerr << "Analyzing super_block " << i / super_block_size << "\n";

						smin_p = i;
						smin_v = data[i]; // init super-block min and it's position

						////std::cerr << "sblock init:" << smin_p << " " << smin_v << "\n";

						it = i;

						for( DTidx j = 0; j < alpha; ++j, ++blk_cnt )
						{
							////std::cerr << "Analyzing block " << blk_cnt << "\n";

							end = it + block_size;
							if( end > size )
							{
								end = size; // compute block end
								j = alpha;  // limit reached, finish the computation
							}

							min_p = it;
							min_v = data[it]; // init block_min and skip to the next item

							////std::cerr << "block init:" << min_p << " " << min_v << "\n";

							while( ++it < end ) // iterate through block
							{
								////std::cerr << "block content " << it << " " << data[it] << "\n";

								if( data[it] < min_v )
									min_p = it, min_v = data[it]; // update min

								////std::cerr << min_p << " " << min_v << "\n";
							}

							////std::cerr << "Final values (block)" << min_p << " " << min_v << "\n";

							M[0][blk_cnt] = min_p;

							if( min_v < smin_v )
								smin_v = min_v, smin_p = min_p; // update super block

						}

						////std::cerr << "Final values (sblock)" << smin_p << " " << smin_v << "\n";

						sM[0][ i/super_block_size ] = smin_p;

						if( end == size )
							break;

					} // end DP initialization


					// //std::cerr << "Blk count " << block_count << "\n";

			//		//std::cerr << "Computed block values:\n";
			//		for( uint i = 0; i < block_count; ++i )
			//			//std::cerr << M[0][i] << "\n";
			//
			//		//std::cerr << "Computed super block values:\n";
			//		for( uint i = 0; i < super_block_count; ++i )
			//			//std::cerr << sM[0][i] << "  " << data[ sM[0][i] ] << "\n";

				// TODO - parallel

				/**
				 *	M[i][j] -> stores the position of the minimum in the subarray
				 *			   data[ j * block_size, (j + 2^i) * block_size - 1 ]
				 */

				for( DTidx i = 1; (unsigned) (1<<i) <= block_count; ++i )
				{
					DTidx pw = 1 << (i-1);
					DTidx sz = block_count - (1<<i);

					M[i] = (DTidx*) calloc( sz + 1, sizeof(DTidx) );

					for( DTidx j = 0; j <= sz; ++j )
					{
						if( data[ M[i-1][j] ] <= data[ M[i-1][j+pw] ] ) // strict order
						{
							M[i][j] = M[i-1][j];
						}
						else
						{
							M[i][j] = M[i-1][j+pw];
						}
					}
				}

				/**
				 *	sM[i][j] -> stores the position of the minimum in the subarray
				 *			    data[ j * block_size, (j + 2^i) * block_size - 1 ]
				 */

				for( DTidx i = 1; (unsigned) (1<<i) <= super_block_count; ++i )
				{
					DTidx pw = 1 << (i-1);
					DTidx sz = super_block_count - (1<<i);

					sM[i] = (DTidx*) calloc( sz+1, sizeof(DTidx) );

					for( DTidx j = 0; j <= sz; ++j )
					{
						if( data[ sM[i-1][j] ] <= data[ sM[i-1][j+pw] ] ) // strict order
						{
							sM[i][j] = sM[i-1][j];
						}
						else
						{
							sM[i][j] = sM[i-1][j+pw];
						}
					}
				}
	}
};

#endif /* RMQ_HPP_ */
