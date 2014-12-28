/*
 * rmq_brute.h
 *
 *  Created on: Oct 13, 2012
 *      Author: alexei
 */

#ifndef RMQ_BRUTE_H_
#define RMQ_BRUTE_H_

#include <stdlib.h>
#include "RangeMinimumQuery.hpp"

template< class DT = int >
class rmq_brute : public RangeMinimumQuery {

private:

	uint size;
	uint** Matrix;

	/**
	 * If the instance case is small enough, use the fastest ( dumbest ) solution available.
	 *
	 * This one takes O(N^2) space, where N -> instance size and query is O(1) ( simple
	 * table interrogation )
	 *
	 */

	uint** Four_Russian_Trick( std::vector<DT>& data )
	{
		size = data.size();

		uint** Matrix = (uint**) calloc( size, sizeof(uint*) );

		DT cmin;
		register uint prevPos;

		for( uint i = 0; i < size; ++i )
		{
			Matrix[i] = (uint*) calloc( size, sizeof(uint) );

			cmin = data[i];
			Matrix[i][i] = prevPos = i;

			for( uint j = i + 1; j < size; ++j )
			{
				if( cmin <= data[j] )
					Matrix[i][j] = prevPos;
				else
				{
					Matrix[i][j] = j;
					cmin = data[j];
				}

				prevPos = Matrix[i][j];
			}
		}

		for( uint i = 0; i < size; ++ i )
			for( uint j = i+1; j < size; ++j )
				Matrix[j][i] = Matrix[i][j];

		return Matrix;
	}



public:

	rmq_brute( std::vector<DT>& data )
	{
		Matrix = this->Four_Russian_Trick( data );
	}

	uint query( uint a, uint b )
	{
		return Matrix[a][b];
	}

	virtual ~rmq_brute()
	{
		for( uint i = 0; i < size; ++i ){
			free( Matrix[i] );
		}

		free( Matrix );
	}
};

#endif /* RMQ_BRUTE_H_ */
