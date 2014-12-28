/*
 * test.cpp
 *
 *  Created on: Oct 13, 2012
 *      Author: alexei
 */
#include<vector>
#include<stdio.h>
#include<time.h>
#include<iostream>

#include "rmq_brute.h"
#include "RMQ.hpp"

//#define TEST_INIT
#define TEST_QUERY

int N = 100000, M = 10000;

int main()
{
	// std::cout << V.size() << "\n";

	//rmq_brute<int> rmq(V);

	//std::cout << rmq.query(1,4) << "\n";

	//int data = new int[15];
	/*
	int* data = new int[19];
	data[0] = 1; data[1] =  2, data[2] =  3, data[3] =  1,data[4] =  3, data[5] =  2,data[6] =  -2048,data[7] = -2048,data[8] = 2,data[9] = 2,data[10]= 1,data[11]= 3,data[12]= 3,data[13]= 2,data[14]= 1;
	data[15] = 3; data[16] = 22; data[17] = 21, data[18] = 24;

	RMQ<int,uint> R( data, 18 );

	//std::cerr << R.query(0,2)  << " <- ans\n";
	//std::cerr << R.query(0,18) << " <- ans\n";
	*/


#ifdef TEST_INIT
	freopen("rmq2.in","w",stdout);

	srand( time(NULL) );
	// std::cout << N << " " << M << "\n";

	for( int i = 0; i < N; ++i ){
		std::cout << rand() % N << "\n";
	}

	for( int i = 0; i < M; ++i ){
		std::cout << rand() % N << "  " << rand() % N << "  " << i << "\n";
	}
#endif

#ifdef TEST_QUERY

	freopen("rmq2.in","r",stdin);
	freopen("rmq2_reference.out","w",stdout);

	int* data_p = (int*) calloc( N+1, sizeof(int) );

	for( int i = 0; i < N; ++i ){
		scanf("%d",&data_p[i]);
	}

	std::cout << "Hello world!" << "\n";

	RMQ<int,uint> R(data_p,N,true);

	bool flag = true;

	for( int i = 0; i < M; ++i )
	{
		int a, b, c, min, p;
		scanf("%d%d%d",&a,&b,&c);

		if( a > b ){
			a^=b;b^=a;a^=b;
		}

		min = data_p[a], p = a;

		for( int j = a+1; j <= b; ++j )
			if( data_p[j] < min )
			{
				min = data_p[j], p = j;
			}

		int pp = R.query(a,b);

		if( p != pp )
		{
			std::cerr << "Test: " << a << "  " << b << "\n";
			// //std::cerr << "Correct answer:" << min << " with pos " << p << "\n";

			std::cout << p << " " << min << " <- brute || ->" << pp << " with value " << data_p[pp] << " " << c << "\n";
			flag = false;
			break;
		}

		//break;
	}

	if( flag )
		std::cout << "Mission success!\n";

	free( data_p );
#endif

	return 0;
}
