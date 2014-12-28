/*
 * RangeMinimumQuery.hpp
 *
 *  Created on: Oct 13, 2012
 *      Author: alexei
 */

#ifndef RANGEMINIMUMQUERY_HPP_
#define RANGEMINIMUMQUERY_HPP_

#include "reference.hpp"

class RangeMinimumQuery{
public:
	virtual uint query( uint pos1 , uint po2 ) = 0;
};


#endif /* RANGEMINIMUMQUERY_HPP_ */
