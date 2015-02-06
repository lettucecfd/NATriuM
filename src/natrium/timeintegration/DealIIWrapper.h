/*
 * DealIIWrapper.h
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#ifndef DEALIIWRAPPER_H_
#define DEALIIWRAPPER_H_

namespace natrium {

/**
 * @short not yet implemented
 */
template <class MATRIX, class VECTOR>
class DealIIWrapper {
private:

public:
	DealIIWrapper();
	virtual ~DealIIWrapper();
	virtual void step(VECTOR& vector, const MATRIX& systemMatrix, const VECTOR& systemVector){
		// TODO Not yet implemented.
	}
};

} /* namespace natrium */

#endif /* DEALIIWRAPPER_H_ */
