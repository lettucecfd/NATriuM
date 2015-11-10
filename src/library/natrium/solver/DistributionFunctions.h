/*
 * DistributionFunctions.h
 *
 *  Created on: 18.03.2014
 *      Author: kraemer
 */

#ifndef DISTRIBUTIONFUNCTIONS_H_
#define DISTRIBUTIONFUNCTIONS_H_

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short This class contains the distribution functions.
 *        As the zero-velocity particles (f0) can be ignored for streaming,
 *        f0 is split from the other distribution functions (fStream) in the implementation.
 *        fStream is a block vector, which has the advantage that it can be multiplied with the SystemMatrix.
 *        The DistributionFunctions class is designed to mime the behaviour of a vector<distributed_vector> .
 *
 */
class DistributionFunctions {
private:

	/// number of discrete velocities
	size_t m_Q;

	/// zero-velocity particle distribution functions
	distributed_vector m_f0;

	/// other distributions functions
	distributed_block_vector m_fStream;

public:
	/**
	 * @short empty constructor. Construction is done through reinit.
	 */
	DistributionFunctions() :
			m_Q(0) {
		// initialization through reinit!
	}

	/**
	 * @short Copy constructor
	 */
	DistributionFunctions(const DistributionFunctions& f) :
			m_Q(f.getQ()), m_f0(f.getF0()), m_fStream(f.getFStream()) {
	}

	/**
	 * @short Copy constructor. Conversion from vector<distributed_vector>.
	 */
	DistributionFunctions(const vector<distributed_vector>& f) ;

	/// Destructor
	virtual ~DistributionFunctions() {

	}

	/**
	 * @short mimes std::vector.at(i)
	 */
	distributed_vector& at(size_t i);

	/**
	 * @short mimes std::vector.at(i)
	 */
	const distributed_vector& at(size_t i) const ;

	/**
	 * @short F0 denotes the vector \f$ f_0 \f$ (zero-velocity particles)
	 */
	const distributed_vector& getF0() const {
		return m_f0;
	}

	/**
	 * @short F0 denotes the vector \f$ f_0 \f$ (zero-velocity particles)
	 */
	distributed_vector& getF0() {
		return m_f0;
	}

	/**
	 * @short F0 denotes the vector \f$ f_0 \f$ (zero-velocity particles)
	 */
	void setF0(const distributed_vector& f0) {
		m_f0 = f0;
	}

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	distributed_block_vector& getFStream() {
		return m_fStream;
	}

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	const distributed_block_vector& getFStream() const {
		return m_fStream;
	}

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	void setFStream(const distributed_block_vector& fStream) {
		m_fStream = fStream;
	}

	/**
	 * @short the number of discrete velocities
	 */
	size_t getQ() const {
		return m_Q;
	}

#ifdef WITH_TRILINOS_MPI
	/**
	 * @short reinitialize the sizes of the distribution functions - without ghost elements
	 */
	void reinit(size_t Q, const dealii::IndexSet &local, const dealii::IndexSet &relevant,
			const MPI_Comm &communicator = MPI_COMM_WORLD);

	/**
	 * @short reinitialize the sizes of the distribution functions - with ghost elements
	 */
	void reinit(size_t Q, const dealii::IndexSet &local,
			const MPI_Comm &communicator = MPI_COMM_WORLD);

#else
	void reinit(size_t Q, size_t size);
#endif

	/**
	 * @short the number of discrete velocities, including zero
	 */
	size_t size() const {
		return m_Q;
	}


	/**
	 * @short call dealii's compress function to all distributed_vectors stored herein. Compress has to
	 * be called, whenever the elements of a vector have been changed by hand. It distributes the local
	 * information to the other processors, if required.
	 * @param[in] operation specifies, if element was inserted or added. Operation has to be
	 * dealii::VectorOperation::add or dealii::VectorOperation::compress
	 */
	void compress(dealii::VectorOperation::values operation);

	/**
	 * @short call operator= for all distribution functions
	 * @note is only allowed for DistributionFunctions of equal size
	 */
	void operator=(const DistributionFunctions& other);

};
/* class DistributionFunctions */

} /* namespace natrium */

#endif /* DISTRIBUTIONFUNCTIONS_H_ */
