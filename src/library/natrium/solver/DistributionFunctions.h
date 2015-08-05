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
	DistributionFunctions(const vector<distributed_vector>& f) :
			m_Q(f.size()), m_f0(f.at(0)) {
#ifdef WITH_TRILINOS
		m_fStream.reinit(m_Q);
#else
		m_fStream.reinit(m_Q, m_f0.size());
#endif
		for (size_t i = 1; i < m_Q; i++) {
			m_fStream.block(i - 1).reinit(f.at(i));
			// reinit does only change the size but not the content
		}
		m_fStream.collect_sizes();
		for (size_t i = 1; i < m_Q; i++) {
			for (size_t j = 0; j < f.at(i).size(); j++){
				m_fStream.block(i-1)(j) = f.at(i)(j);
			}
		}
	}

	/// Destructor
	virtual ~DistributionFunctions() {

	}

	/**
	 * @short mimes std::vector.at(i)
	 */
	distributed_vector& at(size_t i) {
		assert(m_Q > 0);
		assert(i < m_Q);
		if (i == 0) {
			return m_f0;
		} else {
			return m_fStream.block(i - 1);
		}
	}

	/**
	 * @short mimes std::vector.at(i)
	 */
	const distributed_vector& at(size_t i) const {
		assert(m_Q > 0);
		assert(i < m_Q);
		if (i == 0) {
			return m_f0;
		} else {
			return m_fStream.block(i - 1);
		}
	}

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
	const size_t getQ() const {
		return m_Q;
	}

	/**
	 * @short reinitialize the sizes of the distribution functions
	 */
#ifdef WITH_TRILINOS_MPI
void reinit(size_t Q, const dealii::IndexSet &local, const dealii::IndexSet &ghost, const MPI_Comm &communicator=MPI_COMM_WORLD) {
		m_Q = Q;
		m_f0.reinit(local, ghost, communicator);
		m_fStream.reinit(Q-1);
		for (size_t i = 0; i < Q - 1; i++){
			m_fStream.block(i).reinit(m_f0);
		}
		m_fStream.collect_sizes();
	}
#else
	void reinit(size_t Q, size_t size) {
		m_Q = Q;
		m_f0.reinit(size);
#ifdef WITH_TRILINOS
		m_fStream.reinit(Q-1);
		for (size_t i = 0; i < Q - 1; i++){
			m_fStream.block(i).reinit(m_f0);
		}
#else
		m_fStream.reinit(Q - 1, size);
#endif
		m_fStream.collect_sizes();
	}
#endif

	/**
	 * @short the number of discrete velocities, including zero
	 */
	size_t size() const {
		return m_Q;
	}

};

} /* namespace natrium */

#endif /* DISTRIBUTIONFUNCTIONS_H_ */
