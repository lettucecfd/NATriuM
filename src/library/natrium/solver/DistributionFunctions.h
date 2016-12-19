/*
 * DistributionFunctions.h
 *
 *  Created on: 18.03.2014
 *      Author: kraemer
 */

#ifndef DISTRIBUTIONFUNCTIONS_H_
#define DISTRIBUTIONFUNCTIONS_H_

#include "../utilities/BasicNames.h"
#include "../stencils/Stencil.h"
#include "../utilities/NATriuMException.h"

#include "../utilities/Timing.h"

namespace natrium {

/**
 * @short Exception class for AdvectionOperator
 */
class DistributionFunctionsException: public NATriuMException {
private:
	std::string message;
public:
	DistributionFunctionsException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	DistributionFunctionsException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~DistributionFunctionsException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

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

	/// zero-velocity particle distribution functions, including ghost elements
	distributed_vector m_f0Ghosted;

	/// other distributions functions, including ghost elements
	distributed_block_vector m_fStreamGhosted;

	/// flag that indicates whether we have a parallel DG discretization;
	/// in this case, we do not need ghosted elements as all local dof indices
	/// are in fact locally owned on the current MPI process
	bool m_dg;

public:
	/**
	 * @short empty constructor. Construction is done through reinit.
	 * @param dg flag that indicates whether we have a parallel DG discretization;
	 *    	in this case, we do not need ghosted elements as all local dof indices
	 *   	are in fact locally owned on the current MPI process
	 */
	DistributionFunctions(bool dg = true) :
			m_Q(0), m_dg(dg) {
		// initialization through reinit!
	}

	/**
	 * @short Copy constructor
	 */
	DistributionFunctions(const DistributionFunctions& f) :
			m_Q(f.getQ()), m_f0(f.getF0()), m_fStream(f.getFStream()), m_f0Ghosted(f.m_f0Ghosted), m_fStreamGhosted(
					f.m_fStreamGhosted), m_dg(
							f.m_dg) {
	}

	/**
	 * @short Copy constructor. Conversion from vector<distributed_vector>.
	 */
	DistributionFunctions(const vector<distributed_vector>& f, bool dg = true);

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
	const distributed_vector& at(size_t i) const;

	/**
	 * @short mimes std::vector.at(i)
	 */
	distributed_vector& atGhosted(size_t i);

	/**
	 * @short mimes std::vector.at(i)
	 */
	const distributed_vector& atGhosted(size_t i) const;

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
	 * @short F0 denotes the vector \f$ f_0 \f$ (zero-velocity particles)
	 */
	const distributed_vector& getF0Ghosted() {
		if (m_dg)
			return m_f0;
		else {
			// TODO in Debug mode, check at least first element for equality
			return m_f0Ghosted;
		}
	}

	/**
	 * @short F0 denotes the vector \f$ f_0 \f$ (zero-velocity particles)
	 */
	/*distributed_vector& getF0Ghosted() {
	 if (m_dg)
	 return m_f0;
	 else {
	 // TODO in Debug mode, check at least first element for equality
	 return m_f0Ghosted;
	 }
	 }*/

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	/*distributed_block_vector& getFStreamGhosted() {
	 if (m_dg)
	 return m_fStream;
	 else {
	 // TODO in Debug mode, check at least first element for equality
	 return m_fStreamGhosted;
	 }
	 }*/

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	const distributed_block_vector& getFStreamGhosted() const {
		if (m_dg)
			return m_fStream;
		else {
			// TODO in Debug mode, check at least first element for equality
			return m_fStreamGhosted;
		}
	}

	/**
	 * @short FStream denotes the block vector containing the vectors \f$ f_1, ..., f_Q \f$
	 */
	void setFStream(const distributed_block_vector& fStream) {
		m_fStream = fStream;
		if (not m_dg) {
			m_fStreamGhosted = fStream;
		}
	}

	/**
	 * @short the number of discrete velocities
	 */
	size_t getQ() const {
		return m_Q;
	}

	/**
	 * @short reinitialize the sizes of the distribution functions - with ghost elements
	 */
	void reinit(size_t Q, const dealii::IndexSet &local,
			const dealii::IndexSet &relevant, const MPI_Comm &communicator =
			MPI_COMM_WORLD);

	/**
	 * @short reinitialize the sizes of the distribution functions - without ghost elements
	 */
	void reinit(size_t Q, const dealii::IndexSet &local,
			const MPI_Comm &communicator = MPI_COMM_WORLD);

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

	/**
	 * @short checks whether two DistributionFunction objects are equal wrt. a given threshold. Use carefully: vectors are only
	 *        considered equal if they have the same distributions among processors.
	 * @param other Other DistributionFunction
	 * @param threshold threshold (default: 1e-10)
	 * @note  If you want to compare vectors that have different distributions among processors, you can work around this functions
	 *        using Vector::add(v, allow_different_maps = true ) and check whether the entries are close to zero (requires additional memory)
	 */
	bool equals(const DistributionFunctions& other,
			double threshold = 1e-8) const;

	void transferFromOtherScaling(const Stencil& old_stencil,
			const Stencil& new_stencil,
			const dealii::IndexSet& locally_owned_dofs);

	size_t memory_consumption() const {
		return m_f0.memory_consumption() + m_fStream.memory_consumption();
	}

	/**
	 * @short update ghosted distributions, i.e. communicate over MPI processes
	 * It is very! important! that after each manipulation of the distribution functions,
	 * you call updateGhosted(). Otherwise you might get terrible errors at boundaries,
	 * (or in other parts, e.g. calculating the entropy), due to the fact that your distributions
	 * at ghost nodes are not up to date.
	 */
	void updateGhosted() {
		if (m_dg) {
			return;
		}
		// copy elements to ghosted vector
		// "=" is a collective operation and will communicate the ghost elements
		m_f0Ghosted = m_f0;
		m_fStreamGhosted = m_fStream;
		TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
	}

	bool isDg() const {
		return m_dg;
	}
};
/* class DistributionFunctions */

} /* namespace natrium */

#endif /* DISTRIBUTIONFUNCTIONS_H_ */
