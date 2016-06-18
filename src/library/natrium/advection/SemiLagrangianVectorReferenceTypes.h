/*
 * SemiLagrangianVectorReferenceTypes.h
 *
 *  Created on: 16.06.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_
#define LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "mpi.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

typedef dealii::TrilinosWrappers::internal::VectorReference TrilinosVRef;

struct FunctionDepartureValue {
	bool isBoundary;
	size_t secondaryBoundaryDoF;
	vector<double> internalValues;
	vector<size_t> internalDoFs;
	size_t alpha;
	FunctionDepartureValue(size_t boundary_dof) :
			isBoundary(true), secondaryBoundaryDoF(boundary_dof), alpha(0) {
	}
};

/**
 * @short A generalized degree of freedom that can be used for FE degrees of freedom or boundary hits.
 * The generalized dof is required to describe the propagation of information through Boundary Hits.
 */
class GeneralizedDestinationDoF {
private:
	bool m_secondaryBoundaryHit;
	size_t m_index;
	size_t m_alpha;
public:

	/**
	 * @short Constructor
	 * @param[in] is_primary indicates whether this dof is a primary boundary dof (i.e. it does not depend on other boundary dofs)
	 * @param[in] is_secondary indicates whether this dof is a secondary boundary dof (i.e. it does depend on other boundary dofs)
	 * @param[in] index the index (in locally_owned_dofs) of this dof
	 * @param[in] alpha the id of the streaming direction
	 * @note if is_primary == false and is_secondary == false, the dof is not a boundary dof
	 */
	GeneralizedDestinationDoF(bool is_secondary, size_t index, size_t alpha) :
			m_secondaryBoundaryHit(is_secondary), m_index(index), m_alpha(alpha) {

	}
	/**
	 * @short copy constructor
	 */
	GeneralizedDestinationDoF(const GeneralizedDestinationDoF& other) {
		m_secondaryBoundaryHit = other.isSecondaryBoundaryHit();
		m_index = other.getIndex();
		m_alpha = other.getAlpha();
	}
	virtual ~GeneralizedDestinationDoF() {
	}

	/**
	 * @short get the dof index (i.e. one of the locally owned dofs)
	 */
	size_t getIndex() const {
		return m_index;
	}

	/**
	 * @short set index to a the locally owned dof
	 */
	void setIndex(size_t index) {
		m_index = index;
	}

	/**
	 * @short indicates whether this dof is a secondary boundary dof (i.e. other boundary dof depends on its value)
	 */
	bool isSecondaryBoundaryHit() const {
		return m_secondaryBoundaryHit;
	}

	/**
	 * @short set this dof to be a secondary boundary dof
	 */
	void setSecondaryBoundaryHit(bool secondaryBoundaryHit) {
		m_secondaryBoundaryHit = secondaryBoundaryHit;
	}

	/**
	 * @short get the streaming direction of this dof
	 */
	size_t getAlpha() const {
		return m_alpha;
	}

	/**
	 * @short set the streaming direction of this dof
	 */
	void setAlpha(size_t alpha) {
		m_alpha = alpha;
	}
};

class SemiLagrangianVectorAccessException: public NATriuMException {
private:
	std::string message;
public:
	SemiLagrangianVectorAccessException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	SemiLagrangianVectorAccessException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~SemiLagrangianVectorAccessException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short A generalized version of the distributed_block_vector for the semi-Lagrangian boundary description
 * The semi-Lagrangian advection solver requires to track the distribution functions
 * to their departure points across boundaries. To allow for arbitrary nonlinear boundaries, the information
 * at the boundary has to be stored in addition to the usual degrees of freedom. The present class is container
 * for the boundary values. It distinguishes between primary and secondary boundary values. Secondary boundary values
 * denote those that are required by other boundaries -- which only happens when a Lagrangian path hits more
 * than one boundary during a single time step. This situation may occur e.g. in the corners of the computational domain.
 * Each GeneralizedDoFVector is built upon a distributed_block_vector, i.e. a dealii::TrilinosWrappers::MPI::BlockVector.
 */
class SecondaryBoundaryDoFVector {
private:

	/**
	 * @short underlying distributed_block_vector
	 */
	const dealii::IndexSet& m_locallyOwnedDoFs;

	/**
	 * @short the vector that stores the secondary boundary values, i.e. those that are required to calculate other boundary values
	 */
	distributed_vector m_secondaryBoundaryValues;

	/**
	 * @short indices of secondary boundary values (a subset of the locally owned dofs).
	 * The boundary dof indices are arbitrary and have nothing to do with the usual dof indices.
	 * They are only defined as a subset of the locally owned dofs to facilitate the
	 * distribution among processors with a unique indexing. This procedure assumes that
	 * there are more locally owned dofs than boundary hits, which will hopefully be the case.
	 */
	dealii::IndexSet m_secondaryBoundaryIndices;

public:

	/**
	 * Constructor.
	 * @param[in] the underlying distributed_block_vector which should at least contain one block
	 * The underlying distributed_block_vector does not have to be fully filled or initialized.
	 */
	SecondaryBoundaryDoFVector(const dealii::IndexSet& locally_owned_dofs) :
			m_locallyOwnedDoFs(locally_owned_dofs), m_secondaryBoundaryIndices(
					locally_owned_dofs.size()) {

	}

	/**
	 * @short reinitialize the boundary indices. Has to be called after the distributed_block_vector has the right size.
	 */
	/*void reinit() {
	 assert(m_secondaryBoundaryIndices.n_elements() == 0);
	 m_secondaryBoundaryIndices.set_size(
	 m_vector.block(0).locally_owned_elements().size());

	 }*/

	/**
	 * @short add a boundary dof and return its generalized dof index
	 */
	GeneralizedDestinationDoF appendSecondaryBoundaryDoF() {
		size_t index;
		assert(
				m_secondaryBoundaryIndices.n_elements()
						< m_locallyOwnedDoFs.n_elements());
		index = m_locallyOwnedDoFs.nth_index_in_set(
				m_secondaryBoundaryIndices.n_elements());
		m_secondaryBoundaryIndices.add_index(index);

		GeneralizedDestinationDoF result(true, index, 0);
		return result;
	}

	/**
	 * Allocate memory for the primary and secondary boundary values
	 */
	void compress() {
		m_secondaryBoundaryValues.reinit(m_secondaryBoundaryIndices);
	}

	TrilinosVRef operator()(size_t index) {
		TrilinosVRef ref = m_secondaryBoundaryValues(index);
		return ref;
	}

	/**
	 *  @return secondary boundary indices, see class description
	 */
	const dealii::IndexSet& getSecondaryBoundaryIndices() const {
		return m_secondaryBoundaryIndices;
	}

	/**
	 * @return secondary boundary values, see class description
	 */
	const distributed_vector& getSecondaryBoundaryValues() const {
		return m_secondaryBoundaryValues;
	}

};

class SemiLagrangianVectorAccess {
private:
	const DistributionFunctions& m_fOld;
	DistributionFunctions& m_fNew;
	SecondaryBoundaryDoFVector& m_fSecondaryBoundary;
public:
	SemiLagrangianVectorAccess(const DistributionFunctions& f_old,
			DistributionFunctions& f_new,
			SecondaryBoundaryDoFVector& f_secondary_boundary) :
			m_fOld(f_old), m_fNew(f_new), m_fSecondaryBoundary(
					f_secondary_boundary) {
	}

	/**
	 * @short write access to an element of the generalized vector.
	 * @param[in] dof may denote a boundary hit or a usual dof
	 * @note Primary boundary hits are not explicitly stored and can thus not be accessed.
	 */
	TrilinosVRef operator[](const GeneralizedDestinationDoF& dof) {
		if (dof.isSecondaryBoundaryHit()) {
			return m_fSecondaryBoundary(dof.getIndex());
		}
		assert(dof.getAlpha() > 0);
		return m_fNew.at(dof.getAlpha())(dof.getIndex());
		/*TrilinosVRef ref = m_fNew.at(dof.getAlpha())(dof.getIndex());
		 return ref;*/
	}

	/**
	 * @short read only access (to f_old)
	 */
	double operator()(const FunctionDepartureValue& fv) const {
		if (fv.isBoundary) {
			return m_fSecondaryBoundary(fv.secondaryBoundaryDoF);
		}
		assert(fv.internalDoFs.size() == fv.internalValues.size());
		assert(fv.alpha < m_fOld.getQ());
		double result = 0.0;
		for (size_t i = 0; i < fv.internalDoFs.size(); i++) {
			result += m_fOld.at(fv.alpha)(fv.internalDoFs[i])
					* fv.internalValues[i];
		}
		return result;
	}

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_ */
