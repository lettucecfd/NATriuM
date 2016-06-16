/*
 * SemiLagrangianVectorReferenceTypes.h
 *
 *  Created on: 16.06.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_
#define LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_

#include "../utilities/BasicNames.h"

namespace natrium {

typedef dealii::TrilinosWrappers::internal::VectorReference TrilinosVRef;

/**
 * @short A generalized degree of freedom that can be used for FE degrees of freedom or boundary hits.
 * The generalized dof is required to describe the propagation of information through Boundary Hits.
 */
class GeneralizedDoF {
private:
	bool m_primaryBoundaryHit;
	bool m_secondaryBoundaryHit;
	size_t m_index;
	size_t m_Q;
public:
	GeneralizedDoF(bool is_primary, bool is_secondary, size_t index,
			size_t q) :
			m_primaryBoundaryHit(is_primary), m_secondaryBoundaryHit(
					is_secondary), m_index(index), m_Q(q) {

	}
	GeneralizedDoF(const GeneralizedDoF& other) {
		m_primaryBoundaryHit = other.isPrimaryBoundaryHit();
		m_secondaryBoundaryHit = other.isSecondaryBoundaryHit();
		m_index = other.getIndex();
		m_Q = other.getQ();
	}
	virtual ~GeneralizedDoF() {
	}

	size_t getIndex() const {
		return m_index;
	}

	void setIndex(size_t index) {
		m_index = index;
	}

	bool isPrimaryBoundaryHit() const {
		return m_primaryBoundaryHit;
	}

	void setPrimaryBoundaryHit(bool primaryBoundaryHit) {
		m_primaryBoundaryHit = primaryBoundaryHit;
	}

	bool isSecondaryBoundaryHit() const {
		return m_secondaryBoundaryHit;
	}

	void setSecondaryBoundaryHit(bool secondaryBoundaryHit) {
		m_secondaryBoundaryHit = secondaryBoundaryHit;
	}

	size_t getQ() const {
		return m_Q;
	}

	void setQ(size_t q) {
		m_Q = q;
	}
};

class GeneralizedDoFVector {
private:

	const size_t m_numberOfBlocks;

	const dealii::IndexSet& m_locallyOwned;

	const distributed_block_vector& m_vector;

	distributed_vector m_primaryBoundaryValues;

	distributed_vector m_secondaryBoundaryValues;

	dealii::IndexSet m_primaryBoundaryIndices;

	dealii::IndexSet m_secondaryBoundaryIndices;

public:

	GeneralizedDoFVector(const dealii::IndexSet& locally_owned_dofs,
			const distributed_block_vector& vector) :
			m_numberOfBlocks(vector.n_blocks()), m_locallyOwned(
					locally_owned_dofs), m_vector(vector), m_primaryBoundaryIndices(
					locally_owned_dofs.size()), m_secondaryBoundaryIndices(
					locally_owned_dofs.size()) {

	}
	TrilinosVRef operator[](const GeneralizedDoF& dof) {
		if (dof.isPrimaryBoundaryHit()) {
			return m_primaryBoundaryValues(dof.getIndex());
		}
		if (dof.isSecondaryBoundaryHit()) {
			return m_secondaryBoundaryValues(dof.getIndex());
		}
		return m_vector.block(dof.getBlock())(dof.getIndex());
	}

	GeneralizedDoF appendBoundaryDoF(bool primary) {
		size_t index;
		if (primary) {
			assert(m_primaryBoundaryIndices.size() < m_locallyOwned.size());
			index = m_locallyOwned.nth_index_in_set(
					m_primaryBoundaryIndices.size());
			m_primaryBoundaryIndices.add_index(index);
		} else {
			assert(m_secondaryBoundaryIndices.size() < m_locallyOwned.size());
			index = m_locallyOwned.nth_index_in_set(
					m_secondaryBoundaryIndices.size());
			m_secondaryBoundaryIndices.add_index(index);
		}
		GeneralizedDoF result(primary, !primary, index, 0);
		return result;
	}

	void compress() {
		m_primaryBoundaryValues.reinit(m_primaryBoundaryValues, MPI_COMM_WORLD);
		m_secondaryBoundaryValues.reinit(m_secondaryBoundaryValues,
				MPI_COMM_WORLD);
	}

	const dealii::IndexSet&& getLocallyOwned() const {
		return m_locallyOwned;
	}

	const size_t getNumberOfBlocks() const {
		return m_numberOfBlocks;
	}

	const dealii::IndexSet& getPrimaryBoundaryIndices() const {
		return m_primaryBoundaryIndices;
	}

	const distributed_vector& getPrimaryBoundaryValues() const {
		return m_primaryBoundaryValues;
	}

	const dealii::IndexSet& getSecondaryBoundaryIndices() const {
		return m_secondaryBoundaryIndices;
	}

	const distributed_vector& getSecondaryBoundaryValues() const {
		return m_secondaryBoundaryValues;
	}

	const distributed_block_vector&& getVector() const {
		return m_vector;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANVECTORREFERENCETYPES_H_ */
