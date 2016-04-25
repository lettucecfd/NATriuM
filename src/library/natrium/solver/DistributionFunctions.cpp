/**
 * @file DistributionFunctions.cpp
 * @short
 * @date 09.11.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DistributionFunctions.h"

namespace natrium {

DistributionFunctions::DistributionFunctions(
		const vector<distributed_vector>& f) :
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
		m_fStream.block(i - 1) = f.at(i);
	}
}

distributed_vector& DistributionFunctions::at(size_t i) {
	assert(m_Q > 0);
	assert(i < m_Q);
	if (i == 0) {
		return m_f0;
	} else {
		return m_fStream.block(i - 1);
	}
}

const distributed_vector& DistributionFunctions::at(size_t i) const {
	assert(m_Q > 0);
	assert(i < m_Q);
	if (i == 0) {
		return m_f0;
	} else {
		return m_fStream.block(i - 1);
	}
}

void DistributionFunctions::reinit(size_t Q, const dealii::IndexSet &local,
		const dealii::IndexSet &relevant, const MPI_Comm &communicator) {
	m_Q = Q;
	m_f0.reinit(local, relevant, communicator);
	m_fStream.reinit(Q - 1);
	for (size_t i = 0; i < Q - 1; i++) {
		m_fStream.block(i).reinit(m_f0);
	}
	m_fStream.collect_sizes();
}
void DistributionFunctions::reinit(size_t Q, const dealii::IndexSet &local,
		const MPI_Comm &communicator) {
	m_Q = Q;
	m_f0.reinit(local, communicator);
	m_fStream.reinit(Q - 1);
	for (size_t i = 0; i < Q - 1; i++) {
		m_fStream.block(i).reinit(m_f0);
	}
	m_fStream.collect_sizes();
}

void DistributionFunctions::compress(
		dealii::VectorOperation::values operation) {
	m_f0.compress(operation);
	m_fStream.compress(operation);

}

void DistributionFunctions::operator=(const DistributionFunctions& other) {
	m_f0 = other.getF0();
	m_fStream = other.getFStream();
}

bool DistributionFunctions::equals(const DistributionFunctions& other,
		double threshold) const {
	if (size() != other.size()) {
		return false;
	}
	if (size() == 0) {
		// empty vectors are defined equal
		return true;
	}

	// check elements
	bool result = true;
	for (size_t i = 0; i < size(); i++) {
		if (result == false){
			break;
		}
		const distributed_vector& fi = at(i);
		const distributed_vector& gi = other.at(i);
		const dealii::IndexSet& indices = fi.locally_owned_elements();
		const dealii::IndexSet& other_indices = gi.locally_owned_elements();
		dealii::IndexSet::ElementIterator it = indices.begin();
		dealii::IndexSet::ElementIterator end = indices.end();
		for (it = indices.begin(); it != end; it++) {
			size_t j = *it;
			if (not other_indices.is_element(j)) {
				result = false;
				break;
			} else if (fabs(fi(j) - gi(j)) > threshold) {
				result = false;
				break;
			}
		}
	}
	return dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD).min;

}


void DistributionFunctions::transferFromOtherScaling(const Stencil& old_stencil,
		const Stencil& new_stencil,
		const dealii::IndexSet& locally_owned_dofs) {


	assert(new_stencil.getQ() == m_Q);
	assert(old_stencil.getQ() == m_Q);
	assert(old_stencil.getD() == new_stencil.getD());

	if (abs(old_stencil.getSpeedOfSound() - new_stencil.getSpeedOfSound())
			< 1e-5) {
		return;
	}

	// vectors and matrices for transformations
	numeric_vector old_f(m_Q);
	numeric_vector new_f(m_Q);
	numeric_matrix old_f_to_M(m_Q);
	numeric_matrix M_to_new_f(m_Q);
	numeric_matrix T(m_Q); // trafo matrix
	// avoid calls to block() // avoid calls to getDirections
	std::vector<distributed_vector*> f;
	for (size_t i = 0; i < m_Q; i++) {
		distributed_vector* fi = &at(i);
		f.push_back(fi);
	}
	// get transformation matrices
	old_stencil.getMomentBasis(old_f_to_M);
	new_stencil.getInverseMomentBasis(M_to_new_f);
	M_to_new_f.mmult(T, old_f_to_M);

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;

		// fill old_f
		for (size_t j = 0; j < m_Q; j++) {
			old_f(j) = (*f[j])(i);
		}
		// Trafo
		T.vmult(new_f, old_f);
		// assign back to global dofs
		for (size_t j = 0; j < m_Q; j++) {
			(*f[j])(i) = new_f(j);
		}
	}
}

} /* namespace natrium */
