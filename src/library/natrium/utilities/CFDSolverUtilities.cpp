/**
 * @file CFDSolverUtilities.cpp
 * @short
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CFDSolverUtilities.h"

#include <fstream>

#include "deal.II/grid/grid_tools.h"
#include "deal.II/grid/grid_out.h"
#include "deal.II/base/geometry_info.h"
#include "deal.II/base/geometry_info.h"

namespace natrium {

template<size_t dim>
double CFDSolverUtilities::getMinimumDoFDistanceGLL(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement) {
	assert(orderOfFiniteElement >= 1);
	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = CFDSolverUtilities::getMinimumVertexDistance<
			dim>(tria);

	// calculate distance between closest quadrature nodes on a line
	dealii::QGaussLobatto<1> quadrature(orderOfFiniteElement + 1);
	double min_dof_distance = 10000;
	for (size_t i = 0; i < orderOfFiniteElement + 1; i++) {
		for (size_t j = i + 1; j < orderOfFiniteElement + 1; j++) {
			double pointDist = quadrature.get_points().at(i).distance(
					quadrature.get_points().at(j));
			if (pointDist < min_dof_distance) {
				min_dof_distance = pointDist;
			}
		}
	}
	return min_vertex_distance * min_dof_distance;
}
template double CFDSolverUtilities::getMinimumDoFDistanceGLL<2>(
		const Mesh<2>& tria, const size_t orderOfFiniteElement);
template double CFDSolverUtilities::getMinimumDoFDistanceGLL<3>(
		const Mesh<3>& tria, const size_t orderOfFiniteElement);

template<size_t dim>
double CFDSolverUtilities::getMinimumDoFDistance(const Mesh<dim>& tria,
		const dealii::FiniteElement<dim,dim>& fe){
	assert (fe.has_support_points());

	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = CFDSolverUtilities::getMinimumVertexDistance<
			dim>(tria);
	double unit_min = 100.0;
	const std::vector< dealii::Point< dim > > & support_points = fe.get_unit_support_points ();
	for (size_t i = 0; i < support_points.size(); i++){
		for (size_t j = i+1; j < support_points.size(); j++){
			double distance_ij = support_points.at(i).distance(support_points.at(j));
			if ((distance_ij < unit_min) and (distance_ij > 1e-10)){
				unit_min = distance_ij;
			}
		}
	}

	return min_vertex_distance * unit_min;
}

template<size_t dim>
double CFDSolverUtilities::getMinimumVertexDistance(const Mesh<dim>& tria) {
	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = 100000000000.0;
	double distance = 0.0;
	for (typename Mesh<dim>::active_cell_iterator cell = tria.begin_active();
			cell != tria.end(); ++cell) {
		if (cell->is_locally_owned()) {
			distance = cell->minimum_vertex_distance();
			if (distance < min_vertex_distance) {
				min_vertex_distance = distance;
			}
		}
	}
	// sync over all MPI processes
	return dealii::Utilities::MPI::min_max_avg(min_vertex_distance,
	MPI_COMM_WORLD).min;
}
template double CFDSolverUtilities::getMinimumVertexDistance<2>(
		const Mesh<2>& tria);
template double CFDSolverUtilities::getMinimumVertexDistance<3>(
		const Mesh<3>& tria);

template<size_t dim>
double CFDSolverUtilities::calculateTimestep(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement, const Stencil& stencil, double cFL) {
	assert(orderOfFiniteElement >= 1);
	double dx = CFDSolverUtilities::getMinimumVertexDistance<dim>(tria);
	double u = stencil.getMaxParticleVelocityMagnitude();
	// according to Hesthaven, dt ~ p^{-2}
	double dt = cFL * dx / (u * orderOfFiniteElement * orderOfFiniteElement);
	return dt;
}
template double CFDSolverUtilities::calculateTimestep<2>(const Mesh<2>& tria,
		const size_t orderOfFiniteElement, const Stencil& stencil, double cFL);
template double CFDSolverUtilities::calculateTimestep<3>(const Mesh<3>& tria,
		const size_t orderOfFiniteElement, const Stencil& stencil, double cFL);

template<int dim>
void CFDSolverUtilities::mesh_info(const Mesh<dim> &tria,
		const std::string &filename) {
	pout << "Mesh info:" << std::endl << " dimension: " << dim << std::endl
			<< " no. of cells: " << tria.n_active_cells() << std::endl;
	{
		std::map<unsigned int, unsigned int> boundary_count;
		typename Mesh<dim>::active_cell_iterator cell = tria.begin_active(),
				endc = tria.end();
		for (; cell != endc; ++cell) {
			for (unsigned int face = 0;
					face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
				if (cell->face(face)->at_boundary())
					boundary_count[cell->face(face)->boundary_id()]++;
			}
		}
		pout << " boundary indicators: ";
		for (std::map<unsigned int, unsigned int>::iterator it =
				boundary_count.begin(); it != boundary_count.end(); ++it) {
			pout << it->first << "(" << it->second << " times) ";
		}
		pout << endl;
	}
	std::ofstream out(filename.c_str());
	dealii::GridOut grid_out;
	grid_out.write_eps(tria, out);
	pout << " written to " << filename << std::endl << std::endl;
} /*mesh_info */
template void CFDSolverUtilities::mesh_info<2>(const Mesh<2> &tria,
		const std::string &filename);
template void CFDSolverUtilities::mesh_info<3>(const Mesh<3> &tria,
		const std::string &filename);

void CFDSolverUtilities::get_integrator_by_id(size_t id,
		TimeIntegratorName& time_integrator,
		DealIntegratorName& deal_integrator, std::string& integrator_name) {

	deal_integrator = NONE;
	switch (id) {
	case 1: {
		time_integrator = RUNGE_KUTTA_5STAGE;
		break;
	}

	case 2: {
		time_integrator = THETA_METHOD;
		break;
	}

	case 3:
		time_integrator = EXPONENTIAL;
		break;

	case 4: {
		time_integrator = OTHER;
		deal_integrator = FORWARD_EULER;
		break;
	}

	case 5: {
		time_integrator = OTHER;
		deal_integrator = RK_THIRD_ORDER;
		break;
	}

	case 6: {
		time_integrator = OTHER;
		deal_integrator = RK_CLASSIC_FOURTH_ORDER;
		break;
	}
	case 7: {
		time_integrator = OTHER;
		deal_integrator = BACKWARD_EULER;
		break;
	}
	case 8: {
		time_integrator = OTHER;
		deal_integrator = IMPLICIT_MIDPOINT;
		break;
	}
	case 9: {
		time_integrator = OTHER;
		deal_integrator = CRANK_NICOLSON;
		break;
	}
	case 10: {
		time_integrator = OTHER;
		deal_integrator = SDIRK_TWO_STAGES;
		break;
	}
	case 11: {
		time_integrator = OTHER;
		deal_integrator = HEUN_EULER;
		break;
	}
	case 12: {
		time_integrator = OTHER;
		deal_integrator = BOGACKI_SHAMPINE;
		break;
	}
	case 13: {
		time_integrator = OTHER;
		deal_integrator = DOPRI;
		break;
	}
	case 14: {
		time_integrator = OTHER;
		deal_integrator = FEHLBERG;
		break;
	}
	case 15: {
		time_integrator = OTHER;
		deal_integrator = CASH_KARP;
		break;
	}
	default: {
		LOG(ERROR)
				<< "Time integrator not set properly in CFDSolverUtilities::get_integrator_by_id()."
				<< endl;
		break;
	}
	}
	integrator_name.assign(get_integrator_name(time_integrator, deal_integrator));
} /* get_integrator_by_id */

std::string CFDSolverUtilities::get_integrator_name(
		const TimeIntegratorName& time_integrator,
		const DealIntegratorName& deal_integrator) {

	switch (time_integrator) {
	case RUNGE_KUTTA_5STAGE: {
		return "RUNGE_KUTTA_5STAGE";
	}
	case THETA_METHOD: {
		return "THETA_METHOD";
	}
	case EXPONENTIAL: {
		return "EXPONENTIAL";
	}
	case OTHER: {
		switch (deal_integrator) {
		case FORWARD_EULER: {
			return "FORWARD_EULER";
		}
		case RK_THIRD_ORDER: {
			return "RK_THIRD_ORDER";
		}
		case RK_CLASSIC_FOURTH_ORDER: {
			return "RK_CLASSIC_FOURTH_ORDER";
		}
		case BACKWARD_EULER: {
			return "BACKWARD_EULER";
		}
		case IMPLICIT_MIDPOINT: {
			return "IMPLICIT_MIDPOINT";
		}
		case CRANK_NICOLSON: {
			return "CRANK_NICOLSON";
		}
		case SDIRK_TWO_STAGES: {
			return "SDIRK_TWO_STAGES";
		}
		case HEUN_EULER: {
			return "HEUN_EULER";
		}
		case BOGACKI_SHAMPINE: {
			return "BOGACKI_SHAMPINE";
		}
		case DOPRI: {
			return "DOPRI";
		}
		case FEHLBERG: {
			return "FEHLBERG";
		}
		case CASH_KARP: {
			return "CASH_KARP";
		}
		case NONE: {
			return "";
		}
		}
		break;
	}
	default: {
		LOG(ERROR)
				<< "Time integrator not set properly in CFDSolverUtilities::get_integrator_by_id()."
				<< endl;
		break;
	}
	}
	return "";
} /* get_integrator_by_id */

boost::shared_ptr<Stencil> CFDSolverUtilities::make_stencil(size_t d, size_t q,
		size_t scaling) {
	{
		std::stringstream msg;
		msg << "Could not create stencil with d=" << d << ", q=" << q << ".";
		if (2 == d) {
			if (9 == q) {
				return boost::make_shared<D2Q9>(scaling);
			}
            if (25 == q) {
                return boost::make_shared<D2Q25H>(scaling);
            }
		} else if (3 == d) {
			if (15 == q) {
				return boost::make_shared<D3Q15>(scaling);
			}
			if (19 == q) {
				return boost::make_shared<D3Q19>(scaling);
			}
			if (27 == q) {
				return boost::make_shared<D3Q27>(scaling);
			}
            if (45 == q) {
                return boost::make_shared<D3Q45>(scaling);
            }
            if (13 == q) {
                return boost::make_shared<D3Q13>(scaling);
            }
            if (21 == q) {
                return boost::make_shared<D3Q21>(scaling);
            }
		}
		throw CFDSolverUtilitiesException(msg.str());
		return NULL;
	}
}

} /* namespace natrium */

