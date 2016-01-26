/*
 * elasticity_3D_test.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ELASTICITY_3D_ELASTICITY_3D_TEST_H_
#define ELASTICITY_3D_ELASTICITY_3D_TEST_H_

#include "../common_headers_libmesh.h"
#include "assemble_functions_elasticity_3D.h"

// Classes

// ---> Border displacememnt
class border_displacement : public FunctionBase<Number>
{
	private:
		const unsigned int _u_var, _v_var, _w_var;
		const Real _x_displ, _y_displ, _z_displ;

	public:
		border_displacement (	unsigned int u_var,
								unsigned int v_var,
								unsigned int w_var,
								Real x_displ = 0,
								Real y_displ = 0,
								Real z_displ = 0)
								: _u_var(u_var),
								  _v_var(v_var),
								  _w_var(w_var),
								  _x_displ(x_displ),
								  _y_displ(y_displ),
								  _z_displ(z_displ)
		{
			this->_initialized = true;
		}

		virtual Number operator() (const Point&, const Real = 0)
		{
			libmesh_not_implemented();
		}

		virtual void operator() (	const Point& p,
									const Real,
									DenseVector<Number>& output)
		{
			output.resize(3);
			output.zero();

			output(_u_var) = _x_displ;
			output(_v_var) = _y_displ;
			output(_w_var) = _z_displ;
		}

		virtual AutoPtr<FunctionBase<Number> > clone() const
		{
			return AutoPtr<FunctionBase<Number> >
				(new border_displacement(	_u_var, _v_var, _w_var,
											_x_displ, _y_displ, _z_displ));
		}
};



#endif /* ELASTICITY_3D_ELASTICITY_3D_TEST_H_ */
