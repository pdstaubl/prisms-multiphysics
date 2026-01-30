#include "customPDE.h"

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0 - Order Parameter
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, dndt, strain_df");
    set_dependencies_gradient_term_RHS(0, "grad(n)");
	
	// Variable 1 - Time Derivative of Order Parameter
	set_variable_name				(1,"dndt");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

    set_dependencies_value_term_RHS(1, "n, dndt, strain_df");
    set_dependencies_gradient_term_RHS(1, "grad(n)");

	// Variable 2 - Strain driving force
	set_variable_name				(2,"strain_df");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

    set_dependencies_value_term_RHS(2, "strain_df");
    set_dependencies_gradient_term_RHS(2, "");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double>> & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {

	// --- Getting the values and derivatives of the model variables ---

	// The order parameter and its derivatives
	scalarvalueType_pf n = variable_list.get_scalar_value(0);
	scalargradType_pf nx = variable_list.get_scalar_gradient(0);

	// The time derivative of the order parameter
	scalarvalueType_pf dndt = variable_list.get_scalar_value(1);

	// The strain contribiution to the driving force
	scalarvalueType_pf strain_df= variable_list.get_scalar_value(2);

	scalarvalueType_pf mu_twV = constV(delf_tw*1.5)*(4.0*n*(n-1.0)*(n-0.5));
	scalargradType_pf kappagradn;
	scalarvalueType_pf L = constV(0.0);
	
	//Outward Normal vector
	scalargradType_pf nvec = -nx/(std::sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2])+constV(regval));


	// q_point_loc is vectorized, but CPFE functions are not.
	// TODO: create vectorized versions of the relevant CPFE functions, to fix this

	// For now, unroll the vectorization.
	for (unsigned int v = 0; v < q_point_loc[0].size(); v++) {

		// Get the materialID (a.k.a. grainID) from CPFE for this mesh point
		double coords[3] = {q_point_loc[0][v], q_point_loc[1][v], q_point_loc[2][v]};
		unsigned int materialID = this->cpfe_orientations->getMaterialID(coords);

		// Get the crystal orientation as a Rodrigues vector from CPFE
		// (note: although the variable is named euelrAngles, it's actually Rodrigues vectors)
		// TODO: rename the CPFE orientations variable and possibly refactor that code
		dealii::Tensor<1, dim> rot;
		rot.clear();
		rot[0] = this->cpfe_orientations->eulerAngles[materialID][0];
		rot[1] = this->cpfe_orientations->eulerAngles[materialID][1];
		rot[2] = this->cpfe_orientations->eulerAngles[materialID][2];
		dealii::Tensor<2, dim> rotmat;
		rotmat.clear();
		rodrigues_to_rotmat(rotmat, rot);

		// Get the twin direction and twin normal in the crystal frame (for this variant)
		// TODO: get td and tn from the CPFE input files instead of parameters_pf.prm
		//       This will prevent accidental mismatch and reduce the effort of
		//       adding additional twin variants
		dealii::Tensor<1, dim> e_X = td;
		dealii::Tensor<1, dim> e_Y = tn;
		dealii::Tensor<1, dim> e_Z;

		// Compute e_Z = e_X cross e_Y
		e_Z[0] = e_X[1]*e_Y[2] - e_X[2]*e_Y[1];
		e_Z[1] = e_X[2]*e_Y[0] - e_X[0]*e_Y[2];
		e_Z[2] = e_X[0]*e_Y[1] - e_X[1]*e_Y[0];

		// Normalize e_Z
		double norm_eZ = std::sqrt(e_Z*e_Z);
		for (unsigned int i = 0; i < dim; ++i)
				e_Z[i] /= norm_eZ;
		
		// Compute Q, the rotation from twin frame to crystal frame
		dealii::Tensor<2, dim> Q;
		for (unsigned int i = 0; i < dim; ++i) {
				Q[i][0] = e_X[i];
				Q[i][1] = e_Y[i];
				Q[i][2] = e_Z[i];
		}

		// Compute Q^T
		dealii::Tensor<2, dim> Q_T;
		for (unsigned int i = 0; i < dim; ++i)
				for (unsigned int j = 0; j < dim; ++j)
						Q_T[i][j] = Q[j][i];

		// Rotate Kij and Lij from the twin frame to the crystal frame
		//  Ltens_ccref = Q * Lij_tp * Q^T
		dealii::Tensor<2, dim> temp1;
		dealii::Tensor<2, dim> Ltens_ccref;
		dealii::Tensor<2, dim> temp2;
		dealii::Tensor<2, dim> K_ccref;
		temp1.clear();
		temp2.clear();
		Ltens_ccref.clear();
		K_ccref.clear();
		for (unsigned int i = 0; i < dim; ++i)
				for (unsigned int j = 0; j < dim; ++j)
						for (unsigned int k = 0; k < dim; ++k)
							{
								temp1[i][j] += Q[i][k] * Lij_tp[k][j];
								temp2[i][j] += Q[i][k] * Kij_tp[k][j];
							}

		for (unsigned int i = 0; i < dim; ++i)
				for (unsigned int j = 0; j < dim; ++j)
						for (unsigned int k = 0; k < dim; ++k)
							{
								Ltens_ccref[i][j] += temp1[i][k] * Q_T[k][j];
								K_ccref[i][j]     += temp2[i][k] * Q_T[k][j];
							}

		// Rotate Kij and Lij from the crystal frame to the sample frame
		dealii::Tensor<2, dim> Ltens, K;
		Ltens.clear();
		K.clear();

		for (unsigned int i = 0; i < dim; i++)
				for (unsigned int j = 0; j < dim; j++)
						for (unsigned int k = 0; k < dim; k++)
								for (unsigned int a = 0; a < dim; a++)
									{
										// K_ij = R * K' * R^T = R_ik K'_ka R_ja
										K[i][j]     += rotmat[i][k]*K_ccref[k][a]*rotmat[j][a];
										Ltens[i][j] += rotmat[i][k]*Ltens_ccref[k][a]*rotmat[j][a];
									}

		// Kij and Lij are now in the sample frame at the current quad point.
		// TODO: move some of the above calculation to initialization, caching the
		//       values in order to save computation time at the expense of memory

		// --- Setting the expressions for the terms in the governing equations ---

		kappagradn[0][v] = K[0][0]*nx[0][v] + K[0][1]*nx[1][v] + K[0][2]*nx[2][v];
		kappagradn[1][v] = K[1][0]*nx[0][v] + K[1][1]*nx[1][v] + K[1][2]*nx[2][v];
		kappagradn[2][v] = K[2][0]*nx[0][v] + K[2][1]*nx[1][v] + K[2][2]*nx[2][v];

		//Computing the outward mobility (L = grad(nvec) dot Ltens dot grad(nvec))
		for(unsigned int i=0; i < dim; i++) {
			for(unsigned int j=0; j < dim; j++) {
				//Mobility tensor (rotated)
				L[v] += nvec[i][v]*nvec[j][v]*Ltens[i][j];
			}
    }

} // End vectorization unroll

//Applying a filter to localize driving force to the twin boundary 
scalarvalueType_pf strain_df_filter = 1.5*(1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

//Defining the value and gradient terms
scalarvalueType_pf eq_n = (n-constV(userInputs_pf.dtValue)*L*(mu_twV-strain_df_filter));
scalargradType_pf eqx_n = -(constV(userInputs_pf.dtValue)*L*kappagradn);

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(0,eq_n);
variable_list.set_scalar_gradient_term_RHS(0,eqx_n);

}

// Copied from crystalPlasticity<dim>::odfpoint()
template <int dim, int degree>
void customPDE<dim,degree>::rodrigues_to_rotmat(dealii::Tensor<2, dim> &OrientationMatrix, dealii::Tensor<1, dim> r) const
{
    double rdotr = 0.0;

    for(unsigned int i = 0; i < dim; i++){
        rdotr = rdotr + r[i]*r[i];
    }

    double term1 = 1.0 - rdotr;
    double term2 = 1.0 + rdotr;

    OrientationMatrix.clear();
		OrientationMatrix[0][0] = 1.0;
		OrientationMatrix[1][1] = 1.0;
		OrientationMatrix[2][2] = 1.0;

    for(unsigned int i = 0; i < dim; i++)
		  {
        OrientationMatrix[i][i] = OrientationMatrix[i][i]*term1;
      }

    for(unsigned int i = 0; i < dim; i++)
		  {
        for(unsigned int j = 0; j < dim; j++)
				  {
            OrientationMatrix[i][j] = OrientationMatrix[i][j] + 2.0*r[i]*r[j];
          }
      }

    OrientationMatrix[0][1] = OrientationMatrix[0][1] - 2.0*r[2];
    OrientationMatrix[0][2] = OrientationMatrix[0][2] + 2.0*r[1];
    OrientationMatrix[1][2] = OrientationMatrix[1][2] - 2.0*r[0];
    OrientationMatrix[1][0] = OrientationMatrix[1][0] + 2.0*r[2];
    OrientationMatrix[2][0] = OrientationMatrix[2][0] - 2.0*r[1];
    OrientationMatrix[2][1] = OrientationMatrix[2][1] + 2.0*r[0];

    for(unsigned int i = 0; i < dim; i++)
		  {
        for(unsigned int j = 0; j < dim; j++)
				  {
            OrientationMatrix[i][j] = OrientationMatrix[i][j]*1.0/term2;
          }
      }
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
// --- Getting the values and derivatives of the model variables ---

// The order parameter and its derivatives
scalarvalueType_pf n = variable_list.get_scalar_value(0);
scalargradType_pf nx = variable_list.get_scalar_gradient(0);

// The time derivative of the order parameter
scalarvalueType_pf dndt = variable_list.get_scalar_value(1);

// The strain contribiution to the driving force
scalarvalueType_pf strain_df= variable_list.get_scalar_value(2);

scalarvalueType_pf mu_twV = constV(delf_tw)*(4.0*n*(n-1.0)*(n-0.5));

//Outward Normal vector
scalargradType_pf nvec = -nx/(std::sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2])+constV(regval));
scalargradType_pf kappagradn;
scalarvalueType_pf L = constV(0.0);

// q_point_loc is vectorized, but CPFE functions are not.
// TODO: create vectorized versions of the relevant CPFE functions, to fix this

// For now, unroll the vectorization.
for (unsigned int v = 0; v < q_point_loc[0].size(); v++) {

	// Get the materialID (a.k.a. grainID) from CPFE for this mesh point
	double coords[3] = {q_point_loc[0][v], q_point_loc[1][v], q_point_loc[2][v]};
	unsigned int materialID = this->cpfe_orientations->getMaterialID(coords);

	// Get the crystal orientation as a Rodrigues vector from CPFE
	// (note: although the variable is named euelrAngles, it's actually Rodrigues vectors)
	// TODO: rename the CPFE orientations variable and possibly refactor that code
	dealii::Tensor<1, dim> rot;
	rot.clear();
	rot[0] = this->cpfe_orientations->eulerAngles[materialID][0];
	rot[1] = this->cpfe_orientations->eulerAngles[materialID][1];
	rot[2] = this->cpfe_orientations->eulerAngles[materialID][2];
	dealii::Tensor<2, dim> rotmat;
	rotmat.clear();
	rodrigues_to_rotmat(rotmat, rot);

	// Get the twin direction and twin normal in the crystal frame (for this variant)
	// TODO: get td and tn from the CPFE input files instead of parameters_pf.prm
	//       This will prevent accidental mismatch and reduce the effort of
	//       adding additional twin variants
	dealii::Tensor<1, dim> e_X = td;
	dealii::Tensor<1, dim> e_Y = tn;
	dealii::Tensor<1, dim> e_Z;

	// Compute e_Z = e_X cross e_Y
	e_Z[0] = e_X[1]*e_Y[2] - e_X[2]*e_Y[1];
	e_Z[1] = e_X[2]*e_Y[0] - e_X[0]*e_Y[2];
	e_Z[2] = e_X[0]*e_Y[1] - e_X[1]*e_Y[0];

	// Normalize e_Z
	double norm_eZ = std::sqrt(e_Z*e_Z);
	for (unsigned int i = 0; i < dim; ++i)
			e_Z[i] /= norm_eZ;
	
	// Compute Q, the rotation from twin frame to crystal frame
	dealii::Tensor<2, dim> Q;
	for (unsigned int i = 0; i < dim; ++i) {
			Q[i][0] = e_X[i];
			Q[i][1] = e_Y[i];
			Q[i][2] = e_Z[i];
	}

	// Compute Q^T
	dealii::Tensor<2, dim> Q_T;
	for (unsigned int i = 0; i < dim; ++i)
			for (unsigned int j = 0; j < dim; ++j)
					Q_T[i][j] = Q[j][i];

	// Rotate Kij and Lij from the twin frame to the crystal frame
	//  Ltens_ccref = Q * Lij_tp * Q^T
	dealii::Tensor<2, dim> temp1;
	dealii::Tensor<2, dim> Ltens_ccref;
	dealii::Tensor<2, dim> temp2;
	dealii::Tensor<2, dim> K_ccref;
	temp1.clear();
	temp2.clear();
	Ltens_ccref.clear();
	K_ccref.clear();
	for (unsigned int i = 0; i < dim; ++i)
			for (unsigned int j = 0; j < dim; ++j)
					for (unsigned int k = 0; k < dim; ++k)
					  {
							temp1[i][j] += Q[i][k] * Lij_tp[k][j];
							temp2[i][j] += Q[i][k] * Kij_tp[k][j];
						}

	for (unsigned int i = 0; i < dim; ++i)
			for (unsigned int j = 0; j < dim; ++j)
					for (unsigned int k = 0; k < dim; ++k)
					  {
							Ltens_ccref[i][j] += temp1[i][k] * Q_T[k][j];
							K_ccref[i][j]     += temp2[i][k] * Q_T[k][j];
						}

	// Rotate Kij and Lij from the crystal frame to the sample frame
	dealii::Tensor<2, dim> Ltens, K;
	Ltens.clear();
	K.clear();

	for (unsigned int i = 0; i < dim; i++)
	    for (unsigned int j = 0; j < dim; j++)
		      for (unsigned int k = 0; k < dim; k++)
				      for (unsigned int a = 0; a < dim; a++)
								{
									// K_ij = R * K' * R^T = R_ik K'_ka R_ja
									K[i][j]     += rotmat[i][k]*K_ccref[k][a]*rotmat[j][a];
									Ltens[i][j] += rotmat[i][k]*Ltens_ccref[k][a]*rotmat[j][a];
								}

	// --- Setting the expressions for the terms in the governing equations ---

	kappagradn[0][v] = K[0][0]*nx[0][v] + K[0][1]*nx[1][v] + K[0][2]*nx[2][v];
	kappagradn[1][v] = K[1][0]*nx[0][v] + K[1][1]*nx[1][v] + K[1][2]*nx[2][v];
	kappagradn[2][v] = K[2][0]*nx[0][v] + K[2][1]*nx[1][v] + K[2][2]*nx[2][v];

	//Computing the outward mobility (L = grad(nvec) dot Ltens dot grad(nvec))
	for(unsigned int i = 0; i < dim; i++) {
		for(unsigned int j = 0; j < dim; j++) {
			//Mobility tensor (rotated)
			L[v] += nvec[i][v]*nvec[j][v]*Ltens[i][j];;
		}
	}

} // end vectorization unroll

//Applying a filter to localize driving force to the twin boundary 
scalarvalueType_pf strain_df_filter = 1.5*(1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

scalarvalueType_pf eq_dndt = -L*(mu_twV-strain_df_filter);
scalargradType_pf eqx_dndt = -L*kappagradn;

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(1,eq_dndt);
variable_list.set_scalar_gradient_term_RHS(1,eqx_dndt);

variable_list.set_scalar_value_term_RHS(2,strain_df);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
