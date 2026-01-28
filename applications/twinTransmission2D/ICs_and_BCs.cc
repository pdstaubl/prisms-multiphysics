#include "customPDE.h"

// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE 
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

    // The initial condition is a set of overlapping circles/spheres defined
    // by a hyperbolic tangent function. The center of each circle/sphere is
    // given by "center" and its radius is given by "radius".

  double center[3] = {0.5,0.75,0.5};
  dealii::Tensor<1, dim> nc;
  dealii::Tensor<1, dim> n;
  dealii::Tensor<1, dim> n_int;
  dealii::Tensor<1, dim> td_int;
  dealii::Tensor<1, dim> tn_int;
  nc.clear(); n.clear(); n_int.clear(); td_int.clear(); tn_int.clear();
  double dist, edist;
  double b0=a0*std::sqrt((1.0-ecc*ecc));
  double nX, nY, nZ, nZ_reg;
  scalar_IC = 0;

  if (index==0){
    //Calculating distance from center of the system
    dist = 0.0;
    for (unsigned int dir = 0; dir < dim; dir++){
      dist += (p[dir]-center[dir]*userInputs_pf.domain_size[dir])*(p[dir]-center[dir]*userInputs_pf.domain_size[dir]);
    }
    dist = std::sqrt(dist);
    
    //Elliptical seed
    //Calculating distance from center to perimeter of the ellipse
    //Components of the normal vector from the center to point p
    nc[0] = (p[0]-center[0]*userInputs_pf.domain_size[0])/(dist + 1.0e-7);
    nc[1] = (p[1]-center[1]*userInputs_pf.domain_size[1])/(dist + 1.0e-7);
    nc[2] = (p[2]-center[2]*userInputs_pf.domain_size[2])/(dist + 1.0e-7);
	
    // Get the materialID (a.k.a. grainID) from CPFE for this mesh point
	  double coords[3] = {p[0], p[1], p[2]};
	  int materialID = this->cpfe_orientations->getMaterialID(coords);

    //Rotating td and tn according to parent grain orientation
    dealii::Tensor<1, dim> rot;
	  rot.clear();
	  rot[0] = this->cpfe_orientations->eulerAngles[materialID][0];
	  rot[1] = this->cpfe_orientations->eulerAngles[materialID][1];
	  rot[2] = this->cpfe_orientations->eulerAngles[materialID][2];
	  dealii::Tensor<2, dim> rotmat;
	  rotmat.clear();
	  rodrigues_to_rotmat(rotmat, rot);

    //Rotated normal vector, tn_int
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            tn_int[i] += rotmat[i][j] * tn[j];

    //Rotated direction vector, td_int
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            td_int[i] += rotmat[i][j] * td[j];

    //Rotated twin direction and twin normal vectors
    dealii::Tensor<1, dim> e_X = td_int;
    dealii::Tensor<1, dim> e_Y = tn_int;
    dealii::Tensor<1, dim> e_Z;

    // Compute e_Z = e_X cross e_Y
    e_Z[0] = e_X[1]*e_Y[2] - e_X[2]*e_Y[1];
    e_Z[1] = e_X[2]*e_Y[0] - e_X[0]*e_Y[2];
    e_Z[2] = e_X[0]*e_Y[1] - e_X[1]*e_Y[0];

    // Normalize e_Z
    double norm_eZ = std::sqrt(e_Z*e_Z);
    for (unsigned int i = 0; i < dim; ++i)
        e_Z[i] /= norm_eZ;

    // Construct rotation matrix Q as 3 column vectors
    dealii::Tensor<2, dim> Q;
    for (unsigned int i = 0; i < dim; ++i) {
        Q[i][0] = e_X[i];
        Q[i][1] = e_Y[i];
        Q[i][2] = e_Z[i];
    }
    //Intermediate vector, n_int
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            n[i] += Q[j][i] * nc[j];

    //Rotated unit vector with respect to the twin plane    
    nX = n[0];
    nY = n[1];
    nZ = n[2];

    //Distance from center of the ellipse to a point in the ellipse that intersects the line defined by (nX,nY,nZ)
    //Regularized nZ such that nX^2 + nY^2 + nZ^2 = 1
    nZ_reg = std::sqrt(1.0-nX*nX-nY*nY);
    edist = 1.0/std::sqrt((nX/a0)*(nX/a0) + (nY/b0)*(nY/b0) + (nZ_reg/a0)*(nZ_reg/a0));
                     
    scalar_IC =  0.5*(1.0-std::tanh((dist-edist)/(1.0*del0*std::sqrt(edist/a0))));

    if (scalar_IC > 1.0) scalar_IC = 1.0;
                     
  } else {
    scalar_IC = 0.0;
  }
  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
