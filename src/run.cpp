#include "header.hpp"
#include "model.hpp"
#include "derivatives.hpp"
#include "tools.hpp"

using namespace std;

void Model::Pre()
{
  // we make the system relax (without activity)
  if(relax_time>0)
  {

    double save_alpha  = 0; swap(alpha,  save_alpha);
    double save_zetaS  = 0; swap(zetaS,  save_zetaS);
    double save_zetaQ  = 0; swap(zetaQ,  save_zetaQ);
    double save_Dnem  = 0; swap(Dnem,  save_Dnem);
    double save_Dpol  = 0; swap(Dpol,  save_Dpol);
    double save_Jnem  = 0; swap(Jnem,  save_Jnem);
    double save_Jpol  = 0; swap(Jpol,  save_Jpol);
    double save_Kpol  = 0; swap(Kpol,  save_Kpol);
    double save_Knem  = 0; swap(Knem,  save_Knem);
    double save_Wnem  = 0; swap(Wnem,  save_Wnem);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned j=0; j<=npc; ++j) Update(j==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(alpha, save_alpha);
    swap(zetaS, save_zetaS);
    swap(zetaQ, save_zetaQ);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Kpol, save_Kpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);

  }

  if(BC==5 || BC==7) ConfigureWalls(1);
  if(BC==6) ConfigureWalls(0);
}

void Model::Post()
{}

void Model::PreRunStats()
{}

void Model::RuntimeStats()
{}

void Model::RuntimeChecks()
{}

void Model::UpdateSumsAtNode(unsigned n, unsigned q)
{
  const auto k = GetIndexFromPatch(n, q);
  const auto p = phi[n][q];

  sum[k]     += p;
  square[k]  += p*p;
  
  sumS00[k]  += p*S00[n];
  sumS01[k]  += p*S01[n];
  sumS02[k]  += p*S02[n];
  sumS12[k]  += p*S12[n];
  sumS11[k]  += p*S11[n];
  sumS22[k]  += p*S22[n];
  
  //sumQ00[k]  += p*Q00[n];
  //sumQ01[k]  += p*Q01[n];
  sumQ00[k] += 1;
  sumQ01[k] += 0;
    //if (walls[k] == 1.){
    //const double z = GetZPosition(k);
    //if(z < wall_thickness + 10){
    //sumQ00[k]  += 1;
    //sumQ01[k]  += 0;
    //}
  
  P0[k]      += p*polarization[n][0];
  P1[k]      += p*polarization[n][1];
  P2[k]      += p*polarization[n][2];
  
  U0[k]      += p*velocity[n][0];
  U1[k]      += p*velocity[n][1];
  U2[k]      += p*velocity[n][2];
}







void Model::UpdatePotAtNode(unsigned n, unsigned q)
{
  const auto  k  = GetIndexFromPatch(n, q);
  const auto& s  = neighbors[k];
  const auto& sq = neighbors_patch[q];
  const auto p  = phi[n][q];
  const auto a  = vol[n];
  const auto ll = laplacian(phi[n], sq);
  const auto ls = laplacian(sum, s);

  const double internal = (
      // CH term
      + gam*(8*p*(1-p)*(1-2*p)/lambda - 2*lambda*ll)
      // vol conservation term
      - 4*mu/V0*(1-a/V0)*p
    );
  
  // kappa = kappas[n];
  const double interactions = (
      // repulsion term
      
      + 2*kappa_cc/lambda*p*(square[k]-p*p)
      // adhesion term
      - 2*omega_cc*lambda*(ls-ll)
      // repulsion with walls
      + 2*kappa_cs/lambda*p*walls[k]*walls[k]
      // adhesion with walls
      // - 2*SubAdhesion_field(x,y)*lambda*walls_laplace[k]
      - 2*omega_cs*lambda*walls_laplace[k]

    );
  
  // delta F / delta phi_i
  V[n][q] = internal + interactions;
  // pressure
  pressure[k] += p*interactions;

}













void Model::UpdateForcesAtNode(unsigned n, unsigned q)
{
  const auto  k  = GetIndexFromPatch(n, q);
  const auto& s  = neighbors[k];
  const auto& sq = neighbors_patch[q];

  const auto dx  = derivX(phi[n], sq);
  const auto dy  = derivY(phi[n], sq);
  const auto dz  = derivZ(phi[n], sq);
  
  const auto dxs = derivX(sum, s);
  const auto dys = derivY(sum, s);
  const auto dzs = derivZ(sum, s);
    
  Fpressure[n] += { pressure[k]*dx, pressure[k]*dy, pressure[k]*dz };
  
  Fshape[n]    += { zetaS*sumS00[k]*dx + zetaS*sumS01[k]*dy + zetaS*sumS02[k]*dz,
                    zetaS*sumS01[k]*dx + zetaS*sumS11[k]*dy + zetaS*sumS12[k]*dz,
                    zetaS*sumS02[k]*dx + zetaS*sumS12[k]*dy + zetaS*sumS22[k]*dz };
                    
  Fnem[n]  	 += { walls[k]*zetaQ*sumQ00[k]*dx + walls[k]*zetaQ*sumQ01[k]*dy , 	      
  		      walls[k]*zetaQ*sumQ01[k]*dx- walls[k]*zetaQ*sumQ00[k]*dy, 
  		     -(0.5)*walls[k]*zetaQ*dz };  
  		     
  // store derivatives
  phi_dx[n][q] = dx;
  phi_dy[n][q] = dy;
  phi_dz[n][q] = dz;

  // nematic torques
  tau[n]       += phi[n][q] * (sumQ00[k]*Q01[n] - sumQ01[k]*Q00[n]);
  vorticity[n] -= { U2[k]*dy-U1[k]*dz, U0[k]*dz-U2[k]*dx, U2[k]*dx-U0[k]*dy };//got to check the sign
  
  // polarisation torques (not super nice)
  const double ovlap = -( dx*(dxs-dx) + dy*(dys-dy) + dz*(dzs-dz)  );
  const vec<double, 3> P = { P0[k]-phi[n][q]*polarization[n][0], P1[k]-phi[n][q]*polarization[n][1], P2[k]-phi[n][q]*polarization[n][2] };
  
  delta_theta_pol[n] += ovlap*atan2( 
    sqrt(pow( (P[1]*polarization[n][0]-P[0]*polarization[n][1]) ,2) + pow( (P[2]*polarization[n][0]-P[0]*polarization[n][2]) ,2) + pow( (P[2]*polarization[n][1]-P[1]*polarization[n][2]) ,2) ),
    P[0]*polarization[n][0]+P[1]*polarization[n][1]+P[2]*polarization[n][2]                               
                                   );
  
}












void Model::UpdatePhaseFieldAtNode(unsigned n, unsigned q, bool store)
{
  const auto k = GetIndexFromPatch(n, q);

  // compute dphi
  dphi[n][q] =
    // free energy
    - V[n][q]
    // advection term
    - velocity[n][0]*phi_dx[n][q] - velocity[n][1]*phi_dy[n][q] - velocity[n][2]*phi_dz[n][q];
    ;

  // store values
  if(store)
  {
    dphi_old[n][q] = dphi[n][q];
    phi_old[n][q]  = phi[n][q];
  }
  
  // predictor-corrector
  {
    double p = phi_old[n][q]
               + time_step*.5*(dphi[n][q] + dphi_old[n][q]);

    // update for next call
    phi[n][q]    = p;

    com_x[n] += com_x_table[GetXPosition(k)]*p;
    com_y[n] += com_y_table[GetYPosition(k)]*p;
    com_z[n] += com_z_table[GetZPosition(k)]*p;
    vol[n]  += p*p;
 
  }
  
  // reinit values: we do reinit values here for the simple reason that it is
  // faster than having a supplementary loop afterwards. There is a race
  // condition in principle here but since we are setting evth back to 0 it
  // should be fine. Note that this should be done before the patches are updated
  ReinitSumsAtNode(k);
}

void Model::UpdateNematic(unsigned n, bool store)
{
  // euler-marijuana update
  if(store)
    theta_nem_old[n] = theta_nem[n] + sqrt_time_step*Dnem*random_normal();

  double F00 = 0, F01 = 0;
  switch(align_nematic_to)
  {
    case 0:
    {
      const auto ff = velocity[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 1:
    {
      const auto ff = Fpressure[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 2:
      F00 = S00[n];
      F01 = S01[n];
      break;
  }
  const auto strength = pow(F01*F01 + F00*F00, 0.25);

  theta_nem[n] = theta_nem_old[n] - time_step*(
      + Knem*tau[n]
      + Jnem*strength*atan2( F00*Q01[n]-F01*Q00[n], F00*Q00[n]+F01*Q01[n]) );
      // + Wnem*vorticity[n]; need to update for 3D model 
  
  Q00[n] = Snem*cos(2*theta_nem[n]);
  Q01[n] = Snem*sin(2*theta_nem[n]);
}

void Model::UpdatePolarization(unsigned n, bool store)
{
  // euler-marijuana update
  if(store)
    theta_pol_old[n] = theta_pol[n] + sqrt_time_step*Dpol*random_normal();

  vec<double, 3> ff = {0, 0, 0};
  switch(align_polarization_to)
  {
    case 0:
      ff = velocity[n];
      break;
    case 1:
      ff = Fpressure[n];
      break;
  }

  theta_pol[n] = theta_pol_old[n] - time_step*(
      + Kpol*delta_theta_pol[n]
      + Jpol*ff.abs() * atan2( 
    sqrt(pow( (ff[1]*polarization[n][0]-ff[0]*polarization[n][1]) ,2) + pow( (ff[2]*polarization[n][0]-ff[0]*polarization[n][2]) ,2) + pow( (ff[2]*polarization[n][1]-ff[1]*polarization[n][2]) ,2) ),
    ff[0]*polarization[n][0]+ff[1]*polarization[n][1]+ff[2]*polarization[n][2]                               
                             ));
  
  polarization[n] = { Spol*cos(theta_pol[n]), Spol*sin(theta_pol[n]) };
}






void Model::ComputeCoM(unsigned n)
{
  // the strategy to deal with the periodic boundary conditions is to compute
  // all the integrals in Fourier space and come back at the end. This way the
  // periodicity of the domain is automatically taken into account.
  const auto mx = arg(com_x[n]/static_cast<double>(N)) + Pi;
  const auto my = arg(com_y[n]/static_cast<double>(N)) + Pi;
  const auto mz = arg(com_z[n]/static_cast<double>(N)) + Pi;
  com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] , mz/2./Pi*Size[2] };
}

void Model::UpdatePatch(unsigned n)
{
  // obtain the new location of the patch min and max
  const coord com_grd { unsigned(round(com[n][0])), unsigned(round(com[n][1])), unsigned(round(com[n][2])) };
  const coord new_min = ( com_grd + Size - patch_margin ) % Size;
  const coord new_max = ( com_grd + patch_margin - coord {1u, 1u} ) % Size;
  coord displacement  = ( Size + new_min - patch_min[n] ) % Size;

  // I guess there is somehthing better than this...
  if(displacement[0]==Size[0]-1u) displacement[0] = patch_size[0]-1u;
  if(displacement[1]==Size[1]-1u) displacement[1] = patch_size[1]-1u;
  if(displacement[2]==Size[2]-1u) displacement[2] = patch_size[2]-1u;  

  // update offset and patch location
  offset[n]    = ( offset[n] + patch_size - displacement ) % patch_size;
  patch_min[n] = new_min;
  patch_max[n] = new_max;
  
}

void Model::UpdateStructureTensorAtNode(unsigned n, unsigned q)
{
  const auto  dx = phi_dx[n][q];
  const auto  dy = phi_dy[n][q];
  const auto  dz = phi_dz[n][q];

  S00[n] += -(dx*dx - (1./3.)*(dx*dx+dy*dy+dz*dz));
  S01[n] += -dx*dy;
  S02[n] += -dx*dz;
  S12[n] += -dy*dz;
  S11[n] += -(dy*dy - (1./3.)*(dx*dx+dy*dy+dz*dz));
  S22[n] += -(dz*dz - (1./3.)*(dx*dx+dy*dy+dz*dz));

}



void Model::ReinitSumsAtNode(unsigned k)
{
  sum[k] = 0;
  square[k] = 0;

  sumS00[k] = 0;
  sumS01[k] = 0;
  sumS02[k] = 0;
  sumS12[k] = 0;
  sumS11[k] = 0;
  sumS22[k] = 0;
  
  sumQ00[k] = 0;
  sumQ01[k] = 0;
  pressure[k] = 0;
  
  U0[k] = 0;
  U1[k] = 0;
  U2[k] = 0;
  
}

void Model::Update(bool store, unsigned nstart)
{
  // Compute all global sums
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdateSumsAtNode(n, q);
  }

  // Compute stresses
  // We need another loop because the pressure involves a double sum over all
  // the cells.
  
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdatePotAtNode(n, q);
  }
  
  // Compute induced force and passive velocity
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    Fpol[n] = Fshape[n] = Fpressure[n] = Fnem[n] = vorticity[n] = {0, 0, 0};
    // Fpol[n] = Fshape[n] = Fpressure[n] = vorticity[n] = {0, 0, 0};    
    delta_theta_pol[n] = tau[n] = 0;

    // update in restricted patch only
    for(unsigned q=0; q<patch_N; ++q){  
      UpdateForcesAtNode(n, q);
    }

    // normalise and compute total forces and vel
    tau[n]     /= lambda;
    Fpol[n]     = alpha*polarization[n];
    velocity[n] = ( Fpressure[n] + Fnem[n] + Fshape[n] + Fpol[n] )/xi;
  }
  
  // Predictor-corrector function for updating the phase fields
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    com_x[n] = com_y[n] = com_z[n] = S00[n] = S01[n] = S02[n] = S12[n] = S11[n] = S22[n] = vol[n] = 0;

    // only update fields in the restricted patch of field n
    for(unsigned q=0; q<patch_N; ++q)
    {
     UpdatePhaseFieldAtNode(n, q, store);
     UpdateStructureTensorAtNode(n, q);
    }

    // update polarisation
    UpdatePolarization(n, store);
    // update Q-tensor
    UpdateNematic(n, store);
    // update center of mass
    ComputeCoM(n);
    // update patch boundaries
    UpdatePatch(n);
  }
  
}
