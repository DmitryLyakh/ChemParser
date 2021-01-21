!Molecular descriptor generator
!Author: Dmitry I. Lyakh
       module mol_descr
       use chem_parser
       implicit none

       public compute_transition_density

       contains

       subroutine compute_transition_density(mol_params,overlap,cis_a,cis_b,hole_dens,particle_dens)
        type(mol_params_t), intent(in):: mol_params
        real(8), intent(in):: overlap(:,:)
        real(8), intent(in):: cis_a(:,:),cis_b(:,:)
        real(8), allocatable, intent(out):: hole_dens(:,:),particle_dens(:,:)

        return
       end subroutine compute_transition_density

       end module mol_descr
