       program main
        use mol_descr
        use chem_parser
        use stsubs
        implicit none
        logical:: parsed
        type(mol_params_t):: mol_params
        real(8), allocatable:: moa(:,:),mob(:,:),sao(:,:),cisa(:,:,:),cisb(:,:,:)
        real(8), allocatable:: hole_dens(:,:,:),particle_dens(:,:,:)
        type(basis_func_info_t), allocatable:: basis_info(:)

        call test_chem_parser()

        parsed=orca_extract_mol_params('orca.dat',mol_params); if(.not.parsed) stop
        parsed=orca_extract_overlap('orca.dat',mol_params,sao); if(.not.parsed) stop
        parsed=orca_extract_mo_coef('orca.molden',mol_params,moa,mob); if(.not.parsed) stop
        parsed=orca_extract_cis_coef('orca.dat',mol_params,cisa,cisb); if(.not.parsed) stop
        parsed=orca_extract_basis_info('orca.molden',mol_params,basis_info); if(.not.parsed) stop

        call compute_transition_density(mol_params,basis_info,sao,moa,mob,cisa,cisb,hole_dens,particle_dens)
        if(allocated(hole_dens).and.allocated(particle_dens)) then
         write(*,'("Computed AO hole density for excited state 1:")')
         call wr_mat_dp(size(hole_dens,1),size(hole_dens,2),hole_dens(:,:,1))
         write(*,'("Computed AO particle density for excited state 1:")')
         call wr_mat_dp(size(particle_dens,1),size(particle_dens,2),particle_dens(:,:,1))
        else
         write(*,'("#ERROR: Undefined atomic state vectors!")'); stop
        endif

       end program main
