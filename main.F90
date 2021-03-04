       program main
        use mol_descr
        use chem_parser
        use stsubs
        implicit none
        logical:: parsed
        type(mol_params_t):: mol_params
        real(8), allocatable:: moa(:,:),mob(:,:),sao(:,:),cisa(:,:,:),cisb(:,:,:),cis_energy(:)
        real(8), allocatable:: hole_dens(:,:,:),particle_dens(:,:,:),hdiag(:),pdiag(:)
        type(atom_state_t), allocatable:: asv(:,:)
        type(basis_func_info_t), allocatable:: basis_info(:)
        integer:: i,nh,np,state
        real(8):: dnorm

        call test_chem_parser()

        !Extract necessary quantum-chemical data:
        parsed=orca_extract_mol_params('orca.dat',mol_params); if(.not.parsed) stop
        parsed=orca_extract_overlap('orca.dat',mol_params,sao); if(.not.parsed) stop
        parsed=orca_extract_mo_coef('orca.molden',mol_params,moa,mob); if(.not.parsed) stop
        parsed=orca_extract_cis_coef('orca.dat',mol_params,cis_energy,cisa,cisb); if(.not.parsed) stop
        parsed=orca_extract_basis_info('orca.molden',mol_params,basis_info); if(.not.parsed) stop

        !Compute excited state descriptors (atomic state vectors):
        call compute_transition_density(mol_params,basis_info,sao,moa,mob,cisa,cisb,hole_dens,particle_dens)
        if(allocated(hole_dens).and.allocated(particle_dens)) then
         state=2
         nh=size(hole_dens,1); np=size(particle_dens,1)
         allocate(hdiag(nh)); allocate(pdiag(np))
         do i=1,nh; hdiag(i)=hole_dens(i,i,state); enddo
         do i=1,np; pdiag(i)=particle_dens(i,i,state); enddo
         !write(*,'("Computed AO hole density matrix for excited state:")')
         !call wr_mat_dp(nh,nh,hole_dens(:,:,state))
         !write(*,'("Computed AO particle density matrix for excited state:")')
         !call wr_mat_dp(np,np,particle_dens(:,:,state))
         write(*,'("The AO hole density for excited state:")')
         dnorm=0d0
         do i=1,nh
          write(*,*) i,basis_info(i)%atom_id,basis_info(i)%shell_id,basis_info(i)%ao_id,hdiag(i)
          dnorm=dnorm+hdiag(i)
         enddo
         write(*,'("The AO hole density trace = ",D25.14)') dnorm
         write(*,'("The AO particle density for excited state:")')
         dnorm=0d0
         do i=1,np
          write(*,*) i,basis_info(i)%atom_id,basis_info(i)%shell_id,basis_info(i)%ao_id,pdiag(i)
          dnorm=dnorm+pdiag(i)
         enddo
         write(*,'("The AO particle density trace = ",D25.14)') dnorm
         call compute_atomic_state_vectors(mol_params,basis_info,hole_dens,particle_dens,asv)
         !call print_atomic_state_vectors(asv)
        else
         write(*,'("#ERROR: Undefined hole and/or particle density matrices!")'); stop
        endif

       end program main
