!Molecular descriptor generator (atomic state vector for ANI model)
!Author: Dmitry I. Lyakh
       module mol_descr
       use chem_parser
       use stsubs
       implicit none

       type, public:: atom_state_t
        integer:: atom_id=0                           !atom id: [1..M]
        integer:: num_shells=0                        !number of shells: L
        real(8):: hole_density(0:MAX_AO_SHELLS-1)     !hole density [0..L-1]
        real(8):: particle_density(0:MAX_AO_SHELLS-1) !particle density [0..L-1]
       end type atom_state_t

       public compute_transition_density
       public compute_atomic_state_vectors
       public print_atomic_state_vectors

       private filter_small_elems

       contains

       subroutine compute_transition_density(mol_params,basis_info,overlap,mo_a,mo_b,cis_a,cis_b,hole_dens,particle_dens)
        type(mol_params_t), intent(in):: mol_params
        type(basis_func_info_t), intent(in):: basis_info(:) ![1:AO]
        real(8), intent(in):: overlap(:,:)                  ![1:AO,1:AO]
        real(8), intent(in):: mo_a(:,:),mo_b(:,:)           ![1:AO,1:MO]
        real(8), intent(inout):: cis_a(:,:,:),cis_b(:,:,:)  ![1:OCC,1:VIRT,1:STATES]
        real(8), allocatable, intent(out):: hole_dens(:,:,:),particle_dens(:,:,:) ![1:AO,1:AO,1:STATES]
        real(8), allocatable:: hdm(:,:),pdm(:,:),hao(:,:),pav(:,:),maa(:,:),mtx(:,:),mmm(:,:)
        integer:: num_occ,num_virt,num_cis_states,i,j,k,l,m,n
        real(8):: norm_a,norm_b,metrics_ao,metrics_mo,norm

        write(*,'("Computing hole/particle density matrices:")')
 !Check input:
        write(*,'(" Checking input ... ")',ADVANCE='NO')
        if(mol_params%num_ao_orbitals.le.0.or.&
          &mol_params%num_mo_orbitals.le.0.or.&
          &mol_params%num_electrons_a.le.0.or.&
          &mol_params%num_electrons_b.le.0) then
         write(*,'("#ERROR(compute_transition_density): Invalid molecular parameters!")'); stop
        endif
        if(mol_params%num_electrons_a.ne.mol_params%num_electrons_b) then
         write(*,'("#ERROR(compute_transition_density): Spin-polarized molecules are not supported yet!")'); stop
        endif
        num_occ=mol_params%num_electrons_a; num_virt=mol_params%num_mo_orbitals-num_occ
        if(size(basis_info,1).ne.mol_params%num_ao_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid basis info size!")'); stop
        endif
        if(size(overlap,1).ne.mol_params%num_ao_orbitals.or.size(overlap,2).ne.mol_params%num_ao_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid overlap matrix shape!")'); stop
        endif
        if(size(mo_a,1).ne.mol_params%num_ao_orbitals.or.size(mo_a,2).ne.mol_params%num_mo_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid shape of the MO coefficients alpha matrix!")'); stop
        endif
        if(size(mo_b,1).ne.mol_params%num_ao_orbitals.or.size(mo_b,2).ne.mol_params%num_mo_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid shape of the MO coefficients beta matrix!")'); stop
        endif
        num_cis_states=size(cis_a,3)
        if(size(cis_a,1).ne.num_occ.or.size(cis_a,2).ne.num_virt) then
         write(*,'("#ERROR(compute_transition_density): Invalid shape of the CIS coefficients alpha matrix!")'); stop
        endif
        if(size(cis_b,1).ne.num_occ.or.size(cis_b,2).ne.num_virt.or.size(cis_b,3).ne.num_cis_states) then
         write(*,'("#ERROR(compute_transition_density): Invalid shape of the CIS coefficients beta matrix!")'); stop
        endif
        write(*,'("Ok")')
 !Check CIS coefficients:
        write(*,'(" Checking the norm of the CIS coefficients ... ")')
        do n=1,num_cis_states
         norm_a=0d0
         do j=1,num_virt
          do i=1,num_occ
           norm_a = norm_a + cis_a(i,j,n)**2
          enddo
         enddo
         norm_b=0d0
         do j=1,num_virt
          do i=1,num_occ
           norm_b = norm_b + cis_b(i,j,n)**2
          enddo
         enddo
         write(*,'("  State ",i3," norm alpha/beta = ",D25.14,1x,D25.14)') n,norm_a,norm_b
         cis_a(:,:,n) = cis_a(:,:,n) * sqrt(1/norm_a)
         cis_b(:,:,n) = cis_b(:,:,n) * sqrt(1/norm_b)
        enddo
        write(*,'(" Ok: All states renormalized to unity")')
 !Compute CIS hole/particle density matrices in AO basis:
        write(*,'(" Allocating arrays ... ")',ADVANCE='NO')
        allocate(hole_dens(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_cis_states))
        allocate(particle_dens(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_cis_states))
        allocate(hdm(num_occ,num_occ))
        allocate(pdm(num_virt,num_virt))
        allocate(hao(mol_params%num_ao_orbitals,num_occ))
        allocate(pav(mol_params%num_ao_orbitals,num_virt))
        allocate(maa(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        allocate(mtx(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        allocate(mmm(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals))
        write(*,'("Ok")')
        do n=1,num_cis_states
         write(*,'(" Processing electronic state ",i2," ...")') n
         hole_dens(:,:,n)=0d0; particle_dens(:,:,n)=0d0

         write(*,'("  Checking metrics ... ")',ADVANCE='NO')
         maa(:,:)=0d0; mtx(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_mo_orbitals,&
                     &mo_a(:,:),mo_a(:,:),maa(:,:)) !inverse overlap
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &maa(:,:),overlap(:,:),mtx(:,:)) !delta
         !call wr_mat_dp(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mtx(:,:))
         do i=1,mol_params%num_ao_orbitals; mtx(i,i)=mtx(i,i)-1d0; enddo
         metrics_ao=max(abs(maxval(mtx)),abs(minval(mtx)))
         write(*,'(" AO metrics error = ",D25.14,";")',ADVANCE='NO') metrics_ao

         mmm(:,:)=0d0; mtx(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_mo_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),mo_a(:,:),mtx(:,1:mol_params%num_mo_orbitals))
         call mattmat(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals,mol_params%num_ao_orbitals,&
                     &mo_a(:,:),mtx(:,1:mol_params%num_mo_orbitals),mmm(:,:))
         !call wr_mat_dp(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals,mmm(:,:))
         do i=1,mol_params%num_mo_orbitals; mmm(i,i)=mmm(i,i)-1d0; enddo
         metrics_mo=max(abs(maxval(mmm)),abs(minval(mmm)))
         write(*,'(" MO metrics error = ",D25.14)',ADVANCE='NO') metrics_mo

         if(metrics_ao.ge.1d-4.or.metrics_mo.ge.1d-4) then
          write(*,'("  Failed")'); stop
         endif
         write(*,'(" Ok")')

         write(*,'("  Computing hole density matrix ... ")',ADVANCE='NO')
         maa(:,:)=0d0
         hdm(:,:)=0d0; mtx(:,:)=0d0
         call matmatt(num_occ,num_occ,num_virt,cis_a(:,:,n),cis_a(:,:,n),hdm(:,:))
         call matmat(mol_params%num_ao_orbitals,num_occ,num_occ,&
                    &mo_a(:,1:num_occ),hdm(:,:),mtx(:,1:num_occ))
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                     &mtx(:,1:num_occ),mo_a(:,1:num_occ),maa(:,:))
         hdm(:,:)=0d0; mtx(:,:)=0d0
         call matmatt(num_occ,num_occ,num_virt,cis_b(:,:,n),cis_b(:,:,n),hdm(:,:))
         call matmat(mol_params%num_ao_orbitals,num_occ,num_occ,&
                    &mo_b(:,1:num_occ),hdm(:,:),mtx(:,1:num_occ))
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                     &mtx(:,1:num_occ),mo_b(:,1:num_occ),maa(:,:))
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),maa(:,:),hole_dens(:,:,n))
         hole_dens(:,:,n)=hole_dens(:,:,n)*(5d-1) !normalize to unity
         norm=0d0; do i=1,mol_params%num_ao_orbitals; norm=norm+hole_dens(i,i,n); enddo
         write(*,'("Trace = ",D25.14)',ADVANCE='NO') norm
         write(*,'(" Ok")')

         write(*,'("  Computing particle density matrix ... ")',ADVANCE='NO')
         maa(:,:)=0d0
         pdm(:,:)=0d0; mtx(:,:)=0d0
         call mattmat(num_virt,num_virt,num_occ,cis_a(:,:,n),cis_a(:,:,n),pdm(:,:))
         call matmat(mol_params%num_ao_orbitals,num_virt,num_virt,&
                    &mo_a(:,num_occ+1:num_occ+num_virt),pdm(:,:),mtx(:,num_occ+1:num_occ+num_virt))
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_virt,&
                     &mtx(:,num_occ+1:num_occ+num_virt),mo_a(:,num_occ+1:num_occ+num_virt),maa(:,:))
         pdm(:,:)=0d0; mtx(:,:)=0d0
         call mattmat(num_virt,num_virt,num_occ,cis_b(:,:,n),cis_b(:,:,n),pdm(:,:))
         call matmat(mol_params%num_ao_orbitals,num_virt,num_virt,&
                    &mo_b(:,num_occ+1:num_occ+num_virt),pdm(:,:),mtx(:,num_occ+1:num_occ+num_virt))
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_virt,&
                     &mtx(:,num_occ+1:num_occ+num_virt),mo_b(:,num_occ+1:num_occ+num_virt),maa(:,:))
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),maa(:,:),particle_dens(:,:,n))
         particle_dens(:,:,n)=particle_dens(:,:,n)*(5d-1) !normalize to unity
         norm=0d0; do i=1,mol_params%num_ao_orbitals; norm=norm+particle_dens(i,i,n); enddo
         write(*,'("Trace = ",D25.14)',ADVANCE='NO') norm
         write(*,'(" Ok")')

         write(*,'(" Done")')
        enddo
        write(*,'(" Cleaning temporaries ... ")',ADVANCE='NO')
        if(allocated(mmm)) deallocate(mmm)
        if(allocated(mtx)) deallocate(mtx)
        if(allocated(maa)) deallocate(maa)
        if(allocated(pav)) deallocate(pav)
        if(allocated(hao)) deallocate(hao)
        if(allocated(pdm)) deallocate(pdm)
        if(allocated(hdm)) deallocate(hdm)
        write(*,'("Ok")')
        write(*,'("Success: Hole/particle density matrices computed successfully!")')
        return
       end subroutine compute_transition_density

       subroutine compute_atomic_state_vectors(mol_params,basis_info,hole_dens,particle_dens,asv)
        type(mol_params_t), intent(in):: mol_params
        type(basis_func_info_t), intent(in):: basis_info(:) ![1:AO]
        real(8), intent(in):: hole_dens(:,:,:),particle_dens(:,:,:) ![1:AO,1:AO,1:STATES]
        type(atom_state_t), allocatable, intent(out):: asv(:,:) ![1:ATOMS,1:STATES]
        integer:: num_states,num_basis_func,atom,shell,n,m,state
        real(8):: hdens,pdens

        write(*,'("Computing atomic state vectors:")')
        num_states=size(particle_dens,3)
        write(*,'(" Allocating atomic state vector array ",i6," x ",i3," ... ")',ADVANCE='NO')&
        &mol_params%num_atoms,num_states
        allocate(asv(1:mol_params%num_atoms,1:num_states))
        asv(:,:)=atom_state_t(0,0,0,0)
        write(*,'("Ok")')
        num_basis_func=size(basis_info,1)
        do n=1,num_states
         write(*,'(" Computing atomic state vectors for state ",i2," ... ")',ADVANCE='NO') n
         atom=0; shell=0; hdens=0d0; pdens=0d0
         do m=1,num_basis_func
          if(basis_info(m)%atom_id.gt.atom) then !new atom
           if(atom.gt.0) then
            asv(atom,n)%atom_id=atom
            asv(atom,n)%hole_density(shell)=hdens
            asv(atom,n)%particle_density(shell)=pdens
            asv(atom,n)%num_shells=shell+1
           endif
           atom=basis_info(m)%atom_id; shell=basis_info(m)%shell_id; hdens=0d0; pdens=0d0
          elseif(basis_info(m)%atom_id.eq.atom) then !same atom
           if(basis_info(m)%shell_id.gt.shell) then !next shell
            asv(atom,n)%hole_density(shell)=hdens
            asv(atom,n)%particle_density(shell)=pdens
            shell=basis_info(m)%shell_id; hdens=0d0; pdens=0d0
           endif
          endif
          hdens=hdens+hole_dens(m,m,n)
          pdens=pdens+particle_dens(m,m,n)
         enddo
         asv(atom,n)%atom_id=atom
         asv(atom,n)%hole_density(shell)=hdens
         asv(atom,n)%particle_density(shell)=pdens
         asv(atom,n)%num_shells=shell+1
         write(*,'("Ok")')
        enddo
        write(*,'("Success: Atomic state vectors computed successfully!")')
        return
       end subroutine compute_atomic_state_vectors

       subroutine print_atomic_state_vectors(asv)
        type(atom_state_t), intent(in):: asv(1:,1:) ![1:ATOMS,1:STATES]
        integer:: state,atom,shell

        write(*,'("Printing atomic state vectors:")')
        do state=1,size(asv,2)
         write(*,'(1x,"Electronic state ",i3)') state
         do atom=1,size(asv,1)
          write(*,'(2x,"Atom ",i6)') asv(atom,state)%atom_id
          do shell=0,asv(atom,state)%num_shells-1
           write(*,'(3x,i2,3x,D15.7,1x,D15.7)') shell,asv(atom,state)%hole_density(shell),asv(atom,state)%particle_density(shell)
          enddo
         enddo
        enddo
        write(*,'("Success: Atomic state vectors printed successfully!")')
        return
       end subroutine print_atomic_state_vectors

       subroutine save_atomic_state_vectors(asv,cis_energy)
        type(atom_state_t), intent(in):: asv(1:,1:) ![1:ATOMS,1:STATES]
        real(8), intent(in):: cis_energy(1:) ![1:STATES]
        integer:: state,atom,shell

        do state=1,size(asv,2)
         write(*,'("Electronic state ",i3,": Energy = ",D20.7)') state,cis_energy(state)
         do atom=1,size(asv,1)
          write(*,'(32(2x,D15.7,1x,D15.7))')&
          &(/(asv(atom,state)%hole_density(shell),asv(atom,state)%particle_density(shell),shell=0,MAX_AO_SHELLS-1)/)
         enddo
        enddo
        return
       end subroutine save_atomic_state_vectors

       subroutine filter_small_elems(matrix,thresh)
        real(8), intent(inout):: matrix(:,:)
        real(8), intent(in):: thresh
        integer:: i,j

        do j=lbound(matrix,2),ubound(matrix,2)
         do i=lbound(matrix,1),ubound(matrix,1)
          if(abs(matrix(i,j)).le.thresh) matrix(i,j)=0d0
         enddo
        enddo
        return
       end subroutine filter_small_elems

       end module mol_descr
