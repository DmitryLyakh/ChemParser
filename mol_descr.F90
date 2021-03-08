!Molecular descriptor generator
!Author: Dmitry I. Lyakh
       module mol_descr
       use chem_parser
       use stsubs
       implicit none

       logical:: BIORTH_DENSITY=.FALSE.

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
        real(8), intent(in):: cis_a(:,:,:),cis_b(:,:,:)     ![1:OCC,1:VIRT,1:STATES]
        real(8), allocatable, intent(out):: hole_dens(:,:,:),particle_dens(:,:,:) ![1:AO,1:AO,1:STATES]
        real(8), allocatable:: cis_ul(:,:),t0(:,:),t1(:,:),t2(:,:),     t3(:,:)
        integer:: num_occ,num_virt,num_cis_states,i,j,k,l,m,n
        real(8):: norm_a,norm_b

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
         write(*,'("#ERROR(compute_transition_density): Spin-polarized molecules are not supported!")'); stop
        endif
        num_occ=mol_params%num_electrons_a; num_virt=mol_params%num_mo_orbitals-num_occ
        if(size(basis_info,1).ne.mol_params%num_ao_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid basis info!")'); stop
        endif
        if(size(overlap,1).ne.mol_params%num_ao_orbitals.or.size(overlap,2).ne.mol_params%num_ao_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid overlap matrix!")'); stop
        endif
        if(size(mo_a,1).ne.mol_params%num_ao_orbitals.or.size(mo_a,2).ne.mol_params%num_mo_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid MO coefficients alpha matrix!")'); stop
        endif
        if(size(mo_b,1).ne.mol_params%num_ao_orbitals.or.size(mo_b,2).ne.mol_params%num_mo_orbitals) then
         write(*,'("#ERROR(compute_transition_density): Invalid MO coefficients beta matrix!")'); stop
        endif
        num_cis_states=size(cis_a,3)
        if(size(cis_a,1).ne.num_occ.or.size(cis_a,2).ne.num_virt) then
         write(*,'("#ERROR(compute_transition_density): Invalid CIS coefficients alpha matrix!")'); stop
        endif
        if(size(cis_b,1).ne.num_occ.or.size(cis_b,2).ne.num_virt.or.size(cis_b,3).ne.num_cis_states) then
         write(*,'("#ERROR(compute_transition_density): Invalid CIS coefficients beta matrix!")'); stop
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
        enddo
        write(*,'("Ok")')
 !Compute CIS hole/particle density matrices in AO basis:
        write(*,'(" Allocating arrays ... ")',ADVANCE='NO')
        allocate(hole_dens(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_cis_states))
        allocate(particle_dens(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_cis_states))
        allocate(cis_ul(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        allocate(t0(num_occ,mol_params%num_ao_orbitals))
        allocate(t1(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        if(BIORTH_DENSITY) allocate(t2(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))

        allocate(t3(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals))

        write(*,'("Ok")')
        do n=1,num_cis_states
         write(*,'(" Processing electronic state ",i4," ... ")',ADVANCE='NO') n
         hole_dens(:,:,n)=0d0; particle_dens(:,:,n)=0d0

         write(*,'("Cmp*Cnp:")')
         t1(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_mo_orbitals,&
                     &mo_a(:,:),mo_a(:,:),t1(:,:))
         call filter_small_elems(t1,1d-10)
         call wr_mat_dp(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,t1(:,:))

         write(*,'("Cmp*Cmq:")')
         t3(:,:)=0d0
         call mattmat(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals,mol_params%num_ao_orbitals,&
                     &mo_a(:,:),mo_a(:,:),t3(:,:))
         call filter_small_elems(t3,1d-10)
         call wr_mat_dp(mol_params%num_mo_orbitals,mol_params%num_mo_orbitals,t3(:,:))
         stop


  !Alpha contribution:
   !t0(i,n) = cis_a(i,a) * mo_a(n,a):
         t0(:,:)=0d0
         call matmatt(num_occ,mol_params%num_ao_orbitals,num_virt,&
                     &cis_a(:,:,n),mo_a(:,num_occ+1:num_occ+num_virt),t0(:,:))
   !t1(m,n) = mo_a(m,i) * t0(i,n):
         t1(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                    &mo_a(:,1:num_occ),t0(:,:),t1(:,:))
   !cis_ul(n,r) =  t1(n,s) * overlap(r,s):
         cis_ul(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t1(:,:),overlap(:,:),cis_ul(:,:))
         if(BIORTH_DENSITY) then
   !t2(m,n) = t1(m,r) * cis_ul(n,r):
          t2(:,:)=0d0
          call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),t2(:,:))
   !hole_dens(m,n) += t2(m,r) * overlap(r,n):
          call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t2(:,:),overlap(:,:),hole_dens(:,:,n))
         else
   !hole_dens(m,n) = t1(m,r) * cis_ul(n,r):
          call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),hole_dens(:,:,n))
         endif
   !cis_ul(r,n) = overlap(r,s) * t1(s,n):
         cis_ul(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),t1(:,:),cis_ul(:,:))
         if(BIORTH_DENSITY) then
   !t2(m,n) = t1(r,m) * cis_ul(r,n):
          t2(:,:)=0d0
          call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),t2(:,:))
   !particle_dens(m,n) += t2(m,r) * overlap(r,n):
          call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t2(:,:),overlap(:,:),particle_dens(:,:,n))
         else
   !particle_dens(m,n) = t1(r,m) * cis_ul(r,n):
          call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),particle_dens(:,:,n))
         endif
  !Beta contribution:
   !t0(i,n) = cis_a(i,a) * mo_a(n,a):
         t0(:,:)=0d0
         call matmatt(num_occ,mol_params%num_ao_orbitals,num_virt,&
                     &cis_b(:,:,n),mo_b(:,num_occ+1:num_occ+num_virt),t0(:,:))
   !t1(m,n) = mo_a(m,i) * t0(i,n):
         t1(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                    &mo_b(:,1:num_occ),t0(:,:),t1(:,:))
   !cis_ul(n,r) = t1(n,s) * overlap(r,s):
         cis_ul(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t1(:,:),overlap(:,:),cis_ul(:,:))
         if(BIORTH_DENSITY) then
   !t2(m,n) = t1(m,r) * cis_ul(n,r):
          t2(:,:)=0d0
          call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),t2(:,:))
   !hole_dens(m,n) += t2(m,r) * overlap(r,n):
          call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t2(:,:),overlap(:,:),hole_dens(:,:,n))
         else
   !hole_dens(m,n) = t1(m,r) * cis_ul(n,r):
          call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),hole_dens(:,:,n))
         endif
   !cis_ul(r,n) = overlap(r,s) * t1(s,n):
         cis_ul(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),t1(:,:),cis_ul(:,:))
         if(BIORTH_DENSITY) then
   !t2(m,n) = t1(r,m) * cis_ul(r,n):
          t2(:,:)=0d0
          call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),t2(:,:))
   !particle_dens(m,n) += t2(m,r) * overlap(r,n)
          call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t2(:,:),overlap(:,:),particle_dens(:,:,n))
         else
   !particle_dens(m,n) = t1(r,m) * cis_ul(r,n):
          call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                      &t1(:,:),cis_ul(:,:),particle_dens(:,:,n))
         endif
         write(*,'("Ok")')
        enddo
        write(*,'(" Cleaning temporaries ... ")',ADVANCE='NO')
        if(allocated(t2)) deallocate(t2)
        if(allocated(t1)) deallocate(t1)
        if(allocated(t0)) deallocate(t0)
        if(allocated(cis_ul)) deallocate(cis_ul)
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
        write(*,'(" Allocating atomic state vectors array ",i6," x ",i3," ... ")',ADVANCE='NO')&
        &mol_params%num_atoms,num_states
        allocate(asv(1:mol_params%num_atoms,1:num_states))
        write(*,'("Ok")')
        num_basis_func=size(basis_info,1)
        do n=1,num_states
         write(*,'(" Computing atomic state vectors for state ",i3," ... ")',ADVANCE='NO') n
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
           write(*,'(3x,i2,3x,D25.14,1x,D25.14)') shell,asv(atom,state)%hole_density(shell),asv(atom,state)%particle_density(shell)
          enddo
         enddo
        enddo
        write(*,'("Success: Atomic state vectors printed successfully!")')
        return
       end subroutine print_atomic_state_vectors

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
