!Molecular descriptor generator
!Author: Dmitry I. Lyakh
       module mol_descr
       use chem_parser
       use stsubs
       implicit none

       public compute_transition_density

       contains

       subroutine compute_transition_density(mol_params,basis_info,overlap,mo_a,mo_b,cis_a,cis_b,hole_dens,particle_dens)
        type(mol_params_t), intent(in):: mol_params
        type(basis_func_info_t), intent(in):: basis_info(:) ![1:AO]
        real(8), intent(in):: overlap(:,:)                  ![1:AO,1:AO]
        real(8), intent(in):: mo_a(:,:),mo_b(:,:)           ![1:AO,1:MO]
        real(8), intent(in):: cis_a(:,:,:),cis_b(:,:,:)     ![1:OCC,1:VIRT,1:STATE]
        real(8), allocatable, intent(out):: hole_dens(:,:,:),particle_dens(:,:,:) ![1:AO,1:AO,1:STATE]
        real(8), allocatable:: cis_lu(:,:),t0(:,:),t1(:,:)
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
        allocate(cis_lu(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        allocate(t0(num_occ,mol_params%num_ao_orbitals))
        allocate(t1(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals))
        write(*,'("Ok")')
        do n=1,num_cis_states
         write(*,'(" Processing electronic state ",i4," ... ")',ADVANCE='NO') n
         hole_dens(:,:,n)=0d0; particle_dens(:,:,n)=0d0
  !Alpha contribution:
   !t0(i,n) = cis_a(i,a) * mo_a(n,a):
         t0(:,:)=0d0
         call matmatt(num_occ,mol_params%num_ao_orbitals,num_virt,&
                     &cis_a(:,:,n),mo_a(:,num_occ+1:num_occ+num_virt),t0(:,:))
   !t1(m,n) = mo_a(m,i) * t0(i,n):
         t1(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                    &mo_a(:,1:num_occ),t0(:,:),t1(:,:))
   !cis_lu(r,n) = overlap(r,s) * t1(n,s):
         cis_lu(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &overlap(:,:),t1(:,:),cis_lu(:,:))
   !hole_dens(m,n) += t1(m,r) * cis_lu(r,n)
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &t1(:,:),cis_lu(:,:),hole_dens(:,:,n))
   !cis_lu(r,n) = overlap(r,s) * t1(s,n):
         cis_lu(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),t1(:,:),cis_lu(:,:))
   !particle_dens(m,n) += t1(r,m) * cis_lu(r,n):
         call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t1(:,:),cis_lu(:,:),particle_dens(:,:,n))
  !Beta contribution:
   !t0(i,n) = cis_a(i,a) * mo_a(n,a):
         t0(:,:)=0d0
         call matmatt(num_occ,mol_params%num_ao_orbitals,num_virt,&
                     &cis_b(:,:,n),mo_b(:,num_occ+1:num_occ+num_virt),t0(:,:))
   !t1(m,n) = mo_a(m,i) * t0(i,n):
         t1(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,num_occ,&
                    &mo_b(:,1:num_occ),t0(:,:),t1(:,:))
   !cis_lu(r,n) = overlap(r,s) * t1(n,s):
         cis_lu(:,:)=0d0
         call matmatt(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &overlap(:,:),t1(:,:),cis_lu(:,:))
   !hole_dens(m,n) += t1(m,r) * cis_lu(r,n)
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &t1(:,:),cis_lu(:,:),hole_dens(:,:,n))
   !cis_lu(r,n) = overlap(r,s) * t1(s,n):
         cis_lu(:,:)=0d0
         call matmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                    &overlap(:,:),t1(:,:),cis_lu(:,:))
   !particle_dens(m,n) += t1(r,m) * cis_lu(r,n):
         call mattmat(mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,mol_params%num_ao_orbitals,&
                     &t1(:,:),cis_lu(:,:),particle_dens(:,:,n))
         write(*,'("Ok")')
        enddo
        write(*,'(" Cleaning temporaries ... ")',ADVANCE='NO')
        deallocate(t1)
        deallocate(t0)
        deallocate(cis_lu)
        write(*,'("Ok")')
        write(*,'("Success: Hole/particle density matrices computed successfully!")')
        return
       end subroutine compute_transition_density

       end module mol_descr
