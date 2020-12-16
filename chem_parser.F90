!Parser of the Psi4 output.
!Author: Dmitry I. Lyakh
       module chem_parser
       use parse_prim
       use combinatoric
       use stsubs
       implicit none

       logical, private:: VERBOSE=.TRUE.

       type, public:: mol_params_t
        integer:: num_ao_orbitals
        integer:: num_mo_orbitals
        integer:: num_electrons_a
        integer:: num_electrons_b
        integer:: multiplicity
        integer:: charge
        real(8):: nuclear_repulsion
       end type mol_params_t

       contains

       function psi4_extract_mol_params(filename,mol_params) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename          !in: file name
        type(mol_params_t), intent(out):: mol_params !out: molecular parameters
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,spin,ierr,i,j,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val

        parsed=.FALSE.
        open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
        if(VERBOSE) then
         write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
        endif
        eloop: do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          call skip_blanks(str(1:l),first,last)
          if(str(first:first+len_trim('Nuclear repulsion = ')-1).eq.'Nuclear repulsion = ') then
           call charnum(str(first+len_trim('Nuclear repulsion = '):last),mol_params%nuclear_repulsion,i)
           if(VERBOSE) write(*,'("Extracted nuclear repulsion = ",F25.6)') mol_params%nuclear_repulsion
           str=' '; read(10,'(A1024)',end=100) str
           str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
           matched=match_symb_pattern(str(1:l),'Charge       = `',num_pred,pred_offset,pred_length,ierr)
           if(matched) then
            call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,mol_params%charge)
            if(VERBOSE) write(*,'("Extracted molecule charge = ",i6)') mol_params%charge
            str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
            matched=match_symb_pattern(str(1:l),'Multiplicity = `',num_pred,pred_offset,pred_length,ierr)
            if(matched) then
             call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,mol_params%multiplicity)
             if(VERBOSE) write(*,'("Extracted multiplicity = ",i6)') mol_params%multiplicity
             str=' '; read(10,'(A1024)',end=100) str
             str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
             matched=match_symb_pattern(str(1:l),'Nalpha       = `',num_pred,pred_offset,pred_length,ierr)
             if(matched) then
              call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,mol_params%num_electrons_a)
              if(VERBOSE) write(*,'("Extracted number of alpha electrons = ",i6)') mol_params%num_electrons_a
              str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
              matched=match_symb_pattern(str(1:l),'Nbeta        = `',num_pred,pred_offset,pred_length,ierr)
              if(matched) then
               call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,mol_params%num_electrons_b)
               if(VERBOSE) write(*,'("Extracted number of beta electrons = ",i6)') mol_params%num_electrons_b
               do
                str=' '; read(10,'(A1024)') str; l=len_trim(str)
                if(l.gt.0) then
                 call skip_blanks(str(1:l),first,last)
                 if(str(first:first+len_trim('Number of basis function: ')-1).eq.'Number of basis function: ') then
                  call charnum(str(first+len_trim('Number of basis function: '):last),val,mol_params%num_mo_orbitals)
                  if(VERBOSE) write(*,'("Extracted number of basis functions = ",i6)') mol_params%num_mo_orbitals
                  exit
                 endif
                endif
               enddo
              else
               write(*,'("#ERROR: Unable to find number of beta electrons!")'); stop
              endif
             else
              write(*,'("#ERROR: Unable to find number of alpha electrons!")'); stop
             endif
            else
             write(*,'("#ERROR: Unable to find multiplcity!")'); stop
            endif
           else
            write(*,'("#ERROR: Unable to find molecule charge!")'); stop
           endif
           parsed=.TRUE.
           exit eloop
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function psi4_extract_mol_params

       function psi4_extract_overlap(filename,overlap) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename              !in: file name
        real(8), allocatable, intent(out):: overlap(:,:) !out: AO overlap matrix
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,ierr,i,j,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val

        parsed=.FALSE.
        open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
        if(VERBOSE) then
         write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
        endif
        do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          call skip_blanks(str(1:l),first,last)
          if(str(first:first+len_trim('## S')-1).eq.'## S') then
           if(VERBOSE) write(*,'("Detected AO overlap matrix")')
           str=' '; read(10,'(A1024)') str; l=len_trim(str)
           matched=match_symb_pattern(str(1:l),'Irrep: ` Size: ` x `',num_pred,pred_offset,pred_length,ierr)
           if(ierr.eq.0) then
            call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,m)
            call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,n)
            if(VERBOSE) write(*,'("Dimensions = ",i9," x ",i9)') m,n
            if(.not.allocated(overlap)) then
             allocate(overlap(m,n))
            else
             write(*,'("#ERROR: Repeated AO overlap matrix!")'); stop
            endif
            if(VERBOSE) write(*,'("Allocated AO overlap matrix array")')
            do j=1,n,5
             str=' '; read(10,'(A1024)') str
             str=' '; read(10,'(A1024)') str
             str=' '; read(10,'(A1024)') str
             do i=1,m
              read(10,*) l,overlap(i,j:min(n,j+4))
             enddo
            enddo
            if(VERBOSE) write(*,'("Extracted the data")')
            parsed=.TRUE.
           endif
          endif
         endif
        enddo
100     close(10)
        return
       end function psi4_extract_overlap

       function psi4_extract_mo_coef(filename,mo_coef_a,mo_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                        !in: file name
        real(8), allocatable, target, intent(out):: mo_coef_a(:,:) !out: Alpha MO coefficients
        real(8), allocatable, target, intent(out):: mo_coef_b(:,:) !out: Alpha MO coefficients
        real(8), pointer:: mo_coef(:,:)
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,spin,ierr,i,j,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val

        parsed=.FALSE.
        open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
        if(VERBOSE) then
         write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
        endif
        do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          spin=0
          call skip_blanks(str(1:l),first,last)
          if(str(first:first+len_trim('## Alpha MO coefficients')-1).eq.'## Alpha MO coefficients') then
           if(VERBOSE) write(*,'("Detected Alpha MO coefficients")')
           spin=+1
          elseif(str(first:first+len_trim('## Beta MO coefficients')-1).eq.'## Beta MO coefficients') then
           if(VERBOSE) write(*,'("Detected Beta MO coefficients")')
           spin=-1
          endif
          if(spin.ne.0) then
           str=' '; read(10,'(A1024)') str; l=len_trim(str)
           matched=match_symb_pattern(str(1:l),'Irrep: ` Size: ` x `',num_pred,pred_offset,pred_length,ierr)
           if(ierr.eq.0) then
            call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,m)
            call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,n)
            if(VERBOSE) write(*,'("Dimensions = ",i9," x ",i9)') m,n
            if(spin.gt.0) then
             if(.not.allocated(mo_coef_a)) then
              allocate(mo_coef_a(m,n))
              mo_coef=>mo_coef_a
             else
              write(*,'("#ERROR: Repeated Alpha MO coefficients!")'); stop
             endif
            elseif(spin.lt.0) then
             if(.not.allocated(mo_coef_b)) then
              allocate(mo_coef_b(m,n))
              mo_coef=>mo_coef_b
             else
              write(*,'("#ERROR: Repeated Beta MO coefficients!")'); stop
             endif
            endif
            if(VERBOSE) write(*,'("Allocated MO coefficients array")')
            do j=1,n,5
             str=' '; read(10,'(A1024)') str
             str=' '; read(10,'(A1024)') str
             str=' '; read(10,'(A1024)') str
             do i=1,m
              read(10,*) l,mo_coef(i,j:min(n,j+4))
             enddo
            enddo
            if(VERBOSE) write(*,'("Extracted the data")')
            parsed=.TRUE.
           endif
          endif
         endif
        enddo
100     close(10)
        return
       end function psi4_extract_mo_coef

       function psi4_extract_cis_coef(filename,mol_params,cis_coef_a,cis_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                   !in: file name
        type(mol_params_t), intent(in):: mol_params           !in: molecular parameters
        real(8), allocatable, intent(out):: cis_coef_a(:,:,:) !out: CIS coefficients for all states
        real(8), allocatable, intent(out):: cis_coef_b(:,:,:) !out: CIS coefficients for all states
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,inda,indb,ierr,i,j,k,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val,coef

        parsed=.FALSE.
        open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
        if(VERBOSE) then
         write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
        endif
        eloop: do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          call skip_blanks(str(1:l),first,last)
          if(last-first+1.eq.len('Configuration Interaction').and.str(first:last).eq.'Configuration Interaction') then
           if(VERBOSE) write(*,'("Detected configuration interaction section")')
           do
            str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
            if(l.gt.0) then
             matched=match_symb_pattern(str(1:l),'NUM ROOTS      =        ` `',num_pred,pred_offset,pred_length,ierr)
             if(matched) then
              call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,n)
              if(VERBOSE) write(*,'("Detected number of roots = ",i6)') n
              exit
             endif
            endif
           enddo
           allocate(cis_coef_a(mol_params%num_electrons_a,mol_params%num_mo_orbitals-mol_params%num_electrons_a,1:n-1))
           allocate(cis_coef_b(mol_params%num_electrons_b,mol_params%num_mo_orbitals-mol_params%num_electrons_b,1:n-1))
           if(VERBOSE) write(*,'("Allocated CIS coefficient arrays")')
           do m=0,n-1 !number of roots
            do
             str=' '; read(10,'(A1024)') str; l=len_trim(str)
             if(l.gt.0) then
              matched=match_symb_pattern(str(1:l),'CI Root ` energy = `',num_pred,pred_offset,pred_length,ierr)
              if(matched) then
               call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,i)
               call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,j)
               if(VERBOSE) write(*,'("Found CI root ",i6," with energy ",F25.6)') i,val
               if(i.eq.m) then
                exit
               else
                write(*,'("#ERROR: Missing CI root info!")'); stop
               endif
              endif
             endif
            enddo
            str=' '; read(10,'(A1024)') str
            str=' '; read(10,'(A1024)') str
            str=' '; read(10,'(A1024)') str
            k=0
            do
             str=' '; read(10,'(A1024)') str; l=len_trim(str)
             if(l.gt.0) then
              matched=match_symb_pattern(str(1:l),'*   `    `  (    `,    `)  `',num_pred,pred_offset,pred_length,ierr)
              if(matched) then
               k=k+1
               call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),coef,i)
               call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,inda)
               call charnum(str(pred_offset(4):pred_offset(4)+pred_length(4)-1),val,indb)
               if(m.gt.0) then !skip the ground state
                if(inda.gt.0) then
                 i=mol_params%num_electrons_a-mod(inda-1,mol_params%num_electrons_a)
                 j=(inda-1)/mol_params%num_electrons_a+1
                 cis_coef_a(i,j,m)=coef
                elseif(indb.gt.0) then
                 i=mol_params%num_electrons_b-mod(indb-1,mol_params%num_electrons_b)
                 j=(indb-1)/mol_params%num_electrons_b+1
                 cis_coef_b(i,j,m)=coef
                endif
               endif
              else
               write(*,'("#ERROR: Invalid CI vector output!")'); stop
              endif
             else
              if(VERBOSE) write(*,'("Extracted CI root ",i6," with ",i9," coefficients")') m,k
              exit
             endif
            enddo
           enddo
           if(VERBOSE) write(*,'("Extracted all CI roots")')
           parsed=.TRUE.
           exit eloop
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function psi4_extract_cis_coef

       end module chem_parser


       program test_chem_parser
        use chem_parser
        logical:: parsed
        type(mol_params_t):: mol_params
        real(8), allocatable:: moa(:,:),mob(:,:),sao(:,:),cisa(:,:,:),cisb(:,:,:)

        parsed=psi4_extract_mol_params('output.dat',mol_params)
        parsed=psi4_extract_overlap('output.dat',sao)
        !call wr_mat_dp(size(sao,1),size(sao,2),sao) !debug
        parsed=psi4_extract_mo_coef('output.dat',moa,mob)
        !call wr_mat_dp(size(moa,1),size(moa,2),moa) !debug
        !call wr_mat_dp(size(mob,1),size(mob,2),mob) !debug
        parsed=psi4_extract_cis_coef('output.dat',mol_params,cisa,cisb)
        !call wr_mat_dp(size(cisa,1),size(cisa,2),cisa(:,:,1)) !debug
        !call wr_mat_dp(size(cisb,1),size(cisb,2),cisb(:,:,1)) !debug
       end program test_chem_parser
