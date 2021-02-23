!Parser of output files of quantum-chemical software: Psi4, ORCA
!Author: Dmitry I. Lyakh
       module chem_parser
       use parse_prim
       use combinatoric
       use stsubs
       implicit none

       logical, private:: VERBOSE=.TRUE.

       type, public:: mol_params_t
        integer:: num_ao_orbitals=0
        integer:: num_mo_orbitals=0
        integer:: num_electrons=0
        integer:: num_electrons_a=0
        integer:: num_electrons_b=0
        integer:: multiplicity=0
        integer:: charge=0
        real(8):: nuclear_repulsion=0d0
       end type mol_params_t

       type, public:: basis_func_info_t
        integer:: atom_id   !atom id {1..M}
        integer:: shell_id  !shell id within atom {0..L}
        integer:: ao_id     !AO id within atom shell {0..S}
       end type basis_func_info_t

       public psi4_extract_mol_params
       public psi4_extract_overlap
       public psi4_extract_mo_coef
       public psi4_extract_cis_coef

       public orca_extract_mol_params
       public orca_extract_overlap
       public orca_extract_mo_coef
       public orca_extract_cis_coef
       public orca_extract_basis_info

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
               mol_params%num_electrons=mol_params%num_electrons_a+mol_params%num_electrons_b
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

       function psi4_extract_overlap(filename,mol_params,overlap) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename              !in: file name
        type(mol_params_t), intent(in):: mol_params      !in: molecular parameters
        real(8), allocatable, intent(out):: overlap(:,:) !out: AO overlap matrix (AO,AO)
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,ierr,i,j,l,m,n
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
          if(str(first:first+len_trim('## S')-1).eq.'## S') then
           if(VERBOSE) write(*,'("Detected AO overlap matrix")')
           str=' '; read(10,'(A1024)') str; l=len_trim(str)
           matched=match_symb_pattern(str(1:l),'Irrep: ` Size: ` x `',num_pred,pred_offset,pred_length,ierr)
           if(matched) then
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
            if(VERBOSE) write(*,'("Extracted the data successfully")')
            parsed=.TRUE.
            exit eloop
           endif
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function psi4_extract_overlap

       function psi4_extract_mo_coef(filename,mol_params,mo_coef_a,mo_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                        !in: file name
        type(mol_params_t), intent(in):: mol_params                !in: molecular parameters
        real(8), allocatable, target, intent(out):: mo_coef_a(:,:) !out: alpha MO coefficients (AO,MO)
        real(8), allocatable, target, intent(out):: mo_coef_b(:,:) !out: beta MO coefficients (AO,MO)
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
        eloop: do
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
           if(matched) then
            call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,m)
            call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,n)
            if(VERBOSE) write(*,'("Dimensions = ",i9," x ",i9)') m,n
            if(spin.gt.0) then
             if(.not.allocated(mo_coef_a)) then
              allocate(mo_coef_a(m,n))
              mo_coef_a=0d0
              mo_coef=>mo_coef_a
             else
              write(*,'("#ERROR: Repeated Alpha MO coefficients!")'); stop
             endif
            elseif(spin.lt.0) then
             if(.not.allocated(mo_coef_b)) then
              allocate(mo_coef_b(m,n))
              mo_coef_b=0d0
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
            if(VERBOSE) write(*,'("Extracted the data successfully")')
            parsed=.TRUE.
           endif
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function psi4_extract_mo_coef

       function psi4_extract_cis_coef(filename,mol_params,cis_coef_a,cis_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                   !in: file name
        type(mol_params_t), intent(in):: mol_params           !in: molecular parameters
        real(8), allocatable, intent(out):: cis_coef_a(:,:,:) !out: alpha CIS coefficients for all states (occ,virt,state)
        real(8), allocatable, intent(out):: cis_coef_b(:,:,:) !out: beta CIS coefficients for all states (occ,virt,state)
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
           cis_coef_a=0d0; cis_coef_b=0d0
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

       function orca_extract_mol_params(filename,mol_params) result(parsed)
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
          if(str(first:first+len_trim('ORCA SCF')-1).eq.'ORCA SCF') then
           sloop: do
            str=' '; read(10,'(A1024)') str; l=len_trim(str)
            if(l.gt.0) then
             call skip_blanks(str(1:l),first,last)
             if(str(first:first+len_trim('General Settings:')-1).eq.'General Settings:') then
              str=' '; read(10,'(A1024)') str
              str=' '; read(10,'(A1024)') str
              str=' '; read(10,'(A1024)') str; l=len_trim(str)
              matched=match_symb_pattern(str(1:l),'Total Charge`Charge`.... `',num_pred,pred_offset,pred_length,ierr)
              if(matched) then
               call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,mol_params%charge)
               if(VERBOSE) write(*,'("Extracted molecule charge = ",i6)') mol_params%charge
               str=' '; read(10,'(A1024)') str; l=len_trim(str)
               matched=match_symb_pattern(str(1:l),'Multiplicity`Mult`.... `',num_pred,pred_offset,pred_length,ierr)
               if(matched) then
                call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,mol_params%multiplicity)
                if(VERBOSE) write(*,'("Extracted spin multiplicity = ",i6)') mol_params%multiplicity
                str=' '; read(10,'(A1024)') str; l=len_trim(str)
                matched=match_symb_pattern(str(1:l),'Number of Electrons`NEL`.... `',num_pred,pred_offset,pred_length,ierr)
                if(matched) then
                 call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,mol_params%num_electrons)
                 if(VERBOSE) write(*,'("Extracted total number of electrons = ",i6)') mol_params%num_electrons
                 mol_params%num_electrons_a=(mol_params%num_electrons/2)+mod(mol_params%num_electrons,2)
                 mol_params%num_electrons_b=(mol_params%num_electrons/2)
                 if(VERBOSE) write(*,'("Assuming number of alpha/beta electrons = ",i6,1x,i6)')&
                                  &mol_params%num_electrons_a,mol_params%num_electrons_b
                 str=' '; read(10,'(A1024)') str; l=len_trim(str)
                 matched=match_symb_pattern(str(1:l),'Basis Dimension`Dim`.... `',num_pred,pred_offset,pred_length,ierr)
                 if(matched) then
                  call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),val,mol_params%num_ao_orbitals)
                  if(VERBOSE) write(*,'("Extracted number of basis functions = ",i6)') mol_params%num_ao_orbitals
                  mol_params%num_mo_orbitals=mol_params%num_ao_orbitals
                  if(VERBOSE) write(*,'("Assuming number of molecular orbitals = ",i6)') mol_params%num_mo_orbitals
                  str=' '; read(10,'(A1024)') str; l=len_trim(str)
                  matched=match_symb_pattern(str(1:l),'Nuclear Repulsion`ENuc`.... ` Eh',num_pred,pred_offset,pred_length,ierr)
                  if(matched) then
                   call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),mol_params%nuclear_repulsion,i)
                   if(VERBOSE) write(*,'("Extracted nuclear repulsion energy = ",D25.14)') mol_params%nuclear_repulsion
                   exit sloop
                  else
                   write(*,'("#ERROR: Unable to find nuclear repulsion energy!")'); stop
                  endif
                 else
                  write(*,'("#ERROR: Unable to find number of basis functions!")'); stop
                 endif
                else
                 write(*,'("#ERROR: Unable to find number of electrons!")'); stop
                endif
               else
                write(*,'("#ERROR: Unable to find spin multiplicity!")'); stop
               endif
              else
               write(*,'("#ERROR: Unable to find molecule charge!")'); stop
              endif
             endif
            endif
           enddo sloop
           parsed=.TRUE.
           exit eloop
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function orca_extract_mol_params

       function orca_extract_overlap(filename,mol_params,overlap) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename              !in: file name
        type(mol_params_t), intent(in):: mol_params      !in: molecular parameters
        real(8), allocatable, intent(out):: overlap(:,:) !out: AO overlap matrix (AO,AO)
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,ierr,i,j,l,m
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
          if(str(first:first+len_trim('OVERLAP MATRIX')-1).eq.'OVERLAP MATRIX') then
           if(VERBOSE) write(*,'("Detected AO overlap matrix")')
           str=' '; read(10,'(A1024)') str; l=len_trim(str)
           m=mol_params%num_ao_orbitals
           if(m.gt.0) then
            if(VERBOSE) write(*,'("Assumed dimensions = ",i9," x ",i9)') m,m
            if(.not.allocated(overlap)) then
             allocate(overlap(m,m))
            else
             write(*,'("#ERROR: Repeated AO overlap matrix!")'); stop
            endif
            if(VERBOSE) write(*,'("Allocated AO overlap matrix array")')
            do j=1,m,6
             str=' '; read(10,'(A1024)') str
             do i=1,m
              read(10,*) l,overlap(i,j:min(m,j+5))
             enddo
            enddo
            if(VERBOSE) write(*,'("Extracted the data successfully")')
            parsed=.TRUE.
            exit eloop
           else
            write(*,'("#ERROR: Molecular parameters are missing the number of AO!")'); stop
           endif
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function orca_extract_overlap

       function orca_extract_mo_coef(filename,mol_params,mo_coef_a,mo_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                        !in: file name
        type(mol_params_t), intent(in):: mol_params                !in: molecular parameters
        real(8), allocatable, target, intent(out):: mo_coef_a(:,:) !out: alpha MO coefficients (AO,MO)
        real(8), allocatable, target, intent(out):: mo_coef_b(:,:) !out: beta MO coefficients (AO,MO)
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
        eloop: do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          call skip_blanks(str(1:l),first,last)
          if(str(first:first+len_trim('[MO]')-1).eq.'[MO]') then
           if(VERBOSE) write(*,'("Detected MO coefficients")')
           allocate(mo_coef_a(mol_params%num_ao_orbitals,mol_params%num_mo_orbitals))
           allocate(mo_coef_b(mol_params%num_ao_orbitals,mol_params%num_mo_orbitals))
           mo_coef_a=0d0; mo_coef_b=0d0
           if(VERBOSE) write(*,'("Allocated MO coefficients: Dimensions = ",i9," x ",i9)')&
                            &mol_params%num_ao_orbitals,mol_params%num_mo_orbitals
           do j=1,mol_params%num_mo_orbitals
            str=' '; read(10,'(A1024)') str
            str=' '; read(10,'(A1024)') str; l=len_trim(str)
            str=' '; read(10,'(A1024)') str; l=len_trim(str)
            matched=match_symb_pattern(str(1:l),'Spin= `',num_pred,pred_offset,pred_length,ierr)
            if(matched) then
             mo_coef=>NULL()
             if(str(pred_offset(1):pred_offset(1)+pred_length(1)-1).eq.'Alpha') then
              mo_coef=>mo_coef_a
             elseif(str(pred_offset(1):pred_offset(1)+pred_length(1)-1).eq.'Beta') then
              mo_coef=>mo_coef_b
             endif
             if(associated(mo_coef)) then
              str=' '; read(10,'(A1024)') str; l=len_trim(str)
              do i=1,mol_params%num_ao_orbitals
               read(10,*) l,mo_coef(i,j)
              enddo
             else
              write(*,'("#ERROR: Invalid format of MO coefficients: Invalid spin!")'); stop
             endif
            else
             write(*,'("#ERROR: Invalid format of MO coefficients: Spin not found!")'); stop
            endif
           enddo
           if(VERBOSE) write(*,'("Extracted the data successfully")')
           parsed=.TRUE.
           exit eloop
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function orca_extract_mo_coef

       function orca_extract_cis_coef(filename,mol_params,cis_coef_a,cis_coef_b) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                           !in: file name
        type(mol_params_t), intent(in):: mol_params                   !in: molecular parameters
        real(8), allocatable, target, intent(out):: cis_coef_a(:,:,:) !out: alpha CIS coefficients for all states (occ,virt,state)
        real(8), allocatable, target, intent(out):: cis_coef_b(:,:,:) !out: beta CIS coefficients for all states (occ,virt,state)
        real(8), pointer:: cis_coef(:,:,:)
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,nocc,ierr,i,j,k,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val

        parsed=.FALSE.
        open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
        if(VERBOSE) then
         write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
        endif
        n=0
        eloop: do
         str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
         if(l.gt.0) then
          matched=match_symb_pattern(str(1:l),'CIS-EXCITED STATES (`)',num_pred,pred_offset,pred_length,ierr)
          if(matched) then
           if(VERBOSE) then
            write(*,'("Detected CIS excited state vectors:")',ADVANCE='NO')
            write(*,*) str(pred_offset(1):pred_offset(1)+pred_length(1)-1)
           endif
           do m=1,n
            do
             str=' '; read(10,'(A1024)') str; l=len_trim(str)
             if(l.ge.5) then
              if(str(1:5).eq.'STATE') exit
             endif
            enddo
            do
             str=' '; read(10,'(A1024)') str; l=len_trim(str)
             if(l.gt.0) then
              matched=match_symb_pattern(str(1:l),'` -> ` : ` (c= `)',num_pred,pred_offset,pred_length,ierr)
              if(matched) then
               if(str(pred_offset(1)+pred_length(1)-1:pred_offset(1)+pred_length(1)-1).eq.'a') then
                cis_coef=>cis_coef_a; nocc=mol_params%num_electrons_a
               elseif(str(pred_offset(1)+pred_length(1)-1:pred_offset(1)+pred_length(1)-1).eq.'b') then
                cis_coef=>cis_coef_b; nocc=mol_params%num_electrons_b
               else
                write(*,'("#ERROR: Invalid format of CIS coefficients: Invalid spin label")'); stop
               endif
               call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-2),val,i)
               call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-2),val,j)
               call charnum(str(pred_offset(4):pred_offset(4)+pred_length(4)-1),val,l)
               cis_coef(i+1,j+1-nocc,m)=val
              else
               write(*,'("#ERROR: Invalid format of CIS coefficients!")'); stop
              endif
             else
              exit
             endif
            enddo
            if(VERBOSE) write(*,'("Extracted CIS state ",i4)') m
           enddo
           if(VERBOSE) write(*,'("Extracted the data successfully")')
           parsed=.TRUE.
           exit eloop
          else
           matched=match_symb_pattern(str(1:l),'Number of roots to be determined`... `',num_pred,pred_offset,pred_length,ierr)
           if(matched) then
            call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,n)
            if(VERBOSE) write(*,'("Detected number of roots = ",i6)') n
            allocate(cis_coef_a(mol_params%num_electrons_a,mol_params%num_mo_orbitals-mol_params%num_electrons_a,1:n))
            allocate(cis_coef_b(mol_params%num_electrons_b,mol_params%num_mo_orbitals-mol_params%num_electrons_b,1:n))
            cis_coef_a=0d0; cis_coef_b=0d0
            if(VERBOSE) write(*,'("Allocated CIS coefficient arrays")')
           endif
          endif
         endif
        enddo eloop
100     close(10)
        return
       end function orca_extract_cis_coef

       function orca_extract_basis_info(filename,mol_params,basis_info) result(parsed)
        logical:: parsed
        character(*), intent(in):: filename                               !in: file name
        type(mol_params_t), intent(in):: mol_params                       !in: molecular parameters
        type(basis_func_info_t), allocatable, intent(out):: basis_info(:) !out: information of each basis function {1..N}
        integer:: pred_offset(1024),pred_length(1024),num_pred,first,last,nocc,ierr,i,j,l,m,n
        character(1024):: str
        logical:: matched
        real(8):: val

        parsed=.FALSE.
        if(mol_params%num_ao_orbitals.gt.0) then
         open(10,file=filename(1:len_trim(filename)),form='FORMATTED',status='OLD')
         if(VERBOSE) then
          write(*,'("Processing file")',advance='NO'); write(*,*) filename(1:len_trim(filename))
         endif
         eloop: do
          str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
          if(l.gt.0) then
           if(str(1:l).eq.'[GTO]') then
            if(VERBOSE) write(*,'("Detected Gaussian type orbital information:")')
            allocate(basis_info(mol_params%num_ao_orbitals))
            n=0
            do while(n.lt.mol_params%num_ao_orbitals) !loop over atoms
             read(10,*) i,j !i = atom id: [1..N]
             m=0
             do !loop over atomic shells
              str=' '; read(10,'(A1024)',end=100) str; l=len_trim(str)
              if(l.gt.0) then
               call skip_blanks(str(1:l),first,last)
               select case(str(first:first))
               case('s','S')
                do j=1,1
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               case('p','P')
                do j=1,3
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               case('d','D')
                do j=1,5
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               case('f','F')
                do j=1,7
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               case('g','G')
                do j=1,9
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               case('h','H')
                do j=1,11
                 n=n+1
                 basis_info(n)%atom_id=i
                 basis_info(n)%shell_id=m
                 basis_info(n)%ao_id=j-1
                enddo
                m=m+1
               end select
              else
               exit
              endif
             enddo
            enddo
            if(VERBOSE) then
             write(*,'("Extracted ",i7," basis functions successfully")') n
            endif
            parsed=.TRUE.
            exit eloop
           endif
          endif
         enddo eloop
100      close(10)
        else
         write(*,'("#ERROR: Unknown number of AO orbitals in Molecular Parameters!")'); stop
        endif
        return
       end function orca_extract_basis_info

       subroutine test_chem_parser()
        logical:: parsed
        type(mol_params_t):: mol_params
        real(8), allocatable:: moa(:,:),mob(:,:),sao(:,:),cisa(:,:,:),cisb(:,:,:)
        type(basis_func_info_t), allocatable:: basis_info(:)

 !Test Psi4 output parsing:
        parsed=psi4_extract_mol_params('psi4.dat',mol_params)
        parsed=psi4_extract_overlap('psi4.dat',mol_params,sao)
        !call wr_mat_dp(size(sao,1),size(sao,2),sao) !debug
        parsed=psi4_extract_mo_coef('psi4.dat',mol_params,moa,mob)
        !call wr_mat_dp(size(moa,1),size(moa,2),moa) !debug
        !call wr_mat_dp(size(mob,1),size(mob,2),mob) !debug
        parsed=psi4_extract_cis_coef('psi4.dat',mol_params,cisa,cisb)
        !call wr_mat_dp(size(cisa,1),size(cisa,2),cisa(:,:,1)) !debug
        !call wr_mat_dp(size(cisb,1),size(cisb,2),cisb(:,:,1)) !debug
        if(allocated(cisa)) deallocate(cisa)
        if(allocated(cisb)) deallocate(cisb)
        if(allocated(moa)) deallocate(moa)
        if(allocated(mob)) deallocate(mob)
        if(allocated(sao)) deallocate(sao)
        write(*,'()')

 !Test ORCA output parsing:
        parsed=orca_extract_mol_params('orca.dat',mol_params)
        parsed=orca_extract_overlap('orca.dat',mol_params,sao)
        parsed=orca_extract_mo_coef('orca.molden',mol_params,moa,mob)
        !call wr_mat_dp(size(moa,1),size(moa,2),moa) !debug
        !call wr_mat_dp(size(mob,1),size(mob,2),mob) !debug
        parsed=orca_extract_cis_coef('orca.dat',mol_params,cisa,cisb)
        !call wr_mat_dp(size(cisa,1),size(cisa,2),cisa(:,:,1)) !debug
        !call wr_mat_dp(size(cisb,1),size(cisb,2),cisb(:,:,1)) !debug
        parsed=orca_extract_basis_info('orca.molden',mol_params,basis_info)
        if(allocated(basis_info)) deallocate(basis_info)
        if(allocated(cisa)) deallocate(cisa)
        if(allocated(cisb)) deallocate(cisb)
        if(allocated(moa)) deallocate(moa)
        if(allocated(mob)) deallocate(mob)
        if(allocated(sao)) deallocate(sao)
        write(*,'()')
        return
       end subroutine test_chem_parser

       end module chem_parser
