                      !
!> Subroutine to read a .fil file and write a correspondent .mm file or to  !
!> to modify a certain mesh based on the displacement results               !
!                                                                           !
! REVISION HISTORY:                                                         !
!> @date                                                                    !
!> 16/05/2016   .   Added keywords, Miguel Matos                            !
!                                                                           !
!---------------------------------------------------------------------------!   
subroutine ABQMAIN
       
    ! Abaqus include with variable size declarations
    include 'aba_param.inc'
    ! Arrays from abaqus files
    dimension  ARRAY(513), JRRAY(NPRECD,513), LRUNIT(2,1)
    equivalence (ARRAY(1), JRRAY(1,1))
    
    ! Parameters
    integer, parameter          ::  MAX_DOFS = 30
    double precision, parameter ::  NEAR_ZERO = 1.0d-35
    
    ! Filename arrays
    character(len=256)  ::  input_file, fil_file, output_file, nodes_file
    character(len=256)  ::  argument, extra_argument
    character(len=5)    ::  file_ext = '.mm'
    character(len=5)    ::  fil_ext  = '.fil'
    character(len=5)    ::  node_ext = '.inp'
    character(len=50)   ::  node_pre = '_def'
    character(len=1024) ::  node_file_header
    character(len=50)   ::  header_string
    character           ::  y_or_n
    
    ! Local variables
    integer             ::  node_number = 0, elem_number = 0, itrans = 0, count = 0, count_w = 0
    integer             ::  step_number = -1, incr_number = -1
    integer             ::  act_step = -1, act_incr = -1
    integer             ::  i, j, k, li, ientry, numargs, prev_key=-1
    integer             ::  IFILE = 6, input_key = -1
    integer             ::  next_element_number = -1
    logical             ::  to_screen = .False., limit_keys = .False., mod_mesh = .False., missing_last = .True.
    logical             ::  auxl = .False., alloc_done = .False., all_steps = .False., all_incrs = .False.
    double precision    ::  t0,t1
    integer             ::  aloc_stat = -1
    integer             ::  ACTIVE_DOF(MAX_DOFS) = 0
    ! Allocatable
    integer, allocatable            ::  INP_NNUMBER(:)
    double precision, allocatable   ::  DISPLC(:,:)
    double precision, allocatable   ::  INP_COORD(:,:)
    
    ! Analyse input data
    write(*,*) "GetResults - Reading input"
    numargs = iargc()
    do i=1, numargs
        call getarg(i,argument)
        if (argument(1:1) == "-") then
            !write(*,'(I0," ",A)') i,trim(argument) 
            select case (trim(argument(2:)))
            case ("job")
                auxl = .True.
                if (i+1>numargs) goto 9991
                call getarg(i+1,input_file)
                ! Remove extension if needed                
                j = scan(trim(input_file),".", BACK= .true.)             
                if (j .ne. 0) input_file = input_file(1:j-1)
                ! Print jobname
                write(*,'(A,A)') "   job    = ",trim(input_file)
            case ("key")
                if (i+1>numargs) goto 9991
                call getarg(i+1,extra_argument)
                read(extra_argument,*) input_key
                write(*,'(A,I0)') "   key    = ",input_key
                limit_keys = .True.
            case("deform")
                if (.not. auxl) goto 9991   
                write(*,'(A,A)') "   nodes  > ",trim(input_file)//trim(node_pre)//trim(node_ext)  
                mod_mesh = .True.
            case("step")
                if (i+1>numargs) goto 9991
                call getarg(i+1,extra_argument)
                if(extra_argument(1:3) == "all") then
                    all_steps = .True.
                    write(*,'(A,A)') "   step   = ","all"
                else
                    read(extra_argument,*) act_step
                    write(*,'(A,I0)') "   step   = ",act_step
                end if
            case("increment")
                if (i+1>numargs) goto 9991
                call getarg(i+1,extra_argument)
                if(extra_argument(1:3) == "all") then
                    all_incrs = .True.
                    write(*,'(A,A)')    "   increm = ","all"
                else
                    read(extra_argument,*) act_incr
                    write(*,'(A,I0)')   "   increm = ",act_incr
                end if
            case("h")
                goto 9990
            case ("screen")
                to_screen = .True.
            end select 
        end if
    end do 
    
    ! Check if there was a job input
    if(.not. auxl) goto 9991
    
    ! Check if results are present
    fil_file = trim(input_file)//trim(fil_ext)
    inquire(FILE=fil_file, EXIST=auxl)
    if(.not. auxl) goto 9992
    
    ! Initialize abaqus files
    NRU = 1
    LRUNIT(1,NRU) = 8   ! 8 for .fil file
    LRUNIT(2,NRU) = 2   ! 2 for binary, 1 for ASCII
    LOUTF = 0
    call INITPF (input_file, NRU, LRUNIT, LOUTF)
    
    ! Open the results file
    JUNIT = LRUNIT(1,NRU)   ! 8 for .fil file
    call DBRNU (JUNIT)
    
    ! Open output file
    if (.not. to_screen) then
        IFILE = 37
        output_file = trim(input_file)//trim(file_ext)
        write(*,'(A,A)') "   output = ",trim(output_file)
        open(UNIT=IFILE, FILE=output_file, ACTION='WRITE')
    end if
    
    ! Debug only
    !read (*,*) y_or_n
    
    ! Read file
    call CPU_TIME(t0)
    write(*,*) ""
    write(*,*) " Reading file"
    JRCD = 0
    j = 0
    do while (JRCD==0)
        
        ! Read line
        CALL DBFILE (0, ARRAY, JRCD)
        if(JRCD .ne. 0) then
            write(*,'(A,I0,A)') "          ",j," lines"   
            exit
        end if
        

        ! Get key
        KEY = JRRAY(1,2)
        
        ! Read if the key number is 1 - header only
        if(KEY .EQ. 1) then
            next_element_number = JRRAY(1,3)
        end if
        if(KEY .EQ. 1) go to 9902
        
        !START OF EDIT
        if(limit_keys .and. (KEY .EQ. 2000)) then
            prev_key = KEY
            write(IFILE,'(" -",A,G14.6,",",G14.6,",",G14.6)') "Step time, Total time, Increment, ",ARRAY(4),ARRAY(3),ARRAY(13)
        end if
        !END OF EDIT
                
        ! Limit output
        if(limit_keys .and. (KEY .ne. input_key)) go to 9901
            
        if((prev_key .ne. KEY)) then
            if(prev_key .ne. -1) write(*,'(A,I0,A)') "          ",j," lines"           
            write(*,'(A,I0)')   "   > key: ", KEY
            prev_key = KEY
            j = 0
        end if
                
        ! Record type: Model description
        if (KEY == 1921) then
            node_number = JRRAY(1,8)
            elem_number = JRRAY(1,7)
            write(IFILE,*)      " -Model"
            write(IFILE,2040)   "Element number", elem_number
            write(IFILE,2040)   "Node number   ", node_number
            write(IFILE,2041)   "Typical length", ARRAY(9)
            write(IFILE,*)      ""
            
            if (mod_mesh .and. .not. alloc_done) then
                allocate( INP_COORD(3,node_number), STAT=aloc_stat)
                if(aloc_stat .ne. 0) goto 9993
                INP_COORD(1:3,1:node_number) = 0.0d0
                allocate( INP_NNUMBER(node_number), STAT=aloc_stat)
                if(aloc_stat .ne. 0) goto 9993
                allocate( DISPLC(3,node_number),    STAT=aloc_stat)
                if(aloc_stat .ne. 0) goto 9993
                DISPLC(1:3,1:node_number) = 0.0d0
                alloc_done = True                
            end if
            
        ! Record type: Node definitions
        else if (KEY == 1901 .and. mod_mesh) then
            ientry = j+1
            if (ientry>node_number .or. ientry<1) goto 9994            
            ! Store node number and coordinates
            INP_NNUMBER(ientry) = JRRAY(1,3)
            INP_COORD(1,ientry) = dble(ARRAY(4))
            ! write the coordinates only if exist
            if (JRRAY(1,1) .GE. 5) INP_COORD(2,ientry) = dble(ARRAY(5))
            if (JRRAY(1,1) .GE. 6) INP_COORD(3,ientry) = dble(ARRAY(6))
        
        ! Record type: Active degrees of freedom
        else if (KEY == 1902) then
            if (j==0) write(IFILE,*) " -Active degrees of freedom"
            do i=1,30
                if(JRRAY(1,2+i).ne.0) then
                    ACTIVE_DOF(i) = JRRAY(1,2+i)
                    if(i<=3) itrans = itrans+1  ! number of active translational DOF
                    write (IFILE,'(A,I0,A,I0)') "     DOF ",i," active @ location ", ACTIVE_DOF(i)
                end if
            end do
            !write(IFILE,*) ""
        
        ! Record type: Increment start record
        else if (KEY == 2000) then
            !write(*,*) "count=",count,",step_number=",step_number,",act_step=",act_step,",incr_number=",incr_number,"JRRAY(8)",JRRAY(1,8)
            !read (*,*) y_or_n
            if( (mod_mesh .and. count>0) .and. ( &
              & ( all_steps) .or. &   ! read all or read specific or read last incr of step
              & ( step_number == act_step .and. incr_number == act_incr) .or. &
              & ( step_number == act_step .and. incr_number == -1 .and. step_number .ne. JRRAY(1,8)) .or. &
              & ( step_number == act_step .and. all_incrs ) &
              & ) ) then
              
                if( (.not. all_incrs) .and. (.not. all_steps) ) node_pre = "_def"
              
                nodes_file =  trim(input_file)//trim(node_pre)//trim(node_ext)
                call print_deformed_nodes(  nodes_file,         &
                        &                   node_number,        &
                        &                   itrans,             &   
                        &                   INP_COORD,          & 
                        &                   INP_NNUMBER,        &
                        &                   DISPLC,             &
                        &                   node_file_header    )
                count_w = count_w+1
                
                ! the desired step+increment was written
                if( (.not. all_steps) .and. (step_number .ne. JRRAY(1,8)) ) missing_last = .False.
                ! the desired step+increment was written
                if( (.not. all_incrs) .and. (step_number == act_step) .and. (incr_number== act_incr) ) missing_last = .False.
                
            end if
            
            step_number = JRRAY(1,8)
            incr_number = JRRAY(1,9) 
            
            write(IFILE,*)      " -Increment"
            write(IFILE,2041)   "Step time  ",ARRAY(4)
            write(IFILE,2041)   "Total time ",ARRAY(3)
            write(IFILE,2041)   "Increment  ",ARRAY(13)
            write(IFILE,2040)   "Step number",step_number
            write(IFILE,2040)   "Incr number",incr_number
            !write(IFILE,*)      ""
            
            write(node_pre,'(A,I0,A,I0)') "_def_step",step_number,"_incr",incr_number
            write(node_file_header, '(A,A,A,I0,A,I0)') &
                &   "deformed nodes for ", trim(input_file), ", step = ", step_number, ", increment = ", incr_number
            
            count = count+1
            !*
            
        !   Output variable identifier: SDV
        else if (KEY == 5) then
            k = JRRAY(1,1)-2
            if(j==0) then
                write(IFILE,*)      " -Solution-dependent state variables"
                write(IFILE,'(" ",A10)',advance="no") "ELEMENT"
                li=1
                do while (li .LE. k)
                    write(header_string,'("SDV(",I0,")")') (li)
                    write(IFILE,'(",",A24)',advance="no") trim(header_string)
                    li = li+1
                end do 
                write(IFILE,'(A)') " " !this newline is needed
            end if
            
            ! write element number
            write(IFILE,'(" ",I10)',advance='no') next_element_number
            ! write each of the SDVs
            li=1
            do while (li .LE. k)
                write(IFILE,'(",",E24.14)',advance="no") ARRAY(3+li-1)
                li = li+1
            end do
            write(IFILE,'(A)') " "
            
            ! Limit screen output
            if (j==20 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
            
        !	Output variable identifier: U
        else if (KEY == 101) then
            k = JRRAY(1,1)-2
            ientry = j+1
            if (ientry>node_number .or. ientry<1) goto 9994
            if(j==0) then
                write(IFILE,*)      " -Displacement"
                if(k==3) write(IFILE,2020) "NODE","U1","U2","UR3"
                if(k==6) write(IFILE,2030) "NODE","U1","U2","U3","UR1","UR2","UR3"
                if(k .ne. 3 .and. k .ne. 6) write(IFILE,'(A,I0,A)') "  ERROR: Displacement w/ ",k," entries unexpected"
            end if
            
            !Store DOF
            if(mod_mesh) then
                do i=1,itrans
                    if(ARRAY(3+i) > NEAR_ZERO) DISPLC(i,ientry) = ARRAY(3+i)
                end do
                !write(*,1010)   JRRAY(1,3), (DISPLC(i,itrans),i = 1, itrans)
            end if
            
            !Print out
            if (k==3) write(IFILE,2021) JRRAY(1,3),ARRAY(4),ARRAY(5),ARRAY(6)
            if (k==6) write(IFILE,2031) JRRAY(1,3),ARRAY(4),ARRAY(5),ARRAY(6),ARRAY(7),ARRAY(8),ARRAY(9)
            
            ! Limit screen output
            if (j==10 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
        !   Output variable identifier: EPOT
        else if (KEY == 105) then
            k = JRRAY(1,1)-2
            if(j==0) then
                write(IFILE,*)      " -Electrical Potential"
                write(IFILE,2010) "NODE","EPOT"
            end if
            
            write(IFILE,2011) JRRAY(1,3),ARRAY(4)
            
            ! Limit screen output
            if (j==20 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
        !   Output variable identifier: RCHG
        else if (KEY == 119) then
            k = JRRAY(1,1)-2
            if(j==0) then
                write(IFILE,*)      " -Electrical reaction charges"
                write(IFILE,2010) "NODE","RCHG"
            end if
            
            write(IFILE,2011) JRRAY(1,3),ARRAY(4)
            
            ! Limit screen output
            if (j==20 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
        !   Output variable identifier: CECHG
        else if (KEY == 120) then
            k = JRRAY(1,1)-2
            if(j==0) then
                write(IFILE,*)      " -Concentrated electrical nodal charges"
                write(IFILE,2010) "NODE","CECHG"
            end if
            
            write(IFILE,2011) JRRAY(1,3),ARRAY(4)
            
            
            
            
            
            ! Limit screen output
            if (j==20 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
        !   Output variable identifier: RFL
        else if (KEY == 204) then
            k = JRRAY(1,1)-2
            if(j==0) then
                write(IFILE,*)      " -Reaction fluxes"
                write(IFILE,2010) "NODE","RFL"
            end if
            
            write(IFILE,2011) JRRAY(1,3),ARRAY(4)
            
            ! Limit screen output
            if (j==20 .and. to_screen) then
                write(IFILE,'(A)',ADVANCE="NO") " Large output expected. Continue?"
                read (*,*) y_or_n
                if(y_or_n == "N" .or. y_or_n == "n") exit
            end if
            
        end if

               
9901    continue
        j=j+1
        
9902    continue        
        
    end do
    
    
    ! Close file
    if (.not. to_screen) then
        close(IFILE)
    end if
    
    
    ! Deformed mesh checks
    if( mod_mesh ) then
        if( act_step == -1) then
            continue
        elseif( step_number == act_step .and. incr_number == act_incr) then
            continue
        elseif( step_number .ne. act_step) then
            missing_last = .False.
        elseif( step_number == act_step .and. act_incr .ne. -1) then
            missing_last = .False.
        end if
    end if
    ! Last deformed mesh        
    if( mod_mesh .and. (all_steps .or. missing_last) ) then
        if(count_w==0) node_pre = "_def"    !none written
        nodes_file =  trim(input_file)//trim(node_pre)//trim(node_ext)
        call print_deformed_nodes(  nodes_file,         &
                &                   node_number,        &
                &                   itrans,             &   
                &                   INP_COORD,          & 
                &                   INP_NNUMBER,        &
                &                   DISPLC,             &
                &                   node_file_header    )
        count_w = count_w+1
    end if    
       
    ! Deallocate
    if ( mod_mesh ) then
        deallocate(INP_COORD,  STAT=aloc_stat)
        if(aloc_stat .ne. 0) goto 9993
        deallocate(INP_NNUMBER,  STAT=aloc_stat)
        if(aloc_stat .ne. 0) goto 9993
        deallocate(DISPLC,        STAT=aloc_stat)
        if(aloc_stat .ne. 0) goto 9993
        ! check if something was written
        if( count_w==0) goto 9995
    end if
    
    
    ! Present timer
    call CPU_TIME(t1)
    write(*,*)  ""
    write(*,'(A,F8.5)') "  File successfully read. Elapsed CPU time:",t1-t0
    
    return

9990 write(*,'(A,A)')   ACHAR(10),"  Usage: abaqus getresults -job <jobname> [-key <#> -step <all/#> -increment <all/#> -deform]"
     return
9991 write(*,'(A,A)')   ACHAR(10),"  Input error. Exiting... Try -h for help"
     return
9992 write(*,'(A,A,A)') ACHAR(10),"  Input error: Result file not found: ", fil_file     
     return
9993 write(*,'(A,A)')   ACHAR(10),"  Memory error, (de)allocation failed!"
     return
9994 write(*,'(A,A)')   ACHAR(10),"  Memory error, indexes exceeded!"
     return
9995 write(*,'(A,A)')   ACHAR(10),"  File error, desired step/increment not found!"
     return  
     
     ! Node Displacement
1011 format(I7,3(',',G14.6))
1010 format(I6,3(',',1PE14.6))
     
     ! Nodal/1D
2010 format(' ',A10,',',A24)
2011 format(' ',I10,',',E24.14)    
     ! 2D formats
2020 format(' ',A10,',',A24,',',A24,',',A24)
2021 format(' ',I10,',',E24.14,',',E24.14,',',E24.14)     
     ! 3D formats
2030 format(' ',A10,',',A24,',',A24,',',A24',',A24,',',A24,',',A24)
2031 format(' ',I10,',',E24.14,',',E24.14,',',E24.14,',',E24.14,',',E24.14,',',E24.14)       
     ! Integer simple
2040 format(' ','    ',A,' :  ',I0)     
2041 format(' ','    ',A,' :', G11.4)    
     
    end subroutine
    
    
!---------------------------------------------------------------------------!
!                                                                           !
!---------------------------------------------------------------------------!
subroutine print_deformed_nodes(filename,                                   &
    &                           node_number,                                &
    &                           trans_number,                               &   
    &                           INP_COORDINATES,                            & 
    &                           INP_NUMBER,                                 &
    &                           READ_DISPLACEMENTS,                         &
    &                           header_comment                              )

    implicit none

    ! input
    character(len=*)                ::  filename
    integer, intent(in)             ::  node_number, trans_number
    integer, intent(in)             ::  INP_NUMBER(node_number)
    double precision, intent(in)    ::  INP_COORDINATES(3,node_number)
    double precision, intent(in)    ::  READ_DISPLACEMENTS(3,node_number)
    character(len=*), optional      ::  header_comment

    ! local variables
    integer                         ::  i, j
    double precision                ::  deformed_coordinates(3)

    ! Open file
    open(UNIT=192, FILE=trim(filename), ACTION='WRITE') 
    write(*,'(A,A)') "   writing nodes to ",trim(filename)
    
    ! Write header if desired
    if(present(header_comment)) then
        write(192,'(A,A)') "** ", trim(header_comment)
    end if
    
    ! Deform mesh
    do i=1,node_number
        deformed_coordinates(1:3) = 0.0d0
        do j=1,trans_number
            deformed_coordinates(j) = INP_COORDINATES(j,i) + READ_DISPLACEMENTS(j,i)
        end do
        write(192,1010) INP_NUMBER(i), (deformed_coordinates(j),j = 1, trans_number)
    end do
    
    ! Close file
    close(192)
    
1010 format(I6,3(',',1PE14.6))

end subroutine print_deformed_nodes
!---------------------------------------------------------------------------! 
    
    
