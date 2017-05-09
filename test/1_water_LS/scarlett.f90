	module scarlett_mod
	implicit none
        real*8 :: Force_cut, Energy_cut, minimzation_steep !energy and force convergence crit and initial steep
        integer :: n_points ! number of points scaned for lineal search
        integer :: n_min_steeps !number of optimization steps
        logical :: lineal_search !enable lineal search
	real*8, dimension(:,:), allocatable :: r
	integer, dimension(:), allocatable :: Iz
	integer :: natom
	integer :: verbose
	character(len=500) :: do_energy_1, do_energy_2, do_energy_3
	character(len=500) :: do_forces_1, do_forces_2, do_forces_3
	character(len=500) :: make_coords
	end module

	program scarlett
	USE scarlett_mod, only : natom
	implicit none
	integer :: narg, i
	character(len=20) :: argument, inpfile
	real*8 :: E
	E=0.d0
	call defaults_scarlet()

	narg=command_argument_count()
	do i=1, narg
	   call get_command_argument(i,argument)
           select case(adjustl(argument))
           case("-i")
             call get_command_argument(i+1,inpfile)
           case default
           endselect
	enddo

	call read_input(inpFile)

	if (natom.eq.0) stop "0 atoms in your input file"
	call init_scarlet()

	call do_steep(E)
	write(*,*) "final energy ", E
	call finalize_scarlet()
	end program 

	subroutine init_scarlet()
	USE scarlett_mod, only : r, Iz, natom
	implicit none
	allocate (r(natom,3), Iz(natom))
!para prueba
	Iz=1
	IZ(1)=8
	end subroutine init_scarlet

	subroutine finalize_scarlet()
        USE scarlett_mod, only : r, Iz
        implicit none
        deallocate (r, Iz)
        end subroutine finalize_scarlet

	subroutine defaults_scarlet()
	USE scarlett_mod, only :Force_cut, Energy_cut, minimzation_steep, verbose, n_min_steeps, lineal_search, n_points, natom, do_energy_1, do_energy_2, do_energy_3,do_forces_1, do_forces_2, do_forces_3, make_coords


	implicit none
	verbose=0                      ; Force_cut=1D-5                     ;
	Energy_cut= 1D-4               ; minimzation_steep=5D-2             ;
	n_min_steeps = 10              ; lineal_search=.true.               ;
	n_points = 5                   ; natom=0                            ;
	do_energy_1 = ""               ; do_energy_2 = ""                   ;
	do_energy_3 = ""               ; do_forces_1 = ""                   ;
	do_forces_2 = ""               ; do_forces_3 = ""                   ;
	make_coords = ""
	end subroutine defaults_scarlet

	subroutine read_input(inputFile)
	USE scarlett_mod, only :Force_cut, Energy_cut, minimzation_steep, n_points,n_min_steeps,lineal_search,natom,verbose, do_energy_1, do_energy_2, do_energy_3,do_forces_1, do_forces_2, do_forces_3, make_coords
	implicit none
	integer :: ios, iErr
	logical :: fileExists
	character(len=20), intent(in)  :: inputFile
	namelist /jolie/ Force_cut, Energy_cut, minimzation_steep,n_points,n_min_steeps,lineal_search,natom,verbose, natom, do_energy_1, do_energy_2, do_energy_3,do_forces_1, do_forces_2, do_forces_3, make_coords
	Inquire(file = inputFile, exist = fileExists)

	write(*,*) "input file ", inputFile, "2"

	if(fileExists) then
	   write(*,*) "leo"
	   write(*,*) "natom", natom
	   open(unit = 100, file = inputFile, iostat = ios)
	   read(100, nml = jolie, iostat = iErr)
	   if(ierr.gt.0) stop 'Input error in scarlet namelist.'
	   close(100)
	   write(*,*) "natom", natom
	else
	   write(*,*)'File ', adjustl(inputFile), ' not found, using defaults.'
	endif
	write(*,nml=jolie)
	return
	end subroutine read_input



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Steepest Decend algorithm for geometry optimization using a generic
! program that calculate energy and forces
!
! Nicolas Foglia, 2017
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	subroutine do_steep(E)
	USE scarlett_mod, only : Force_cut, Energy_cut, minimzation_steep, n_min_steeps, natom, r, lineal_search, n_points, do_energy_1, do_energy_2, do_energy_3, do_forces_1, do_forces_2, do_forces_3, make_coords


	IMPLICIT NONE
	real*8, intent(inout) :: E !energy
	real*8 :: Emin, step_size, Fmax, d_E !Energy of previus steep, max displacement (Bohrs) of an atom in each steep, max |force| on an atoms, E(steep i)-E(steep i-1)
	integer :: n !steep number
	integer :: i,j !auxiliars

	double precision, allocatable, dimension(:) :: Energy !array of energy for lineal search
	double precision :: lambda !optimal displacement for minimice energy
	logical :: require_forces ! control force calculations in NO lineal search case
	double precision, dimension(natom, 3) :: r_scrach !positions scratch
	double precision, dimension(3, natom) :: gradient 
	logical :: stop_cicle, converged !for stop geometry optimization cicle

	gradient=0.d0

	if ( lineal_search ) write(*,*) "starting Steepest Descend with linear search"
	if ( .not. lineal_search ) write(*,*) "starting Steepest Descend without linear search"

	open (unit=12, file='traj.xyz') ! trajectory file
	open (unit=13, file='optimization.out') ! 
	write(13,4800)

	n=0
	stop_cicle=.false.
	converged=.false.
	step_size= minimzation_steep
	Fmax=0.d0
	lambda=0.d0
	require_forces=.true.


!first energy calculation
	CALL system(do_energy_1)
	CALL system(do_energy_2)
	CALL system(do_energy_3)
	open(unit=50,file="energy-buffer",form='unformatted',access='stream',status='old')
        read(50) E
        close(50)

	open(unit=52,file="position-buffer",form='unformatted',access='stream',status='old')
	do i=1,natom
	  do j=1,3
	    read(52) r(i,j)
	  enddo
	enddo
	close(52)



	Emin=E

	if ( lineal_search ) allocate(Energy(n_points))


	DO WHILE (n .lt. n_min_steeps .and. .not. stop_cicle)
	  write(13,4801) n, E, Fmax, step_size, lambda
	  call save_traj()

	  if ( require_forces) then
	    CALL system(do_forces_1)
	    CALL system(do_forces_2)
	    CALL system(do_forces_3)

	    open(unit=51,file="force-buffer",form='unformatted',access='stream',status='old')
            do i=1,natom
              do j=1,3
                read(51) gradient(j,i)
              enddo
            enddo
            close(51)
	    call max_force(gradient, Fmax) !find max force
	  end if

	  gradient=-gradient

	  if ( lineal_search ) then
	    call make_E_array(gradient, n_points,Energy, step_size, E,  Fmax) !obtain E(i) with r(i) = r_old - gradient/|gradient| * step_size *i
	    call line_search(n_points,Energy, step_size, lambda ) !predict best lambda that minimice E moving r(lambda) = r_old - gradient/|gradient| * lambda
	  else
	    r_scrach=r
	    lambda=step_size
	  end if


!	  write(*,*) "r"
!	  write(*,*) r
!	  write(*,*) "f"
!	  write(*,*) gradient

	  call move(lambda, Fmax,gradient)

!	  write(*,*)  "rnew"
!	  write(*,*) r

	 open(unit=53,file="position-buffer2")
	 DO i=1,NATOM
	    write(53,*) r(i,1), r(i,2), r(i,3)
	 END DO
	 close(53)
	call system(make_coords)


! energy calculation
	CALL system(do_energy_1)
	CALL system(do_energy_2)
	CALL system(do_energy_3)
	open(unit=50,file="energy-buffer",form='unformatted',access='stream',status='old')
	read(50) E
	close(50)



	  d_E=abs(E-Emin)
	  if ( lineal_search ) then
	    Emin=E
	    if (lambda .lt. 0.1d0 * step_size) step_size=step_size*0.1d0
	    if (lambda .gt. (dble(n_points)-0.1d0) * step_size) then
	      step_size=step_size*1.5d0
	      require_forces=.false.
	    else 
	      require_forces=.true.
	    end if
	  else
	    if (E .le. Emin) then
	      Emin=E
	      require_forces=.true.
	      step_size=step_size*1.2d0
	    else
	      r=r_scrach
	      require_forces=.false.
	      step_size=step_size*0.5d0
	    end if
	  end if

!convergence criterium
	  if ( step_size .lt. 1D-12) stop_cicle=.true.
	  if ( Fmax .lt. Force_cut .and. d_E .lt. Energy_cut .and. require_forces) then
	    stop_cicle=.true.
	    converged=.true.
	  end if
	  n=n+1
	END DO

        write(13,4801) n, E, Fmax, step_size, lambda
        call save_traj()

	IF (converged) THEN
	  WRITE(*,*) "system converged to required presition in ", n, " steeps"
	ELSE
	  WRITE(*,*) "system has not convergetd to required presition in ", n, " steeps"
	  WRITE(*,*) "Fmax component is: ", Fmax, "and energy diferences is: ", d_E
	END IF

	close(12)
	close(13)

	if ( lineal_search ) deallocate(Energy)
 4800   FORMAT(6x, "steep", 7x, "energy", 14x, "Fmax ", 9x, "steep size", 10x, "lambda")
 4801   FORMAT(2x,"!",2x, i4, 2x, 4(f16.10,2x))
 4802   FORMAT(2x,3(f16.10,2x))
	END subroutine


	subroutine max_force(gradient, Fmax)
	use scarlett_mod, only : natom
	implicit none
	double precision, intent(out) :: Fmax
	double precision, intent(in) :: gradient(3,natom)
	double precision :: F_i
	integer :: i
	Fmax=0.d0
        do i=1, natom
          F_i=gradient(1,i)**2d0 + gradient(2,i)**2d0 + gradient(3,i)**2d0
          F_i=sqrt(F_i)
          if (F_i .gt. Fmax) Fmax=F_i
        end do
	return
	end subroutine max_force


	subroutine make_E_array(gradient, n_points,Energy, step_size, E, Fmax)
!generate E(i) moving system r_new=r_old - i*step_size*gradient
	use scarlett_mod, only : r, natom, verbose, do_energy_1, do_energy_2, do_energy_3, make_coords
	implicit none
	double precision, intent(in) :: gradient(3,natom), step_size, Fmax
        double precision, intent(inout) :: E
	integer, intent(in) :: n_points
	double precision, intent(out) :: Energy(n_points)
	double precision, dimension(natom, 3) :: r_ini 
	double precision :: a
	double precision :: max_move
	integer :: i,j,k,l

!define step that normalice max force
	max_move=step_size/Fmax

!generates Energy(i)
	r_ini=r
	Energy(1)=E
	do i=2, n_points
	  do j=1, natom  !move positions along gradient direction
	    do k=1, 3
	      r(j,k)=r_ini(j,k)-max_move*dble(i)*gradient(k,j)
	    end do
	  end do
	 
	  open(unit=53,file="position-buffer2")
	  DO l=1,NATOM
	    write(53,*) r(l,1), r(l,2), r(l,3)
	  END DO
	  close(53)
	  call system(make_coords)
! energy calculation
          CALL system(do_energy_1)
          CALL system(do_energy_2)
          CALL system(do_energy_3)
	  open(unit=50,file="energy-buffer",form='unformatted',access='stream',status='old')
          read(50) E
          close(50)
	  Energy(i)=E
	end do

	if(verbose) then !ponerle luego un nivel alto de verbose
	  do i=1, n_points
	    write(13,5008) Energy(i)
	  end do
	end if

	r=r_ini

	return
 5008 FORMAT(2x,"Lineal search Energy ",2x,f16.10)
	end subroutine make_E_array



        subroutine line_search(n_points,Energy, step_size, lambda )
!if minimun value of E is E(n_points) return lambda = step_size * n_points
!if minimun value of E is E(1) return lambda = 0.d0
!else returns an expected lambda that minimice E(lambda) using a parabolic interpolation
        implicit none
        integer, intent(in) :: n_points
        double precision, dimension(n_points), intent(in) :: Energy
        double precision, intent(in) :: step_size
        double precision, intent(out) :: lambda
        integer :: i
        integer :: min_Energy_position
        double precision :: dE1, dE2, modif_fac

        if (n_points .le. 2) then
	  write(*,*) "wrong n_points in lineal search, n_points need to be > 2"
	  stop
	end if
!find min value in Energy elements
        min_Energy_position=1
        do i=2, n_points
           if (Energy(i) .lt. Energy(min_Energy_position)) min_Energy_position=i
        end do
        
        if (min_Energy_position .eq. 1) then
          lambda=0.d0
          return
        elseif (min_Energy_position .eq. n_points) then
          lambda=step_size*dble(n_points)
          return
        end if

        dE2=abs(Energy(min_Energy_position) - Energy(min_Energy_position+1))
        dE1=abs(Energy(min_Energy_position) - Energy(min_Energy_position-1))

        modif_fac=step_size*(dE2-dE1)/(dE1+dE2)
        lambda=step_size*dble(min_Energy_position) - 0.5d0 * modif_fac
        return
        end subroutine line_search


        SUBROUTINE move(lambda, Fmax, gradient) !calculate new positions
	USE scarlett_mod, only : r, natom
        IMPLICIT NONE
        INTEGER :: i,j
        DOUBLE PRECISION, INTENT(IN) :: lambda, Fmax, gradient(3,natom)
	DOUBLE PRECISION :: a
        a=lambda/Fmax
        DO i=1,natom
          DO j=1,3
             r(i,j)= r(i,j)-a*gradient(j,i)
          END DO
        END DO
	END SUBROUTINE move

	SUBROUTINE save_traj()
	USE scarlett_mod, only : IZ, r, natom
	IMPLICIT NONE
	integer :: i
	REAL*8 :: b_2_ang
	b_2_ang=0.52917725D0
        write(12,*) NATOM
        write(12,*)
        DO i=1,NATOM
          write(12,5000) IZ(i), r(i,1)*b_2_ang, r(i,2)*b_2_ang, r(i,3)*b_2_ang
	  write(*,5000) IZ(i), r(i,1)*b_2_ang, r(i,2)*b_2_ang, r(i,3)*b_2_ang
        END DO
 5000 FORMAT(2x,i2,2x,3(f16.10,2x))
	END SUBROUTINE save_traj
