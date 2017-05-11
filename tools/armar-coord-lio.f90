	program armar_coords_lio
	implicit none
	integer :: natom, i
	real*8, dimension(:,:), allocatable :: r
	integer, dimension(:), allocatable :: Iz
	real*8 :: b_2_ang
	b_2_ang = 0.52917725D0
	natom = 3
	allocate(r(natom,3), IZ(natom))
	
         open(unit=50,file="agua.xyz")
         DO i=1,NATOM
            read(50,*) IZ(i),r(i,1), r(i,2), r(i,3)
         END DO
	 close(50)

         open(unit=51,file="position-buffer2")
         DO i=1,NATOM
            read(51,*) r(i,1), r(i,2), r(i,3)
         END DO
	 close(51)

         open(unit=52,file="agua.xyz")
         DO i=1,NATOM
            write(52,4803) IZ(i),r(i,1)*b_2_ang, r(i,2)*b_2_ang, r(i,3)*b_2_ang
         END DO
         close(52)
 4802   FORMAT(2x,3(f16.10,2x))
 4803   FORMAT(2x,I2,2x,3(f30.15,2x))
	end program


