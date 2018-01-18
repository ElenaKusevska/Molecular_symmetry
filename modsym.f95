!
!------------------------------------------------------------------
!Module containing the subroutines for determining
!each of the symmetry elements:
!-----------------------------------------------------------------
!
module symmetry_elements
implicit none
contains
!
!------------------------------------------------------------------
!Subroutine for inversion center:
!-----------------------------------------------------------------
!
subroutine inversion (mat1, n, logic, stdev)
implicit none
!
!-----------------------------------------------------------------
!This subroutine checks if the molecule has an inversion center
!The only input that it requires is the specification of the molecular
!structure in matrix form (mat1) and the number of atoms
!in the molecule
!----------------------------------------------------------------------
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), allocatable, dimension(:,:), intent(in) :: mat1 !input matrix
real(dp), allocatable, dimension(:,:) :: mat2 !inverted matrix
integer,intent(in) :: n !dimension of matrix
integer, intent(out) :: logic !1 or 0.
integer :: i, j !counters
real(dp), intent(out) :: stdev  !standard deviation from perfect symmetry
!
allocate (mat2(4,n)) 
!
do i = 1, n
    mat2(1,i) = mat1(1,i) !the first column contains the atomic masses
end do
!
do i = 2, 4 !it starts from 2 because the first column contains atomic mass
    do j = 1, n
        mat2(i,j) = -(mat1(i,j))
    end do
end do
!
call compare (mat1, mat2, n, logic, stdev)
!
end subroutine inversion 
!
!--------------------------------
!subroutine to check for the existence of a horisontal 
!reflection plane
!--------------------------------
!
subroutine horisontal(mat1, n, logic, stdev)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), allocatable, dimension(:,:), intent(in) :: mat1 !input matrix
real(dp), allocatable, dimension(:,:) :: mat2 !reflected matrix
integer,intent(in) :: n !dimension of matrix
integer, intent(out) :: logic !1 or 0.
integer :: i, j !counters
real(dp), intent(out) :: stdev
!
allocate (mat2(4,n))
!
do i = 1, n
    mat2(1,i) = mat1(1,i)
    mat2(2,i) = mat1(2,i)
    mat2(3,i) = mat1(3,i)
    mat2(4,i) = -mat1(4,i) !reflection in the xy - plane
end do
!
!I am assuming that when this subroutine is run, the molecule is already
!aligned along the main rotation axis
!
call compare (mat1, mat2, n, logic, stdev)
!
end subroutine horisontal
!
!-------------------------------------------
!subroutine for vertical reflection plane
!----------------------------------------------
!
subroutine vertical(mol1, ref_point, n, logic, stdev)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), dimension(3,1), intent(in) :: ref_point !SEA or midpoint
real(dp), dimension(3,3) :: Rz
real(dp) :: sinz, cosz
real(dp), allocatable, dimension(:,:), intent(in) :: mol1
real(dp), allocatable, dimension(:,:) :: mol2, mol3, coordinates, rotated
integer :: n !dimension
integer, intent(out) :: logic !0 or 1
integer :: i, j
real(dp), intent(out) :: stdev
!
allocate( coordinates(3,n), rotated(3,n), mol2(4,n), mol3(4,n) )
!
!----------------------------------------------------
!aligning refference point with x - axis:
!-----------------------------------------------------
!
do i = 1, 3
    do j = 1, n
        coordinates(i,j) = mol1(i+1,j) !transfer coordinates from molecule
                                        !matrix
    end do
end do
!
cosz = ref_point(1,1)/sqrt( (ref_point(1,1))**(2.)+(ref_point(2,1))**(2.) )
sinz = ref_point(2,1)/sqrt( (ref_point(1,1))**(2.)+(ref_point(2,1))**(2.) )
!
call z_rotation_matrix(cosz, sinz, Rz)
!
rotated = matmul(Rz, coordinates)
!        
do i = 1, 3
    do j = 1, n
        mol2(i+1,j) = rotated(i,j) !transfer rotated coordinates to
                                    !to molecule matrix. 
    end do
end do
!
do i = 1, n
    mol2(1,i) = mol1(1,i) !atomic masses
end do
!
!---------------------------------------------
!molecule after reflection along the xz axis:
!---------------------------------------------
!
do i = 1, n
    mol3(1,i) = mol2(1,i)
    mol3(2,i) = mol2(2,i)
    mol3(3,i) = - mol2(3,i) ! y - coordinate (xz plane)
    mol3(4,i) = mol2(4,i)
end do 
!
call compare(mol2, mol3, n, logic, stdev)
!
end subroutine vertical
!
!----------------------------------------------------------
!subroutine to check for secondary rotation axis
!-----------------------------------------------------------
!
subroutine secondary_rotation (mol, ref_point, n, logic, stdev)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), dimension(3,1), intent(in) :: ref_point !the point through which
                             !the axis should pass
real(dp), allocatable, dimension(:,:), intent(in) ::  mol
integer, intent(in) :: n
real(dp), allocatable, dimension(:,:) :: mol1, mol2
real(dp), dimension(3,1) :: vecz
integer, intent(out) :: logic !1 or 2
real(dp), dimension(3,3) :: Rz, Ry
integer, parameter :: r = 2 !second order rotation matrix
real(dp), intent(out) :: stdev 
!
allocate (mol1(4,n), mol2(4,n))
!
call align (ref_point, mol, n, vecz, mol1)
!
call rotation(mol1, n, r, mol2)
!
call compare (mol1, mol2, n, logic, stdev)
!
end subroutine secondary_rotation
!
!-----------------------------------------
!Subroutine for the S2n point groups:
!-----------------------------------------
!
subroutine rotoinversion (mol1, n, ord, logic, stdev)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
real(dp), allocatable, dimension(:,:), intent(in) :: mol1
real(dp), allocatable, dimension(:,:) :: mol2, mol3
integer, intent(in) :: ord !order of rotation
integer, intent(in) :: n !dimension of matrix
integer, intent(out) :: logic ! 1 or 0
integer :: i !counter
real(dp), intent(out) :: stdev
!
allocate(mol2(4,n), mol3(4,n))
!
call rotation (mol1, n, ord, mol2)
!
do i = 1, n
    mol3(1,i) = mol2(1,i)
    mol3(2,i) = mol2(2,i)
    mol3(3,i) = mol2(3,i)
    mol3(4,i) = -mol2(4,i) !reflection in the xy - plane
end do
!
call compare (mol1, mol3, n, logic, stdev)
!
end subroutine rotoinversion
!
!-----------------------------------------------------------------
!Subroutine to compare two matrices (molecular structures), to see
!if they are trully identical.
!-----------------------------------------------------------------
!
subroutine compare (matrix1, matrix2, n, prod, stdev)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
integer :: i, j, k !counters
!
real(dp), allocatable, dimension(:,:), intent(in) :: matrix1, matrix2
real(dp), allocatable, dimension(:,:) :: coordinates, coordinates1
integer, intent(in) :: n !dimension
integer, allocatable, dimension (:) :: comparison, B, BB ! comparison array
integer, intent(out) :: prod ! for the product for the comparison
real, allocatable, dimension(:) :: var !variance array
real(dp) :: variance,  numb, dist, dist1 !variance value
real(dp), intent(out) :: stdev
                                    ! deviation from perfect symmetry
!
allocate (comparison(n), B(n), BB(n), var(n))
allocate (coordinates(3,n), coordinates1(3,n))
!
do i = 1, 3
    do j = 1, n
        coordinates(i,j) = matrix1(i+1,j)
        coordinates1(i,j) = matrix2(i+1,j)
    end do
end do
do i = 1, n !set them to 0 at the start of the comparison.
    B(i) = 0
    comparison(i) = 0
    BB(i) = 0
end do
!
!--------------------------------------------------------------------
!check if for every atom in matrix1, there is an atom in matrix2 that
!has is of the same element and has the same set of coordinates:
!--------------------------------------------------------------------
!
call norm (n, coordinates)
call norm (n, coordinates1)
do i = 1, n
    do j = 1, n
        if (j == B(j)) cycle
        if (i == BB(i)) cycle
        if (abs(matrix1(1,i) - matrix2(1,j)) .LE. 0.015) then
            if (abs(matrix1(2,i) - matrix2(2,j)) .LE. 0.015) then
                if (abs(matrix1(3,i) - matrix2(3,j)) .LE. 0.015) then
                    if (abs(matrix1(4,i) - matrix2(4,j)) .LE. 0.015) then
                        comparison(i) = 1
                        B(j) = j
                        BB(i) = i
        dist = sqrt( (coordinates(i,1))**(2.) + coordinates(i,2)**(2.) + &
            coordinates(i,3)**(2.) )
        dist1 = sqrt( (coordinates1(j,1))**(2.) + coordinates1(j,2)**(2.) &
            + coordinates1(j,3)**(2.) )
                    var(i) = (dist - dist1)**(2.)
                    end if
                end if
            end if 
        end if
    end do
end do
!
!---------------------------------------------------------------------
!and now to see if all the matrix elements have been found to be equal
!and if so, what is the deviation from perfect symmetry:
!----------------------------------------------------------------------
!
prod = 1
do i = 1, n
    prod = prod * comparison(i)
end do 
!
if (prod == 1) then

    variance = 0
    numb = real(n) !number of atoms
    do i = 1, n
        variance = variance + var(i)
    end do
    variance = variance / numb
    stdev = sqrt(variance)
end if
!
end subroutine compare 
!
!----------------------------------------------------------
!Align the molecule so that the selcted moment of inertia 
!or refference point is along
!the z - axis of the cartesian coordinate system:
!-----------------------------------------------------------
!
subroutine align (vector, mol, n, vecz, mol1)
implicit none 
!
!-----------------------------------------------------------------
!This subroutine aligns the axis of interest with the z - axis.
!it does this by determining the sine and cosine for the
!axis in question based on its position relative to the coordinate system, 
!building the rotation matrices for rotation around the y - and x - 
!matrices and then applying them to every atom in the molecule. 
!
!This is equivalent to either rotating the entire
!molecule, or "moving space", i.e. rotating the coordinate system,
!to allow the moment of inertia to be aligned with the z - axis
!--------------------------------------------------------------
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), dimension(3,1), intent(in) :: vector
real(dp), dimension(3,1) :: vecxz !rotation around x-axis to 
                !move the vector
                    !to the xz plane
real(dp), dimension(3,1), intent(out) :: vecz ! the vector is aligned
                    !with z - axis after rotation around y - axis
integer, intent(in) :: n !number of atoms in the molecule
real(dp), allocatable, dimension(:,:), intent(in) :: mol !input molecule
real(dp), allocatable, dimension(:,:) :: coordinates, rotated1, rotated2
                    !working functions
real(dp) :: cosx, sinx, cosy, siny ! to rotate space about the x & y - axes
real(dp), dimension(3,3) :: Rx, Ry!Rotation about x & y - axes
real(dp), allocatable, dimension(:,:), intent(out) :: mol1
integer :: i, j !counter
!
allocate (coordinates(3,n))
allocate (rotated1(3,n), rotated2(3,n))
allocate (mol1(4,n))
!
do i = 1, 3
    do j = 1, n
        coordinates(i,j) = mol(i+1,j)
    end do
end do
!
if (vector(1,1) == 0) then
    if (vector(2,1) /= 0) then
        !
        !----------------------------
        ! [0,y,z]
        !--------------------------
        !
            !----------------------------------------------
            !find matrix for rotation about the x-axis
            !-----------------------------------------------
!
            cosx=vector(3,1)/sqrt( vector(2,1)**(2.) + vector(3,1)**(2.) )
            sinx=vector(2,1)/sqrt( vector(2,1)**(2.) + vector(3,1)**(2.) )
!
            call x_rotation_matrix (cosx, sinx, Rx)
!
            !--------------------------------------------------------
            !rotate space about the x - axis to bring the inertia vector in
            !the xz-plane:
            !----------------------------------------------------------
!
            rotated1 = matmul(Rx,coordinates)
            vecz = matmul(Rx,vector)
!
            do i = 1, n
                mol1(1,i) = mol(1,i) ! just transferring the atomic masses
            end do
!
            do i = 1, 3
                do j = 1, n
                    mol1(i+1,j) = rotated1(i,j)
                end do
            end do
!
    else if (vector(2,1) == 0) then
    !
    !-------------------------------------------
    ! [0,0,z]
    !-----------------------------------------------
    !
        vecz = vector
        mol1 = mol
    end if
end if
!
!
if (vector(2,1) == 0) then
    if (vector(1,1) /= 0) then
    !
    !----------------------------------------------
    ! [x, 0, z]
    !----------------------------------------------------
    !
        !---------------------------------------------------
        !find matrix for rotation about the y-axis
        !---------------------------------------------------
!
        cosy = vector(3,1)/sqrt(vector(1,1)**(2.) + vector(3,1)**(2.))
        siny = vector(1,1)/sqrt(vector(1,1)**(2.) + vector(3,1)**(2.))
!
        call y_rotation_matrix (cosy, siny, Ry)
!
        !-------------------------------------------------------------
        !rotate space around the y - axis, so that the rotation axis lies
        !along the positive z - axis
        !----------------------------------------------------------------
!
        rotated2 = matmul(Ry,coordinates)
        vecz = matmul(Ry,vector) !should have the form [0,0,z]
!
        do i = 1, n
            mol1(1,i) = mol(1,i) ! just transferring the atomic masses
        end do
!
        do i = 1, 3
            do j = 1, n
                mol1(i+1,j) = rotated2(i,j)
            end do
        end do
    end if 
end if
!
if (vector(1,1) /= 0 ) then
    if (vector(2,1) /= 0) then 
    !
    !--------------------------------------
    ! [x, y, z]
    !---------------------------------------
    !
        !----------------------------------------------
        !find matrix for rotation about the x-axis
        !-----------------------------------------------
!
        cosx = vector(3,1)/sqrt( vector(2,1)**(2.) + vector(3,1)**(2.) )
        sinx = vector(2,1)/sqrt( vector(2,1)**(2.) + vector(3,1)**(2.) )
!
        call x_rotation_matrix (cosx, sinx, Rx)
!        
        !--------------------------------------------------------
        !rotate space about the x - axis to bring the inertia vector in
        !the xz-plane:
        !----------------------------------------------------------
!
        rotated1 = matmul(Rx,coordinates)
        vecxz = matmul(Rx,vector)
!
        !---------------------------------------------------
        !find matrix for rotation about the y-axis
        !---------------------------------------------------
!
        cosy = vecxz(3,1)/sqrt( vecxz(1,1)**(2.) + vecxz(3,1)**(2.) )
        siny = vecxz(1,1)/sqrt( vecxz(1,1)**(2.) + vecxz(3,1)**(2.) )
!
        call y_rotation_matrix (cosy, siny, Ry)
!
        !-------------------------------------------------------------
        !rotate space around the y - axis, so that the rotation axis lies
        !along the positive z - axis
        !----------------------------------------------------------------
!
        rotated2 = matmul(Ry,rotated1)
        vecz = matmul(Ry,vecxz) !should have the form [0,0,z]
!
        do i = 1, n
            mol1(1,i) = mol(1,i) ! just transferring the atomic masses
        end do
!
        do i = 1, 3
            do j = 1, n
                mol1(i+1,j) = rotated2(i,j)
            end do
        end do
    end if
end if
!
end subroutine align
!
!---------------------------------------------------------
!Subroutine to perform rotations around the z - axis:
!---------------------------------------------------------
!
subroutine rotation(mol, n, r, molz)
implicit none
!
!------------------------------------------------------------
!This subroutine rotates a molecule around the z - axis of the coordinate
!system by an angle defined in the call of the subroutine. 
!-------------------------------------------------------------
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), allocatable, dimension(:,:), intent(in) :: mol
real(dp), allocatable, dimension(:,:) :: coordinates, rotated
real(dp), dimension(3,3) :: Rz ! rotation axis
integer, intent(in) :: n !dimension of molecule
integer, intent(in) :: r !order of rotation axis
real(dp), allocatable, dimension(:,:), intent(out) :: molz
real(dp), parameter :: pi = 3.141592653589793
real(dp) :: cosz, sinz, k ! = real(n)
integer :: i, j !counters
!
allocate (molz(4,n), coordinates(3,n), rotated(3,n))
!
do i = 1, 3
    do j = 1, n
        coordinates(i,j) = mol(i+1,j)
    end do
end do
!
k = real(r)
cosz = cos(((2.)*pi)/k)
sinz = sin(((2.)*pi)/k)
!
call z_rotation_matrix (cosz, sinz, Rz)
!
rotated = matmul(Rz,coordinates)
!
do i = 1, n
    molz(1,i) = mol(1,i) ! just transferring the atomic masses
end do
! 
do i = 1, 3
    do j = 1, n
        molz(i+1,j) = rotated(i,j)
    end do
end do
!
end subroutine rotation
!
!-----------------------------------------------------------
!subroutines to generate the rotation matrices:
!---------------------------------------------------------
!
subroutine x_rotation_matrix (cosx, sinx, Rx)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: cosx, sinx
real(dp), dimension(3,3), intent(out) :: Rx
!
Rx(1,1)=1
Rx(1,2)=0
Rx(1,3)=0
Rx(2,1)=0
Rx(2,2)=cosx
Rx(2,3)=-sinx
Rx(3,1)=0
Rx(3,2)=sinx
Rx(3,3)=cosx
!
end subroutine x_rotation_matrix
!
subroutine y_rotation_matrix (cosy, siny, Ry)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: cosy, siny
real(dp), dimension(3,3), intent(out) :: Ry
!
Ry(1,1)=cosy
Ry(1,2)=0
Ry(1,3)=-siny
Ry(2,1)=0
Ry(2,2)=1
Ry(2,3)=0
Ry(3,1)=siny
Ry(3,2)=0
Ry(3,3)=cosy
!
end subroutine y_rotation_matrix
!
subroutine z_rotation_matrix(cosz, sinz, Rz)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: sinz, cosz
real(dp), dimension(3,3), intent(out) :: Rz
!
Rz(1,1)=cosz
Rz(1,2)=-sinz
Rz(1,3)=0
Rz(2,1)=sinz
Rz(2,2)=cosz
Rz(2,3)=0
Rz(3,1)=0
Rz(3,2)=0
Rz(3,3)=1
!
end subroutine z_rotation_matrix
!
!-----------------------------------------
!subroutine to scale the longest distance to the CM
!in the molecule to 1:
!-----------------------------------------
!
subroutine norm (n, coord)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), allocatable, dimension(:,:), intent(inout) :: coord
real(dp), allocatable, dimension(:) :: distance
integer, intent(in) :: n
integer :: i, j 
real(dp) :: maxdist, scalefac !maximum distance from center of mass and
                        !scale factor
allocate( distance(n) )
!
do j = 1, n 
    distance(j) = sqrt( coord(1,j)**(2.) + coord(2,j)**(2.) &
        + coord(3,j)**(2.) )
end do
!
maxdist = maxval(distance)
!
scalefac = (1.)/maxdist
!
do i = 1, 3
    do j = 1, n
        coord(i,j) = coord(i,j)*scalefac
    end do
end do
!
end subroutine norm
!
end module symmetry_elements
