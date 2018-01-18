!----------------------------------------------------
!Program to determine the symmetry point 
!group of a molecule.
!
!June 2014
!----------------------------------------------------
!
program symmetry    
use symmetry_elements
use other_routines
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
integer :: i, j, k, l, h, p, gg, f !counters
!
!-------------------------------------------------------------
!Variables for specifying the molecular structure:
!-------------------------------------------------------------
!
real(dp), allocatable, dimension(:) :: mass, x, y, z
character(len=3), allocatable, dimension(:) :: labels
integer :: m !number of atoms in the molecule
character(len=20) :: title
!
!---------------------
!Periodic table:
!---------------------
!
character(len=3), dimension(118) :: atoms
real(dp), dimension(118) :: masses
!
!---------------------------------
!Center of mass calculation:
!---------------------------------
!
real(dp), dimension(3) :: center !center of mass
real(dp) :: sumxa, sumxb, sumya, sumyb, sumza, sumzb
            !sums for the calculation
!
!-------------------------------------------------------------
!Calculation of the inertia tensor and principal moments
!of inertia: 
!-------------------------------------------------------------
!
real(dp), dimension(3,3) :: inert !inertia tensor
integer :: n, lda, info, lwork !for lapack dysyev routine
real(dp), dimension(3) :: w !eigenvalues
character :: jobz, uplo !for lapack dysyev routine
real(dp), allocatable, dimension(:,:) :: work !for lapack routine
real(dp) :: lwork1 !because the output for lwork is a real number
real(dp), dimension(3,1) :: eigenv1, eigenv2, eigenv3 !eigenvectors
character(len=50) :: workstring !for converting the numerical value
                                  !of the variable work(1,1) to
                                  !a string, for neater printing
character(len=25) :: infostring ! same as workstring, but for info
!
!----------------------------------------------------------------
!For the determination of the symmetry of the entire molecule:
!------------------------------------------------------------------
!
real(dp), allocatable, dimension(:,:) :: molecule
real(dp), allocatable, dimension(:,:) :: mola, molb !molecule after alignment
real(dp), dimension(3) :: vec1 !eigenvectors after alignemnt
real(dp), dimension(3) :: main_axis !direction of the main rotational axis.
integer :: r, rr ! order of rotation axis
integer :: test ! 0 or 1
integer, dimension(3) :: rot_order !vector that contains the maximum 
                                !order of the rotations around each of 
                                !the axes corresponding to the three 
                                !principal moments of inertia
integer :: order ! order of rotations
integer :: main_direction !main rotation order
reaL(dp) :: stdv !standard deviation
real(dp), dimension(3) :: standard_dev_array
character(len=26) :: horisontal1, vertical1, rotation1, rotation2, rotoinv
character(len=26) :: inversion1
!
!-----------------------------
!distance matrix and SEAs:
!-----------------------------
!
real(dp), allocatable, dimension(:,:) :: distance_matrix
!character, allocatable, dimension(:,:) :: dmatrix
real(dp), dimension(3,1) :: SEA1, SEA2, midpoint, SEA11, SEA21, midpoint1
real(dp), allocatable, dimension(:) :: D1, D2 !columns of distance matrix
integer, allocatable, dimension(:) :: C, D, T, AA !working f-tions
integer :: prod1 !comparison
!
!Cpu_time:
!
real(dp) :: start, finish, cputime
!
call CPU_time(start)
!
horisontal1 = 'horisontal reflection plane'
vertical1 = 'vertical reflection plane'
rotation1 = 'secondary rotation axis'
inversion1 = 'inversion center'
rotation2 = 'main rotation axis'
rotoinv = 'improper rotation axis'
!
!
!------------------------------------------------
!Construct a periodic table - two arrays, one
!containing the atomic symbols, the other,
!the atomic mass of the corresponding element:
!------------------------------------------------
!
call periodic (atoms, masses)
!
!------------------------
!!open output file:
!------------------------
!
open (unit=3, file='input', status='old', action='read')
open (unit=2, file='output', status='new', action='write')
!
!-----------------------------------------------------------------
!assigning the atoms and coordinates - user generated input:
!--------------------------------------------------------------
!
read(3,*) m
read(3,'(A20)') title
!
allocate (labels(m))
allocate (x(m))
allocate (y(m))
allocate (z(m))
allocate (mass(m))
!
do i = 1, m
    read(3,*) labels(i), x(i), y(i), z(i)
end do
!
write(2,*) 'input geometry for' 
write(2,*) trim(title) //  ':'
write(2,*) '--------------------------------------'
do i = 1, m
    write(2,*) labels(i), x(i), y(i), z(i)
end do
write(2,*)
!
!--------------------------------
!If the molecule is a diatiomic:
!--------------------------------
!
if (m==2) then
    if ( labels(1) == labels(2)) then
        write(2,*) 'the point group is Dh-inf'
        call CPU (Start)
        stop
    else if ( labels(1) /= labels(2)) then
        write(2,*) 'the point group is Cv-inf'
        call CPU (start)
        stop
    end if
end if
!
!-------------------------------------
!Assigning the masses:
!
!(I.e. read the masses of the atoms
!from the periodic table)
!-------------------------------------
!
do i = 1, m
    do k = 1, 118
        if (labels(i) == atoms(k)) then
            mass(i) = masses(k)
        end if
    end do 
end do 
! 
!-------------------------------
!Determine center of mass:
!-------------------------------
!
!x-coordinate:
!
sumxa = 0.0
do i = 1, m 
    sumxa = sumxa + mass(i)*x(i)
end do
!
sumxb = 0.0
do i =1, m
    sumxb = sumxb + mass(i)
end do
!
center(1) = sumxa/sumxb
!
!y - coordinate:
!
sumya = 0.0
do i = 1, m
    sumya = sumya + mass(i)*y(i)
end do
!
sumyb = 0.0
do i = 1, m
    sumyb = sumyb + mass(i)
end do
!
center (2) = sumya/sumyb
!
!z - coordinate:
!
sumza = 0.0
do i = 1, m 
    sumza = sumza + mass(i)*z(i)
end do
!
sumzb = 0.0
do i = 1, m
    sumzb = sumzb + mass(i)
end do 
!
center(3) = sumza/sumzb
!
!---------------------------------------
!place the center of mass at (0,0,0):
!---------------------------------------
!
do i = 1, m
    x(i) = x(i) - center(1)
end do
!
do i = 1, m 
    y(i) = y(i) - center(2)
end do
!
do i = 1, m
    z(i) = z(i) - center(3)
end do
!
!-----------------------------------
!Calculate the inertia tensor: 
!-----------------------------------
!
inert(1,1) = 0.0
do i = 1, m 
    inert(1,1) = inert(1,1) + mass(i)*( (y(i))**(2.) + (z(i))**(2.) )
end do
!
inert(1,2) = 0.0
do i = 1, m 
    inert(1,2) = inert(1,2) - mass(i)*y(i)*x(i)
end do
!
inert(2,1) = inert(1,2)
!
inert(1,3) = 0.0
do i = 1, m
    inert(1,3) = inert(1,3) - mass(i)*x(i)*z(i)
end do
!
inert(3,1) = inert(1,3)
!
inert(2,2) = 0.0
do i = 1,m 
    inert(2,2) = inert(2,2) + mass(i)*( (x(i))**(2.) + (z(i))**(2.) )
end do
!
inert(2,3) = 0.0
do i = 1, m
    inert(2,3) = inert(2,3) - mass(i)*y(i)*z(i)
end do
!
inert(3,2) = inert(2,3)
!
inert(3,3) = 0.0
do i = 1, m
    inert(3,3) = inert(3,3) + mass(i)*( (x(i))**(2.) + (y(i))**(2.) )
end do 
!
!-----------------------------------------------------------------
!Diagonalize the inertia tensor - determine the three principal
!components of the inertia tensor. 
!-----------------------------------------------------------------
!
n = 3
lda = 3
jobz = 'N'
uplo = 'U'
lwork = -1
!
allocate (work(1,3))
!
    call DSYEV(jobz, uplo, n, inert, lda, w, work, lwork, info)
    write(infostring,*)info
    write(workstring,*) work(1,1)
    write(2,*) 'info: ' // adjustl(trim(infostring)), ' lwork1: ' // adjustl(trim(workstring))
    lwork = int(work(1,1))
!
deallocate (work)
allocate (work(1,lwork))
!
    jobz = 'V'
    call DSYEV(jobz, uplo, n, inert, lda, w, work, lwork, info)
    write(2,*) 'eigenvalues of the inertia tensor:', w 
    write(2,*)  'info', info
    write(2,*)
!
deallocate (work)
!
!on input, the DSYEV subroutine takes the matrix inert, which represents
!the inertia tensor. On output it gives the same matrix (inert), but now
!it contains the eigenvectors of the inertia tensor in its columns. 
!
do i = 1, 3
    eigenv1(i,1) = inert(i,1)
    eigenv2(i,1) = inert(i,2)
    eigenv3(i,1) = inert(i,3)
end do

write(2,*) 'eigenvectors of the inertia tensor:'
write(2,*) '---------------------------------------'
write(2,*) 'eigenvector 1:', eigenv1
write(2,*) 'eigenvector 2:', eigenv2
write(2,*) 'eigenvector 3:', eigenv3
write(2,*)
!
!----------------------------
!cubic groups:
!----------------------------
!
if (dabs(w(1) - w(2)) .le. 0.5) then
    if (dabs(w(2) - w(3)) .le. 0.5) then
        if (dabs(w(1) - w(3)) .le. 0.5) then
            write (2,*) 'the molecule belongs to one of &
            the point groups with higher order (cubic) symmetry' 
            write (2,*)
            call CPU (start)
            Stop
        end if
    end if
end if
!
!----------------------------
!create molecule matrix:
!----------------------------
!
!the n - columns of this matrix represent the n - atoms
!
allocate (molecule(4,m))
!
do i = 1, m
    molecule(1,i) = mass(i) !row 1 contains the masses
    molecule(2,i) = x(i) !row 2 contains the x- coordinates
    molecule(3,i) = y(i) !row 3 contains the y - coordinates
    molecule(4,i) = z(i) !row 4 contains the z - coordinates
end do
!
!---------------------
!for linear molecules:
!---------------------
!
If ( (w(1) .le. 0.01) .or. (w(2) .le. 0.01) .or. (w(3) .le. 0.01) ) then
    call inversion (molecule, m, test, stdv)
    if (test == 1) then
        write(2,*) 'the point group of the molecule is: D-inf-h'
        call standard_dev (stdv, inversion1)
    else if (test /= 1) then
        write(2,*) 'the point group of the molecule is C-inf-v'
    end if
    call CPU (start)
    stop
end if
!
!-----------------------------------------------
!If the molecule is neither cubic, nor linear:
!-----------------------------------------------
!----------------------------------------------------------------
!Find the maximum order of rotation around every inertia vector
!----------------------------------------------------------------
!
allocate (mola(4,m), molb(4,m))
!
call align(eigenv1, molecule, m, vec1, mola) !align molecule with
                                             !inertia vector eigenv1
!
do i = 6, 1, -1
    r = i
    rot_order(1) = i
    call rotation(mola, m, r, molb)
    call compare(mola, molb, m, test, stdv)
    if (test == 1) exit !the maximum rotation order is recorded in
                        !the array rot_order
    standard_dev_array(1) = stdv
end do
!
call align(eigenv2, molecule, m, vec1, mola)
!
do i = 6, 1, -1
    r = i
    rot_order(2) = i
    call rotation(mola, m, r, molb)
    call compare(mola, molb, m, test, stdv)
    if (test == 1) exit
    standard_dev_array(2) = stdv
end do
!
call align(eigenv3, molecule, m, vec1, mola)
!
do i = 6, 1, -1
    r = i
    rot_order(3) = i
    call rotation(mola, m, r, molb)
    call compare(mola, molb, m, test, stdv)
    if (test == 1) exit
    standard_dev_array(3) = stdv
end do
!
order = maxval(rot_order)
!
!-----------------------------------------------
!in case there is no rotational symmetry:
!------------------------------------------------
!
!search for a horisontal plane:
if (order == 1) then
    call align(eigenv1, molecule, m, vec1, mola)
    call horisontal (mola, m, test, stdv)
    if (test == 1) then 
        write(2,*) 'the molecule belongs to the point group: Cs'
        call standard_dev (stdv, horisontal1)
        call CPU (start)
        stop
    else if (test == 0) then
        call align (eigenv2, molecule, m, vec1, mola)
        call horisontal (mola, m, test, stdv)
        if (test == 1) then
            write(2,*) 'the molecule belongs to the point group: Cs'
            call standard_dev (stdv, horisontal1)
            call CPU (start)
            stop
        else if(test == 0) then
            call align(eigenv3, molecule, m, vec1, mola)
            call horisontal (mola, m, test, stdv)
            if (test == 1) then
                write(2,*) 'the molecule belongs to the point group: Cs'
                call standard_dev (stdv, horisontal1)
                call CPU (start)
                stop
            else if (test == 0) then 
                !if there is no horisontal reflection plane, search for
                        !an inversion center:
                call inversion (molecule, m, test, stdv)
                if (test == 1) then
                    write(2,*) 'the molecule belongs to the point group: Ci'
                    call standard_dev (stdv, inversion1)
                    call CPU (start)
                    stop
                    if (test == 0) then
                        write(2,*) 'the molecule belongs to the point &
                            group: C1'
                        call CPU (start)
                        stop
                    end if
                end if
            end if
        end if
    end if
end if
!
!---------------------------------------------------------------
!when there is rotational symmetry, determine the direction of the
!main axis of rotation.
!---------------------------------------------------------------
!
write(2,*) 'the order of the main rotation axis is:', order
!
main_direction=maxloc(rot_order,1)!column of the inertia eigenvectors matrix
            !that corresponds to the main rotation axis. i.e. eigenvector
            !of the inertia matrix that corresponds to the main rotation
            !axis.
!
do i = 1, 3
    main_axis(i) = inert(i,main_direction)!the direction of the 
        !main rotational axis and of the moment of inertia eigenvector
        !corresponding to it
end do
!
stdv = standard_dev_array(main_direction) 
call standard_dev (stdv, rotation2)
!
!-------------------------------------------------------------------
!Build distance matrix in order to find two symetrically equivalent atoms
!-------------------------------------------------------------
!
call align(main_axis, molecule, m, vec1, mola) !From now on the molecule 
                                    !will be aligned with the z-axis
!
allocate (distance_matrix(m,m))
!
do i = 1, m
    do j = 1, m
        distance_matrix(i,j) = sqrt( (mola(2,i)-mola(2,j))**(2.) + &
            (mola(3,i)-mola(3,j))**(2.) + (mola(4,i)-mola(4,j))**(2.) )
    end do
end do
!
write(2,*)'Distance matrix:'
write(2,*)'----------------------------------------'
write(2,'(A7,20A7)') 'DM', (labels(i), i = 1, m)
!
do i = 1, m
        write(2,'(A7,20F7.3)') labels(i), (distance_matrix(i,j), j=1,m)
end do
!
!-----------------------------------------
!find two symetrically equivalent atoms:
!-----------------------------------------
!
do i = 1, 3
    SEA1(i,1) = 0.0
    SEA2(i,1) = 0.0
    midpoint(i,1) = 0.0
    SEA11(i,1) = 0.0
    SEA21(i,1) = 0.0
    midpoint1(i,1) = 0.0
end do

allocate (D1(m), D2(m), C(m), D(m), T(m), AA(m))

do i = 1, m
    D(i) = 0
    C(i) = 0
    T(i) = 0
    AA(i) = 0
end do

!First, select two columns (i and k) from the distance matrix:
determine_SEA: do i = 1, m
    AA(i) = i
    do k = 1, m
        if (i == k) cycle !because its the same atom
        if (k == AA(k)) cycle !so that I don't do the checking 
            !twice (1 with 2 and then 2 with 1 again).
            do h = 1, m
            D(h) = 0 !return the D variable used for comparison to 0.
            T(h) = 0
            C(h) = 0
        end do
    !Now, fill the D1 and D2 vectors used for the comparison with elements
    !from the i and k columns:
        do j = 1, m 
            do l = 1, m  
                D1(j) = distance_matrix(i,j)
                D2(l) = distance_matrix(k,l)
            end do
        end do
    !Compare these two columns of the distance matrix
    !to see if they have the same elements
        do gg = 1, m
            do p = 1, m
                if (p == C(p)) cycle
                if (gg == T(gg)) cycle
                if ( dabs(D1(gg)-D2(p)) .LE. 0.01 ) then
                    D(gg) = 1
                    C(p) = p
                    T(gg) = gg
                end if
            end do
        end do
        prod1 = 1 !to see if for every element in D1 a corresponding equal 
                    !element in D2 has been found.
        do f = 1, m
            prod1 = prod1 * D(f)
        end do 
        if (prod1 == 1) then
            !meaning all the distances for atom i and atom k
            !are equal
            do h = 1, 3
                SEA1(h,1) = mola(h+1,i)
                SEA2(h,1) = mola(h+1,k)
                midpoint(h,1) = (SEA1(h,1)+SEA2(h,1))/(2.)
            end do                        
            !to check if these SEA can be used to determine the secondary
            !rotation axes:
            if ( (dabs(SEA1(3,1)) .le. 0.01) .or. (dabs(SEA2(3,1)) .le. &
                0.01) .or. (dabs(midpoint(3,1)) .le. 0.01) ) then 
                                                     !meaining they are 
                                                     !lying on the xy plane
                if ( (dabs(midpoint(2,1)) .ge.0.01) .or. &
                    (dabs(midpoint(1,1)) .ge. 0.01) ) then
                                        !meaning the midpoint is not
                                        !lying on the main rotational axis.
                                        !or in the center of
                                        !the coordinate system
                    midpoint1 = midpoint
                    SEA11=SEA1
                    SEA21=SEA2
                    exit determine_SEA
                end if
            end if
        end if
    end do
end do determine_SEA
!
write(2,*) 'refference points'
write(2,*) '--------------------------'
write(2,*) 'symmetrically equivalent atoms:'
write(2,*) SEA1
write(2,*) SEA2
write(2,*) 'midpoint between the two symmetrically equivalent atoms:', &
    midpoint1
!
!---------------------------------------------------------------------
!when it is certain that the point group must be cyclic because we haven't
!found SEAs that fulfill the conditions at line 570:
!---------------------------------------------------------------------
!
if (midpoint1(1,1) == 0) then         
    if (midpoint1(2,1) == 0) then    
        if (midpoint1(3,1) == 0) then 
            call horisontal(mola, m, test, stdv)
            if (test == 1) then
                write(2,*) 'The symmetry group of the molecule is C', &
                    order, 'h'
                call standard_dev (stdv, horisontal1)
                call CPU(start)
                stop
            else if (test == 0) then
                call vertical(mola, midpoint, m, test, stdv)
                if (test == 1) then
                    write(2,*) 'The symmetry group is C', order, 'v'
                    call standard_dev (stdv, vertical1)
                    call CPU (start)
                    stop
                else if (test == 0) then
                    call vertical(mola, SEA1, m, test, stdv)
                    if (test == 1) then
                        write(2,*) 'The symmetry group is C', order, 'v'
                        call standard_dev (stdv, vertical1)
                        call CPU (start)
                        stop
                    else if (test == 0) then
                        call vertical(mola, SEA2, m, test, stdv)
                        if (test == 1) then
                            write(2,*) 'The symmetry group is C', order, 'v'
                            call standard_dev (stdv, vertical1)
                            call CPU (start)
                            stop
                        else if (test == 0) then
                            rr = 2*order
                            call rotoinversion (mola, m, rr, test, stdv)
                            if (test == 1) then
                                write(2,*) 'the symmettry group is S', rr
                                call standard_dev (stdv, rotoinv)
                                call CPU (start)
                            else if (test == 0) then
                                write(2,*) 'The symmetry group is C', order
                                call CPU(start)
                                stop
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end if
end if
!
!----------------------------------------------
!Determine symmetry when there is a possibility that the point group is
!dihedral:
!-------------------------------------------
!
call secondary_rotation(mola, SEA11, m, test, stdv)
if (test == 1) then
    write(2,*) 'the molecule has a secondary rotation axis.'
    call standard_dev (stdv, rotation1)
    call horisontal(mola, m, test, stdv)
    if (test == 1) then
        write(2,*) 'The symmetry group of the molecule is D', order, 'h'
        write(2,*) 'for the horisontal reflection plane:'
        call standard_dev (stdv, horisontal1)
        call CPU (start)
        stop
    else if (test == 0) then
        call vertical(mola, midpoint1, m, test, stdv)
        if (test == 1) then
            write(2,*) 'The symmetry group is D', order, 'd'
            write(2,*) 'for the vertical reflection plane:'
            call standard_dev (stdv, vertical1)
            call CPU(start)
            stop
        else if (test == 0) then
            write(2,*) 'The symmetry group is D', order
            call CPU (start)
            stop
        end if
    end if
else if (test == 0) then
    call secondary_rotation(mola, midpoint1, m, test, stdv)
    write(2,*) 'the molecule has a secondary rotation axis.'
    call standard_dev (stdv, rotation1)
    if (test == 1) then
        call horisontal(mola, m, test, stdv)
        if (test == 1) then
            write(2,*) 'The symmetry group of the molecule is D', order, 'h'
            call standard_dev (stdv, horisontal1)
            call CPU (start)
            stop
        else if (test == 0) then
            call vertical(mola, SEA11, m, test, stdv)
            if (test == 1) then
                write(2,*) 'The symmetry group is D', order, 'd'
                call standard_dev (stdv, vertical1)
                call CPU (start)
                stop
            else if (test == 0) then
                call vertical(mola, SEA2, m, test, stdv)
                if (test == 1) then
                    write(2,*) 'The symmetry group is D', order, 'd'
                    call standard_dev (stdv, vertical1)
                    call CPU (start)
                    stop
                else if (test == 0) then
                    write(2,*) 'The symmetry group is D', order
                    call CPU(start)
                    stop
                end if
            end if
        end if
    else if (test == 0) then
      call horisontal(mola, m, test, stdv)
         if (test == 1) then
             write(2,*) 'The symmetry group of the molecule is C', &
                 order, 'h'
             call standard_dev (stdv, horisontal1)
             call CPU (start)
             stop
         else if (test == 0) then
            call vertical(mola, midpoint, m, test, stdv)
            if (test == 1) then
                write(2,*) 'The symmetry group is C', order, 'v'
                call standard_dev (stdv, vertical1)
                call CPU (start)
            else if (test == 0) then
                call vertical(mola, SEA1, m, test, stdv)
                if (test == 1) then
                    write(2,*) 'The symmetry group is C', order, 'v'
                    call standard_dev (stdv, vertical1)
                    call CPU (start)
                else if (test == 0) then
                    call vertical(mola, SEA2, m, test, stdv)
                    if (test == 1) then
                        write(2,*) 'The symmetry group is C', order, 'v'
                        call standard_dev (stdv, vertical1)
                        call CPU (start)
                        stop
                    else if (test == 0) then
                        rr = 2*order
                        call rotoinversion (mola, m, rr, test, stdv)
                        if (test == 1) then 
                            write(2,*) 'The symmetry group is S', rr
                            call CPU (start)
                            stop
                        else if (test == 0) then
                            write(2,*) 'The symmetry group is C', order
                            call CPU (start)
                            stop
                        end if
                    end if
                 end if
             end if
         end if
     end if
 end if


end program symmetry
