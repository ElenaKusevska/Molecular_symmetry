!-----------------------------------------------------
!Module containing a few other subroutines:
!-----------------------------------------------------
!
module other_routines
implicit none
contains
!
!---------------------------------------------
!Periodic system:
!---------------------------------------------
!
subroutine periodic (atoms, masses)
Implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(15)
!
    real(dp), dimension(118), intent(inout) :: masses
    character(len=3), dimension(118), intent(inout) :: atoms 
!
    atoms = (/'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ', 'O  ', &
'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ', 'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', &
'K  ', 'Ca ', 'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', &
'Cu ', 'Zn ', 'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', 'Rb ', 'Sr ', &
'Y  ', 'Zr ', 'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', &
'In ', 'Sn ', 'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ', 'La ', 'Ce ', &
'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', &
'Tm ', 'Yb ', 'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ', 'Pt ', &
'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', 'Fr ', 'Ra ', &
'Ac ', 'Th ', 'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ', &
'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', 'Rf ', 'Db ', 'Sg ', 'Bh ', 'Hs ', &
'Mt ', 'Ds ', 'Rg ', 'Cn ', 'Uut', 'Fl ', 'Uup', 'Lv ', 'Uus', 'Uuo'/)
!
    masses = (/1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, &
18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, &
39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, &
58.933, 58.693, 63.546, 65.38, 69.723, 72.63, 74.922, 78.96, 79.904, &
83.798, 85.468, 87.62, 89.906, 91.224, 92.906, 95.96, 97.91, 101.07, &
102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, &
131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 144.91, 150.36, &
151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, &
178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, &
204.38, 207.2, 208.98, 208.98, 209.99, 222.02, 223.02, 226.03, 227.03, &
232.04, 231.04, 238.03, 237.05, 244.06, 243.06, 247.07, 247.07, 251.08, &
252.08, 257.10, 258.10, 259.10, 262.11, 265.12, 268.13, 271.13, 270.0, &
277.15, 276.15, 281.16, 280.16, 285.17, 284.18, 289.19, 288.19, 293.0, &
294.0, 294.0/)
!
end subroutine periodic
!
!--------------------------------
!CPU time:
!--------------------------------
!
subroutine CPU (start)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) :: start
real(dp) :: finish, cputime
!
call CPU_time(finish)
cputime = finish - start
write(2,*) 'CPU_TIME:', cputime
!
end subroutine CPU
!
!--------------------------------
!Print standard deviation:
!--------------------------------
!
subroutine average_dev (st, element)
implicit none
!
integer, parameter :: dp = SELECTED_REAL_KIND(15)
real(dp), intent(in) ::  st
character(len=26), intent(in) :: element
character(len=50) :: doublechar
!
write(doublechar,*) st
write(2,*) 'average deviation of the atomic positions &
    from the symmetrically perfect case for the ' // &
    adjustl(trim(element)) // ': ' // adjustl(trim(doublechar))
!
end subroutine average_dev
!
end module other_routines
