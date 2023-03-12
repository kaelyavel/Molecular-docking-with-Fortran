program main

    use vdw_type
    use molecule_type

    implicit none

    type(atomvdw) :: atomvd

    type(molecule) :: molecula

    type(molecule) :: ligand

    type(atomvdw) :: readvdw

    type(atomvdw), dimension(53) :: vdw_array

    character(100) :: fileName

    integer :: i, j

    real :: x1, x2, y1, y2, z1, z2,d

    !type(atomvdw), dimension(52):: vdw_array

    !character(len=2):: na
    !integer :: i
    !real :: rad, radinser

    !na = "H "
    !rad = 2.5

    !call atomvd%init_vdw(na, rad)

    !print '(f4.2)',atomvd%radius
    
    !call atomvd%read_vdw_file(vdw_array)

    !print '(f4.2)', vdw_array(1)%radius

    !call atomvd%get_atom_radius("H ", vdw_array, radinser)


    !print '(f4.2)', radinser


    call getarg(1,fileName)

    call readvdw%read_vdw_file(fileName, vdw_array)

    call getarg(2,fileName)

    call ligand%read(fileName)

    

    call getarg(3,fileName)

    call molecula%read(fileName)

    
    call molecula%set_radius(vdw_array)

    print'(f4.2)', molecula%atoms(1)%radius


    !call molecula%print_mol2()

    !call ligand%print_mol2()



    ! x1 = molecula%atoms(1)%coordinates(1)
    ! y1 = molecula%atoms(1)%coordinates(2)
    ! z1 = molecula%atoms(1)%coordinates(3)
  
    ! x2 = ligand%atoms(1)%coordinates(1)
    ! y2 = ligand%atoms(1)%coordinates(2)
    ! z2 = ligand%atoms(1)%coordinates(3)

    !call molecula%atoms(1)%print_atom2()
    !call ligand%atoms(1)%print_atom2()
    
    !d = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    !print '(f9.5)', d
    ! print '(i4)', ligand%nb_atoms
    ! print '(i4)', molecula%nb_atoms

    ! do j=1, molecula%nb_atoms

    !     !print '(a2)', molecula%atoms(j)%element
    ! enddo

    ! do i=1, ligand%nb_atoms

    !     x1 = molecula%atoms(i)%coordinates(1)
    !     y1 = molecula%atoms(i)%coordinates(2)
    !     z1 = molecula%atoms(i)%coordinates(3)
        
    !     do j=1, molecula%nb_atoms
            
          
    !         x2 = ligand%atoms(j)%coordinates(1)
    !         y2 = ligand%atoms(j)%coordinates(2)
    !         z2 = ligand%atoms(j)%coordinates(3)

    !         d = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    !         print '(f9.5)', d
    !         !print '(i4, i4)', i, j
    !     enddo
    ! enddo


    deallocate(molecula%atoms)
    deallocate(ligand%atoms)

end program main