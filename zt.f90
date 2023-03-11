program zt

    use vdw_type
    use molecule_type
    implicit none

    type(atomvdw) :: atomvd

    type(molecule) :: molecula

    type(atomvdw), dimension(52):: vdw_array

    character(len=2):: na
    integer :: i
    real :: rad, radinser

    na = "H "
    rad = 2.5

    !call atomvd%init_vdw(na, rad)

    !print '(f4.2)',atomvd%radius
    
    !call atomvd%read_vdw_file(vdw_array)

    !print '(f4.2)', vdw_array(1)%radius

    !call atomvd%get_atom_radius("H ", vdw_array, radinser)


    !print '(f4.2)', radinser

    call molecula%read()

    call molecula%print_mol2()



end program zt