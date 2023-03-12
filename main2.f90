program main

    use vdw_type
    use molecule_type

    implicit none

    type(atomvdw) :: atomvd

    type(molecule) :: molecula

    type(molecule) :: ligand

    type(molecule), dimension(10) :: ligand_random

    type(atomvdw) :: readvdw

    type(atomvdw), dimension(53) :: vdw_array

    type(atom) :: temp

    character(100) :: fileName

    integer :: i, j, ok, nbligandrandom

    real :: x1, x2, y1, y2, z1, z2,d, teta, xx, yy, zz, xt, yt, zt

    real :: deltax, deltay, deltaz


    nbligandrandom = 10


    call getarg(1,fileName)

    call readvdw%read_vdw_file(fileName, vdw_array)

    call getarg(2,fileName)

    call ligand%read(fileName)    

    call getarg(3,fileName)

    call molecula%read(fileName)

    call molecula%set_radius(vdw_array)

    print'(f4.2)', molecula%atoms(1)%radius

    teta = 1

    !borne max
    ! 10

    deltax = 0.1
    deltay = 0.25
    deltaz = 0.75

    do i=1, nbligandrandom
        do j = 1, ligand%nb_atoms

            xt = ligand%atoms(j)%coordinates(1) + deltax*j
            yt = ligand%atoms(j)%coordinates(2) + deltay*i
            zt = ligand%atoms(j)%coordinates(3) + deltaz*(j/i)


            xx = xt * cos(teta)
            yy = yt * sin(teta)
            zz = zt
            

            temp%element = ligand%atoms(j)%element
            temp%coordinates(1) = xx
            temp%coordinates(2) = yy
            temp%coordinates(3) = zz


            allocate(ligand_random(i)%atoms(ligand%nb_atoms), stat=ok)
            
            ligand_random(i)%atoms(j) = temp

            teta = teta + 1


        enddo
    enddo
    do i=1, nbligandrandom
        do j = 1, ligand%nb_atoms
            print '(f9.5)', ligand_random(i)%atoms(j)%coordinates(1)
            print '(f9.5)', ligand_random(i)%atoms(j)%coordinates(2)
            print '(f9.5)', ligand_random(i)%atoms(j)%coordinates(3)

        enddo
    enddo


    do i=1, nbligandrandom
        deallocate(ligand_random(i)%atoms)  
    enddo



    deallocate(molecula%atoms)
    deallocate(ligand%atoms)

end program main