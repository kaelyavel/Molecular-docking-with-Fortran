program main

    use vdw_type
    use molecule_type

    implicit none

    
    ! Initialisation of all the variables used
    type(atomvdw) :: atomvd

    type(molecule) :: molecula

    type(molecule) :: ligand

    type(molecule) :: best_solution

    type(molecule), dimension(100000) :: ligand_random, ligand_no_collision

    type(atomvdw) :: readvdw

    type(atomvdw), dimension(53) :: vdw_array

    type(molecule), dimension(50) :: tab_ligand

    character(100) :: fileName, fileNameOutput

    integer :: i, j, n


    logical :: inside = .false.
    type(atom) :: temp

    integer :: ok, nbligandrandom, a, size_no_collision, vdIndex

    real :: x1, x2, y1, y2, z1, z2,d,xx, yy, zz, xt, yt, zt,volume, r1,r2, best_volume

    real :: deltax, deltay, deltaz

    real :: start_time, end_time, elapsed_time

    call cpu_time(start_time)

    a = 1

    ! Number of Rotation and Translation applied to the original ligand
    nbligandrandom = 100000


    call getarg(1,fileName)

    call readvdw%read_vdw_file(fileName, vdw_array)

    call getarg(2,fileName)

    call ligand%read(fileName)    

    call getarg(3,fileName)

    call molecula%read(fileName)

    call molecula%set_radius(vdw_array)
    call ligand%set_radius(vdw_array) 

    ! Double Loop to verify if the original ligand touches the site
    do i=1, ligand%nb_atoms

        x1 = molecula%atoms(i)%coordinates(1)
        y1 = molecula%atoms(i)%coordinates(2)
        z1 = molecula%atoms(i)%coordinates(3)
        
        do j=1, molecula%nb_atoms
            
          
            x2 = ligand%atoms(j)%coordinates(1)
            y2 = ligand%atoms(j)%coordinates(2)
            z2 = ligand%atoms(j)%coordinates(3)

            ! Distance between the atoms
            d = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            r1 = molecula%atoms(i)%radius
            r2 = ligand%atoms(j)%radius

            if(d <= r1 + r2) then
                inside = .true.
                exit
            endif

            !print '(f9.5)', d
            !print '(i4, i4)', i, j
        enddo
    enddo
    
    if(inside .eqv. .false.) then
        print *,"Le ligand est accepte"
    else
        print *,"Le ligand est rejete"
    endif

    vdIndex = 1

    ! Creation of the ligands resulting of the rotations and translations
    do i=1, nbligandrandom
        do j = 1, ligand%nb_atoms
            deltax = rand()*35
            deltay = rand()*35
            deltaz = rand()*35
            xt = ligand%atoms(j)%coordinates(1) + deltax
            yt = ligand%atoms(j)%coordinates(2) + deltay
            zt = ligand%atoms(j)%coordinates(3) + deltaz



            xx = xt * cos(rand()*360)
            yy = yt * sin(rand()*360)
            zz = zt


            temp%element = ligand%atoms(vdIndex)%element
            temp%coordinates(1) = xx
            temp%coordinates(2) = yy
            temp%coordinates(3) = zz

            temp%radius= ligand%atoms(vdIndex)%radius



            allocate(ligand_random(i)%atoms(ligand%nb_atoms), stat=ok)

            ligand_random(i)%atoms(j) = temp

            vdIndex = vdIndex + 1

            if(vdIndex == ligand%nb_atoms + 1) then
                vdIndex = 1
            endif


        enddo
        ligand_random(i)%nb_atoms = ligand%nb_atoms
    enddo

    
    inside = .false.

    ! Triple loop to verify if the ligands generated collide with the site
    do n=1, nbligandrandom
        inside = .false.
        do i=1, ligand_random(n)%nb_atoms

            x1 = ligand_random(n)%atoms(i)%coordinates(1)
            y1 = ligand_random(n)%atoms(i)%coordinates(2)
            z1 = ligand_random(n)%atoms(i)%coordinates(3)
            
            do j=1, molecula%nb_atoms
                
              
                x2 = molecula%atoms(j)%coordinates(1)
                y2 = molecula%atoms(j)%coordinates(2)
                z2 = molecula%atoms(j)%coordinates(3)

                
    
                d = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                r1 = molecula%atoms(j)%radius
                
                r2 = ligand_random(n)%atoms(i)%radius
               


                if(d <= r1 + r2) then
                    inside = .true.
                    exit
                endif
                
            enddo

            if(inside .eqv. .true.)then
                exit
            endif
        enddo
        
        if(inside .eqv. .false.) then
            ! If the ligand does not collide we keep it in a new tab that contains lignads that not collide
            ligand_no_collision(a) = ligand_random(n)
            a = a + 1
        endif    
    enddo

    size_no_collision = a - 1

    best_solution = ligand_no_collision(1)
    
    call molecula%box(ligand_no_collision(1),volume)
    best_volume = volume

    ! Iterations through the ligand_no_collision tab to extract the best_solution
    do i=2, size_no_collision
        call molecula%box(ligand_no_collision(i),volume)
        if(volume < best_volume) then 
            best_solution = ligand_no_collision(i)
            best_volume = volume
        endif
    enddo

    !print'(f12.4)', best_volume

    !call best_solution%print_mol2()


    fileNameOutput = "best_ligand.xyz"
    
    ! Writing of the best ligand into a xyz file
    call best_solution%write(fileNameOutput)

    call cpu_time(end_time)
    elapsed_time = end_time - start_time

    print'("Temps ecoule : ", f10.5,"s")',elapsed_time 

    do i=1, nbligandrandom
        deallocate(ligand_random(i)%atoms)  
    enddo

    

    deallocate(molecula%atoms)
    deallocate(ligand%atoms)

end program main