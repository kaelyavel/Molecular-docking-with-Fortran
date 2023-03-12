!******************************************************************
!
! CHPS 1002 PROJECT : MOLECULAR DOCKING WITH FORTRAN
!
! MADE BY : KHADRAOUI Mohamed El Bachir & EL BAHRAOUI IMADE
!
! TO COMPILE : gfortran .\vdw_type.f90 .\atom_type.f90 .\molecule_type.f90 .\main.f90 -o main
! 
! TO RUN : ./main .\rayons.vdw .\ligand.xyz .\molecule.xyz
!
!******************************************************************

program main

    use vdw_type
    use molecule_type

    implicit none

    
    ! Initializing all used variables

    !! Structures

    type(atomvdw) :: readvdw ! to read .vdw file

    type(molecule) :: molecula ! to read and manipulate molecule

    type(molecule) :: ligand ! to read and manipulate ligand

    type(molecule) :: best_solution ! to store best ligand variation

    type(atom) :: temp ! atom from temporary storing

    !! Array of structures

    type(molecule), dimension(100000) :: ligand_random, ligand_no_collision


    type(atomvdw), dimension(53) :: vdw_array


    !! variables

    integer :: i, j, n, a, vdIndex ! for loops

    integer :: nbligandrandom, size_no_collision ! for sizes

    integer :: ok ! for confirming

    real :: x1, x2, y1, y2, z1, z2,d,xx, yy, zz, xt, yt, zt ! for positions & rotations

    real :: volume, r1,r2, best_volume ! for storing

    real :: deltax, deltay, deltaz ! for random positioning

    real :: start_time, end_time, elapsed_time ! for time

    character(100) :: fileName, fileNameOutput ! file input/output

    logical :: inside = .false.
   
    ! Initializing end


    ! Start time
    call cpu_time(start_time)


    ! Setting variables

    !! Number of Rotation and Translation applied to the original ligand

    nbligandrandom = 100000

    ! Reading files

    call getarg(1,fileName)

    call readvdw%read_vdw_file(fileName, vdw_array)

    call getarg(2,fileName)

    call ligand%read(fileName)    

    call getarg(3,fileName)

    call molecula%read(fileName)

    ! Setting the radius for the molecule and ligand 
    
    call molecula%set_radius(vdw_array)

    call ligand%set_radius(vdw_array) 



    ! Double Loop to verify if the original ligand touches the site
    ! if its does 

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

            ! condition if one of ligand's atom is inside the molecule

            if(d <= r1 + r2) then 
                inside = .true.
                exit
            endif

        enddo
    enddo
    
    if(inside .eqv. .false.) then
        print '(a20)',"Le Ligand est accepte"
    else
        print  '(a20)',"Le ligand est rejete"
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
    a = 1

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

    print '(a60, i3)', "Nombre de ligands translates n'ayant pas de collisions avec la molecule : ", size_no_collision
    print '(a48 , i6)', "Nombre de ligands translates testes au totale : ", nbligandrandom
    
    call molecula%box(ligand_no_collision(1),volume)
    
    ! setting best volume as the first one
    best_volume = volume

    ! Iterations through the ligand_no_collision array to extract the best_solution
    do i=2, size_no_collision
        call molecula%box(ligand_no_collision(i),volume)
        if(volume < best_volume) then 
            best_solution = ligand_no_collision(i)
            best_volume = volume
        endif
    enddo

    print'(a20, f12.4)', "Le meilleur volume est : " ,best_volume

    !call best_solution%print_mol2()

    ! file to write te
    fileNameOutput = "best_ligand.xyz"
    
    ! Writing of the best ligand into a xyz file
    call best_solution%write(fileNameOutput)

    call cpu_time(end_time)

    elapsed_time = end_time - start_time

    print'("Temps ecoule : ", f10.5,"s")',elapsed_time 


    do i=1, nbligandrandom
        deallocate(ligand_random(i)%atoms)  
    enddo

    print'("Fichier avec la meilleur strcuture : ", a15)',fileNameOutput

    !Deallocation (otherwise crashes are expected)
    

    deallocate(molecula%atoms)
    deallocate(ligand%atoms)

end program main