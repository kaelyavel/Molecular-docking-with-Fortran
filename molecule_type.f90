module molecule_type

    use atom_type
        implicit none
    
        type molecule
        !array of molecule atoms
        type(atom), dimension(:), allocatable, public :: atoms
        !number of atoms in molecule (size of atoms array)
        integer , public :: nb_atoms
    
        contains
        procedure :: print_mol, box, write, set_radius, read

    
        end type molecule
    
    contains
    
        

        ! function that prints a molecule
        subroutine print_mol(m)
            class (molecule),intent(in) :: m
            integer :: i

            do i=1, m%nb_atoms
                print '(a2)', m%atoms(i)%element
                print '(f9.5, f9.5, f9.5)', m%atoms(i)%coordinates(1), m%atoms(i)%coordinates(2), m%atoms(i)%coordinates(3)
            enddo 
        end subroutine
    

        ! function that reads an .xyz file defining a molecule

        subroutine read(m, fileName)
            class(molecule), intent(inout) :: m
            character(len=100), intent(in) :: filename
            
            character(len=100):: ligne
            character(len=3) :: element
            integer :: i, nb_atoms, ok, ok2
            
            real, dimension(3) :: coord

           
            print '(/,a,a)', "XYZ file to read = ", trim(fileName)
            
            
            open(unit=10,file=fileName,iostat=ok,status='old')
    

            if(ok/=0) then
             print '(a,4x,a)', "Error during opening", fileName
             stop 20
            end if

            
            read(10, '(a)', iostat=ok) ligne
            read(ligne(1:5),'(i4)') nb_atoms
            
            m%nb_atoms = nb_atoms
            allocate(m%atoms(nb_atoms))
           
            read(10, '(a)', iostat=ok) ligne
            
            
            do i=1, nb_atoms
                
                read(10, '(a80)', iostat=ok) ligne
                
                read(ligne(3:4), '(a80)') element
                read(ligne(10:18), '(f9.5)') coord(1)
                read(ligne(27:35), '(f9.5)') coord(2)
                read(ligne(41:50), '(f9.5)') coord(3)

                m%atoms(i)%element = element
                m%atoms(i)%coordinates = coord
                
                
            enddo 
        end subroutine

        ! set each atom of the molecule with the appropriate radius gotten for van der walls array
        subroutine set_radius(m, vdw_array)
            class(molecule), intent(inout) :: m
            type(atomvdw), dimension(53), intent(in) :: vdw_array
            character(2) :: buff,buff2
            integer :: i,j


            do j=1, m%nb_atoms
                buff2 = m%atoms(j)%element
                do i=1, 52
                    buff = vdw_array(i)%atom_name
                
    
                    if(buff == buff2) then
                        m%atoms(j)%radius = vdw_array(i)%radius
                        exit
                    endif
        
                    
                enddo
            enddo

        end subroutine set_radius



        ! takes ligand and molecule and forms a minimal box containing the both
        subroutine box(m, ligand, volume)
            class(molecule), intent(inout) :: m
            type(molecule), intent(in) :: ligand
            real, intent(inout) :: volume
    
            type(atom) :: a, a1
            integer :: i
            real :: max_x, min_x, max_y, min_y, max_z, min_z
            real :: x_side, y_side, z_side
            a1 = m%atoms(1)
    
            max_x = a1%coordinates(1)
            min_x = a1%coordinates(1)
            max_y = a1%coordinates(2)
            min_y = a1%coordinates(2)
            max_z = a1%coordinates(3)
            min_z = a1%coordinates(3)
    
            do i=2, m%nb_atoms
                a = m%atoms(i)
                if (max_x < a%coordinates(1)) then
                    max_x = a%coordinates(1)
                end if
                if (min_x > a%coordinates(1)) then
                    min_x = a%coordinates(1)
                end if
                if (max_y < a%coordinates(2)) then
                    max_y = a%coordinates(2)
                end if
                if (min_y > a%coordinates(2)) then
                    min_y = a%coordinates(2)
                end if
                if (max_z < a%coordinates(3)) then
                    max_z = a%coordinates(3)
                end if
                if (min_z > a%coordinates(3)) then
                    min_z = a%coordinates(3)
                end if
            end do
    
            do i=2, ligand%nb_atoms
                a = ligand%atoms(i)
                if (max_x < a%coordinates(1)) then
                    max_x = a%coordinates(1)
                end if
                if (min_x > a%coordinates(1)) then
                    min_x = a%coordinates(1)
                end if
                if (max_y < a%coordinates(2)) then
                    max_y = a%coordinates(2)
                end if
                if (min_y > a%coordinates(2)) then
                    min_y = a%coordinates(2)
                end if
                if (max_z < a%coordinates(3)) then
                    max_z = a%coordinates(3)
                end if
                if (min_z > a%coordinates(3)) then
                    min_z = a%coordinates(3)
                end if
            end do
    
            x_side = max_x - min_x
            y_side = max_y - min_y
            z_side = max_z - min_z
    
            volume = x_side * y_side * z_side
        end subroutine box

        ! write a file containing a molecule (/ligand) in the .xyz format
        subroutine write(m, filename)
            class(molecule), intent(in) :: m
            character(100), intent(in) :: fileName
            integer :: i, ok
            type(atom) :: a
            
            ! Open File provided
            !print '(/,a,a)', "File to write = ",trim(fileName)
            open(unit=10,file=fileName,iostat=ok,status='new')
            if(ok/=0) then
             print '(a,4x,a)', "Error during opening", fileName
             stop 20
            end if
    
            
            write(10, '(i2)', iostat=ok) m%nb_atoms
    
            
            write(10, '(a)', iostat=ok) "ligand"
    
            
            do i=1, m%nb_atoms
                a = m%atoms(i)
                write(10, '(a4,3f16.5)', iostat=ok) a%element, a%coordinates(1), a%coordinates(2), a%coordinates(3)
            enddo 
        endsubroutine
    
    end module molecule_type