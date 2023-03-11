module vdw_type
implicit none

    type atomvdw
        character(len=2), public :: atom_name
        real, public :: radius
            
        contains
        procedure,public :: init_vdw, print_vdw, read_vdw_file,  get_atom_radius
        
    end type atomvdw 

    contains
        subroutine init_vdw(atvdw, name, rad)

            class(atomvdw), intent(inout) :: atvdw 
            character(2), intent(in) :: name 
            real, intent(in) :: rad

            atvdw%atom_name = name
            atvdw%radius = rad
    
    
        end subroutine init_vdw
    
        subroutine print_vdw(atomvdwa)
            class(atomvdw), intent(in) :: atomvdwa

            print '(a8, f4.2)', atomvdwa%atom_name, atomvdwa%radius
            
        end subroutine print_vdw

        subroutine read_vdw_file(atomvdwa,  vdw_atom_arrays)
            class(atomvdw), intent(in) :: atomvdwa
            character(80) :: fileName, buffer, line
            integer :: i, ok
            class(atomvdw), dimension(53), intent(inout) :: vdw_atom_arrays
            
            if(iargc() /= 1) then
                print '(a)', "Please provide a VDW file"
                stop 10
            end if
        
            call getarg(1,fileName)
            print '(/,a,a)', "File to read = ",trim(fileName)

            open(unit=10,file=fileName,iostat=ok,status='old')
    
            if(ok/=0) then
            print '(a,4x,a)', "Error during opening", fileName
            stop 20
            end if

            ! Skip first 6 lines
            do i=1, 5
                read(10,'(a80)') line
                print '(a)',line
            enddo

            ! Read the content
            do i=1, 53
                read(10, '(A8)', iostat=ok) buffer
                if(ok/=0) then
                    exit
                endif
                
                vdw_atom_arrays(i)%atom_name = buffer(1:2)
                read(buffer(3:6), '(f4.2)')  vdw_atom_arrays(i)%radius
            enddo

            do i=1, 52
                print '(a,x,f4.2)',  vdw_atom_arrays(i)%atom_name, vdw_atom_arrays(i)%radius
            enddo
            
            
    
        end subroutine read_vdw_file

        subroutine get_atom_radius(atomvdwa, atom_n, vdw_atom_arrays, actual_radius)
            class(atomvdw), intent(in) :: atomvdwa
            integer :: i
            real, intent(out) :: actual_radius
            character(2), intent(in) :: atom_n
            character(2) :: buff
            class(atomvdw), dimension(52), intent(in) :: vdw_atom_arrays
           
            
            do i=1, 52
                buff = vdw_atom_arrays(i)%atom_name
                print '(a2, a2)', vdw_atom_arrays(i)%atom_name, atom_n

                if(buff == atom_n) then
                    print *, "equal"
                    actual_radius = vdw_atom_arrays(i)%radius
                    exit
                endif
    
                
            enddo
           
            
            
        end subroutine get_atom_radius

end module vdw_type