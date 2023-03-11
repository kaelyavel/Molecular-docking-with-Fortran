program LectureVdw
    implicit none

    type atom_vdw

        character(2) :: atom_name
        real  :: radius

    end type atom_vdw

    integer :: ok, i
    character(80) :: line
    character(10) :: buffer
    character(len=128) :: fileName
    type(atom_vdw), dimension(52) :: atom_array

    ! Verify args
    if(iargc() /= 1) then
        print '(a)', "Please provide a WFN file"
        stop 10
    end if

    call getarg(1,fileName)
    print '(/,a,a)', "File to read = ",trim(fileName)

    fileName = radius_array(fileName)
    

end program LectureVdw

function radius_array(fileName)
    implicit none
 
    type atom_vdw

        character(2) :: atom_name
        real  :: radius

    end type atom_vdw

    character(50), intent(in) :: fileName
    integer :: ok, i
    character(80) :: line
    character(10) :: buffer
    type(atom_vdw), dimension(52):: atom_array
    type(atom_vdw), dimension(52):: radius_array

    open(unit=10,file=fileName,iostat=ok,status='old')
    
    if(ok/=0) then
     print '(a,4x,a)', "Error during opening", fileName
     stop 20
    end if

   ! Skip first 6 lines
    do i=1, 6
        read(10,'(a80)') line
        print '(a)',line
    enddo

! Read the content
    do i=1, 53
        read(10, '(A8)', iostat=ok) buffer
        if(ok/=0) then
            exit
        endif
        atom_array(i)%atom_name = buffer(1:2)
        read(buffer(3:6), '(f4.2)') atom_array(i)%radius
    enddo

    do i=1, 53
        print '(a,x,f4.2)',  atom_array(i)%atom_name, atom_array(i)%radius
    enddo
 
   radius_array = atom_array

end function radius_array