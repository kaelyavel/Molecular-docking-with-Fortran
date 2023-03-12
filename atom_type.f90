module atom_type

  use vdw_type
  implicit none

    type atom 
       character(len=2), public :: element
       real, dimension(3), public :: coordinates
       real, public :: radius

     contains
       procedure :: print_atom, init_atom, print_atom2
  
      
  
    end type atom
  
  contains
  
    subroutine init_atom(a, elt, coordinates)
      class (atom),intent(inout) :: a
      character(len=2), intent(in)::elt
      real, dimension(3), intent(in) :: coordinates
  
      a%element=elt
      a%coordinates=coordinates
  
    end subroutine init_atom
  
    subroutine print_atom(a, unit, iotype, v_list, iostat, iomsg)
      class (atom),intent(in) :: a
  
      integer, intent(in) :: unit
      character(len=*), intent(in) :: iotype
      integer, dimension(:), intent(in):: v_list
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg
  
      ! the '/' character add a carriage return
      write(unit=unit, fmt='(a2,3f8.3/)', iostat=iostat, iomsg=iomsg) trim(a%element),a%coordinates
  
    end subroutine print_atom

    subroutine print_atom2(a)
      class (atom),intent(inout) :: a

      print '(3f9.5)', a%coordinates(1), a%coordinates(2), a%coordinates(3)

    end subroutine print_atom2
    

  

  end module atom_type