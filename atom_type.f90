module atom_type

  use vdw_type
  implicit none

    type atom 
       character(len=2), public :: element
       real, dimension(3), public :: coordinates
       real, public :: radius

     contains
       procedure :: print_atom, init_atom
  
      
  
    end type atom
  
  contains
  
    subroutine init_atom(a, elt, coordinates)
      class (atom),intent(inout) :: a
      character(len=2), intent(in)::elt
      real, dimension(3), intent(in) :: coordinates
  
      a%element=elt
      a%coordinates=coordinates
  
    end subroutine init_atom

    subroutine print_atom(a)
      class (atom),intent(inout) :: a

      print '(3f9.5)', a%coordinates(1), a%coordinates(2), a%coordinates(3)

    end subroutine print_atom
    

  

  end module atom_type