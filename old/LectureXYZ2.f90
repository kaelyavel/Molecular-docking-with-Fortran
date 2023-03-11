program LectureXYZ2

    
    use molecule_type
    use atom_type

    implicit none

    character(len=100):: ligne, fileName, line, buffer
    character(len=2):: atomname
    integer :: i, size, ok, cache, nbatoms
    type(molecule) :: mol, ligand
    type(atom) :: atomone
    real :: x,y,z

    
    ! Verify args
    if(iargc() /= 1) then
        print '(a)', "Please provide a WFN file"
        stop 10
    end if

    ! Open File provided
    call getarg(1,fileName)
    print '(/,a,a)', "File to read = ",trim(fileName)

    open(unit=10,file=fileName,iostat=ok,status='old')

    if(ok/=0) then
        print '(a,4x,a)', "Error during opening", fileName
        stop 20
    end if

    !call mol%write(fileName)   

    read(10, '(a8)', iostat=ok) buffer
    read(buffer, '(i4)') nbatoms
    print '(i4)', nbatoms

    read(10, '(a8)', iostat=ok) buffer
    print '(a8)', buffer

    do i=1, nbatoms
        read(10, '(a3, 3(f15.5))', iostat=ok) atomname, x,y,z
        print '(a2,f15.5, f15.5, f15.5)',atomname,x,y,z
    
    enddo






    ! Lire nombre d'atomes du ligand
    !read(10, '(a)', iostat=ok) ligne
    !read(ligne(1:5),'(i20)') size

    !print '(i20)', size


    !call init_mol(ligand, size)

    ! Ligand fourni en entrée ?
    !read(10, '(a)', iostat=ok) ligne
    !if (ligne /= 'ligand') then
        !print "Ce n'est pas un fichier de ligand : ", '(a,4x,a)', fileName
        !stop 30
    

    ! Read the content
    do i=1, 256
        !read(10, '(A8)', iostat=ok) ligne
        !read(ligne(3:4), '(A8)') element
        !read(ligne(10:18), '(f9,5)') x
        !read(ligne(27:35), '(f9,5)') y
        !read(ligne(41:50), '(f9,5)') z
        !init_atom(atom, element, (x,y,z))
        !add_atom(ligand, atom)
    enddo
    
    
    
    
    
end program LectureXYZ2
