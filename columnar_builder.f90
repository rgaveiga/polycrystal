! Module containing global variables
module global
	integer*8		:: seed
	real*8			:: r 
end module

! Main program
program columnar_builder
	use global
	
	implicit none

	character*10            :: fileformat
	character*4,allocatable :: t(:),ti(:,:),tcurr(:),tucell(:)
	integer                 :: atom_index=0,index0,index1,cellindex
	integer                 :: i,j,h,k,l,m,n,nrepz
	integer                 :: ngrains,natucell,natoms,natomsf,ncells,npartitions,cell,maxatcell=0
	integer,allocatable     :: natomscell(:),index(:,:),neighborscell(:,:),gindex(:)
	integer                 :: incx,incy,incz
	integer,parameter       :: MIN_NGRAINS=4,NQUADRANTS=4,NNCELLS=9
	real*4                  :: xtrans,ytrans,ztrans
	real*4                  :: latpar,a0(3),dr0,dr1,dx0,dy0,dz0,dx1,dy1,dmin,gdist
	real*4                  :: box_length,column_length,dispx,dispy,dispz
	real*4,allocatable      :: gbx(:),gby(:)
	real*4,allocatable      :: x(:),y(:),z(:),xi(:,:),yi(:,:),zi(:,:),xcurr(:),ycurr(:),zcurr(:),xcell(:),ycell(:)
	real*4,allocatable      :: xucell(:),yucell(:),zucell(:)
	real*4                  :: dcell,max_rot_angle
	real*4,allocatable      :: theta(:)
	real*4,parameter        :: PARTITION_SIZE=10.0
	real*8,parameter        :: pi=3.14159265,degree2radian=0.0174532925
	logical                 :: rnd_grain_nuclei,rnd_rot_angle
	logical,allocatable     :: writable(:)
	
	open(10,file='polycrystal.tmp')
	
	natoms=0

	write(6,*) '-----------------------------------------------------------'
	write(6,*) '                     COLUMNAR_BUILDER                      '
	write(6,*) 'A program for generating columnar polycrystals with tilt   '
	write(6,*) 'boundaries developed by Dr. Roberto Gomes de Aguiar Veiga  '
	write(6,*) '      at Universidade de São Paulo, Brazil (2012)          '
	write(6,*) '-----------------------------------------------------------'
	write(6,*)

	! User input
1	write(6,*) '===> Randomly determine the position of grain nuclei (true/false) and/or rotation angle (true/false)?'
	read(*,*) rnd_grain_nuclei,rnd_rot_angle

2	write(6,*) '===> Number of atoms in the unit cell:'
	read(*,*) natucell

	if(natucell<=0)then
		write(6,*) '   ERROR: Number of atoms in the unit cell must be a positive number greater than zero!'
		
		goto 2
	endif

	allocate(xucell(natucell))
	allocate(yucell(natucell))
	allocate(zucell(natucell))
	allocate(tucell(natucell))

3	write(6,*) '===> Lattice parameter, in Angstroms:'
	read(*,*) latpar

	if(latpar<=0.0)then
        	write(6,*) '   ERROR: The lattice parameter must be greater than zero!'

        	goto 3
    	endif

4	write(6,*) '===> Crystallographic parameters A, B, and C, in units of the lattice parameter:'
	read(*,*) a0(1),a0(2),a0(3)

	if(a0(1)<1.0.or.a0(2)<1.0.or.a0(3)<1.0)then
		write(6,*) '   ERROR: The crystallographic parameters must be greater than or equal to one!'

		goto 4
	endif

	a0(1)=a0(1)*latpar
	a0(2)=a0(2)*latpar
	a0(3)=a0(3)*latpar

5	write(6,*) '===> Provide the atom type (string, 4 chars maximum) and the coordinates (x,y,z), in units of'
	write(6,*) '     the lattice parameter, for each atom in the unit cell:'

	do i=1,natucell
		read(*,*) tucell(i),xucell(i),yucell(i),zucell(i)

		if(xucell(1)/=0.0.or.yucell(1)/=0.0.or.zucell(1)/=0.0)then
			write(6,*) '   ERROR: First atom of the unit cell must be at (0,0,0)!'
			
			goto 5
		endif

		xucell(i)=xucell(i)*latpar
		yucell(i)=yucell(i)*latpar
		zucell(i)=zucell(i)*latpar
	enddo

6	write(6,*) '===> Box length in X and Y, in Angstroms:'
	read(*,*) box_length
	
	if(box_length<a0(1).or.box_length<a0(2))then
		write(6,*) '   ERROR: Box length must be greater than the unit cell length in X and Y!'
		
		goto 6
	endif

7	write(6,*) '===> Number of unit cells repeated in the Z (textured) direction:'
	read(*,*) nrepz

	if(nrepz<=0)then
		write(6,*) '   ERROR: The number of cells in Z must be greater than zero!'

		goto 7
	endif

	column_length=real(nrepz)*a0(3)

	write(6,*) '   Box volume is ',box_length*box_length*column_length,' Angs^3.'

8	write(6,*) '===> Number of grains:'
	read(*,*) ngrains

	if(ngrains<MIN_NGRAINS)then	
		write(6,*) '   ERROR: At least ',MIN_NGRAINS,' grains are required!'
		
		goto 8
	endif

	allocate(gbx(ngrains))
	allocate(gby(ngrains))
	allocate(theta(ngrains))
	allocate(xi(natucell,ngrains))
	allocate(yi(natucell,ngrains))
	allocate(zi(natucell,ngrains))
	allocate(ti(natucell,ngrains))
	allocate(xcurr(natucell))
	allocate(ycurr(natucell))
	allocate(zcurr(natucell))
	allocate(tcurr(natucell))
	
9	write(6,*) '===> Format of the output file (lammps or xyz):'
	read(*,*) fileformat
	
	if(.not.(fileformat=='lammps'.or.fileformat=='xyz'))then
		write(6,*) '   ERROR: Wrong file format!'
		
		goto 9
	endif

10	if(rnd_grain_nuclei.or.rnd_rot_angle)then
		write(6,*) '===> This code uses a random number generator. Please provide a seed for the generator (an integer):'
		read(*,*) seed

		if(seed<=0)then
			write(6,*) '   ERROR: The seed must be a positive integer!'

			goto 10
		endif

		if(rnd_grain_nuclei)then
			write(6,*) '===> Minimum distance between two grain nuclei, in Angstroms:'
			read(*,*) gdist
			
			if(gdist<=0.0)then
				write(6,*) '   ERROR: The distance between grain nuclei must be grater than zero!'
				
				goto 10
			endif
		endif

		if(rnd_rot_angle)then
			write(6,*) '===> You have to provide the maximum rotation angle, in degrees:'
			read(*,*) max_rot_angle
			
			if(max_rot_angle<0.0.or.max_rot_angle>180.0)then
				write(6,*) '   ERROR: The maximum rotation angle must be within the range 0-180!'
				
				goto 10
			endif

			max_rot_angle=max_rot_angle*degree2radian
		endif
	endif

	! Manually define the position of grain nuclei and/or the rotation angle
	if((.not.rnd_grain_nuclei).and.(.not.rnd_rot_angle))then
		write(6,*) '===> You have to provide the positions of grain nuclei (x,y), in Angstroms, and the rotation angle about z,'
		write(6,*) '     in degrees, for each grain:'

		do i=1,ngrains
			read(*,*) gbx(i),gby(i),theta(i)

			theta(i)=theta(i)*degree2radian
		enddo
	elseif(.not.rnd_grain_nuclei.and.rnd_rot_angle)then
		write(6,*) '===> You have to provide the positions of grain nuclei (x,y) for each grain, in Angstroms:'

		do i=1,ngrains
			read(*,*) gbx(i),gby(i)
		enddo
	elseif(.not.rnd_rot_angle.and.rnd_grain_nuclei)then
		write(6,*) '===> You have to provide the rotation angle for each grain, in degrees:'

		do i=1,ngrains
			read(*,*) theta(i)

			theta(i)=theta(i)*degree2radian
		enddo
	endif

	! Randomly generate the positions of grain nuclei
	do i=1,ngrains
		if(rnd_grain_nuclei)then
			if(i==1)then
				gbx(i)=box_length/2.0
				gby(i)=box_length/2.0
			else
50				call rand
				gbx(i)=r*box_length

				call rand
				gby(i)=r*box_length
		
				do j=1,i-1
					dx0=gbx(j)-gbx(i)
					dy0=gby(j)-gby(i)

					if(dx0>box_length/2.0)then
						dx0=dx0-box_length
					elseif(dx0<-box_length/2.0)then
						dx0=dx0+box_length
					endif

					if(dy0>box_length/2.0)then
						dy0=dy0-box_length
					elseif(dy0<-box_length/2.0)then
						dy0=dy0+box_length
					endif

					dr0=sqrt(dx0*dx0+dy0*dy0)

					if(dr0<gdist)then
						goto 50
					endif
				enddo
			endif
		endif

		! Create the first cell, the one that will be repeated in 3 dimensions	
		do j=1,natucell
			ti(j,i)=tucell(j)
			xi(j,i)=gbx(i)+xucell(j)
			yi(j,i)=gby(i)+yucell(j)
			zi(j,i)=zucell(j)
		end do
	enddo
	
!!! A ser modificado	
	
	! In order to compare the distance between atoms in large simulations boxes efficiently, the simulation box should be partitioned
    npartitions=nint(box_length/PARTITION_SIZE)	
	ncells=npartitions*npartitions

	write(6,*) '   The system will be divided into ',ncells,' cells for comparing interatomic distances more efficiently.'
	
	allocate(neighborscell(ncells,NNCELLS))
	allocate(natomscell(ncells))
	allocate(xcell(ncells))
	allocate(ycell(ncells))

	natomscell(:)=0
	
	dcell=box_length/npartitions
	k=0

	! Here the center of each cell is assigned
	do i=1,npartitions
		do j=1,npartitions
			k=k+1
			
			xcell(k)=i*dcell-dcell/2.0
			ycell(k)=j*dcell-dcell/2.0
		enddo
	enddo

	dr1=sqrt(2.0)*dcell

	! In the next loop, the first and second nearest neighbors of each cell are determined
	do i=1,ncells
		l=0

		do j=1,ncells
			dx0=xcell(j)-xcell(i)
			dy0=ycell(j)-ycell(i)

			if(dx0>box_length/2.0)then
				dx0=dx0-box_length
			elseif(dx0<-box_length/2.0)then
				dx0=dx0+box_length
			endif

			if(dy0>box_length/2.0)then
				dy0=dy0-box_length
			elseif(dy0<-box_length/2.0)then
				dy0=dy0+box_length
			endif

			dr0=sqrt(dx0*dx0+dy0*dy0)
			
			if(dr0<=dr1)then
				l=l+1
				
				neighborscell(i,l)=j
			endif
		enddo
	enddo
	
	! Grains grow from their respective nuclei
	do i=1,ngrains
		write(6,*) '   Generating grain',i
		write(6,*) '   From the grain nucleus at (',gbx(i),gby(i),')'
		
		! The rotation angle is randomly defined...
		if(rnd_rot_angle)then
			call rand
			
			theta(i)=r*max_rot_angle
		endif

		write(6,*) '   Rotation angle:',theta(i)/pi,'*pi'

		! Fill a cylinder around the grain nucleus with atoms
		do j=1,NQUADRANTS
			if(j==1)then
				incx=1
				incy=1

				write(6,*) '   25%...'
			elseif(j==2)then
				incx=-1
				incy=1

				write(6,*) '   50%...'
			elseif(j==3)then
				incx=-1
				incy=-1

				write(6,*) '   75%...'
			elseif(j==4)then
				incx=1
				incy=-1

				write(6,*) '   Done!'
			endif

			! Begin: Z direction
			do l=1,nrepz
				dispz=(l-1)*a0(3)

				! Begin: X direction
				do h=1,int(box_length/(2.0*a0(1)))				
					dispx=(h-1)*a0(1)*incx

					! Begin: Y direction
					do k=1,int(box_length/(2.0*a0(2)))
						dispy=(k-1)*a0(2)*incy

						! Begin: Current atom
						do m=1,natucell			
							xtrans=xi(m,i)+dispx-gbx(i)
							ytrans=yi(m,i)+dispy-gby(i)
							ztrans=zi(m,i)+dispz

							xcurr(m)=(xtrans*cos(theta(i))-ytrans*sin(theta(i)))+gbx(i)
							ycurr(m)=(xtrans*sin(theta(i))+ytrans*cos(theta(i)))+gby(i)
							zcurr(m)=ztrans
							tcurr(m)=ti(m,i)

							! Put the atom back into the box, if necessary
							do n=1,2
								if(xcurr(m)>box_length)then
									xcurr(m)=xcurr(m)-box_length
								elseif(xcurr(m)<0.0)then
									xcurr(m)=xcurr(m)+box_length
								endif
							
								if(ycurr(m)>box_length)then
									ycurr(m)=ycurr(m)-box_length
								elseif(ycurr(m)<0.0)then
									ycurr(m)=ycurr(m)+box_length
								endif
								
								if(zcurr(m)>column_length)then
									zcurr(m)=zcurr(m)-column_length
								elseif(zcurr(m)<0.0)then
									zcurr(m)=zcurr(m)+column_length
								endif
							enddo
							
							! Check distances from the current grain				
							dx0=xcurr(m)-gbx(i)
							dy0=ycurr(m)-gby(i)

							if(dx0>box_length/2.0)then
								dx0=dx0-box_length
							elseif(dx0<-box_length/2.0)then
								dx0=dx0+box_length
							endif

							if(dy0>box_length/2.0)then
								dy0=dy0-box_length
							elseif(dy0<-box_length/2.0)then
								dy0=dy0+box_length
							endif

							dr0=sqrt(dx0*dx0+dy0*dy0)
					
							! Begin: Check distance from the other grains
							if(m==1)then
								do n=1,ngrains
									if(n/=i)then																
										dx1=xcurr(m)-gbx(n)
										dy1=ycurr(m)-gby(n)

										if(dx1>box_length/2.0)then
											dx1=dx1-box_length
										elseif(dx1<-box_length/2.0)then
											dx1=dx1+box_length
										endif

										if(dy1>box_length/2.0)then
											dy1=dy1-box_length
										elseif(dy1<-box_length/2.0)then
											dy1=dy1+box_length
										endif

										dr1=sqrt(dx1*dx1+dy1*dy1)
			
										! Check if this atom is closer to another grain nucleus
										if(dr1<dr0)then
											goto 500
										endif
									endif
700								end do          ! End: Check distance to other grains
							endif

							natoms=natoms+1

							! Begin: Assign a cell to the current atom
							dr1=box_length

							do n=1,ncells
								dx0=xcurr(m)-xcell(n)
								dy0=ycurr(m)-ycell(n)

								if(dx0>box_length/2.0)then
									dx0=dx0-box_length
								elseif(dx0<-box_length/2.0)then
									dx0=dx0+box_length
								endif

								if(dy0>box_length/2.0)then
									dy0=dy0-box_length
								elseif(dy0<-box_length/2.0)then
									dy0=dy0+box_length
								endif

								dr0=sqrt(dx0*dx0+dy0*dy0)
                               
								if(dr0<=dcell/2.0)then
									cell=n

									exit
								elseif(dr0<=dr1)then
									dr1=dr0
									cell=n
								endif
							enddo ! End: Assign a cell to the current atom
        
							natomscell(cell)=natomscell(cell)+1

							write(10,*) tcurr(m),xcurr(m),ycurr(m),zcurr(m),cell,i
600						enddo					! End: Current atom
500					enddo						! End: Y direction					
400				enddo							! End: X direction
300			enddo								! End: Z direction
200		enddo									! End: Quadrant
100	enddo										! End: Grain

	write(6,*) '   Until now, ',natoms,' atoms were created. Some may be too close to other atoms.'
	
	deallocate(gbx)
	deallocate(gby)
    deallocate(theta)
	deallocate(xi)
	deallocate(yi)
	deallocate(zi)
	deallocate(ti)
	deallocate(xucell)
	deallocate(yucell)
	deallocate(tucell)
	deallocate(xcurr)
	deallocate(ycurr)
	deallocate(zcurr)
	deallocate(tcurr)

    ! Check the maximum number of atoms in a cell
    do i=1,ncells
        if(natomscell(i)>maxatcell) maxatcell=natomscell(i)
    enddo

    write(6,*) '   The maximum number of atoms per partition cell is ',maxatcell,'.'

	rewind(10)
	
	! Use an arbitrarily defined cutoff distance to remove atoms too close to other atoms
11	write(6,*) '===> Minimum distance between two atoms, in Angstroms:'
	read(*,*) dmin
	
	if(dmin<0.0)then
		write(6,*) 'Minimum distance must be greater than zero!'
	
		goto 11
	endif
		
	allocate(x(natoms))
	allocate(y(natoms))
	allocate(z(natoms))
	allocate(t(natoms))
	allocate(gindex(natoms))
	allocate(writable(natoms))
    allocate(index(ncells,maxatcell))

	writable(:)=.true.
    index(:,:)=0
		
    ! Store in matrices all information about the atoms in the grains
	do i=1,natoms
		read(10,*) t(i),x(i),y(i),z(i),cellindex,gindex(i)

        do j=1,natomscell(cellindex)
            if(index(cellindex,j)==0)then
                index(cellindex,j)=i

                exit
            endif
        enddo
	enddo

	natomsf=natoms
	
    ! Check distances between atoms of adjacent cells
	write(6,*) '   Checking distances between atoms in adjacent cells...'

    ! Begin: Cell
    do i=1,ncells
        ! Begin: Atom in the cell
        do j=1,natomscell(i)
            index0=index(i,j)

            if(.not.writable(index0))then
                goto 900
            endif

            ! Begin: Neighbor cell
            do k=1,NNCELLS
                cellindex=neighborscell(i,k)

                if(cellindex<i)then
                    goto 1000
                endif

                ! Begin: Atom in the neighbor cell
                do l=1,natomscell(cellindex)
                    index1=index(cellindex,l)

                    if((cellindex==i.and.index1<=index0).or.(.not.writable(index1)))then
                        goto 1100
                    endif
					   
                    dx0=x(index1)-x(index0)
                    dy0=y(index1)-y(index0)
                    dz0=z(index1)-z(index0)

                    if(dx0>box_length/2.0)then
                        dx0=dx0-box_length
                    elseif(dx0<-box_length/2.0)then
                        dx0=dx0+box_length
                    endif

                    if(dy0>box_length/2.0)then
                        dy0=dy0-box_length
                    elseif(dy0<-box_length/2.0)then
                        dy0=dy0+box_length
                    endif

                    if(dz0>column_length/2.0)then
                        dz0=dz0-column_length
                    elseif(dz0<-column_length/2.0)then
                        dz0=dz0+column_length
                    endif

                    dr0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)

                    if(dr0<=dmin)then
                        writable(index1)=.false.
                        natomsf=natomsf-1
                    endif
1100            enddo                   ! End: Atom in the neighbor cell
1000        enddo                       ! End: Neighbor cell
900     enddo                           ! End: Atom in the cell
800 end do                              ! End: Cell
	
	! Write the coordinates
    write(6,*) '   Writing coordinates...',natomsf,'atoms left!'

	if(fileformat=='lammps')then
		open(11,file='polycrystal.lammpstrj')
	
		write(11,'(a)') 'ITEM: TIMESTEP'
		write(11,'(a)') '0'
		write(11,'(a)') 'ITEM: NUMBER OF ATOMS'
		write(11,*) natomsf
		write(11,'(a)') 'ITEM: BOX BOUNDS'
		write(11,*) '0.0',box_length
		write(11,*) '0.0',box_length
		write(11,*) '0.0',column_length
		write(11,'(a)') 'ITEM: ATOMS id type xs ys zs grain'
	elseif(fileformat=='xyz')then
		open(11,file='polycrystal.xyz')
		
		write(11,*) natomsf
		write(11,*) 'Box dimensions:',box_length,box_length,column_length
	endif
	
	do i=1,natoms
		if(writable(i))then
			if(fileformat=='lammps')then
				atom_index=atom_index+1
			
				write(11,*) atom_index,t(i),x(i)/box_length,y(i)/box_length,z(i)/column_length,gindex(i)
			elseif(fileformat=='xyz')then
				write(11,*) t(i),x(i),y(i),z(i),gindex(i)
			endif
		endif
	enddo
	
	deallocate(x)
	deallocate(y)
	deallocate(z)
	deallocate(t)
	deallocate(writable)
    deallocate(index)
    deallocate(natomscell)
    deallocate(neighborscell)
		
	close(10)
	close(11)

	stop 'Done!'
end program columnar_builder

! Subroutine that generates random numbers
subroutine rand()
	use global

	implicit none

	integer*8, parameter    :: m=714025
	integer*8, parameter    :: a=150889
	integer*8, parameter    :: c=1366
	real*8, parameter       :: d=30629.0
	real*8                  :: tmp

	seed=mod((a*seed+c),m)
     
	tmp=real(seed)/d

	r=tmp-int(tmp)
end subroutine
