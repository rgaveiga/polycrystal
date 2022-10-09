! Module containing global variables
module global
    integer*8               :: seed
    real*8                  :: r 
end module

! Main program
program grain_builder
    use global
	
    implicit none

    character*10			:: lattice
    character*10			:: fileformat
    character*4,allocatable	:: tcus(:),t(:),ti(:,:),tcurr(:)
    integer					:: atom_index=0,index0,index1,cellindex
    integer					:: i,j,h,k,l,m,n
    integer					:: ngrains,natucell,natoms,natomsf,ncells,npartitions,cell,maxatcell=0
    integer,allocatable		:: natomscell(:),index(:,:),neighborscell(:,:),gindex(:)
    integer					:: incx,incy,incz
    integer,parameter		:: MIN_NGRAINS=8,NQUADRANTS=8,NNCELLS=27
    real*4					:: xtrans,ytrans,ztrans
    real*4					:: a0,dr0,dr1,dx0,dy0,dz0,dx1,dy1,dz1,dmin,gdist,a,b,c
    real*4					:: box_length,dispx,dispy,dispz
    real*4,allocatable		:: gbx(:),gby(:),gbz(:)
    real*4,allocatable		:: x(:),y(:),z(:),xi(:,:),yi(:,:),zi(:,:),xcurr(:),ycurr(:),zcurr(:),xcus(:),ycus(:),zcus(:)
    real*4,allocatable		:: xcell(:),ycell(:),zcell(:)
    real*4					:: dcell,max_rot_angle
    real*4,allocatable		:: phi(:),theta(:),psi(:)
    real*4					:: a11,a12,a13,a21,a22,a23,a31,a32,a33
    real*4,parameter        :: PARTITION_SIZE=10.0
    real*8,parameter		:: pi=3.14159265,degree2radian=0.0174532925
    logical					:: rnd_grain_nuclei,rnd_euler
    logical,allocatable		:: writable(:)
	
    open(10,file='polycrystal.tmp')
	
    natoms=0

    write(6,*) '---------------------------------------------------------------'
    write(6,*) '                       GRAIN_BUILDER                           '
    write(6,*) 'A program for generating polycrystals with three-dimensional   '
    write(6,*) '  structure developed by Dr. Roberto Gomes de Aguiar Veiga,    '
    write(6,*) '       at Universidade de São Paulo, Brazil (2012).            '
    write(6,*) '---------------------------------------------------------------'
    write(6,*)

    ! User input
1   write(6,*) '===> Randomly determine the position of grain seeds (true/false) and/or Euler angles (true/false)?'
    read(*,*) rnd_grain_nuclei,rnd_euler

2   write(6,*) '===> Lattice (bcc, fcc, diamond, hcp, or custom):'
    read(*,*) lattice

    if(.not.(lattice=='bcc'.or.lattice=='fcc'.or.lattice=='diamond'.or.lattice=='hcp'.or.lattice=='custom'))then
        write(6,*) '   ERROR: Lattice not supported!'
		
        goto 2
    endif

3   write(6,*) '===> Lattice parameter, in Angstroms:'
    read(*,*) a0

    if(a0<=0.0)then
        write(6,*) '   ERROR: The lattice parameter must be greater than zero!'

        goto 3        
    endif

4   if(lattice=='custom')then
        write(6,*) 'Provide the crystallographic parameters A, B, and C, in units of the lattice parameter:'
        read(*,*) a,b,c
		
		if(a<1.0.or.b<1.0.or.c<1.0)then
			write(6,*) '   ERROR: The crystallographic parameters must be greater than or equal to one!'
			
			goto 4
		endif

        a=a*a0
        b=b*a0
        c=c*a0
    else
    	if(lattice=='bcc'.or.lattice=='fcc'.or.lattice=='diamond')then
        	a=a0
        	b=a0
        	c=a0
    	else
    		a=a0
    		b=sqrt(3.0)*a0
    		c=sqrt(8.0/3.0)*a0
    	endif
    endif

5	write(6,*) '===> Box length, in Angstroms:'
	read(*,*) box_length
	
	if(box_length<=a.or.box_length<=b.or.box_length<=c)then
		write(6,*) '   ERROR: Box length must be greater than the unit cell length!'
		
		goto 5
	endif

	write(6,*) '   Box volume is ',box_length*box_length*box_length,' Angs^3.'

6	write(6,*) '===> Number of grains:'
	read(*,*) ngrains

	if(ngrains<MIN_NGRAINS)then	
		write(6,*) '   ERROR: At least ',MIN_NGRAINS,' grains are required!'
		
		goto 6
	endif

	allocate(gbx(ngrains))
	allocate(gby(ngrains))
	allocate(gbz(ngrains))
	allocate(phi(ngrains))
	allocate(theta(ngrains))
	allocate(psi(ngrains))

	if(lattice=='bcc')then
		natucell=2
	elseif(lattice=='fcc')then
		natucell=4
	elseif(lattice=='diamond')then
		natucell=8
	elseif(lattice=='hcp')then
		natucell==4
	elseif(lattice=='custom')then
7		write(6,*) '===> Number of atoms in the unit cell:'	
		read(*,*) natucell
		
		if(natucell<=0)then
			write(6,*) '   ERROR! Number of atoms must be greater than zero!'
			
			goto 7
		endif
	endif
	
	allocate(xi(natucell,ngrains))
	allocate(yi(natucell,ngrains))
	allocate(zi(natucell,ngrains))
	allocate(ti(natucell,ngrains))
	allocate(xcurr(natucell))
	allocate(ycurr(natucell))
	allocate(zcurr(natucell))
	allocate(tcurr(natucell))
	allocate(tcus(natucell))
	allocate(xcus(natucell))
	allocate(ycus(natucell))
	allocate(zcus(natucell))

	if(lattice=='custom')then
8		write(6,*) '===> Provide the atom type (string, 4 chars max) followed by the atomic coordinates,'
		write(6,*) '     in units of the lattice parameter, for each atom in the custom unit cell:'
			
		do i=1,natucell
			read(*,*) tcus(i),xcus(i),ycus(i),zcus(i)
				
			if(xcus(i)<0.0.or.ycus(i)<0.0.or.zcus(i)<0.0)then
				write(6,*) '   ERROR: Atomic coordinates of the custom unit cell must be greater than zero!'

				goto 8
			elseif(xcus(1)>0.0.or.ycus(1)>0.0.or.zcus(1)>0.0)then
				write(6,*) '   ERROR: The coordinates of the first atom in the custom unit cell must be (0,0,0)!'

				goto 8
			endif

			xcus(i)=xcus(i)*a0
			ycus(i)=ycus(i)*a0
			zcus(i)=zcus(i)*a0
		enddo
	endif
	
9	write(6,*) '===> Format of the output file (lammps or xyz):'
	read(*,*) fileformat
	
	if(.not.(fileformat=='lammps'.or.fileformat=='xyz'))then
		write(6,*) '   ERROR: Wrong file format!'
		
		goto 9
	endif

	if(rnd_grain_nuclei.or.rnd_euler)then
10      	write(6,*) '===> This code uses a random number generator. Please provide a seed for the generator (an integer):'
        	read(*,*) seed

		if(seed<=0)then
			write(6,*) '   ERROR: The seed must be a positive integer!'

			goto 10
		endif

		if(rnd_grain_nuclei)then
			write(6,*) '===> Minimum distance between two grain nuclei, in Angstroms:'
			read(*,*) gdist
			
            if(gdist<0.0)then
                write(6,*) '   ERROR: Minimum distance must be greater than or equal to zero!'
				
                goto 10
            endif
        endif
		
        if(rnd_euler)then
            write(6,*) '===> Maximum angle that can be randomly generated, in degrees:'
            read(*,*) max_rot_angle
			
            if(max_rot_angle<0.0.or.max_rot_angle>180.0)then
				write(6,*) '   ERROR: Maximum angle must be within the range 0-180!'
				
				goto 10
			endif
			
			max_rot_angle=max_rot_angle*degree2radian
		endif
    endif

    ! Manually define the position of grain nuclei and/or Euler angles
    if((.not.rnd_grain_nuclei).and.(.not.rnd_euler))then
        write(6,*) '===> You have to provide the positions of grain seeds (x,y,z), in Angstroms, and Euler' 
        write(6,*) '     angles (phi,theta,psi), in degrees, for each grain:'

        do i=1,ngrains
            read(*,*) gbx(i),gby(i),gbz(i),phi(i),theta(i),psi(i)
			
            phi(i)=phi(i)*degree2radian
            theta(i)=theta(i)*degree2radian
            psi(i)=psi(i)*degree2radian
        enddo
    elseif(.not.rnd_grain_nuclei.and.rnd_euler)then
        write(6,*) '===> You have to provide the positions of grain seeds (x,y,z) for each grain, in Angstroms:'

        do i=1,ngrains
            read(*,*) gbx(i),gby(i),gbz(i)
        enddo
    elseif(.not.rnd_euler.and.rnd_grain_nuclei)then
        write(6,*) '===> You have to provide the Euler angles (phi,theta,psi) for each grain, in degrees:'

        do i=1,ngrains
            read(*,*) phi(i),theta(i),psi(i)

            phi(i)=phi(i)*degree2radian
            theta(i)=theta(i)*degree2radian
            psi(i)=psi(i)*degree2radian
        enddo
    endif

 	! Use an arbitrarily defined cutoff distance to remove atoms too close to other atoms
11	write(6,*) '===> Minimum distance between two atoms at grain boundaries, in Angstroms:'
	read(*,*) dmin
	
	if(dmin<0.0)then
		write(6,*) 'Minimum distance must be greater than or equal to zero!'
	
		goto 11
	endif

	! Save information about the grain (position of the seed and Euler angles)
	open(12,file='grain_info.dat')

	write(12,*) '# Index X Y Z PHI THETA PSI'
	write(12,*) '# -------------------------'

    ! Randomly generate the positions of grain nuclei and/or grain orientations
    do i=1,ngrains
        if(rnd_grain_nuclei)then
            if(i==1)then
                gbx(i)=box_length/2.0
                gby(i)=box_length/2.0
                gbz(i)=box_length/2.0
            else
50              call rand
                gbx(i)=r*box_length

                call rand
                gby(i)=r*box_length
		
                call rand
                gbz(i)=r*box_length

                do j=1,i-1
                    dx0=gbx(j)-gbx(i)
                    dy0=gby(j)-gby(i)
                    dz0=gbz(j)-gbz(i)

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

                    if(dz0>box_length/2.0)then
                        dz0=dz0-box_length
                    elseif(dz0<-box_length/2.0)then
                        dz0=dz0+box_length
                    endif

                    dr0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)

                    if(dr0<gdist)then
                        goto 50
                    endif
                enddo
            endif
        endif

	if(rnd_euler)then
		call rand
		phi(i)=r*max_rot_angle

		call rand
		theta(i)=r*max_rot_angle

		call rand
		psi(i)=r*max_rot_angle
	endif

        ! Create the first cell, the one that will be repeated in 3 dimensions	
        if(lattice=='custom')then
            do j=1,natucell
				ti(j,i)=tcus(j)
				xi(j,i)=xcus(j)+gbx(i)
				yi(j,i)=ycus(j)+gby(i)
				zi(j,i)=zcus(j)+gbz(i)
			enddo
		else
			do j=1,natucell
				ti(j,i)='1'
			enddo
		
			xi(1,i)=gbx(i)
			yi(1,i)=gby(i)
			zi(1,i)=gbz(i)
	
			if(lattice=='bcc')then
				xi(2,i)=xi(1,i)+a0*0.5
				yi(2,i)=yi(1,i)+a0*0.5
				zi(2,i)=zi(1,i)+a0*0.5
			elseif(lattice=='fcc')then
				xi(2,i)=xi(1,i)+a0*0.5
				yi(2,i)=yi(1,i)+a0*0.5
				zi(2,i)=zi(1,i)			
				xi(3,i)=xi(1,i)+a0*0.5
				yi(3,i)=yi(1,i)
				zi(3,i)=zi(1,i)+a0*0.5
				xi(4,i)=xi(1,i)
				yi(4,i)=yi(1,i)+a0*0.5
				zi(4,i)=zi(1,i)+a0*0.5
			elseif(lattice=='diamond')then
				xi(2,i)=xi(1,i)+a0*0.25
				yi(2,i)=yi(1,i)+a0*0.25
				zi(2,i)=zi(1,i)+a0*0.25
				xi(3,i)=xi(1,i)+a0*0.5
				yi(3,i)=yi(1,i)+a0*0.5
				zi(3,i)=zi(1,i)
				xi(4,i)=xi(1,i)+a0*(3.0/4.0)
				yi(4,i)=yi(1,i)+a0*(3.0/4.0)
				zi(4,i)=zi(1,i)+a0*0.25
				xi(5,i)=xi(1,i)+a0*0.5
				yi(5,i)=yi(1,i)
				zi(5,i)=zi(1,i)+a0*0.5
				xi(6,i)=xi(1,i)
				yi(6,i)=yi(1,i)+a0*0.5
				zi(6,i)=zi(1,i)+a0*0.5
				xi(7,i)=xi(1,i)+a0*(3.0/4.0)
				yi(7,i)=yi(1,i)+a0*0.25
				zi(7,i)=zi(1,i)+a0*(3.0/4.0)
				xi(8,i)=xi(1,i)+a0*0.25
				yi(8,i)=yi(1,i)+a0*(3.0/4.0)
				zi(8,i)=zi(1,i)+a0*(3.0/4.0)
			elseif(lattice=='hcp')then
				xi(2,i)=xi(1,i)+a0*0.5
				yi(2,i)=yi(1,i)+a0*sqrt(3.0)/2.0
				zi(2,i)=zi(1,i)			
				xi(3,i)=xi(1,i)+a0*0.5
				yi(3,i)=yi(1,i)+a0*1.44337577280899203334
				zi(3,i)=zi(1,i)+a0*sqrt(8.0/3.0)/2.0
				xi(4,i)=xi(1,i)
				yi(4,i)=yi(1,i)+a0*sqrt(3.0)/3.0
				zi(4,i)=zi(1,i)+a0*sqrt(8.0/3.0)/2.0
			endif
		endif

		write(12,*) i,gbx(i),gby(i),gbz(i),phi(i),theta(i),psi(i)
	enddo

	close(12)	

	! In order to compare the distance between atoms in large simulations boxes efficiently, the simulation box should be partitioned
	npartitions=nint(box_length/PARTITION_SIZE)
	ncells=npartitions*npartitions*npartitions

	write(6,*) '   The system will be divided into ',ncells,' cells in order to make more efficient distance calculations'
	write(6,*) '   in later steps.'
	
	allocate(neighborscell(ncells,NNCELLS))
	allocate(natomscell(ncells))
	allocate(xcell(ncells))
	allocate(ycell(ncells))
	allocate(zcell(ncells))

	natomscell(:)=0
	
	dcell=box_length/npartitions
	l=0

	! Here the center of each cell is assigned
	do i=1,npartitions
		do j=1,npartitions
			do k=1,npartitions
				l=l+1
			
				xcell(l)=i*dcell-dcell/2.0
				ycell(l)=j*dcell-dcell/2.0
				zcell(l)=k*dcell-dcell/2.0
			enddo
		enddo
	enddo

	dr1=sqrt(3.0)*dcell

	! In the next loop, the first and second nearest neighbors of each cell are determined
	do i=1,ncells
		l=0

		do j=1,ncells
			dx0=xcell(j)-xcell(i)
			dy0=ycell(j)-ycell(i)
			dz0=zcell(j)-zcell(i)

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

			if(dz0>box_length/2.0)then
				dz0=dz0-box_length
			elseif(dz0<-box_length/2.0)then
				dz0=dz0+box_length
			endif

			dr0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)
			
			if(dr0<=dr1)then
				l=l+1
				
				neighborscell(i,l)=j
			endif
		enddo
	enddo

	! Grains grow from their respective nuclei
	do i=1,ngrains
		write(6,*) '   Generating grain',i
		write(6,*) '   From the grain nucleus at (',gbx(i),gby(i),gbz(i),')'
		write(6,*) '   Euler angles:',phi(i)/pi,'*pi,',theta(i)/pi,'*pi,',psi(i)/pi,'*pi'
		
		! ... and we have the corresponding direction cosine matrix (DCM)
		a11=cos(theta(i))*cos(psi(i))
		a12=-cos(phi(i))*sin(psi(i))+sin(phi(i))*sin(theta(i))*cos(psi(i))
		a13=sin(phi(i))*sin(psi(i))+cos(phi(i))*sin(theta(i))*cos(psi(i))
		a21=cos(theta(i))*sin(psi(i))
		a22=cos(phi(i))*cos(psi(i))+sin(phi(i))*sin(theta(i))*sin(psi(i))
		a23=-sin(phi(i))*cos(psi(i))+cos(phi(i))*sin(theta(i))*sin(psi(i))
		a31=-sin(theta(i))
		a32=sin(phi(i))*cos(theta(i))
		a33=cos(phi(i))*cos(theta(i))

        ! Fill a sphere around the grain nucleus with atoms
        do j=1,NQUADRANTS
            if(j==1)then
                incx=1
                incy=1
                incz=1

                write(6,*) '   12.5%...'
            elseif(j==2)then
                incx=-1
                incy=1
                incz=1

                write(6,*) '   25%...'
            elseif(j==3)then
                incx=-1
                incy=-1
                incz=1

                write(6,*) '   37.5%...'
            elseif(j==4)then
                incx=1
                incy=-1
                incz=1

                write(6,*) '   50%...'
            elseif(j==5)then
                incx=1
                incy=1
                incz=-1

                write(6,*) '   62.5%...'
            elseif(j==6)then
                incx=-1
                incy=1
                incz=-1

                write(6,*) '   75%...'
            elseif(j==7)then
                incx=-1
                incy=-1
                incz=-1

                write(6,*) '   87.5%...'
            elseif(j==8)then
                incx=1
                incy=-1
                incz=-1

                write(6,*) '   Done!'
            endif

            ! Begin: X direction
            do h=1,int(box_length/(2.0*a))				
                dispx=(h-1)*a*incx

                ! Begin: Y direction
                do k=1,int(box_length/(2.0*b))
                    dispy=(k-1)*b*incy

                    ! Begin: Z direction
                    do l=1,int(box_length/(2.0*c))
                        dispz=(l-1)*c*incz

                        ! Begin: Current atom
                        do m=1,natucell
							tcurr(m)=ti(m,i)
                            xtrans=xi(m,i)+dispx-gbx(i)
                            ytrans=yi(m,i)+dispy-gby(i)
                            ztrans=zi(m,i)+dispz-gbz(i)

                            xcurr(m)=(xtrans*a11+ytrans*a12+ztrans*a13)+gbx(i)
                            ycurr(m)=(xtrans*a21+ytrans*a22+ztrans*a23)+gby(i)
                            zcurr(m)=(xtrans*a31+ytrans*a32+ztrans*a33)+gbz(i)

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
							
                                if(zcurr(m)>box_length)then
                                    zcurr(m)=zcurr(m)-box_length
                                elseif(zcurr(m)<0.0)then
                                    zcurr(m)=zcurr(m)+box_length
                                endif
                            enddo
							
			! Check distances from the current grain				
			if(m==1)then
				dx0=xcurr(m)-gbx(i)
								dy0=ycurr(m)-gby(i)
								dz0=zcurr(m)-gbz(i)

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

								if(dz0>box_length/2.0)then
									dz0=dz0-box_length
								elseif(dz0<-box_length/2.0)then
									dz0=dz0+box_length
								endif

								dr0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)
					
								! Begin: Check distance from the other grains
								do n=1,ngrains
									if(n/=i)then																
										dx1=xcurr(m)-gbx(n)
										dy1=ycurr(m)-gby(n)
										dz1=zcurr(m)-gbz(n)

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

										if(dz1>box_length/2.0)then
											dz1=dz1-box_length
										elseif(dz1<-box_length/2.0)then
											dz1=dz1+box_length
										endif

										dr1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
			
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
                                dz0=zcurr(m)-zcell(n)

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

                                if(dz0>box_length/2.0)then
                                    dz0=dz0-box_length
                                elseif(dz0<-box_length/2.0)then
                                    dz0=dz0+box_length
                                endif

                                dr0=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)
                                
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
600                     enddo					! End: Current atom
500                 enddo						! End: Z direction					
400             enddo							! End: Y direction
300         enddo								! End: X direction
200     enddo									! End: Quadrant
100	enddo										! End: Grain

    write(6,*) '   Until now, ',natoms,' atoms were created. Some may be too close to other atoms.'
	
	deallocate(gbx)
	deallocate(gby)
	deallocate(gbz)
    deallocate(phi)
    deallocate(theta)
    deallocate(psi)
	deallocate(xi)
	deallocate(yi)
	deallocate(zi)
	deallocate(ti)
	deallocate(xcurr)
	deallocate(ycurr)
	deallocate(zcurr)
	deallocate(tcurr)
	deallocate(tcus)
	deallocate(xcus)
	deallocate(ycus)
	deallocate(zcus)

    ! Check the maximum number of atoms in a cell
    do i=1,ncells
        if(natomscell(i)>maxatcell) maxatcell=natomscell(i)
    enddo

    write(6,*) '   The maximum number of atoms per partition cell is ',maxatcell,'.'

	rewind(10)
		
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

                    if(dz0>box_length/2.0)then
                        dz0=dz0-box_length
                    elseif(dz0<-box_length/2.0)then
                        dz0=dz0+box_length
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
		write(11,*) '0.0',box_length
		write(11,'(a)') 'ITEM: ATOMS id type xs ys zs grain'
	elseif(fileformat=='xyz')then
		open(11,file='polycrystal.xyz')
		
		write(11,*) natomsf
		write(11,*) 'Box dimensions:',box_length,box_length,box_length
	endif
	
	do i=1,natoms
		if(writable(i))then
			if(fileformat=='lammps')then
				atom_index=atom_index+1
			
				write(11,*) atom_index,t(i),x(i)/box_length,y(i)/box_length,z(i)/box_length,gindex(i)
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
end program grain_builder

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
