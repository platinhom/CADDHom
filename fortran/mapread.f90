!* Author Zhixiong Zhao
!* 2014.12.19
!* Read the density in CCP4 Map file and translate to dat/sdat/dx
!* Support assigned polish iteration.
!* 2015.1.9 support the CCP4 format.
!*
!* Todo: 

program mapread
logical*1 :: readFromMap,writedat,writedx,writesdat,writeddat,removeneg
character*80 mapname,basename,extname,tmpstr,Spolish
integer*4 dotPlace,error,ios,i,j,k,doit,nPolish
integer*4 :: datatype,ix,iy,iz,nx,ny,nz
integer*8 :: ixyz,maxvalue,val
integer*4 :: headinfo(256)
integer*2,allocatable:: iv(:,:,:)
integer*4 :: Tpolish(255),Kpolish,ntime,lasttime,iti
real*4,allocatable :: rv(:,:,:)
real*8 dcel,maxv,val_max
real*8 scale
integer*2,allocatable:: bv(:,:,:)

!* Geometric flow
real*8  w1(-1:1),w2(-1:1),wxy(-1:1,-1:1),up(0:1),down(-1:0)
real*8  dx,dt
real*8,allocatable:: phi(:,:,:),phinew(:,:,:)
real*8  phix,phiy,phiz,dphi,phixx,phiyy,phizz,phixy,phixz,phiyz

!* For program interactive 
integer ( kind = 2 ) :: argc    !input number
character ( len = 80 ) :: argstr    !input parameter

writedat  = .true.
writedx   = .true.
writesdat = .false.
writeddat = .true.
removeneg=.true.; !write when removing data<0 in sdat/ddat/dx
dcel=0.01
scale=1000d0
!* for geometric flow
dx = 0.01;
dt = 1.0d-5;!dt<=dx^2 /2

print*, 'Use as:'
print*, 'Extract map information:'
print*, 'Mapread abc.map'
print*, 'Extract map information with ntime polishing:'
print*, 'Mapread abc.map 0 5 10 15 20'
print*, 'Only output the polishing at assigned time'
print*, 'Mapread abc.map 5 10 15 20'
PRINT*, '**********************************************'
!* For input arguments
!* 1st: file; 2nd: grid size; 3rd: sigma; 4th: scale
argc=iargc();
CALL GETARG(1, mapname)
if (argc < 1) then
    print *,'Error, Please assign the filename!'
    print *,'1:filename '
    stop
end if
!mapname='emd_1894.map'
dotPlace=index(mapname,'.',back = .true.);
basename=mapname(:dotPlace-1);
extname=mapname(dotPlace+1:);
if (.not. (extname(1:3) == 'map' .or. extname(1:4) =='ccp4')) then
	print *,'Error, Please assign the filename!'
    print *,'1:filename '
    stop
endif

readFromMap=.false.
Kpolish=0
j=0
if (argc>1) then
	do i=2,argc
		call GETARG(i,tmpstr)
		read(tmpstr,*) j
		if (j == 0 ) then
			readFromMap=.true.
            writedat=.false.
		else
			Kpolish=Kpolish+1
			Tpolish(Kpolish)=j
		end if 
	enddo
	!* sort from small to large
	call bubble_sort(Tpolish,Kpolish);
else !directly read
	readFromMap=.true.
endif


!! CCP4 Map header: 56 4-byte fields and space for ten 80 character(200 4-byte)
!! Reference: http://lsbr.niams.nih.gov/3demc/3demc_maplib.html
!! 1-3,int: total nx ny nz for grip; Dimensions(voxels)
!! 4,int: data type, normally =2 (real*4, 32bit) or =1(integer*2, 16bit)
!! 5-7,int: first position of x,y,z in map ; Origin(voxels)
!! 
!! 11-13,float: Cell Dimensions
!! 14-16,float: Cell Angles(alpha,beta,gamma)
!! 17-19,int: 
!! 20-22,float: Min/Max/Average density value
!! 26-34,float[3][3]: skew matrix

open (77, file = trim(mapname), iostat = ios, status='old', &
    	ACCESS='DIRECT', Form = 'Unformatted', RecL = 4)
	if (ios /= 0) stop
		read(77,rec=1) nx
		read(77,rec=2) ny
		read(77,rec=3) nz
		read(77,rec=4,iostat=ios) datatype
	if (ios /= 0) stop
	ixyz=nx*ny*nz
	print '(a,I5,I5,I5,I15,I2)','NX,NY,NZ,NX*NY*NZ,datatype: ',nx,ny,nz,ixyz,datatype
close (77)

open (78, file = trim(mapname), iostat = ios, &
    	ACCESS='DIRECT', Form = 'Unformatted', RecL = ixyz*4+1024)

allocate(rv(nx,ny,nz),stat=error)
if (error/=0) then
	write(*,*) 'No enough memory!'
	stop
endif
if (datatype==2) then
	read(78,rec=1,iostat = ios) headinfo(:),rv(:,:,:)
elseif (datatype == 1) then
	allocate(iv(nx,ny,nz))
	read(78,rec=1,iostat = ios) headinfo(:),iv(:,:,:)
	rv=dble(iv)
elseif (datatype == 0 ) then
	allocate(bv(nx,ny,nz))
	read(78,rec=1,iostat = ios) headinfo(:),bv(:,:,:)
	rv=dble(bv)
endif
if (ios /= 0) stop
close(78)

!let the  1 < iso <10000
maxv=maxval(rv);
if (maxv<=0.1) then
    scale=100000.0;
elseif (maxv<=1.0) then
    scale=10000.0;
elseif (maxv<=10.0) then
    scale=1000.0;
elseif (maxv<=100.0) then
    scale=100.0;
elseif (maxv<=1000.0) then
    scale=10.0;
! elseif (maxv<=10000.0) then
!     scale=1.0;
elseif (maxv<=100000.0) then
    scale=1.0;
elseif (maxv<=1000000.0) then
    scale=0.01;
else
    scale=1.0;
end if
maxvalue=ceiling(scale*maxv);

if (readFromMap) then
	if(writedat) then
		write(*,*) 'Write down the origin dat file....'
		open(88, file=trim(basename)//'.dat',iostat=ios)
			if (ios /= 0) stop
			write(88,'(I1)') 3
			write(88,'(I8)') nx
			write(88,'(I8)') ny
			write(88,'(I8)') nz
			do k=1,nz
			do j=1,ny
			do i=1,nx
				write(88,*) rv(i,j,k)
			end do
			end do
			end do
		close(88)
	endif
	if (writedx) then
		write(*,*) 'Write down the dx file....'
	    open(89,file=trim(basename)// '.dx', iostat=ios);
	    if (ios /= 0) stop
	    write(89,'(a,F10.7,a,F10.1,a)') '# isovalue of FRI surface. Max values for ',maxv, '*', scale, ' #'
	    write(89,'(a,I5,I5,I5)') 'object 1 class gridpositions counts',nx,ny,nz;
	    write(89,'(a)')            'origin    0.000   0.000    0.000';
	    write(89,'(a,F6.3,a)') 'delta    ',dcel, '   0.000    0.000'
	    write(89,'(a,F6.3,a)') 'delta     0.000  ',dcel,'    0.000'
	    write(89,'(a,F6.3)')   'delta     0.000   0.000   ',dcel
	    write(89,'(a,I5,I5,I5)') 'object 2 class gridpositions counts',nx, ny, nz;
	    write(89,'(a,I15,a)') 'object 3 class array type double rank 0 items ', ixyz,  ' data follows';
		do i=1,nx
		do j=1,ny
		do k=1,nz
			if (rv(i,j,k)<0 .and. removeneg) then
	            val= maxvalue;
	        else
	            val= ceiling(scale*(maxv-rv(i,j,k)));
	        endif
	        write(89,'(I8)') val; 
		end do
		end do
		end do
	    close(89)
	endif

	if(writesdat .or. writeddat) then
		write(*,*) 'Write down the sdat and/or modified dat file....'
		if (writesdat) then 
			open(90, file=trim(basename)//'.sdat',iostat=ios)
			if (ios /= 0) stop
			write(90,'(I1)') 3
		endif
		if (writeddat) then
			open(91, file=trim(basename)//'.ddat',iostat=ios)
			if (ios /= 0) stop
			write(91,'(I1)') 3
			write(91,'(I8)') nx
			write(91,'(I8)') ny
			write(91,'(I8)') nz
		endif
			do k=1,nz
			do j=1,ny
			do i=1,nx
			if (rv(i,j,k)<0 .and. removeneg) then
	            val= maxvalue;
	        else
	            val= ceiling(scale*(maxv-rv(i,j,k)));
	        endif
	        if (writesdat) write(90,'(I8)') i,j,k,val;
	        if (writeddat) write(91,'(I8)') val; 
			end do
			end do
			end do
		close(90)
		close(91)
	endif
endif

!* Do the geometric-flow based polishing
if (Kpolish>0) then
allocate(phi(nx,ny,nz),phinew(nx,ny,nz))
!!!!!!!!!!!!! assign the original isovalues to phi (exclude<0)
do k=1,nz
do j=1,ny
do i=1,nx
    if(rv(i,j,k)<0.0d0) then  !!!!!!!!!!!!!!!!!!!!! get rid of the negative part
        phi(i,j,k) = 0.0d0;
    else
        phi(i,j,k) = rv(i,j,k);
    end if
end do
end do
end do
phinew = phi;
write(*,*) "Before polish,maxval(phinew),minval(phinew)", maxval(phinew), minval(phinew)

!* Start parameter for geometric flow
   w1(-1) =-0.5d0/dx
   w1(0 ) = 0.0d0
   w1(1 ) = 0.5d0/dx

   w2(-1) = 1.0d0/dx/dx
   w2(0 ) =-2.0d0/dx/dx
   w2(1 ) = 1.0d0/dx/dx

   wxy(-1,-1) =  0.25d0/dx/dx
   wxy(-1, 0) =  0.0d0
   wxy(-1, 1) = -0.25d0/dx/dx
   wxy(0, -1) =  0.0d0
   wxy(0,  0) =  0.0d0
   wxy(0,  1) =  0.0d0
   wxy(1, -1) = -0.25d0/dx/dx
   wxy(1,  0) =  0.0d0
   wxy(1,  1) =  0.25d0/dx/dx

   up(0) = -1.0d0/dx
   up(1) =  1.0d0/dx

   down(0 ) = 1.0d0/dx
   down(-1) =-1.0d0/dx

nPolish=0
lasttime=0
do doit=1,Kpolish
	ntime=Tpolish(doit)
	write(Spolish,*) ntime
	if (ntime <= nPolish) CYCLE;
!******************************************************************************
	do iti=lasttime+1,ntime
		write(*,*)"polish iteration times",nPolish+1
		!Geometric flow
		do ix=1,nx
		do iy=1,ny
		do iz=1,nz
	        ixm1 = ix-1;				
			iym1 = iy-1;				
			izm1 = iz-1;

			ixp1 = ix+1;
			iyp1 = iy+1;
			izp1 = iz+1;

	        if(ix==1) ixm1 = nx;			
			if(iy==1) iym1 = ny;
	        if(iz==1) izm1 = nz;
		
	      	if(ix==nx) ixp1 = 1;
	        if(iy==ny) iyp1 = 1;
			if(iz==nz) izp1 = 1;
	                
			phix =w1(-1)*phi(ixm1,iy  ,iz  )+w1(1)*phi(ixp1,iy  ,iz  )
			phiy =w1(-1)*phi(ix  ,iym1,iz  )+w1(1)*phi(ix  ,iyp1,iz  )
			phiz =w1(-1)*phi(ix  ,iy  ,izm1)+w1(1)*phi(ix  ,iy  ,izp1)

			phizz=w2(-1)*phi(ix  ,iy  ,izm1)+w2(0)*phi(ix,iy,iz)    &
			    +w2( 1)*phi(ix  ,iy  ,izp1)
			phiyy=w2(-1)*phi(ix  ,iym1,iz  )+w2(0)*phi(ix,iy,iz)    &
			    +w2( 1)*phi(ix  ,iyp1,iz  )
			phixx=w2(-1)*phi(ixm1,iy  ,iz  )+w2(0)*phi(ix,iy,iz)    &
			    +w2( 1)*phi(ixp1,iy  ,iz  )
							   
			phixz =  wxy(-1,-1)*phi(ixm1,iy  ,izm1)    &
			       +wxy(-1, 1)*phi(ixm1,iy  ,izp1)    &
			       +wxy( 1,-1)*phi(ixp1,iy  ,izm1)    &
			       +wxy( 1, 1)*phi(ixp1,iy  ,izp1)
			phixy =  wxy(-1,-1)*phi(ixm1,iym1,iz  )    &
			       +wxy(-1, 1)*phi(ixm1,iyp1,iz  )    &
			       +wxy( 1,-1)*phi(ixp1,iym1,iz  )    &
			       +wxy( 1, 1)*phi(ixp1,iyp1,iz  )
			phiyz =  wxy(-1,-1)*phi(ix  ,iym1,izm1)    &
			       +wxy(-1, 1)*phi(ix  ,iym1,izp1)    &
			       +wxy( 1,-1)*phi(ix  ,iyp1,izm1)    &
			       +wxy( 1, 1)*phi(ix  ,iyp1,izp1)
							   
			dphi=(1.0d0+phix**2+phiy**2)*phizz+&
			    (1.0d0+phix**2+phiz**2)*phiyy+&
			    (1.0d0+phiy**2+phiz**2)*phixx+&                               
			(-2.0d0)*(phixy*phix*phiy+phixz*phix*phiz+phiyz*phiy*phiz)
					   		
			gram=(1.d0+phix**2+phiy**2+phiz**2)

			phinew(ix,iy,iz) = phi(ix,iy,iz) + dt*dphi/gram	 
	                   
			if(phinew(ix,iy,iz)<0.0d0) then  !!!!!!!!!!!!!!!!!!!!! get rid of the negative part
			   phinew(ix,iy,iz) = 0.0d0;
			end if        
		enddo 
		enddo 
		enddo   
    	phi=phinew
    	nPolish=nPolish+1
    enddo
    lasttime=nPolish;
!*****************************************************************************************************
	write(*,*)"After polish ntime,maxval(phinew),minval(phinew)",ntime, maxval(phinew),minval(phinew)
	val_max=maxval(phinew);
	maxvalue=ceiling(scale*val_max);

	if(writeddat) then
		write(*,'(a,I4)') 'Write down the scaled dat file for Polish: ',ntime
		open(88, file=trim(basename)//'_po'//trim(adjustl(Spolish))//'.ddat',iostat=ios)
		if (ios /= 0) stop
		write(88,'(I1)') 3
		write(88,'(I8)') nx
		write(88,'(I8)') ny
		write(88,'(I8)') nz
		do k=1,nz
		do j=1,ny
		do i=1,nx
            val= ceiling(scale*(val_max-phinew(i,j,k)));
			write(88,'(I8)') val; 
		end do
		end do
		end do
		close(88)
	endif

	if (writedx) then
		write(*,'(a,I4)') 'Write down the dx file for Polish: ',ntime
	    open(89,file=trim(basename)//'_po'//trim(adjustl(Spolish))//'.dx', iostat=ios);
	    if (ios /= 0) stop
	    write(89,'(a,F10.7,a,F10.1,a)') '# isovalue of FRI surface. Max values for ',maxv, '*', scale, ' #'
	    write(89,'(a,I5,I5,I5)') 'object 1 class gridpositions counts',nx,ny,nz;
	    write(89,'(a)')            'origin    0.000   0.000    0.000';
	    write(89,'(a,F6.3,a)') 'delta    ',dcel, '   0.000    0.000'
	    write(89,'(a,F6.3,a)') 'delta     0.000  ',dcel,'    0.000'
	    write(89,'(a,F6.3)')   'delta     0.000   0.000   ',dcel
	    write(89,'(a,I5,I5,I5)') 'object 2 class gridpositions counts',nx, ny, nz;
	    write(89,'(a,I15,a)') 'object 3 class array type double rank 0 items ', ixyz,  ' data follows';
		do i=1,nx
		do j=1,ny
		do k=1,nz
	        val= ceiling(scale*(val_max-phinew(i,j,k)));
	        write(89,'(I8)') val; 
		end do
		end do
		end do
	    close(89)
	endif
enddo
endif!!if (Kpolish>0) then

!!deallocate(rv,phi,phinew)
end program

subroutine bubble_sort(a,n)
  implicit none
  integer :: n,a(n)
  integer i,j, temp
  do i=n-1,1,-1  
    do j=1,i
    ! if a(j) > a(j+1) change place!
      if ( a(j) > a(j+1) ) then
        temp=a(j)
        a(j)=a(j+1)
        a(j+1)=temp
      end if
    end do
  end do
  return
end subroutine


! subroutine read_datfile(fileID, nx,ny,nz,data)
! 	integer*4 n,nx,ny,nz,ios,fileID
! 	real*8 val
! 	real*8 data(:,:,:)
! 	read (fileID, *, iostat=ios) n,nx,ny,nz
! end subroutine

