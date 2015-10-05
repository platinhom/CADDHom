!* Last updated: 2014.12.16
!* Updated by Zhixiong Zhao
!* 2014.12.15: Add : Cell Listing
!* 2014.12.16: Add :support mol2/xyzr
!* 2014.12.26: Add control for using atom_number weight.
!* 
!* Todo: (1) Support pqr; (2) More input options


!* p=sum( wi * exp( -1 * wexp * (dis/sigma)^nexp )
!* In <Structure 16,295–307>, wi=Zi/(2*pi*sigma)^0.5. Zi is atom number
!* Normally, wi=1.0 enough since scale later except that using Zi.

include 'pdb_read.f90'

program main
    implicit none
  
    !!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !            Variables and Parameters
    !!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !* Set up the constant parameter,important for max atom number
    integer ( kind = 4 ), parameter :: atom_max = 300000
    integer ( kind = 4 ), parameter :: maxres = 50000
    integer ( kind = 4 ), parameter :: mxpatm = 18

    !* For program interactive 
    integer ( kind = 2 ) :: argc    !input number
    character ( len = 80 ) :: argstr    !input parameter       
    
    !* Control steps
    logical ( kind = 1 ) scubtop    !whether generate scubtop file
    logical ( kind = 1 ) writedxfile; !whether generate dx file
    logical ( kind = 1 ) celllisting; !whether using cell_listing method
    logical ( kind = 1 ) writelog;  !whether write log file with date
    logical ( kind = 1 ) useAtomNum; !whether use atom num to scale
    integer ( kind = 4 ) useAtomNumInt; 

    !* For Read from the PDB file and save informations
    character ( len = 80 )  ::  fileinput,filepdb,fname,extname
    integer ( kind = 4 ) ios !io state
    integer ( kind = 4 ) :: input_unit !index for input file
    integer ( kind = 4 )  :: atom_num   !total atoms number in the file
    real ( kind = 8 )  :: coord(atom_max,3) !Save all the x,y,z
    integer ( kind = 4 )  :: hiatom
    integer ( kind = 4 )  :: hires
    integer ( kind = 4 ) numchain
    integer ( kind = 4 ) numline
    integer ( kind = 4 ) numlost
    integer ( kind = 4 ) numres
    integer ( kind = 4 ) prtatm(maxres,mxpatm)
    character ( len = 3 ) resnam(atom_max)
    integer ( kind = 4 ) resnum(atom_max)
    character (len = 4 ) atom_names(atom_max)
    character (len = 2 ) elements(atom_max),element
    real ( kind = 8 ) tfactor(atom_max) !save the B-facter
    real ( kind = 8 ) xmax,xmin,ymax,ymin,zmax,zmin !min/max XYZ in origin data
    real*8 xleft,xright,yleft,yright,zleft,zright !right:max left:min value for XYZ with boundary additive
    !* integer ( kind = 4 ) iarg,iargc.ierror,ilen,ipxfargc
    !* integer ( kind = 4 ) lenc,numarg,npow
    !real ( kind = 8 ) xmax1,xmin1,ymax1,ymin1,zmax1,zmin1,POS(3)   !POS(3)for x,y,z tmp
    
    !* important stored data
    real*8, allocatable :: valiso(:,:,:) ;   !Save tje density
    !atoms with its coordinate and other parameters,ATOMS(4) used to save atom index and used in gbox atom
    real*8, allocatable :: ATOMS(:,:) ;   
    !real*8, allocatable :: atom_x(:),atom_y(:),atom_z(:)!To save all the x,y,z coordinates
    real*8 :: mincoor(3),maxcoor(3);    !min and max of XYZ with boundary dismag
    real*8 :: rad(atom_max);      !possible use for radius of atom from pqr

    !* tmp variables
    integer ( kind = 4) :: na; !Num of atoms, use in "for"
    real ( kind = 8) :: dis;    !distance between two atoms
    real ( kind = 8) :: xx,yy,zz !tmp data for x/y/z and xyz
    integer (kind=4):: i,j,k,ii,i1,j1,k1,jj,kk,ig,jg,kg,inum,km; !tmp data control "for"
    real (kind = 8 ) :: t1,t2,t12;   !time cost
    character*80 ::datestr,tmpstr;     !special time
    real*8 :: readArg;  !tmp argv

    !* parameters for grids
    real ( kind = 8) :: dcel ;  !grid size
    real ( kind = 8) :: dismag ; !boundary additive  
    integer nx,ny,nz;  !number of grid in X/Y/Z 

    !* parameters for boxes
    real ( kind = 8) :: hbox;    !Box size
    real ( kind = 8 ) :: vdw_max ;  !max vdw radius, may affect grid number in a box
    real ( kind = 8 ) :: boxsigma; !affect box size ~ boxsigma * sigma
    integer ( kind = 4 ) :: nboxgrid ; !each box containing grid number.
    integer nx1,ny1,nz1;  !number of box in X/Y/Z , maybe used instead by nbox
    integer ( kind = 4 ) :: nbox(3);   !Box number in XYZ
    integer ( kind = 4 ) :: idx (3); !Box index for XYZ, save atom and operate box idx
    integer, allocatable :: natom(:,:,:);    !atom numbers in XYZ box
    integer, allocatable :: m(:,:,:);    !use to save acumulated atom num in a certain box when assign atoms into box
    integer, allocatable :: ijkdx(:,:,:,:); ! grid range (xmin-xmax..) in a certrain box
    integer*4 :: kmin,kmax,jmin,jmax,imin,imax   !index of near boxes of certain box
    integer*4 :: kgmin,kgmax,jgmin,jgmax,igmin,igmax   !near grid range
    real*8 :: atomx,atomy,atomz;   !tmp to save xyz
    !using the structure box
    type box
        real*8,pointer:: atom(:,:);!first box idx, second is xyz and atom index
    end type
    type(box),pointer:: gbox(:,:,:);

    !* parameters for gaussian
    real ( kind = 8 ) :: sigma  ; ! sigma in e(-(dis/sigma)^nexp)
    real ( kind = 8 ) :: weight,wii ;   !weight for adjust the size, wii indeed weight
    integer ( kind = 4 ) :: nexp ; !num for exponent e^(()^nexp)
    integer ( kind = 4 ) :: elemNum ; !num of electron in atom
    real ( kind = 8 ) :: wexp ;   !weight for adjust the value in exp

    !* parameters for adjusting output data
    real ( kind = 8 ) :: scale ; !To scale the output data
    real ( kind = 8 ) :: val_max ;  !Save the max output data
    real ( kind = 8 ) :: scale_iso ;   !min value for output (else=1)
    real ( kind = 8 ) :: stop_iso ; !max value for output (else max)
    integer ( kind = 4 ) :: iso ;   !Save tmp ISO value (integer for Perseus)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !   Parameters setup
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    writelog=.true. !.true.
    writedxfile=.true. !control output dx file
    scubtop=.false. !control output either cubtop file or scubtop file
    celllisting=.true. ;!!control for cell_listing
    
    useAtomNum=.false.  !control by input args
    vdw_max=2.275;      !!!!!! the maximal ver dan Waals raidus of protein atoms
    boxsigma=3.0*vdw_max/sqrt(2.0);
    !boxsigma=3.0;

    dcel=1d0;          !!!!!!!!!!!!!!! at least 0.08
    sigma=1d0;       !!atom radius related
    dismag=3d0;      !boundary
    nexp=1;       !(dis/sigma)^nexp, 
    wexp=1d0;       !exp(wexp *-(dis/sigma)^nexp)
    weight=1d0;     ! the w value for weight function.
    scale=1.0e3;    ! scale for the origin data to fit persuse input
    scale_iso=500.0;    !At least data to use in dat file.
    stop_iso=1100.0;    !Max data to be used in dat file.

    !* For input arguments
    !* 1st: file; 2nd: grid size; 3rd: sigma; 4th: scale
    !* 5th: nexp; 6th:wexp; 7th: use atom_num(>0), 8th: boxsigma.
    print *,'| ************** Arguments input *************** |'    
    print *,'| 1:filename, 2:grid_size, 3:sigma, 4:scale      |'
    print *,'| 5:nexp, 6:wexp, 7:use atom_num(>0),8:boxsigma  |'
    print *,'| ********************************************** |'
    argc=iargc();
    CALL GETARG(1, fileinput)
    if (argc < 1) then
        print *,'Error, Please assign the filename!'
        print *,'1:filename, 2:grid_size, 3:sigma, 4:scale'
        print *,'5:nexp, 6:wexp, 7:use atom_num(>0),8:boxsigma'
        stop
    end if
    CALL GETARG(2, argstr)
    if (argstr /= '') then 
        read(argstr, *) readArg
        dcel=readArg
    endif
    CALL GETARG(3, argstr)
    if (argstr /= '') then 
        read(argstr, *) sigma
    endif
    CALL GETARG(4, argstr)
    if (argstr /= '') then 
        read(argstr, *) scale
    endif
    CALL GETARG(5, argstr)
    if (argstr /= '') then 
        read(argstr, *) nexp
    endif
    CALL GETARG(6, argstr)
    if (argstr /= '') then 
        read(argstr, *) wexp
    endif
    
    CALL GETARG(7, argstr)
    if (argstr /= '') then 
        read(argstr, *) useAtomNumInt;
        if (useAtomNumInt > 0) then 
            useAtomNum=.true.
            weight=weight/sqrt(6.2832*sigma)!to use as Structure formula.
        endif
    endif

    CALL GETARG(8, argstr)
    if (argstr /= '') then 
        read(argstr, *) boxsigma;
    endif


    !!!nboxgrid: grid number in a box;  each box have 3*(sigma*vdw/sqrt(2))/dcel grids.
    nboxgrid=ceiling(boxsigma * sigma /dcel); 
    !!hbox is the length of each box(total box is larger than original box).
    hbox=dcel*real(nboxgrid,kind=8);

    !fileinput='test.mol2'
    call filenameParse(fileinput,fname,extname)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !            Read in protein information
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (extname(1:3) == 'pdb' ) then !large case PDB?
        filepdb=trim(fname)//'.pdb';
        input_unit=345;
        open ( unit = input_unit, file = filepdb, status = 'old', iostat = ios )
        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
            write ( *, '(a)' ) 'Could not open the pdb file.'
            stop
        end if
        
        !*  Read the ATOM information from the file.
        call pdb_read_atom ( coord, hiatom, hires, input_unit, atom_max, maxres, &
            mxpatm, atom_num, numchain, numline, numlost, numres, prtatm, &
            resnam, resnum, xmax, xmin, ymax, ymin, zmax, zmin, tfactor, elements, atom_names )
        close ( unit = input_unit )

        !*  Print a summary of the information.
        call pdb_summary ( hiatom, hires, atom_max, maxres, mxpatm, &
            atom_num, numchain, numline, numlost, numres, xmax, xmin, &
            ymax, ymin, zmax, zmin )

        !*  Reset any data that is out of bounds.
        call pdb_check ( hires, hiatom, atom_max, maxres, atom_num, numres )
        !call pdb_print_coord ( coord, atom_max, atom_num, resnam, resnum,elements, atom_names)  
        !*  Print out the PRTATM array.

    else if (extname(1:4)=='mol2' ) then 
        filepdb=trim(fname)//'.mol2';
        input_unit=346;
        open ( unit = input_unit, file = filepdb, status = 'old', iostat = ios )
        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
            write ( *, '(a)' ) 'Could not open the mol2 file.'
            stop
        end if   
        !* read mol2 file
        call read_one_mol2 (input_unit, atom_max,atom_num,coord,xmax, xmin, ymax, ymin, zmax, zmin,elements)
        close (unit = input_unit )

    else if (extname(1:4)=='xyzr' ) then 
        filepdb=trim(fname)//'.xyzr';
        input_unit=347;
        open ( unit = input_unit, file = filepdb, status = 'old', iostat = ios )
        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
            write ( *, '(a)' ) 'Could not open the xyzr file.'
            stop
        end if   
        !* read xyzr file
        call read_xyzr (input_unit, atom_max,atom_num,coord,xmax, xmin, ymax, ymin, zmax, zmin,rad)
        close (unit = input_unit )

    else 
        write(*,*) 'Error to read the input file:'//trim(fileinput)
        stop
    end if 

    !* modified ATOMS(dim-5) for element number
    allocate(ATOMS(atom_num,5));
    do na=1,atom_num
        ATOMS(NA,1)=coord(na,1);
        ATOMS(NA,2)=coord(na,2);
        ATOMS(NA,3)=coord(na,3);
        ATOMS(NA,4)=NA; !Save the atom index
        if (extname(1:4)=='mol2' .or. extname(1:3)=='pdb') then
            call element2atomNum(elements(NA),inum)
            ATOMS(NA,5)=dble(inum);
        else 
            ATOMS(NA,5) = 1d0;
        endif 
        !write(*,*) elements(na),int(atoms(na,5));
    end do

    !!!just to write out a file for Xia's Program
    !if (.false.) then
    !    open(unit=7676, file='3J6C-A24-2.xyzr')
    !    do na = 1 , atom_num
    !        write(7676,*) ATOMS(na,1),ATOMS(na,2),ATOMS(NA,3),1d0
    !    enddo
    !end if

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !   FRI surface construction
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cpu_time(t1);
    !process output name,dcel*100_sigma*100
    write(tmpstr,'(I4)') ceiling(dcel*100);
    fname=trim(fname)//'_'//trim(adjustl(tmpstr));
    write(tmpstr,'(I4)') ceiling(sigma*100);
    fname=trim(fname)//'_'//trim(adjustL(tmpstr));
    
    !open(12,file='1vii_NH_500.dat')
    open(12,file=trim(fname)//'.dat')
    call GetDateTimeString(datestr);
    
    if (writelog) then
        open(1234,file=trim(fname)//'_'//trim(datestr)//'.log')
        write(1234,*) 'Work started at time: '//datestr;
    endif

    xleft =xmin-dismag;!!min value, basis
    yleft =ymin-dismag;
    zleft =zmin-dismag;
    xright=xmax+dismag;!!max value, change for box!
    yright=ymax+dismag;
    zright=zmax+dismag;

    !!needed origin box number
    nx1=ceiling((xright-xleft)/hbox);
    ny1=ceiling((yright-yleft)/hbox);
    nz1=ceiling((zright-zleft)/hbox);
    nbox=(/nx1,ny1,nz1/);

    if (celllisting) then
        !indeed grid number when using cell_list.
        !nx=nx1*nboxgrid+1; !xia
        !ny=ny1*nboxgrid+1; !xia
        !nz=nz1*nboxgrid+1; !xia
        nx=nx1*nboxgrid;
        ny=ny1*nboxgrid;
        nz=nz1*nboxgrid;
        xright=xleft+real(nx1)*hbox;
        yright=yleft+real(ny1)*hbox;
        zright=zleft+real(nz1)*hbox;
    else 
        !!needed grid, using in non-cell_list case.
        nx=int((xright-xleft)/dcel)+1;
        ny=int((yright-yleft)/dcel)+1;
        nz=int((zright-zleft)/dcel)+1;
        xright=xleft+real(nx-1)*dcel;
        yright=yleft+real(ny-1)*dcel;
        zright=zleft+real(nz-1)*dcel;
    end if

    !!obtain min/max value in x,y,z with 
    !MINCOOR=MINVAL(ATOMS,DIM=1)-dismag 
    !MAXCOOR=MAXVAL(ATOMS,DIM=1)+dismag
    MINCOOR=(/xleft, yleft, zleft/);
    MAXCOOR=[xright, yright, zright];

    allocate(valiso(nx,ny,nz));

    write(*,'(a,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)') " * xmax/min, ymax/min, zmax/min",&
                                                xmax, xmin, ymax, ymin, zmax, zmin
    write(*,'(a,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)') " * x/y/zright,x/y/zleft",&
                                    xright,yright,zright,xleft,yleft,zleft
    write(*,'(a,I6,I6,I6,I6,I6,I6)') " * ngrid:x,y,z; nbox:x,y,z",nx,ny,nz,nx1,ny1,nz1
    write(*,'(a,F8.3,F6.2,I2,F6.3,F8.3,F12.3)') " * dcel,sigma,nexp,wexp,wi,scale",dcel,sigma,nexp,wexp,weight,scale
    if (useAtomNum) write(*,'(a)') " ** Using atom number to scale."
    write(*,'(a,F8.3,I6,F8.3,I6,I6,I6)') " * hbox,nboxgrid,boxsigma",hbox,nboxgrid,boxsigma
    if (writelog) then
        write(1234,*) " xmax/min, ymax/min, zmax/min",&
                        xmax, xmin, ymax, ymin, zmax, zmin
        write(1234,*) " x/y/zright,x/y/zleft",&
                        xright,yright,zright,xleft,yleft,zleft
        write(1234,*) " ngrid:x,y,z; nbox:x,y,z",nx,ny,nz,nx1,ny1,nz1
        write(1234,*) "dcel,sigma,nexp,wexp,wi,scale",dcel,sigma,nexp,wexp,weight,scale
        if (useAtomNum) write(1234,*) "**Using atom number to scale."
        write(1234,*) "hbox,nboxgrid,boxsigma",hbox,nboxgrid,boxsigma
    endif

    !! initiate the iso value.
    valiso=0.0d0;

    if (celllisting) then   !!!cell_listing method.
        !!natom: atoms num in each box
        !!ijkdx: near 6 boxes 
        allocate(natom(nbox(1),nbox(2),nbox(3)),ijkdx(nbox(1),nbox(2),nbox(3),6))
        allocate(m(nbox(1),nbox(2),nbox(3)))
        allocate(gbox(nbox(1),nbox(2),nbox(3)))

        !* Calculate atom number in box
        natom=0
        do i=1,atom_num       
           idx=ceiling((ATOMS(i,1:3)-mincoor)/hbox) !location of the box which the ith atom belongs to 
           natom(idx(1),idx(2),idx(3))=natom(idx(1),idx(2),idx(3))+1 !compute the number of atoms in each box
        end do

        !!!Save atoms
        do k=1,nbox(3)
        do j=1,nbox(2)
        do i=1,nbox(1)
            if (natom(i,j,k)/=0) then  
                !allocate the space in each box to save the atoms 
                allocate(gbox(i,j,k)%atom(1:natom(i,j,k),1:5))      
            end if
            !!!!!!!! grid index in each box
           ijkdx(i,j,k,1)=nboxgrid*(i-1)+1;  
           ijkdx(i,j,k,2)=nboxgrid*(j-1)+1;
           ijkdx(i,j,k,3)=nboxgrid*(k-1)+1;
           ijkdx(i,j,k,4)=nboxgrid*i;         
           ijkdx(i,j,k,5)=nboxgrid*j;
           ijkdx(i,j,k,6)=nboxgrid*k;
           !if(i==nbox(1))  ijkdx(i,j,k,4)=nboxgrid*i+1;
           !if(j==nbox(2))  ijkdx(i,j,k,5)=nboxgrid*j+1;
           !if(k==nbox(3))  ijkdx(i,j,k,6)=nboxgrid*k+1;
        end do
        end do
        end do

        !* save atom into box and re-index atom in gbox%atom
        M=1
        do i=1,atom_num
            idx=ceiling((ATOMS(i,1:3)-mincoor)/hbox)
            !send the ith atom to gbox(idx(1),idx(2),idx(3))
            gbox(idx(1),idx(2),idx(3))%atom(m(idx(1),idx(2),idx(3)),1:5)=ATOMS(i,1:5);
            m(idx(1),idx(2),idx(3))=m(idx(1),idx(2),idx(3))+1;
        end do

        !* Go into a box and know its neigbor
        do k=1,nbox(3)
            if (k>=2 .and. k<=(nbox(3)-1)) then
                kmin=k-1
                kmax=k+1
            else if (k==1) then
                if (nbox(3)>1) then
                    kmin=1
                    kmax=k+1
                end if
                if (nbox(3)==1) then
                    kmin=1
                    kmax=1
                end if
            else
                kmax=k
                kmin=k-1
            end if
        
        do j=1,nbox(2)  
            if(j>=2 .and. j<=(nbox(2)-1)) then
                jmin=j-1
                jmax=j+1
            else if(j==1) then
                if (nbox(2)>1) then
                    jmin=1
                    jmax=j+1
                end if
                if (nbox(2)==1) then
                    jmin=1
                    jmax=1
                end if
            else
                jmax=j
                jmin=j-1
            end if
        
        do i=1,nbox(1)       
            if (i>=2 .and. i<=(nbox(1)-1)) then
                imin=i-1
                imax=i+1
            else if (i==1) then
                if (nbox(1)>1) then
                    imin=1
                    imax=i+1
                end if
                if (nbox(1)==1) then
                    imin=1
                    imax=1
                end if
            else
                imax=i
                imin=i-1
            endif
        !---------------------------
        !now you are in the  gbox(i,j,k)
        !start to loop all the neighborhood using the following loops

            !!grid ranges
            igmin=ijkdx(i,j,k,1);   
            jgmin=ijkdx(i,j,k,2);
            kgmin=ijkdx(i,j,k,3);
            igmax=ijkdx(i,j,k,4);
            jgmax=ijkdx(i,j,k,5);
            kgmax=ijkdx(i,j,k,6);

            !* loop all the grid,xx,yy,zz is grix coord.
            do kg=kgmin,kgmax
                zz=zleft+real(kg-1)*dcel;
            do jg=jgmin,jgmax
                    yy=yleft+real(jg-1)*dcel;
            do ig=igmin,igmax
                        xx=xleft+real(ig-1)*dcel;
        !do ii=1,natom(i,j,k) !loop atoms in the box, Christ
               
                do k1=kmin,kmax
                do j1=jmin,jmax
                do i1=imin,imax
!---------------------------------!
!!now you are in the neighborhood of gbox(i,j,k), gbox(i1,j1,k1)        
!---------------------------------
                    do jj=1,natom(i1,j1,k1)
                        atomx=gbox(i1,j1,k1)%atom(jj,1);
                        atomy=gbox(i1,j1,k1)%atom(jj,2);
                        atomz=gbox(i1,j1,k1)%atom(jj,3);
                        inum=int(gbox(i1,j1,k1)%atom(jj,4));

                        if (useAtomNum) then
                        !according to element
                            wii = weight * gbox(i1,j1,k1)%atom(jj,5);
                        else
                            wii = weight ;
                        endif

                        dis=sqrt((atomx-xx)**2+(atomy-yy)**2+(atomz-zz)**2);
                        !!Calculate the iso value at this grid. following Xia use atom radius.
                    !   valiso(ig,jg,kg) = valiso(ig,jg,kg)+ exp(-(dis)**2/(sigma*rad(inum))**2);  
                    !   valiso(ig,jg,kg) = valiso(ig,jg,kg)+ exp(-dis/(sigma*rad(inum)));  
                    !   valiso(ig,jg,kg) = valiso(ig,jg,kg)+ exp(-(dis**2-rad(inum)**2));  
                        valiso(ig,jg,kg) = valiso(ig,jg,kg) + wii * exp((-(dis/sigma)**nexp) * wexp);
                    enddo 
                enddo 
                enddo
                enddo
            enddo          
            enddo
            enddo
        enddo              
        enddo
        enddo

    else  !!!normal method without cell_listing
        if ( .true. ) then !!all atoms with same sigma/w
            do kk=1,nz
                zz=zleft+real(kk-1)*dcel;
                do jj=1,ny
                    yy=yleft+real(jj-1)*dcel;
                    do ii=1,nx
                        xx=xleft+real(ii-1)*dcel;
                        do inum=1,atom_num
                            dis=sqrt((ATOMS(inum,1)-xx)**2+(ATOMS(inum,2)-yy)**2+(ATOMS(inum,3)-zz)**2);
                                    !p=sum(wi*exp(-((|r-ri|)/sigma)^nexp))
                            valiso(ii,jj,kk) = valiso(ii,jj,kk) + wii * exp(-(dis/sigma)**nexp);
                        end do
                    end do
                end do
            end do
        else        !!special sigma on some atoms.
            do kk=1,nz
                zz=zleft+real(kk-1)*dcel;
                do jj=1,ny
                    yy=yleft+real(jj-1)*dcel;
                    do ii=1,nx
                        xx=xleft+real(ii-1)*dcel;
                        do inum=1,atom_num
                            dis=sqrt((ATOMS(inum,1)-xx)**2+(ATOMS(inum,2)-yy)**2+(ATOMS(inum,3)-zz)**2);
                                    !p=sum(wi*exp(-((|r-ri|)/sigma)^nexp))
                            if (elements(inum) ==' O' .or. elements(inum) ==' N') then
                                valiso(ii,jj,kk) = valiso(ii,jj,kk)+ wii * exp(-(dis/sigma)**nexp);
                            else
                                valiso(ii,jj,kk) = valiso(ii,jj,kk)+ wii * exp(-(dis/(sigma*0.5))**nexp);
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end if

    val_max=maxval(valiso);
    write(*,*)"minimal isovalue",minval(valiso);
    write(*,*)"maximal isovalue",val_max;
    if (writelog) then
        write(1234,*)"minimal isovalue",minval(valiso);
        write(1234,*)"maximal isovalue",val_max;   
    endif
    !pause
    
    do kk=1,nz
        do jj=1,ny
            do ii=1,nx
                km=km+1
                ! scale and take the complement to the max value.
                valiso(ii,jj,kk)= scale*(val_max-valiso(ii,jj,kk));
                ! valiso(ii,jj,kk)= scale*valiso(ii,jj,kk);
            enddo
        enddo
    enddo

    call cpu_time(t2);
    write(*,*)"time elipse", t1,t2,t2-t1;
    if (writelog) write(1234,*)"time elipse", t1,t2,t2-t1;

    if (writedxfile) then 
        OPEN(11,file=trim(fname)//'.dx')
        write(11,"('# isovalue of FRI surface, ',$)")
        write(11,"('#')")
        write(11,"('object 1 class gridpositions counts 'I5, I5, I5)"), nx, ny, nz
        write(11,"('origin ' f9.3, x, f9.3,x, f9.3)"), xleft, yleft, zleft
        write(11,"('delta  ' 3f8.3)"), dcel, 0.0, 0.0
        write(11,"('delta  ' 3f8.3)"), 0.0, dcel, 0.0
        write(11,"('delta  ' 3f8.3)"), 0.0, 0.0, dcel
        write(11,"('object 2 class gridpositions counts 'I5, I5, I5)"), nx, ny, nz
        write(11,"('object 3 class array type double rank 0 items ' I8 '  data follows')"), nx*ny*nz;
        km=0
        do ii=1,nx
            do jj=1,ny
                do kk=1,nz
                    km=km+1
                    write(11,*)valiso(ii,jj,kk);
                enddo
            enddo
        enddo
        CLOSE(11)
    write(*,*)"finish dx file"
    end if

    !control using scubtop at Control variable
    !scubtop=.true.
    write(12,'(I6)'),3;
    if (.not. scubtop) then
        write(12,'(I3)'),NX;
        write(12,'(I3)'),NY;
        write(12,'(I3)'),NZ;
    end if

    km=0
    do kk=1,nz
        do jj=1,ny
            do ii=1,nx
                km=km+1
                if ( scubtop ) then
                    if(valiso(ii,jj,kk).le.scale_iso) then
                        write(12,'(I6,I6,I6,I5)')ii,jj,kk,1
                    elseif(valiso(ii,jj,kk)>scale_iso.and.valiso(ii,jj,kk)<stop_iso) then
                        iso=ceiling(valiso(ii,jj,kk)-scale_iso);
                        write(12,'(I6,I6,I6,I5)')ii,jj,kk,iso;
                    end if
                else
                    write(12,'(I5)') ceiling(valiso(ii,jj,kk));
                end if
            enddo
        enddo
    enddo
    CLOSE(12)
    write(*,*)"finish cubic file";

    if (writelog) CLOSE(1234)
end



!!Get a string for the date and times.
subroutine GetDateTimeString( outstr )
    integer date_time(8)
    character*10 b(3),str;
    character*80 outstr

    call date_and_time(b(1), b(2), b(3), date_time)
    !print *,'date_time    array values:'
    !print *,'year=',date_time(1)
    !print *,'month_of_year=',date_time(2)
    !print *,'day_of_month=',date_time(3)
    !print *,'time difference in minutes=',date_time(4)
    !print *,'hour of day=',date_time(5)
    !print *,'minutes of hour=',date_time(6)
    !print *,'seconds of minute=',date_time(7)
    !print *,'milliseconds of second=',date_time(8)
    !print *, 'DATE=',b(1)
    !print *, 'TIME=',b(2)
    !print *, 'ZONE=',b(3)
    outstr=trim(b(1)(1:6))//'_'//trim(b(2)(1:6))
!     write(str,'(I2)') date_time(5);
!     outstr=trim(outstr)//str
!     write(str,'(I2)') date_time(6);
!     outstr=trim(outstr)//str
!     write(str,'(I2)') date_time(7);
!     outstr=trim(outstr)//trim(str)
end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  normal distribution with mean U and standard deviation G, N is the length of the random number
  !   A(N), R is the intial value for generating random number.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE NGRNS(U,G,R,N,A)
    real*8 R,S,W,V,T
    real*8 u,g
    real*8 A(N)
    S=65536.0
    W=2053.0
    V=13849.0
    DO 20 J=1,N
      T=0.0
      DO 10 I=1,12
        R=W*R+V
        M=R/S
        R=R-M*S
        T=T+R/S
10    CONTINUE
      A(J)=U+G*(T-6.0)
20  CONTINUE
    RETURN
    END

subroutine filenameParse (filename,basename,extname)
    character*80 filename,basename,extname
    integer*4 dotPlace
    dotPlace=index(filename,'.',back = .true.);
    basename=filename(:dotPlace-1);
    extname=filename(dotPlace+1:)
end subroutine