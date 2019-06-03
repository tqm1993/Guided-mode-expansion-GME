
!This program give the results of parameters scan on photonic crystal membrane 
program driver_scan

	use phys_cte
	use geometry_defs
	use geometry_func
	use mode_solver
	use scan_process
	
	implicit none
	integer ::nbs(2),mdim(2),ind(2),refl,maxguide,rang
	real(dp) :: kvect(2),inpls(6),Grat,scan1bds(2),scan2bds(2)
	real(dp) :: step1,step2
	real(dp),allocatable :: scan1(:),scan2(:)
	character(len=120) :: optcell,outputfile
	real(dp),allocatable ::outE(:,:),outloss(:,:)
	real(dp),allocatable ::outQ(:,:)
	real(dp),allocatable ::outflag(:,:)
	
	!local variable
	integer ::i
	
	!local var for symmetric pts (atm for hexagonal latt)
	!should be folded into subroutine 
	real(dp) ::glatt(2,2),g1(2),g2(2),per,phi
	real(dp) ::Gpts(2),Mpts(2),Kpts(2),g2ndir(2),xprop
	real(dp) ::Xptssq(2),Mptssq(2),g12(2)
	!*********************************************************
	!*********************************************************
	!DECLARATION OF INPUT PARAMETERS
	
	!output filename
	outputfile="test_scan.dat"
	
	!number of scan points and index in inpls of each parameters (up to 2 parameters)
	!BEWARE: At the moment, scans on period and the angle ohi between basis vectors are prohibited.
	nbs=(/31,0/) !number of scan points for each parameters. If nbs(2)=0 the scan will be done for only the first parameters
	ind=(/2,6/) !index in inpls of each paramers
	

	!Lower bound and upper bound of each scan parameter
	scan1bds=(/0.15d0,0.45d0/)
	scan2bds=(/0.5d0,0.9d0/)
	!************

	!type of unit cell
	optcell="general"
	inpls=(/1.d0,0.3d0,pi/3.d0,1.d0,12.11d0,0.5d0/)
	!general:(period,radius,angle,eps0,eps1,d) General lattice, with 1 circular hole per unit cell
	!hybrid: inputls(period,rhex,rgraph,eps0,eps1,d) Hybrid hexagonal-graphite lattice, with 2 circular hole per unit cell. 60 deg angle between basis vector
	!kagome: (period,a,rkag,eps0,eps1,d) Kagome lattice, with 3  circular holes per unit cell. 60 deg angle between basis vector
	!trianghole: (period,L,angle,eps0,eps1,d)  General lattice with 1 equilateral triangular hole per unit cell
	
	!extract period and phi
	per=inpls(1)
	phi=pi/3.d0
	

	!calculate the reciprocal basis vector
	call glatt_init(per,phi,glatt)
	g1=glatt(:,1)
	g2=glatt(:,2)
	
	!-----------------------
	!sympts for hexagonal lattice
	Gpts=0.d0*g1+0.d0*g2
	Mpts=0.5d0*g2
	g12=g2-g1
	Kpts=(/g2(1)-(g2(2)-Mpts(2))*(g12(1)/g12(2)),Mpts(2)/)
	!-----------------------

	!-----------------------
	!sympts for square latt
	Mptssq=0.5*g1+0.5*g2
	Xptssq=0.5*g1
	!-----------------------
	
	!the k-points on which the scan will be conducted
	kvect=0.01d0*Kpts
	!the order of energy band on which the scan will be conducted
	rang=3
	
	!Cutoff of the wavevector basis =Grat*|g1|
	Grat=5.5d0
	!*******
	
	!mode symmetry: odd (-1) or even (1)
	refl=1
	
	!map resolution
	mdim=(/128,128/)
	!max order of guided mode in the basis
	maxguide=4


	!END OF DECLARATION OF PARAMETERS
	!*******************************************************
	!*******************************************************
	
	!verbose part and process
	
	allocate(scan1(nbs(1)))
	if (nbs(2) ==0) then
		allocate (scan2(1))
		print *,"Single scan: Number of scan: ",nbs(1)," Index: ",ind(1)
		print *,""
		print *, "Scan1 from ",scan1bds(1), " to ",scan1bds(2) 
		print *,""
		scan2(1)=0.d0
		
	else
		allocate (scan2(nbs(2)))
		print *,"Dual scan: Number of scan: ",nbs," Indexes: ",ind
		print *,""
		print *, "Scan 1 from ",scan1bds(1), " to ",scan1bds(2)
		print *,""
		print *, "Scan 2 from ",scan2bds(1), " to ",scan2bds(2)
		print *,""
	end if
	
	
	step1=(scan1bds(2)-scan1bds(1))/(1.d0*(nbs(1)-1))
	do i=1,nbs(1)
		scan1(i)=scan1bds(1)+1.d0*(i-1)*step1
	end do
	
	if (nbs(2) .ne. 0) then
		step2=(scan2bds(2)-scan2bds(1))/(1.d0*(nbs(2)-1))
		do i=1,nbs(2)
			scan2(i)=scan2bds(1)+1.d0*(i-1)*step2

		end do
	end if
	

	print *, "Original input params: ", inpls
	print *,""
	print *, "k-vector coord: ",kvect
	print *,""
	if (refl==-1) then
		print *,"Scanning odd modes"
		print *,""
	else if (refl==1) then
		print *,"Scanning even modes"
		print *,""
	end if
	
	if( optcell .eq. "general") then
		print *,"Calculation for general unit cell"
		print *,""
	else if( optcell .eq. "hybrid") then
		print *,"Calculation for hybrid unit cell (hex-graphite)"
		print *,""
	else if( optcell .eq. "kagome") then
		print *,"Calculation for kagome unit cell"
		print *,""
	else if( optcell .eq. "trianghole") then
		print *,"Calculation for general cell (triangular hole)"
		print *,""
		
	end if

	print *,"Band range: ",rang
	print *,""
	
	print *,"Map resolution: ",mdim
	print *,""

	
	print *, "Starting scan ......."
	print *,""
	call loss_scan(kvect,nbs,ind,scan1,scan2,optcell,inpls,Grat,&
								&mdim,refl,maxguide,rang,&
								&outE,outloss,outQ,outflag,outputfile)

	
	deallocate (outE,outloss,outQ,outflag,scan1,scan2)
	
	print *,"Scan process is done. Data was saved"
	
end program driver_scan
