
!This program give the band structure of different type of photonic crystal membrane
program driver_bs
	
	use phys_cte
	use geometry_defs
	use geometry_func
	use mode_solver
	use scan_process
	
	
	implicit none
	integer ::reflflag,lossflag,nbsympts,nbdiv,rang(2),maxguide,mdim(2)
	character(len=120) ::output,optcell
	real(dp) :: inpls(6),Grat
	real(dp),allocatable ::ksympts(:,:)
	real(dp),allocatable :: outkls(:,:),outeven(:,:),outodd(:,:)
	real(dp),allocatable :: outlossod(:,:),outQod(:,:),outflagod(:,:)
	real(dp),allocatable :: outlossev(:,:),outQev(:,:),outflagev(:,:)
	
	!local var for symmetric pts (atm for hexagonal latt)
	real(dp) ::glatt(2,2),g1(2),g2(2),per,phi
	real(dp) ::Gpts(2),MMpts(2),Mpts(2),Kpts(2),g12(2),xprop
	real(dp) :: Mptssq(2),Xptssq(2),Xpts1D(2)
	
	!**********************************************
	!***********************************************
	!DECLARE INPUT PARAMETERS HERE
	
	!Output filename
	output="toto.dat" 
	
	!scan option (reflflag for odd and even modes, lossflag for loss calc option)
	reflflag=1 !Calculate for even (1) or odd mode (-1)
	lossflag=0 !Calculate intrinsic loss (1) or skip it (0)
	
	!type of unit cell
	optcell="general"
	!general:(period,radius,angle,eps0,eps1,d) General lattice, with 1 circular hole per unit cell
	!hybrid: inputls(period,rhex,rgraph,eps0,eps1,d) Hybrid hexagonal-graphite lattice, with 2 circular hole per unit cell. 60 deg angle between basis vector
	!kagome: (period,a,rkag,eps0,eps1,d) Kagome lattice, with 3  circular holes per unit cell. 60 deg angle between basis vector
	!bragg_cavity: (period,wdiel,eps0,eps1,hdummy,d) 1D Bragg cavity
	!trianghole: (period,L,angle,eps0,eps1,d)  General lattice with 1 equilateral triangular hole per unit cell
	
	!Input geometrical parameters and dielectric permittivity values
	!Depends on the type of unit cell (see geometry_func module)
	!*******
	inpls=(/1.d0,0.3d0,pi/3.d0,1.d0,12.11d0,0.5d0/) !example for general hexagonal lattice
	
	!High symmetric pts (defined Mpts,Gpts,Kpts here)
	per=inpls(1) !period (depend on inpls declared for each type of unit cell, BE CAREFUL!)
	phi=inpls(3) !phi angle (depend on inpls declared for each type of unit cell, BE CAREFUL!)
	
	if (optcell .eq. "bragg_cavity") then
		call glatt_1D(per,glatt)
	else
		call glatt_init(per,phi,glatt)
	end if
	g1=glatt(:,1)
	g2=glatt(:,2)

	!--------------------
	!sympts for hexagonal lattice
	Gpts=0.d0*g1+0.d0*g2
	Mpts=0.5d0*g2
	!calc for Kpts
	g12=g2-g1
	Kpts=(/g2(1)-(g2(2)-Mpts(2))*(g12(1)/g12(2)),Mpts(2)/)
	!--------------------

	!--------------------
	!sympts for square latt
	Mptssq=0.5*g1+0.5*g2
	Xptssq=0.5*g1
	!--------------------
	
	!--------------------
	!sympts for 1D Bragg mirror
	Xpts1D=0.5*g1
	!--------------------
	
	!Cutoff of the wavevector basis =Grat*|g1|
	Grat=5.5d0
	!max order of guided mode and the number of high symmetric points used to define the k-path of band structure.
	maxguide=4
	nbsympts=4
	
	allocate(ksympts(2,nbsympts))
	ksympts(:,1)=Gpts
	ksympts(:,2)=Mpts
	ksympts(:,3)=Kpts
	ksympts(:,4)=Gpts

	!number of division in each interval of k-points
	nbdiv=12
	
	!number of energy band to calculate
	rang=(/1,16/)
	
	!map resolution 
	mdim=(/128,128/)
	!************************************
	
	!END OF DECLARATION OF PARAMETERS
	!******************************************************
	!******************************************************
	
	
	
	!verbose part and process
	if (reflflag==0) then
		print *,"Scanning both odd and even modes"
		print *,""
	else if (reflflag==-1) then
		print *,"Scanning odd modes"
		print *,""
	else if (reflflag==1) then
		print *,"Scanning even modes"
		print *,""
	end if
	
	if (lossflag==0) then
		print *, "Loss calculation will NOT be included"
		print *,""
	else if (lossflag==1) then
		print *, "Loss calculation will be included"
		print *,""
	end if

	print *, "Input params: ", inpls
	print *,""
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
	
	else if (optcell .eq. "bragg_cavity") then
		print *,"Calculation for 1D Bragg lattice"
		print *,""
	end if

	print *,"Band range: ",rang
	print *,""

	print *,"Map resolution: ",mdim
	print *,""
	

	
	call band_structure(nbsympts,ksympts,nbdiv,optcell,rang,inpls,&
											&Grat,mdim,reflflag,maxguide,lossflag,&
											&outkls,outeven,outodd,outlossod,outQod,outflagod,&
											&outlossev,outQev,outflagev,output)
	

	
	deallocate(ksympts)	

end program driver_bs
