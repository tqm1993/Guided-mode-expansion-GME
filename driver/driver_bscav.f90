
!This program give the band structure of H1, H2 and L3 cavity, in the supercell approach
program driver_bscav
	
	use phys_cte
	use geometry_defs
	use geometry_func
	use mode_solver
	use scan_process
	
	
	implicit none
	integer ::reflflag,lossflag,nbsympts,nbdiv,rang(2),maxguide,mdim(2)
	character(len=120) ::output,optcell
	real(dp) :: inpls(8),Grat
	real(dp),allocatable ::ksympts(:,:)
	real(dp),allocatable :: outkls(:,:),outeven(:,:),outodd(:,:)
	real(dp),allocatable :: outlossod(:,:),outQod(:,:),outflagod(:,:)
	real(dp),allocatable :: outlossev(:,:),outQev(:,:),outflagev(:,:)
	
	!local var for symmetric pts (atm for hexagonal latt)
	real(dp) ::glatt(2,2),g1(2),g2(2),per,phi
	real(dp) ::Gpts(2),MMpts(2),Mpts(2),Kpts(2),g12(2),xprop
	real(dp) :: Mptssq(2),Xptssq(2)
	
	!**********************************************
	!***********************************************
	!DECLARATION OF INPUT PARAMETERS
	
	!output filename
	output="h1test.dat"

	!scan option (reflflag for odd and even modes, lossflag for loss calc option)
	reflflag=1 !Calculate for even (1) or odd mode (-1)
	lossflag=0 !Calculate intrinsic loss (1) or skip it (0)
	!********
	
	!type of unit cell
	optcell="h1_cavity"
	inpls=(/1.d0,0.3d0,pi/3.d0,1.d0,12.11d0,1.d0,0.d0,0.d0/)
	!For full definition of inpls in each case see geometry_func module
	!h1_cavity inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
	!h2_cavity inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
	!l3_cavity inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
	
	!High symmetric pts of supercell(defined Mpts,Gpts,Kpts here)
	acellnb=8 !number of cell used to construct supercell (Default value is 8, defined in the geometry_func module)
	per=inpls(1)*acellnb*1.d0 !supercell period
	phi=inpls(3) !phi angle
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

	!Cutoff of the wavevector basis =Grat*|g1|
	Grat=12.5d0
	
	!max order of guided mode and the number of high symmetric points used to define the k-path of band structure.
	maxguide=2
	nbsympts=2
	allocate(ksympts(2,nbsympts))
	ksympts(:,1)=Gpts
	ksympts(:,2)=Mpts

	
	!number of division in each interval of k-points
	nbdiv=10
	
	!number of energy band to calculate (should be high for the calculation of cavity with the band folding effect resulted from the use of supercell)
	rang=(/2,900/)
	!map resolution
	mdim=(/512,512/)
	
	
	!END OF DECLARATION OF INPUT PARAMETERS
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
	if( optcell .eq. "h1_cavity") then
		print *,"Calculation for H1 cavity"
		print *,""
	else if( optcell .eq. "h2_cavity") then
		print *,"Calculation for H2 cavity"
		print *,""
	else if( optcell .eq. "l3_cavity") then
		print *,"Calculation for L3 cavity"
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
	

end program driver_bscav