program driver_dos
	
	use phys_cte
	use geometry_defs
	use geometry_func
	use mode_solver
	use scan_process
	
	
	implicit none
	integer ::reflflag,lossflag,nbdiv(2),rang(2),maxguide,mdim(2)
	character(len=120) ::inputfile,optcell
	real(dp) :: inpls(6),Grat
	real(dp),allocatable :: outkls(:,:),outeven(:,:),outodd(:,:)
	real(dp),allocatable :: outlossod(:,:),outQod(:,:),outflagod(:,:)
	real(dp),allocatable :: outlossev(:,:),outQev(:,:),outflagev(:,:)
	
	!local var for symmetric pts (atm for hexagonal latt)
	real(dp) ::glatt(2,2),g1(2),g2(2),per,phi
	real(dp) ::Gpts(2),MMpts(2),Mpts(2),Kpts(2),g12(2),xprop
	real(dp) :: Mptssq(2),Xptssq(2),Xpts1D(2)
	
	!for DOS
	character(len=120) ::fmtin,outdos
	real(dp),allocatable ::Mband(:,:),limfreq(:)
	integer ::nbk,nbband
	real(dp) ::broad,dosin,dosout,dosfull
	real(dp) ::Efirst,Elast,Ener
	integer ::ii,step
	logical ::exist
	
	
	!**********************************************
	!***********************************************
	!DECLARATION OF INPUT PARAMETERS
	
	!inputfile to load band structure
	inputfile="output/dosband"
	!output file for DOS calculation results
	outdos="output/DOSfile"
	!data format of the input file
	fmtin='(10000f20.6)'
	
	!mode symmetry: even (1) or odd (-1)
	reflflag=1
	!flag for loss calculation: Proceed (1) or skip (0)
	lossflag=0
	!********
	
	!type of unit cell
	optcell="general"
	inpls=(/1.d0,0.3d0,pi/3.d0,1.d0,12.11d0,0.5d0/)
	!general:(period,radius,angle,eps0,eps1,d) General lattice, with 1 circular hole per unit cell
	!hybrid: inputls(period,rhex,rgraph,eps0,eps1,d) Hybrid hexagonal-graphite lattice, with 2 circular hole per unit cell. 60 deg angle between basis vector
	!kagome: (period,a,rkag,eps0,eps1,d) Kagome lattice, with 3  circular holes per unit cell. 60 deg angle between basis vector
	!trianghole: (period,L,angle,eps0,eps1,d)  General lattice with 1 equilateral triangular hole per unit cell

	!period and phi angle between basis vector
	per=inpls(1)
	phi=inpls(3)

	call glatt_init(per,phi,glatt)

	g1=glatt(:,1)
	g2=glatt(:,2)
	
	!number of division in each basis vector, following Monkhorst-Pack grid definition
	nbdiv=(/60,60/)
	nbk=nbdiv(1)*nbdiv(2)
	
	!Cutoff of the wavevector basis =Grat*|g1|
	Grat=5.5d0
	!Max order of guided modes which form the basis 
	maxguide=4
	
	!Maximum number of band taken account in the calculation 
	nbband=28
	rang=(/1,nbband/)
	
	!broading of the gaussian function, in normalized frequency (a/lambda) unit
	broad=0.003d0 
	
	!map resolution
	mdim=(/128,128/)

	

	!Scan range for photon energy (in normalized frequency (a/lambda) unit)
	Efirst=0.0d0
	Elast=0.7d0
	step=210
	
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
	end if

	print *,"Band range: ",rang
	print *,""

	print *,"Map resolution: ",mdim
	print *,""
	


	call bzscan_band (nbdiv,optcell,rang,inpls,&
										&Grat,mdim,reflflag,maxguide,lossflag,&
										&outkls,outeven,outodd,output,&
										&outlossod,outQod,outflagod,outlossev,outQev,outflagev)

	call load_phoband(inputfile,nbband,nbk,fmtin,Mband,limfreq)
	

	inquire(file=outdos,exist=exist)
	if (exist) then
		open(1, file=outdos, status="replace", action="write",form="formatted")
	else
		open(1, file=outdos, status="new", action="write",form="formatted")
	end if
	
	do ii=0,step
		Ener=Efirst+(Elast-Efirst)*ii/(1.d0*step)
		call dos_pho_2d(Ener,nbk,nbband,Mband,limfreq,broad,g1,g2,&
										&dosin,dosout,dosfull)
										
		write (1,'(4f20.6)')Ener,dosin,dosout,dosfull
	end do
	
	close(1)

end program driver_dos






