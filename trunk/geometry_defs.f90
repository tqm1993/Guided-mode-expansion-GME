
! This module define the geometry of the photonic crystal membrane
! PCstack contains full info of the optical stack configuration (total number of layers, list of layer average index, list of layer thickness

module geometry_defs
	
	
	use phys_cte
	include 'fftw3.f'

	!type PCstack:
	!nbtotal: Total number of layers
	!navg: List of average refractive index
	!eps: List of average dielectric constant
	!thickness: List of layer thickness
	
	type PCstack
	
    integer :: nbtotal
    real(dp),allocatable :: navg(:)
    real(dp),allocatable :: eps(:)
    real(dp),allocatable :: thickness(:)
    real(dp),allocatable :: tj(:)
    
    
  end type PCstack
  
	contains
		!***********************************************************************************************
		!***********************************************************************************************  
		!allocating memory for PCstack object
		
		subroutine alloc_PCstack(nblayer,stackobj)
		
			implicit none
			integer,intent(in) :: nblayer
			type(PCstack),intent(inout) :: stackobj
			
			stackobj%nbtotal=nblayer
			
			allocate(stackobj%navg(nblayer),stackobj%eps(nblayer),stackobj%thickness(nblayer),stackobj%tj(nblayer+1))
			
		end subroutine alloc_PCstack
		
		!***********************************************************************************************
		!***********************************************************************************************  
		!deallocating memory for PCstack object
		
		subroutine dealloc_PCstack(stackobj)
		
			implicit none
			type(PCstack),intent(inout) :: stackobj
			
			deallocate(stackobj%navg,stackobj%eps,stackobj%thickness,stackobj%tj)
			
		end subroutine dealloc_PCstack

		!***********************************************************************************************
		!*********************************************************************************************** 
		
		!layer_stack(nblayer, navgls, thickls,stackobj)
		!nblayer: Total number of layer
		!epsls: List of average dielectric constant
		!thickls: List of layer thickness
		!stackobj: PCstack object (output)	
		subroutine layer_stack(nblayer, epsls, thickls,stackobj)
			
			implicit none
			
			integer,intent(in) :: nblayer
			real(dp), intent (in) :: epsls(*)
			real(dp),intent (in) :: thickls(*)
			type (PCstack), intent(inout) :: stackobj
			
			!local variable
			integer :: i
			
			
			!Build the PC stack
			call alloc_PCstack(nblayer,stackobj)
			
			!upper bound has to be declared for thickls
			stackobj%thickness=thickls(1:nblayer)
			stackobj%tj(1)=-thickls(1)
			
			do i=1, nblayer
				stackobj%eps(i)=epsls(i)
				stackobj%navg(i)=sqrt(epsls(i))
				stackobj%tj(i+1)=stackobj%tj(i)+thickls(i)
			end do
	
		end subroutine layer_stack
		

		subroutine glatt_init(per,phi,glatt)
			
			implicit none
			real(dp),intent(in) :: per,phi
			real(dp),intent(out) ::glatt(2,2)
			
			real(dp) ::b0
			b0=2.d0*pi/(per*sin(phi))
			
			glatt(:,1)=(/b0*sin(phi),-b0*cos(phi)/)
			glatt(:,2)=(/0.d0,b0/)
			!print *,"glatt",glatt(:,1),glatt(:,2)
			
		end subroutine glatt_init
		
		!for 1D stripes, glatt(:,2) is set to (/0.d0,0.d0/)
		subroutine glatt_1D(per,glatt)
			real(dp),intent(in):: per
			real(dp),intent(out):: glatt(2,2)
			
			real(dp)::b0
			
			b0=2.d0*pi/per
			glatt(:,1)=(/b0,0.d0/)
			glatt(:,2)=(/0.d0,0.d0/)
		
		end subroutine glatt_1D
		
		!unit cell generation for 2D periodic photonic crystal
		!period: Pc period length
		!phi: diagonal angle of the rectangular lattice
		!center: center of the unit supercell
		!fmap: function to generate unit cell map (will be defined in a separate module, specified for each kind of geometry)
		!divxy: number of step interval for x and y direction
		!map: unit cell map (output)
		!alatt: real space primitive cell unit vectors (output)
		!asc :r real space super cell unit vectors (output)
		!glatt: primitive reciprocal lattice unit vectors (output)
		!glattxy: rectangular supercell reciprocal lattice (output) 
		subroutine unit_cell2D(period,phi,center,inpls,fmap,divxy,map,alatt,asc,glatt,glattxy)
		
			implicit none
			
			interface
			
				subroutine fmap (x,y,inputls,val)
					
					use phys_cte
					implicit none
					
					real(dp),intent(in) :: x,y
					real(dp),intent(in) :: inputls(*)
					real(dp),intent(inout) :: val
				
				end subroutine fmap
				
			end interface
			
			real(dp),intent(in) :: center(2),inpls(*)
			real(dp),intent(in) :: period, phi
			integer,intent(in) :: divxy(2)
			real(dp),allocatable, intent(out) :: map(:,:)
			real(dp),intent(out) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2)
			
			!local variable
			
			real(dp) ::b0,stepx,stepy,val
			real(dp) ::refpts(2),currpts(2)
			integer ::dx,dy
			
			
			b0=2.d0*pi/(period*sin(phi))
			
			
			!asc will be used to plot field
			asc(:,1)=(/2.d0*period*cos(phi/2.d0),0.d0/)
			asc(:,2)=(/0.d0,2.d0*period*sin(phi/2.d0)/)
			
			
			
			stepx=asc(1,1)/(1.d0*divxy(1))
			stepy=asc(2,2)/(1.d0*divxy(2))
			
			refpts=center-0.5*asc(:,1)-0.5*asc(:,2)
			
			!TODO: check val variable declaration

			allocate(map(divxy(2)+1,divxy(1)+1))

			
			do dx=1,(divxy(1)+1)
				do dy=1,(divxy(2)+1)
				
					currpts=refpts+(/stepx*(dx-1),stepy*(dy-1)/)
					call fmap(currpts(1),currpts(2),inpls,val)
					map(dy,dx)=val
					
				end do
			end do
			
			!redefining the lattice, important for scan
			alatt(:,1)=(/period,0.d0/)
			alatt(:,2)=(/period*cos(phi),period*sin(phi)/)
			
			glatt(:,1)=(/b0*sin(phi),-b0*cos(phi)/)
			glatt(:,2)=(/0.d0,b0/)
			
		end subroutine unit_cell2D

		
		!unit cell generation for 1D photonic crystal cell (Bragg mirror)
		subroutine unit_cell1D(period,center,inpls,fmap,divxy,map,alatt,asc,glatt,glattxy)
		
			implicit none
			
			interface
			
				subroutine fmap (x,y,inputls,val)
					
					use phys_cte
					implicit none
					
					real(dp),intent(in) :: x,y
					real(dp),intent(in) :: inputls(*)
					real(dp),intent(inout) :: val
				
				end subroutine fmap
				
			end interface
			
			real(dp),intent(in) :: center(2),inpls(*)
			real(dp),intent(in) :: period
			integer,intent(in) :: divxy(2)
			real(dp),allocatable, intent(out) :: map(:,:)
			real(dp),intent(out) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2)
			
			!local variable
			
			real(dp) ::b0,stepx,stepy,val
			real(dp) ::refpts(2),currpts(2)
			integer ::dx,dy
			
			
			b0=2.d0*pi/period

			asc(:,1)=(/period,0.d0/)
			asc(:,2)=(/0.d0,inpls(5)/)
			
			
			stepx=asc(1,1)/(1.d0*divxy(1))
			stepy=asc(2,2)/(1.d0*divxy(2))
			
			refpts=center-0.5*asc(:,1)-0.5*asc(:,2)
			

			allocate(map(divxy(2)+1,divxy(1)+1))

			
			do dx=1,(divxy(1)+1)
				do dy=1,(divxy(2)+1)
				
					currpts=refpts+(/stepx*(dx-1),stepy*(dy-1)/)
					call fmap(currpts(1),currpts(2),inpls,val)
					map(dy,dx)=val
					
				end do
			end do
			
			glatt(:,1)=(/b0,0.d0/)
			glatt(:,2)=(/0.d0,0.d0/)

			
		end subroutine unit_cell1D
		
		! supercell generation for 2D periodic photonic crystal
		subroutine unit_sc2D(phi,inpls,fmap,divxy,map,alatt,asc,glatt,glattxy)
		
			implicit none
			
			interface
			
				subroutine fmap (dimm,inputls,epsmap,superper)
					
					use phys_cte
					implicit none
					
					
					real(dp),intent(out) ::superper
					integer,intent(in) :: dimm(2)
					real(dp),intent(in) :: inputls(*)
					real(dp),intent(inout) :: epsmap(dimm(1),dimm(2))
				
				end subroutine fmap
				
			end interface
			
			real(dp),intent(in) :: inpls(*),phi
			integer,intent(in) :: divxy(2)
			real(dp),allocatable, intent(out) :: map(:,:)
			real(dp),intent(out) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2)
			
			!local variable
			
			integer :: dimmap(2)
			real(dp) ::scper,b0
			
			
			dimmap=(/divxy(2)+1,divxy(1)+1/)
			allocate (map(dimmap(1),dimmap(2)))
	
			call fmap(dimmap,inpls,map,scper)

			
			b0=2.d0*pi/(scper*sin(phi))
			
			
			asc(:,1)=(/2.d0*scper*cos(phi/2.d0),0.d0/)
			asc(:,2)=(/0.d0,2.d0*scper*sin(phi/2.d0)/)
			
			alatt(:,1)=(/scper,0.d0/)
			alatt(:,1)=(/scper*cos(phi),scper*sin(phi)/)
			
			glatt(:,1)=(/b0*sin(phi),-b0*cos(phi)/)
			glatt(:,2)=(/0.d0,b0/)
			
		end subroutine unit_sc2D

		
		! convert to supercell coordinate (integer version)
		subroutine iconvert_to_sc(Mtrans,coordin,coordout)
			
			implicit none
			
			real(dp),intent(in) :: Mtrans(2,2)
			real(dp),intent(in) :: coordin(2)
			integer,intent(out) ::coordout(2)
			
			!local variable
			real(dp) ::tmpcoord(2)
			
			call dgemv('N',2,2,1.d0,Mtrans,2,coordin,1,0.d0,tmpcoord,1)
			
			coordout=(/nint(tmpcoord(1)),nint(tmpcoord(2))/)
			
		end subroutine iconvert_to_sc	
		
	
		
		
		!Fourier transform of the dielectric map using fftw
		subroutine map_fft(map,mdim,fftmap)
			
			implicit none
			integer,intent(in) :: mdim(2)
			real(dp),intent (in) :: map(mdim(1),mdim(2))
			complex(dp),allocatable,intent(out) :: fftmap(:,:)
			
			
			!local variable
			integer :: plan,fdim(2),i,j
			complex(dp) ::mapc(mdim(1),mdim(2))
			
			do i=1,mdim(1)
				do j=1,mdim(2)
					mapc(i,j)=complex(map(i,j),0.d0)
				end do
			end do
			
			allocate(fftmap(mdim(1),mdim(2)))
				
			call dfftw_plan_dft_2d(plan,mdim(1),mdim(2),mapc,fftmap,FFTW_FORWARD,FFTW_ESTIMATE)
			call dfftw_execute_dft(plan,mapc,fftmap)
			call dfftw_destroy_plan(plan)
			
			!scale fftmap with 1/n (2D case)
			fftmap=(1.d0/((1.d0*mdim(1))**2.d0))*fftmap
			
			!backward fft (to debug)
			!call dfftw_plan_dft_2d(plan,mdim(1),mdim(2),fftmap,fftmap,FFTW_BACKWARD,FFTW_ESTIMATE)
			!call dfftw_execute_dft(plan,fftmap,fftmap)
			!call dfftw_destroy_plan(plan)
		
		end subroutine map_fft
		
		!***********************************************************************************************
		!***********************************************************************************************
		!extract fft complex coefficient from a defined list of coordinate 
		subroutine extract_coeff(nbpts,mn,mdim,map_fft,val)
		
			implicit none
			integer,intent(in) :: nbpts
			integer,intent(in) ::mn(2,*),mdim(2)
			complex(dp),intent(in) ::map_fft(mdim(1),mdim(2))
			complex(dp),allocatable,intent(out) :: val(:)
			
			!local variables
			integer ::redflag=0
			integer ::mtmp,ntmp,i
			
			allocate(val(nbpts))
			
			!print *,"aaa"
			do i=1,nbpts
				
				mtmp=mn(1,i)+1
				ntmp=mn(2,i)+1
				!print *,"**********"
				!print *,mn(1,i),mn(2,i)
				
				
				
				!TODO: check for negative freq and map indexes
				if ((abs(ntmp)>mdim(2)) .or. (abs(mtmp)>mdim(1))) then
					redflag=1
					val(i)=(0.d0,0.d0)
					print *, "index higher than matrix dimension.Process stopped"
				
				end if
				
				if (mn(2,i) <0) then
					ntmp=mdim(2)+mn(2,i)+1
					!print *,"neg-n"
				end if
				
				if (mn(1,i) <0) then
					mtmp=mdim(1)+mn(1,i)+1
					!print *,"neg-m"
				end if
				
				
				if (redflag==0) then
					val(i)=map_fft(mtmp,ntmp)
				end if
				
			end do
			
		end subroutine extract_coeff
			
		
		!write dielectric map to file (only real part)
		subroutine dwrite_map_tofile(mdim,map,output)
		
			implicit none
			integer,intent(in) :: mdim(2)
			real(dp),intent(in) :: map(mdim(1),mdim(2))
			character(len=120),intent(in) :: output
			
			!local variables
			logical :: exist
			integer :: i,j
			
			inquire(file=output,exist=exist)
			if (exist) then
				open(1, file=output, status="replace", action="write",form="formatted")
			else
				open(1, file=output, status="new", action="write",form="formatted")

			end if
	
			do i=1,mdim(1)
				!y col is written horizontally
				write(1, '(10000e20.8)')( map(i,j) ,j=1,mdim(2)) 
			end do
				close(1)
				
		end subroutine dwrite_map_tofile
		
		!write dielectric map to file (both the real part and complex part, in separate file)
		subroutine zwrite_map_tofile(mdim,map,output)
		
			implicit none
			integer,intent(in) :: mdim(2)
			complex(dp),intent(in) :: map(mdim(1),mdim(2))
			character(len=120),intent(in) :: output
			
			!local variables
			logical :: exist
			integer :: i,j
			character(len=120):: outputr,outputi
			
			outputr=trim(output) // "_r"
			outputi=trim(output) // "_i"
			
			inquire(file=outputr,exist=exist)
			if (exist) then
				open(1, file=outputr, status="replace", action="write",form="formatted")
			else
				open(1, file=outputr, status="new", action="write",form="formatted")

			end if
	
			do i=1,mdim(1)
				write(1, '(10000e20.8)')( real(map(i,j)) ,j=1,mdim(2)) 
			end do
				close(1)
				
			inquire(file=outputi,exist=exist)
			if (exist) then
				open(2, file=outputi, status="replace", action="write",form="formatted")
			else
				open(2, file=outputi, status="new", action="write",form="formatted")

			end if
	
			do i=1,mdim(1)
				write(2, '(10000e20.8)')( imag(map(i,j)) ,j=1,mdim(2)) 
			end do
				close(2)
				
		end subroutine zwrite_map_tofile

end module geometry_defs
		
