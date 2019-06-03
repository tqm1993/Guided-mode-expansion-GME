
!This module solve the hamiltonian for a defined list of k-points to generate band-structure, and contains I/O subroutine to write data into file

module scan_process
	
	use phys_cte
	use geometry_defs
	use geometry_func
	use mode_solver
	
	contains
	
	
	
	
		!Write band structure data into file (without calculation of intrinsic loss for the photonic modes)
		subroutine save_bs_noloss(nbpts,kvect,g1,g2,eigls,dimeig,rang,eigcoeff,epsls,output)
		
			implicit none
			integer,intent(in) :: nbpts,rang(2),dimeig
			real(dp),intent(in) ::kvect(2,nbpts),g1(2),g2(2),epsls(3)
			real(dp),intent(in) :: eigls(dimeig,nbpts),eigcoeff
			character(len=120),intent(in) :: output
			
			!local variable
			integer ::i,j
			real(dp) :: dummy,kdist(nbpts),tmpdist,tmpdk(2),lline(nbpts)
			logical ::exist
			
			!calculating the plot distance for all k-pts
			dummy=0.d0
			kdist(1)=dummy
			!lightline
			lline(1)=eigcoeff*sqrt((kvect(1,1)**2.d0)+(kvect(2,1)**2.d0))/sqrt(epsls(1))
			do i=2,nbpts
				!kvect already in cartesian coords
				!tmpdk=(kvect(1,i)*g1+kvect(2,i)*g2)-(kvect(1,i-1)*g1+kvect(2,i-1)*g2)
				tmpdk=kvect(:,i)-kvect(:,i-1)
				tmpdist=sqrt((tmpdk(1)**2.d0)+(tmpdk(2)**2.d0))
				dummy=dummy+tmpdist
				kdist(i)=dummy
				lline(i)=eigcoeff*sqrt((kvect(1,i)**2.d0)+(kvect(2,i)**2.d0))/sqrt(epsls(1))
			end do
			
			!check carefully
			kdist=(1.d0/dummy)*kdist
			
			inquire(file=output,exist=exist)
			if (exist) then
				open(1, file=output, status="replace", action="write",form="formatted")
			else
				open(1, file=output, status="new", action="write",form="formatted")

			end if
			
			do i=1,nbpts

				write(1, '(10000f20.6)')kdist(i),kvect(1,i),kvect(2,i),lline(i),( eigls(j,i) ,j=rang(1),rang(2))
				
			end do
			
			close(1)

		end subroutine save_bs_noloss
		
		!Write band structure data into file (with calculation of intrinsic loss for the photonic modes)
		subroutine save_bs_loss(nbpts,nbloss,kvect,g1,g2,eigls,dimeig,eigcoeff,epsls,&
												&lossrang,lossval,Qval,flagval,output)
		
			implicit none
			integer,intent(in) :: nbpts,dimeig,nbloss,lossrang(2)
			real(dp),intent(in) ::kvect(2,nbpts),g1(2),g2(2),epsls(3)
			real(dp),intent(in) :: eigls(dimeig,nbpts),eigcoeff
			real(dp),intent(in) :: lossval(nbloss,nbpts),Qval(nbloss,nbpts),flagval(nbloss,nbpts)
			character(len=120),intent(in) :: output
			
			!local variable
			integer ::i,j,marker
			real(dp) :: dummy,kdist(nbpts),tmpdist,tmpdk(2),lline(nbpts)
			logical ::exist
			
			!calculating the plot distance for all k-pts
			dummy=0.d0
			kdist(1)=dummy
			!lightline
			lline(1)=eigcoeff*sqrt((kvect(1,1)**2.d0)+(kvect(2,1)**2.d0))/sqrt(epsls(1))
			do i=2,nbpts
				!g1 and g2 are already in cartesian units
				!tmpdk=(kvect(1,i)*g1+kvect(2,i)*g2)-(kvect(1,i-1)*g1+kvect(2,i-1)*g2)
				tmpdk=kvect(:,i)-kvect(:,i-1)
				tmpdist=sqrt((tmpdk(1)**2.d0)+(tmpdk(2)**2.d0))
				dummy=dummy+tmpdist
				kdist(i)=dummy
				lline(i)=eigcoeff*sqrt((kvect(1,i)**2.d0)+(kvect(2,i)**2.d0))/sqrt(epsls(1))
			end do
			!check carefully
			kdist=(1.d0/dummy)*kdist
			
			inquire(file=output,exist=exist)
			if (exist) then
				open(1, file=output, status="replace", action="write",form="formatted")
			else
				open(1, file=output, status="new", action="write",form="formatted")

			end if
			marker=4
			write (1,'("# "i5,i5,i5,i5,i5,i5,i5)')nbloss,marker,marker+nbloss,marker+2*nbloss,&
																						&marker+3*nbloss,lossrang(1),lossrang(2)
			do i=1,nbpts
			!check carefully
				write(1, '(10000e20.6)')kdist(i),kvect(1,i),kvect(2,i),lline(i),&
														&( eigls(j,i) ,j=lossrang(1),lossrang(2)),&
														&( lossval(j,i) ,j=1,nbloss),&
														&( Qval(j,i) ,j=1,nbloss),&
														&( flagval(j,i) ,j=1,nbloss)
			end do
			
			close(1)

		
		end subroutine save_bs_loss
		
		!Save results (Mode energy, loss, q factor, mode position regard to light cone) from a parameter scan 
		subroutine save_sscan(nbs1,nbs2,scan1,scan2,outE,outloss,outQ,outflag,output)
			
			implicit none
			integer,intent(in) :: nbs1,nbs2
			character(len=120),intent(in) :: output
			real(dp),intent(in) :: scan1(*),scan2(*)
			real(dp),intent(in) :: outE(nbs1,*),outloss(nbs1,*)
			real(dp),intent(in) :: outQ(nbs1,*),outflag(nbs1,*)
			
			!local variable
			integer :: i,j,tt,tmpnbs2
			logical ::exist
			real(dp),allocatable ::tmpmat(:,:)
			character(len=120) :: tmpfile
			
			
			do tt=1,4
				
				if (nbs2==0) then
					tmpnbs2=1
					allocate(tmpmat(nbs1,1))
				
				else 
					tmpnbs2=nbs2
					allocate(tmpmat(nbs1,nbs2))
				end if
				
				!4 output files separated
				if (tt==1) then
					tmpfile=trim(output) // "_e"
					tmpmat=outE(1:nbs1,1:tmpnbs2)
				else if (tt==2) then
					tmpfile=trim(output) // "_l"
					tmpmat=outloss(1:nbs1,1:tmpnbs2)
				else if (tt==3) then
					tmpfile=trim(output) // "_q"
					tmpmat=outQ(1:nbs1,1:tmpnbs2)
				else if (tt==4) then
					tmpfile=trim(output) // "_f"
					tmpmat=outflag(1:nbs1,1:tmpnbs2)
				
				end if
				
				inquire(file=tmpfile,exist=exist)
				if (exist) then
					open(1, file=tmpfile, status="replace", action="write",form="formatted")
				else
					open(1, file=tmpfile, status="new", action="write",form="formatted")

				end if
				
				
				write(1, '(" # ",10000f20.6)')(scan1(i) ,i=1,nbs1)
				if (nbs2 .ne. 0) then
					write(1, '(" # ",10000f20.6)')(scan2(i) ,i=1,nbs2)
				end if
				
				do i =1,nbs1
					write(1, '(10000e20.6)')(tmpmat(i,j) ,j=1,tmpnbs2)
				end do
				
				close (1)
				
				deallocate(tmpmat)
			end do

		end subroutine save_sscan
		
		!solve the hamiltonaian for a list of k-pts. Loss calculation will be activated if flag==1
		subroutine full_h_solver(kvect,nbpts,kcount,maxguide,Gls,g1,g2,d,&
														&refl,epsls,etamap,coeff,eigls,lossrang,lossval,Qval,flagval,flag)
		
			implicit none
			
			integer,intent(in) :: nbpts,kcount,maxguide,refl,flag,lossrang(2)
			real(dp),intent(in) :: kvect(2,nbpts),d,epsls(3),coeff
			complex(dp),intent(in) :: etamap(kcount,kcount)
			real(dp),intent(in) :: Gls(3,kcount),g1(2),g2(2)
			!maximum case
			!real(dp),intent(out) :: eigls(kcount*maxguide,nbpts)
			!lossval contain im leak mode flag,(w/c)*coeff, Q factor
			real(dp),allocatable,intent(out) :: eigls(:,:),lossval(:,:),Qval(:,:),flagval(:,:)
			!integer,allocatable,intent(out) :: nbeig(:)
			
			!local variable
			integer :: i,j,dimH,nbloss,tt,tttt,verbose
			character(len=120) ::outmod="tmpouttotok2"
			real(dp),allocatable :: kmod(:,:)
			integer,allocatable :: Tarr(:,:),Hindex(:,:)
			complex(dp),allocatable :: ABarr(:,:),H(:,:)
			real(dp),allocatable :: earr(:,:),karr(:,:)
			real(dp) :: normloss
			
			!TODO: fix the zheev part
			!real(dp):: tmpeig(dimH),RWORK(3*dimH-2)
			real(dp),allocatable:: tmpeig(:),RWORK(:),tmploss(:),tmpeig1(:)
			complex(dp),allocatable :: WORK(:)
			integer :: LWORK,INFO
			
			allocate (eigls(kcount*maxguide,nbpts))
			
			!output loss,Q and flag mode
			if (flag==1) then
				nbloss=lossrang(2)-lossrang(1)+1
				allocate(lossval(nbloss,nbpts),Qval(nbloss,nbpts),flagval(nbloss,nbpts))
				
			end if
				
			
			!for each k-pts: Solve for gmod, info_array and extract the hamiltonian
			do i=1,nbpts
				verbose=0
				call gmodlist(kvect(:,i),Gls,kcount,g1,g2,0,outmod,kmod) 

				call info_array(kcount,maxguide,kmod,d,refl,epsls,Tarr,earr,ABarr,karr,Hindex,dimH,verbose)

				call h_guided(kcount,maxguide,dimH,kmod,d,refl,epsls,etamap,Hindex,Tarr,earr,ABarr,karr,H,0)
				!end if
				
				!solve the hamiltonian
				allocate (WORK(1),tmpeig(dimH),RWORK(3*dimH-2))
				call zheev('N','U',dimH,H,dimH,tmpeig,WORK,-1,RWORK,INFO)
				LWORK=nint(real(WORK(1)))
				deallocate(WORK)
				allocate(WORK(LWORK))
				call zheev('V','U',dimH,H,dimH,tmpeig,WORK,LWORK,RWORK,INFO)
				
				!convert to normalized frequencies (the eigenvalue is in k**2 unit)
				do j=1,dimH 
					eigls(j,i)=sqrt(tmpeig(j))*coeff
				end do
				
				
				!loss calculation
				
				if (flag==1) then
					!print *,"lossss"
					allocate (tmpeig1(dimH),tmploss(dimH))
					
					!entry eigen value
					do j=1,dimH
						tmpeig1(j)=sqrt(tmpeig(j))
					end do
					
					call loss_calc(kcount,maxguide,dimH,d,epsls,kvect(:,i),&
												&lossrang,Hindex,Tarr,kmod,earr,ABarr,karr,H,etamap,tmpeig1,tmploss,i,refl)
												
					do tt=lossrang(1),lossrang(2)
						normloss=coeff*tmploss(tt)
						lossval(tt-lossrang(1)+1,i)=normloss
						
						!check leak mode
						if (tmploss(tt) < 9998.d0) then
							Qval(tt-lossrang(1)+1,i)=coeff*tmpeig1(tt)/(2.d0*normloss)
							flagval(tt-lossrang(1)+1,i)=1.d0
						else
							Qval(tt-lossrang(1)+1,i)=0.d0
							flagval(tt-lossrang(1)+1,i)=-1.d0
						end if
					
					end do
				
				!deallocating the tmp array of loss
				deallocate (tmploss,tmpeig1)
				end if


				!deallocating local stored matrix
				deallocate(Tarr,Hindex,ABarr,H,earr,karr,kmod)
				!deallocating zheev working array
				deallocate(tmpeig,RWORK,WORK)
			
			end do

		end subroutine full_h_solver
		
		
		!calculate the band structure and loss for a defined structure and save data to ouput file
		subroutine band_structure(nbsympts,ksympts,nbdiv,optcell,rang,inpls,&
															&Grat,mdim,refl,maxguide,lossflag,&
															&outkls,outeven,outodd,outlossod,outQod,outflagod,&
															&outlossev,outQev,outflagev,outputfile)
															
			implicit none
			
			integer,intent(in) :: nbsympts,nbdiv,rang(2),mdim(2)
			integer,intent(in) :: refl,maxguide,lossflag
			real(dp),intent(in) :: ksympts(2,nbsympts),inpls(*)
			real(dp),intent(in) :: Grat
			character(len=120),intent(in) :: optcell,outputfile
			real(dp),allocatable,intent(out) :: outkls(:,:),outeven(:,:),outodd(:,:)
			real(dp),allocatable,intent(out) :: outlossod(:,:),outQod(:,:),outflagod(:,:)
			real(dp),allocatable,intent(out) :: outlossev(:,:),outQev(:,:),outflagev(:,:)
			
			!local variables
			
			integer ::i,j,nbpts,nbrang,kcount,tmprefl,maxdimeig,divxy(2)
			character(len=120) :: outk="tmpoutk",outfileod,outfileev
			real(dp),allocatable ::Gls(:,:),kmod(:,:),map(:,:)
			real(dp) :: per,phi,center(2),epsclad,Mtrans(2,2),coeff,epsls(3),Gcut
			real(dp) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2),g1(2),g2(2),d
			complex(dp),allocatable ::fftmap(:,:),etamap(:,:)
			
			
			
			
			center=(/0.d0,0.d0/)
			divxy=(/mdim(1)-1,mdim(2)-1/)
			
			outfileod=trim(outputfile) // "_odd"
			outfileev=trim(outputfile) // "_even"			
			nbrang=rang(2)-rang(1)+1
			
			
			!preparation step
			
			!build unit cell
			per=inpls(1)
			
			!inpls(period,radius,angle,eps0,eps1)
			if (optcell .eq. "general") then
				phi= inpls(3)
				call unit_cell2D(per,phi,center,inpls,gen_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,rhex,rgraph,eps0,eps1)
			else if (optcell .eq. "hybrid") then
				phi= pi/3.d0
				call unit_cell2D(per,phi,center,inpls,hexgraph_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,a,rkag,eps0,eps1)
			else if (optcell .eq. "kagome") then
				phi=pi/3.d0
				call unit_cell2D(per,phi,center,inpls,kagome_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,L,angle,eps0,eps1,d)
			else if (optcell .eq. "trianghole") then
				phi=inpls(3)
				call unit_cell2D(per,phi,center,inpls,triang_gen_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
				
			else if (optcell .eq. "crosshatch") then
				phi=inpls(3)
				call unit_cell2D(per,phi,center,inpls,crosshatch_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
				
			!for H cavity
			!inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
			else if (optcell .eq. "h1_cavity") then
				phi=inpls(3)
				call unit_sc2D(phi,inpls,h1_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)

			else if (optcell .eq. "h2_cavity") then
				phi=inpls(3)
				call unit_sc2D(phi,inpls,h2_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			
			else if (optcell .eq. "l3_cavity") then
				phi=inpls(3)
				call unit_sc2D(phi,inpls,l3_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
				
			!for 1D cavity
			
			else if (optcell .eq. "bragg_cavity") then
				phi=0.d0
				call unit_cell1D(per,center,inpls,bragg_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(3) !eps_clad index changed !!
			end if
			
			
			d=inpls(6)
			print *,"Period: ", per
			print *,""
			

			print *,"Angle (deg): ", phi*180.d0/pi
			print *,""

			g1=glatt(:,1)
			g2=glatt(:,2)
			
			Gcut=Grat*sqrt(g1(1)*g1(1)+g1(2)*g1(2))
			call listk(nbsympts,ksympts,nbdiv,outkls,nbpts)
			!fft transform for dielectric map. building epsls
			if (optcell .eq. "bragg_cavity") then
				call glist(g1,g2,Gcut,0,Gls,outk,kcount,1)
			else
				call glist(g1,g2,Gcut,0,Gls,outk,kcount,0)
			end if
			
			call map_fft(map,mdim,fftmap)
			epsls=(/epsclad,real(fftmap(1,1)),epsclad/)
			maxdimeig=kcount*maxguide
			print *,"Stack structure: ",epsls
			print *,""
			print *,"Layer thickness: ",d
			print *,""
			
			!fft transform of inverse dielectric map (Ho's method)
			!Mtrans later can be generalized (depend on the cell vector choice)
			Mtrans(:,1)=(/1.d0,-1.d0/)
			Mtrans(:,2)=(/1.d0,1.d0/)
			if (optcell .eq. "bragg_cavity") then
				Mtrans(:,1)=(/1.d0,0.d0/)
				Mtrans(:,2)=(/0.d0,1.d0/)
				
				print *, "WARNING: 1D Mtrans"
			end if
			
			call etamat(kcount,Gls,mdim(1),Mtrans,fftmap,etamap)
			print *,"teteteta",fftmap(29,29)
			coeff=per/(2.d0*pi)
			
			print *, "Starting the band structure scan: "
			print *,""
			!full solver
			!add save data later
			if (refl ==0) then
				tmprefl=-1
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
														&tmprefl,epsls,etamap,coeff,outodd,rang,&
														&outlossod,outQod,outflagod,lossflag)
				
				
				tmprefl=1
				
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
										&tmprefl,epsls,etamap,coeff,outeven,rang,&
										&outlossev,outQev,outflagev,lossflag)
										
				if (lossflag==0) then
					call save_bs_noloss(nbpts,outkls,g1,g2,outodd,maxdimeig,rang,coeff,epsls,outfileod)
					call save_bs_noloss(nbpts,outkls,g1,g2,outeven,maxdimeig,rang,coeff,epsls,outfileev)
				else if (lossflag==1) then
					call save_bs_loss(nbpts,nbrang,outkls,g1,g2,outodd,maxdimeig,coeff,epsls,&
														&rang,outlossod,outQod,outflagod,outfileod)
														
					call save_bs_loss(nbpts,nbrang,outkls,g1,g2,outeven,maxdimeig,coeff,epsls,&
														&rang,outlossev,outQev,outflagev,outfileev)
				end if
				
										
			else if (refl==-1) then
				
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
										&refl,epsls,etamap,coeff,outodd,rang,&
										&outlossod,outQod,outflagod,lossflag)
				
				if (lossflag==0) then
					call save_bs_noloss(nbpts,outkls,g1,g2,outodd,maxdimeig,rang,coeff,epsls,outputfile)
				else if (lossflag==1) then
					call save_bs_loss(nbpts,nbrang,outkls,g1,g2,outodd,maxdimeig,coeff,epsls,&
														&rang,outlossod,outQod,outflagod,outputfile)
				end if
			
			else if (refl==1) then
				
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
													&refl,epsls,etamap,coeff,outeven,rang,&
													&outlossev,outQev,outflagev,lossflag)
													
													
													
				if (lossflag==0) then
					call save_bs_noloss(nbpts,outkls,g1,g2,outeven,maxdimeig,rang,coeff,epsls,outputfile)
				else if (lossflag==1) then
					call save_bs_loss(nbpts,nbrang,outkls,g1,g2,outeven,maxdimeig,coeff,epsls,&
														&rang,outlossev,outQev,outflagev,outputfile)
				end if
													
			end if
				
			
			!deallocating
			deallocate(map,fftmap,etamap,Gls)
			print *,""
			print *, "Calculation done. Output file saved"
		
		end subroutine band_structure
		
		!Scan of intrinsic loss as function of maximum 2 parameters and then save results to outputfile
		subroutine loss_scan(kvect,nbs,ind,scan1,scan2,optcell,inpls,Grat,&
													&mdim,refl,maxguide,rang,&
													&outE,outloss,outQ,outflag,outputfile)
			
			implicit none
			integer,intent(in) :: mdim(2),nbs(2),ind(2),maxguide,refl,rang
			character(len=120),intent(in) :: optcell,outputfile
			real(dp),intent(in) :: Grat
			real(dp),intent(in) :: scan1(*),scan2(*),kvect(2)
			real(dp),intent(inout) :: inpls(*)
			real(dp),allocatable,intent(out) ::outE(:,:),outloss(:,:)
			real(dp),allocatable,intent(out) ::outQ(:,:)
			real(dp),allocatable,intent(out) ::outflag(:,:)
			
			!local var
			integer ::divxy(2),i,j,tt,kcount,tmpnbs2,dimH
			real(dp),allocatable :: Gls(:,:),kmod(:,:),map(:,:)
			complex(dp),allocatable :: fftmap(:,:),etamap(:,:)
			real(dp) :: per,phi,center(2),epsclad,Gcut,epsls(3),d
			real(dp) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2)
			real(dp) :: Mtrans(2,2),coeff,lossval,flagval,g1(2),g2(2)
			character(len=120) :: outk="tmpoutk",outkag="testkagome"
			integer,allocatable :: Tarr(:,:),Hindex(:,:)
			complex(dp),allocatable:: ABarr(:,:),H(:,:)
			real(dp),allocatable :: earr(:,:),karr(:,:)
			real(dp),allocatable ::tmpeig2(:),tmpeig(:),tmploss(:),RWORK(:)
			complex(dp),allocatable :: WORK(:)
			integer :: LWORK,INFO
			
			
			center=(/0.d0,0.d0/)
			divxy=(/mdim(1)-1,mdim(2)-1/)

			! test outcome dimension, activating outsvar2 if 2-dim scan
			if (nbs(2) .ne. 0) then
				allocate (outE(nbs(1),nbs(2)),outloss(nbs(1),nbs(2)))
				allocate(outQ(nbs(1),nbs(2)),outflag(nbs(1),nbs(2)))

				tmpnbs2=nbs(2)
			
			else 
				allocate (outE(nbs(1),1),outloss(nbs(1),1))
				allocate(outQ(nbs(1),1),outflag(nbs(1),1))

				tmpnbs2=1
			end if
			
			!print *,"tmpnbs2",tmpnbs2
			
			do i=1,nbs(1)
				do j=1,tmpnbs2
					
					!rearranging inpls and create unit cell map
					
					inpls(ind(1))=scan1(i)
					if (tmpnbs2 .ne. 1) then
						inpls(ind(2))=scan2(j)
					end if
					per=inpls(1)
					
					!inpls(period,radius,angle,eps0,eps1,d)
					if (optcell .eq. "general") then
						phi= inpls(3)
						
						call unit_cell2D(per,phi,center,inpls,gen_latt,divxy,map,alatt,asc,glatt,glattxy)
						!print *,"tmpnbs2",tmpnbs2,i,j
						epsclad=inpls(4)
						
					!inpls(period,rhex,rgraph,eps0,eps1,d)
					else if (optcell .eq. "hybrid") then
						phi= pi/3.d0
						call unit_cell2D(per,phi,center,inpls,hexgraph_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
					!inpls(period,a,rkag,eps0,eps1,d)
					else if (optcell .eq. "kagome") then
						phi=pi/3.d0
						call unit_cell2D(per,phi,center,inpls,kagome_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
						!call dwrite_map_tofile(mdim,map,outkag)
					!inpls(period,L,angle,eps0,eps1,d)
					else if (optcell .eq. "trianghole") then
						phi=inpls(3)
						call unit_cell2D(per,phi,center,inpls,triang_gen_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
	
					else if (optcell .eq. "crosshatch") then
						phi=inpls(3)
						call unit_cell2D(per,phi,center,inpls,crosshatch_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
					
					!for H cavity
					!inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
					else if (optcell .eq. "h1_cavity") then
						phi=inpls(3)
						call unit_sc2D(phi,inpls,h1_latt,divxy,map,alatt,asc,glatt,glattxy)
						!call dwrite_map_tofile(mdim,map,outkag)
						epsclad=inpls(4)

					else if (optcell .eq. "h2_cavity") then
						phi=inpls(3)
						call unit_sc2D(phi,inpls,h2_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
						
					else if (optcell .eq. "l3_cavity") then
						phi=inpls(3)
						call unit_sc2D(phi,inpls,l3_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(4)
					
					!for 1D cavity
					
					else if (optcell .eq. "bragg_cavity") then
						phi=0.d0
						call unit_cell1D(per,center,inpls,bragg_latt,divxy,map,alatt,asc,glatt,glattxy)
						epsclad=inpls(3)
						
					end if
					d=inpls(6)
					!fft transform for dielectric map. building epsls
					g1=glatt(:,1)
					g2=glatt(:,2)
					
					!common part (selecting cutoff k-vector)
					Gcut=Grat*sqrt(g1(1)*g1(1)+g1(2)*g1(2))
					if (optcell .eq. "bragg_cavity") then
						call glist(g1,g2,Gcut,0,Gls,outk,kcount,1)
					else
						call glist(g1,g2,Gcut,0,Gls,outk,kcount,0)
					end if
					
					call gmodlist(kvect,Gls,kcount,g1,g2,0,outk,kmod)
					!end of g-vector process
					
					call map_fft(map,mdim,fftmap)
					epsls=(/epsclad,real(fftmap(1,1)),epsclad/)
					print *,"*****************************************************"
					print *,"Scan step (1,2): ",i," , ",j
					print *,"Value of 1st variable: ",inpls(ind(1))
					if (tmpnbs2 .ne. 1) then
						print *,"Value of 2nd variable: ",inpls(ind(2))
					end if
					!fft transform of inverse dielectric map (Ho's method)
					!Mtrans later can be generalized (depend on the cell vector choice)
					! put in the variable list ?
					Mtrans(:,1)=(/1.d0,-1.d0/)
					Mtrans(:,2)=(/1.d0,1.d0/)
					
					if (optcell .eq. "bragg_cavity") then
						Mtrans(:,1)=(/1.d0,0.d0/)
						Mtrans(:,2)=(/0.d0,1.d0/)
						print *, "WARNING: 1D Mtrans"
					end if
					call etamat(kcount,Gls,mdim(1),Mtrans,fftmap,etamap)
					
					!infoarray
					
					call info_array(kcount,maxguide,kmod,d,refl,epsls,Tarr,earr,ABarr,karr,Hindex,dimH,0)
					
					!building hamiltonian and solve for eigval and eigvect
					call h_guided(kcount,maxguide,dimH,kmod,d,refl,epsls,etamap,Hindex,Tarr,earr,ABarr,karr,H,0)
					!print *,"tototatat",dimH
					allocate (WORK(1),tmpeig2(dimH),RWORK(3*dimH-2))
					call zheev('N','U',dimH,H,dimH,tmpeig2,WORK,-1,RWORK,INFO)
					LWORK=nint(real(WORK(1)))
					deallocate(WORK)
					allocate(WORK(LWORK))
					
					!all process will store eigenvector from now on (useful to plot field pattern later)
					call zheev('V','U',dimH,H,dimH,tmpeig2,WORK,LWORK,RWORK,INFO)
					allocate (tmpeig(dimH),tmploss(dimH))
					!convert to w/c
					do tt=1,dimH
						tmpeig(tt)=sqrt(tmpeig2(tt))
					end do
					
					!normalized energy (aw/2pic)
					coeff= per/(2.d0*pi)
					outE(i,j)=coeff*tmpeig(rang)
					
					!loss calculation
					
					call loss_calc(kcount,maxguide,dimH,d,epsls,kvect,&
												&(/rang,rang/),Hindex,Tarr,kmod,earr,ABarr,karr,H,etamap,tmpeig,tmploss,0,refl)
												
					lossval=coeff*tmploss(rang)
					
					if (tmploss(rang) <9998.d0) then
						flagval=1.d0
					else
						flagval=-1.d0
					end if
					
					outflag(i,j)=flagval
					outloss(i,j)=lossval
					outQ(i,j)=outE(i,j)/(2.d0*lossval)
					
					!end of calculation in each loop
					
					!deallocating
					deallocate(map,fftmap,etamap)
					deallocate(Tarr,earr,ABarr,karr,Hindex,H)
					deallocate (WORK,RWORK,tmpeig2,tmpeig,tmploss)
					deallocate(Gls,kmod)
				end do
			end do
			
			!deallocating

			call save_sscan(nbs(1),nbs(2),scan1,scan2,outE,outloss,outQ,outflag,outputfile)

		end subroutine loss_scan		


		!calculate the band structure and loss (if lossflag activated=1)
		!refl=0 means both even and odd modes will be calculated
		!outkls will have transposed dimension (due to what was defined in mode_solver code)
		!180428: add save procedures. d moved in inpls
		!input kls in reduced coordinate, remember to convert to cartesian coordinate before 
		!calling fullhsolver subroutine
		
		subroutine save_bs_BZ(nbpts,kvect,g1,g2,eigls,dimeig,rang,eigcoeff,epsls,output)
		
			implicit none
			integer,intent(in) :: nbpts,rang(2),dimeig
			real(dp),intent(in) ::kvect(2,nbpts),g1(2),g2(2),epsls(3)
			real(dp),intent(in) :: eigls(dimeig,nbpts),eigcoeff
			character(len=120),intent(in) :: output
			
			!local variable
			integer ::i,j
			real(dp) :: dummy,kdist(nbpts),tmpdist,tmpdk(2),lline(nbpts)
			logical ::exist
			
			!calculating the plot distance for all k-pts
			dummy=0.d0
			kdist(1)=dummy
			!lightline
			lline(1)=eigcoeff*sqrt((kvect(1,1)**2.d0)+(kvect(2,1)**2.d0))/sqrt(epsls(1))
			do i=2,nbpts
				!g1 and g2 are already in cartesian units
				!tmpdk=(kvect(1,i)*g1+kvect(2,i)*g2)-(kvect(1,i-1)*g1+kvect(2,i-1)*g2)
				tmpdk=kvect(:,i)-kvect(:,i-1)
				tmpdist=sqrt((tmpdk(1)**2.d0)+(tmpdk(2)**2.d0))
				dummy=dummy+tmpdist
				kdist(i)=dummy
				lline(i)=eigcoeff*sqrt((kvect(1,i)**2.d0)+(kvect(2,i)**2.d0))/sqrt(epsls(1))
			end do
			
			!check carefully
			kdist=(1.d0/dummy)*kdist
			
			inquire(file=output,exist=exist)
			if (exist) then
				open(1, file=output, status="replace", action="write",form="formatted")
			else
				open(1, file=output, status="new", action="write",form="formatted")

			end if
			
			do i=1,nbpts

				write(1, '(10000f20.6)')lline(i),( eigls(j,i) ,j=rang(1),rang(2))
				
			end do
			
			close(1)

		
		end subroutine save_bs_BZ
		
		
		subroutine bzscan_band 		(nbdiv,optcell,rang,inpls,&
															&Grat,mdim,refl,maxguide,lossflag,&
															&outkls,outeven,outodd,outputfile,&
															&outlossod,outQod,outflagod,outlossev,outQev,outflagev)
															
			implicit none
			
			integer,intent(in) :: rang(2),mdim(2),nbdiv(2)
			integer,intent(in) :: refl,maxguide,lossflag
			real(dp),intent(in) :: inpls(*)
			real(dp),intent(in) :: Grat
			character(len=120),intent(in) :: optcell,outputfile
			real(dp),allocatable,intent(out) :: outkls(:,:),outeven(:,:),outodd(:,:)
			real(dp),allocatable,intent(out) :: outlossod(:,:),outQod(:,:),outflagod(:,:)
			real(dp),allocatable,intent(out) :: outlossev(:,:),outQev(:,:),outflagev(:,:)
			
			!local variables
			
			integer ::i,j,nbpts,nbrang,kcount,tmprefl,maxdimeig,divxy(2)
			character(len=120) :: outk="tmpoutk",outfileod,outfileev
			real(dp),allocatable ::Gls(:,:),kmod(:,:),map(:,:)
			real(dp) :: per,phi,center(2),epsclad,Mtrans(2,2),coeff,epsls(3),Gcut
			real(dp) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2),g1(2),g2(2),d
			complex(dp),allocatable ::fftmap(:,:),etamap(:,:)
			
			
			
			
			center=(/0.d0,0.d0/)
			divxy=(/mdim(1)-1,mdim(2)-1/)
			
			outfileod=trim(outputfile) // "_odd"
			outfileev=trim(outputfile) // "_even"
			!creating k-list and allocating the outputs
			
			nbrang=rang(2)-rang(1)+1
			
			
			!preparation step
			
			!build unit cell
			per=inpls(1)
			
			!inpls(period,radius,angle,eps0,eps1)
			if (optcell .eq. "general") then
				phi= inpls(3)
				call unit_cell2D(per,phi,center,inpls,gen_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,rhex,rgraph,eps0,eps1)
			else if (optcell .eq. "hybrid") then
				phi= pi/3.d0
				call unit_cell2D(per,phi,center,inpls,hexgraph_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,a,rkag,eps0,eps1)
			else if (optcell .eq. "kagome") then
				phi=pi/3.d0
				call unit_cell2D(per,phi,center,inpls,kagome_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
			!inpls(period,L,angle,eps0,eps1,d)
			else if (optcell .eq. "trianghole") then
				phi=inpls(3)
				call unit_cell2D(per,phi,center,inpls,triang_gen_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)
				
			else if (optcell .eq. "crosshatch") then
				phi=inpls(3)
				call unit_cell2D(per,phi,center,inpls,crosshatch_latt,divxy,map,alatt,asc,glatt,glattxy)
				epsclad=inpls(4)

			end if
			
			
			d=inpls(6)
			print *,"Period: ", per
			print *,""
			

			print *,"Angle (deg): ", phi*180.d0/pi
			print *,""

			g1=glatt(:,1)
			g2=glatt(:,2)
			
			Gcut=Grat*sqrt(g1(1)*g1(1)+g1(2)*g1(2))
			
			call bz_2d_sampling(sqrt(epsclad),per,g1,g2,nbdiv,outkls)
			nbpts=nbdiv(1)*nbdiv(2)
			!Modified for list of kpts
			
			!fft transform for dielectric map. building epsls
			
			!To be modified in the new code with introduction of k into cutoff
			call glist(g1,g2,Gcut,0,Gls,outk,kcount,0)

			
			call map_fft(map,mdim,fftmap)
			epsls=(/epsclad,real(fftmap(1,1)),epsclad/)
			maxdimeig=kcount*maxguide
			print *,"Stack structure: ",epsls
			print *,""
			print *,"Layer thickness: ",d
			print *,""
			
			!fft transform of inverse dielectric map (Ho's method)
			!Mtrans later can be generalized (depend on the cell vector choice)
			! put in the variable list ?
			Mtrans(:,1)=(/1.d0,-1.d0/)
			Mtrans(:,2)=(/1.d0,1.d0/)
			
			call etamat(kcount,Gls,mdim(1),Mtrans,fftmap,etamap)
			coeff=per/(2.d0*pi)
			
			print *, "Starting the band structure scan: "
			print *,""
			!full solver
										
			if (refl==-1) then
				
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
										&refl,epsls,etamap,coeff,outodd,rang,&
										&outlossod,outQod,outflagod,lossflag)
				

				call save_bs_BZ(nbpts,outkls,g1,g2,outodd,maxdimeig,rang,coeff,epsls,outputfile)
			
			else if (refl==1) then
				
				call full_h_solver(outkls,nbpts,kcount,maxguide,Gls,g1,g2,d,&
													&refl,epsls,etamap,coeff,outeven,rang,&
													&outlossev,outQev,outflagev,lossflag)
													
					
				call save_bs_BZ(nbpts,outkls,g1,g2,outeven,maxdimeig,rang,coeff,epsls,outputfile)

													
			end if
				
			
			!deallocating
			deallocate(map,fftmap,etamap,Gls)
			print *,""
			print *, "Calculation done. Output file saved"
		
		end subroutine bzscan_band
		
		
	!load photonic band from file for the calculation of density of states (DOS)
	subroutine load_phoband(filename,nbband,nbk,fmtout,Mband,limfreq)
	
		implicit none
		character(len=120),intent(in) ::filename,fmtout
		integer,intent(in)::nbband,nbk
		real(dp),allocatable,intent(out) ::Mband(:,:),limfreq(:)
		
		!local variables
		integer ::ii,jj,ierr
		
		allocate(Mband(nbk,nbband),limfreq(nbk))
		
		open(1, file=filename, action="read")
		
		do ii=1,nbk
			read(1,fmtout,iostat=ierr) limfreq(ii),(Mband(ii,jj),jj=1,nbband)
		end do
		
		close(1)
		
	end subroutine load_phoband
	
	!Calculation of DOS
	subroutine dos_pho_2d(Ener,nbk,nbband,Mband,limfreq,broad,g1,g2,&
												&dosin,dosout,dosfull)
		
		implicit none
		real(dp),intent(in) ::Ener,broad,g1(2),g2(2)
		integer,intent(in) ::nbk,nbband
		real(dp),intent(in)::Mband(nbk,nbband),limfreq(nbk)
		real(dp),intent(inout)::dosin,dosout,dosfull
		
		!local variables
		integer :: ii,jj
		real(dp) ::tmpval,flim,Ebjj,dVBZ2d
		
		dosout=0.d0
		dosin=0.d0
		dosfull=0.d0
		
		dVBZ2d=abs(g1(1)*g2(2)-g2(1)*g1(2))
		do ii=1,nbk
			flim=limfreq(ii)
			do jj=1,nbband
				Ebjj=Mband(ii,jj)
				tmpval=(dVBZ2d/(4.d0*pi*pi))*(broad/(2.d0*pi))/(((Ener-Ebjj)**2.d0)+&
								&(broad*broad/4.d0))
				
				!upper side
				if (Ebjj >= flim) then
					dosout=dosout+tmpval
				else
					dosin=dosin+tmpval
				end if
			end do
		end do
		
		dosfull=dosin+dosout

	end subroutine dos_pho_2d
	
		
end module scan_process		