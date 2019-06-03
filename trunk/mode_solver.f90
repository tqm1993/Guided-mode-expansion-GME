
!This module contains solver for the guided mode which form the basis, and construct the hamiltonian for photonic crystal
!Details of the theory can be found in L. C. Andreani and D. Gerace work: PRB 73, 235114 (2016)

module mode_solver
	
	use phys_cte
	use algo
	use geometry_defs
	
	
	implicit none
	
	real(dp) :: dk=1.d-5
	real(dp) :: nbsplit=1600
	real(dp) :: pgtol=1.d-5,factr=1.d+5,tol=3.d-3
	integer :: m_BFGS=5
	real(dp) ::tolval=1.d-6
	real(dp) ::failsafe=1.d-6
	real(dp) ::divz=50
	
	contains
		
		!k stands for w/c
		! g is the considered wavevector
		!trans=0,1 for TE,TM
		!refl=-1,1 for odd and even
		!val: Function value
		!in mirror symmetry
		subroutine eq_mirror(d,k,g,epsls,trans,refl,val)
			
			implicit none
			
			real(dp),intent(in) :: d,k,g,epsls(3)
			integer,intent(in) :: trans,refl
			real(dp),intent(out) ::val
			
			!local variable
			real(dp) ::q,khi1
			
			q=sqrt(epsls(2)*(k**2.d0)-(g**2.d0))
			khi1=sqrt((g**2.d0)-epsls(1)*(k**2.d0))
			!TE-1
			if ((trans==0) .and. (refl==-1)) then
				val=abs(q*cos(q*d/2.d0)+khi1*sin(q*d/2.d0))
			
			!TM-1
			else if ((trans==1) .and. (refl==-1)) then
				val=abs((q/epsls(2))*sin(q*d/2.d0)-(khi1/epsls(1))*cos(q*d/2.d0))
			
			!TE+1
			else if ((trans==0) .and. (refl==1)) then
				val=abs(q*sin(q*d/2.d0)-khi1*cos(q*d/2.d0))	
				
			!TM+1
			else if ((trans==1) .and. (refl==1)) then
				val=abs((q/epsls(2))*cos(q*d/2.d0)+(khi1/epsls(1))*sin(q*d/2.d0))
				
			end if

		end subroutine eq_mirror
		
		subroutine profile(d,g,kbds,epsls,trans,refl,output)
			
			implicit none
			real(dp),intent(in) :: d,g,kbds(2),epsls(3)
			integer,intent(in) :: trans,refl
			character(len=120),intent(in) ::output
			
			!local variable
			integer ::nbpts=10000
			integer ::i
			real(dp) ::step,k,val
			real(dp),allocatable :: M(:,:)
			logical ::exist
			
			inquire(file=output,exist=exist)
			if (exist) then
				open(1, file=output, status="replace", action="write",form="formatted")
			else
				open(1, file=output, status="new", action="write",form="formatted")

			end if
			
			step=(kbds(2)-kbds(1))/(nbpts-1)
			
			allocate(M(2,nbpts-2))
			do i=1,nbpts-2
				k=kbds(1)+1.d0*i*step
				call eq_mirror(d,k,g,epsls,trans,refl,val)
				M(:,i)=(/k,val/)
				write(1, '(f20.6,f20.6)')M(1,i),M(2,i)
			
			end do

			close(1)
			deallocate(M)
					
		end subroutine profile

		! Calculate the function value and function derivative to be used in the BFGS solver for the energy of guided modes in the basis
		subroutine func_deriv(k,d,epsls,g,trans,refl,f,df)
		
			implicit none
			
			real(dp),intent(in) :: k(1),epsls(3),g,d
			integer,intent(in) :: trans,refl
			real(dp),intent(inout) ::f,df(1)
			
			!local variable
			real(dp) :: q,khi1,khi3
			real(dp) :: valb,valf
			integer :: i
			
			
				
			call eq_mirror(d,k(1),g,epsls,trans,refl,f)
			call eq_mirror(d,k(1)-dk,g,epsls,trans,refl,valb)
			call eq_mirror(d,k(1)+dk,g,epsls,trans,refl,valf)
				

			df(1)=(valf-valb)/(2.d0*dk)
		
		end subroutine func_deriv
	
		!Solver for guided mode in the basis
		!refl=-1,1 for odd and even mode
		!pola: Contains information on whether the guided mode is TE-polarized or TM-polarized
		!if number of modes found is less than alpha, then nbsol <alpha
		subroutine k_solver(alpha,d,g,epsls,refl,k,pola,nbsol,verbose)
			
			implicit none
			
			real(dp),intent(in) ::epsls(3),g,d
			integer,intent(in) :: refl,alpha,verbose
			real(dp),allocatable,intent(out) :: k(:)
			integer,allocatable,intent(out) :: pola(:)
			integer,intent(inout) ::nbsol
			
				
			!local variable
			real(dp) ::kbds(2),ktmp,kT(2,alpha),tmpval,incstep,valTE,valTM
			real(dp),allocatable :: kall(:),polaall(:)
			integer :: n = 1, iprint = -1
      character(len=60) :: task, csave
      logical :: lsave(4)
      integer :: isave(44),nbmode(2),i,j,trans,nbfound,flag,loccount,iniflag,pair
      integer,allocatable ::idx(:)
      real(dp) :: f,ff
      real(dp) :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(dp), allocatable  :: x(:), l(:), u(:), df(:),dff(:), wa(:)
      
			allocate ( nbd(n), x(n), l(n), u(n), df(n),dff(n) )
      allocate ( iwa(3*n) )
      allocate ( wa(2*m_BFGS*n + 5*n + 11*m_BFGS*m_BFGS + 8*m_BFGS) )
      
      kbds=(/(g/sqrt(epsls(2))),(g/dmax1(sqrt(epsls(1)),sqrt(epsls(3))))/)

			nbfound=0
			loccount=0
			iniflag=0
			pair=0
			allocate (kall(2*alpha),polaall(2*alpha))
			incstep=(kbds(2)-kbds(1))/(1.d0*nbsplit)
			if (refl==-1) then
				trans=1
			else 
				trans=0
			end if
			!normal state
			l(1)=kbds(1)+dk
			u(1)=l(1)+incstep
			flag=0
			nbd(1)=2

			mainloop: do i=1,2*alpha
				!starting with a guess value for the frequency
				task='START'
				!initializing f
				x(1)=l(1)+dk
				call func_deriv(x,d,epsls,g,trans,refl,f,df)
				
				
				do while (f>tol)
					
					if ((loccount>0) .and. (iniflag>0)) then
						l(1)=l(1)+incstep
						u(1)=l(1)+incstep
						x(1)=l(1)+dk
					end if
					
					if (l(1)<kbds(2)) then
						task='START'
						!L-BFGS-B call.
						do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 

							call setulb (n,m_BFGS,x,l,u,nbd,f,df,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)
							
							if (task(1:2) .eq. 'FG') then
								
								call func_deriv(x,d,epsls,g,trans,refl,f,df)

							end if
						end do
						iniflag=iniflag+1

					else 
						flag=1
						exit
						
					end if
					loccount=loccount+1
					
				end do
				
				
				if (flag==0) then
					!reset the solver condition
					l(1)=x(1)
					u(1)=l(1)+incstep
						
					!end if
					if (x(1)< (kbds(2)-failsafe)) then
						loccount=0
						nbfound=nbfound+1
						
						kall(nbfound)=x(1)
						
						polaall(nbfound)=trans
					end if
					
					if (trans==1) then
						trans=0
					else
						trans=1
					end if
				
				else
					!exit the main loop. No more guided modes can be found
					exit
				
				end if
			
			end do mainloop
								

			!sorting states
			allocate(idx(nbfound))
			
			call dheapsort(nbfound,kall(1:nbfound),idx) !increasing order
			
			if (nbfound > alpha) then
				nbsol=alpha
			else
				nbsol=nbfound
			end if
			
			allocate(k(nbsol),pola(nbsol))
			
			k(1:nbsol)=kall(idx(1:nbsol))
			pola(1:nbsol)=polaall(idx(1:nbsol))

			deallocate(kall,polaall,idx)
			
	
		
		end subroutine k_solver
		

		!Extract the wavefunction amplitude of the guided modes
		subroutine coeff_abcd(k,epsls,d,g,trans,A2,A3,B1,B2,rawraw)
		
			implicit none
			
			real(dp),intent(in) ::k,epsls(3),d,g
			integer,intent(in) ::trans
			complex(dp),intent(inout) :: A2,A3,B1,B2
			real(dp),intent(inout) :: rawraw
			
			!local variable
			real(dp) :: khi1,q,khi3,cosqd,sinqd,rawcoeff2,rawcoeff,qd
			complex(dp) :: expqd,expqdm
			
			q=sqrt(epsls(2)*(k**2.d0)-(g**2.d0))
			khi1=sqrt((g**2.d0)-epsls(1)*(k**2.d0))
			khi3=sqrt((g**2.d0)-epsls(3)*(k**2.d0))

			qd=q*d
			cosqd=cos(q*d)
			sinqd=sin(q*d)
			expqd=exp(cmplx(0.d0,qd/2.d0,dp))
			expqdm=exp(cmplx(0.d0,-qd/2.d0,dp))
			

			!Initializing B1
			if (trans==0) then
				B1=cmplx(1.d0,0.d0,dp)
				B2=cmplx(q,khi1,dp)*expqdm*B1/(2.d0*q)
				A2=cmplx(q,-khi1,dp)*expqd*B1/(2.d0*q)
				A3=cmplx((q*(khi3-khi1)*cosqd+(q*q+khi1*khi3)*sinqd)*(1.d0/(2.d0*q*khi3)),0.d0,dp)
			else if (trans==1) then
				B1=cmplx(1.d0,0.d0,dp)
				B2=cmplx(q/epsls(2),khi1/epsls(1),dp)*expqdm*B1/(2.d0*q/epsls(2))
				A2=cmplx(q/epsls(2),-khi1/epsls(1),dp)*expqd*B1/(2.d0*q/epsls(2))
				A3=cmplx(((q/epsls(2))*(khi3/epsls(3)-khi1/epsls(1))*cosqd+&
								&((q/epsls(2))*(q/epsls(2))+(khi1/epsls(1))*(khi3/epsls(3)))*sinqd)*&
								&(1.d0/(2.d0*(q/epsls(2))*(khi3/epsls(3)))),0.d0,dp)
			end if
			

			
			!normalization for TE (trans==0) and TM(trans==1)
			if (trans==0) then
				
				rawcoeff2=((khi1**2.d0)+(g**2.d0))*(abs(B1)**2.d0)/(2.d0*khi1)+&
								&((khi3**2.d0)+(g**2.d0))*(abs(A3)**2.d0)/(2.d0*khi3)+&
								&d*((g**2.d0)+(q**2.d0))*2.d0*(abs(A2)**2.d0)+&
								&((g**2.d0)-(q**2.d0))*(real(dconjg(A2)*B2+dconjg(B2)*A2))*sinqd/q
			
			else if (trans==1) then
				
				rawcoeff2=(abs(B1)**2.d0)/(2.d0*khi1)+(abs(A3)**2.d0)/(2.d0*khi3)+&
								&d*2.d0*(abs(A2)**2.d0)+(real(dconjg(A2)*B2+dconjg(B2)*A2))*sinqd/q
			
			end if
			
			rawcoeff=sqrt(rawcoeff2)
			rawraw=rawcoeff
			
			B1=B1/rawcoeff
			B2=B2/rawcoeff
			A2=A2/rawcoeff
			A3=A3/rawcoeff
			
		end subroutine coeff_abcd
	
	
	
		!calculate the integral appear in the element of the hamiltonian
		subroutine guided_integ(ku,kv,epsls,gu,gv,d,I1,I2m,I2f,I3)
			
			implicit none
			real(dp),intent(in) :: ku,kv,epsls(3),gu,gv,d
			real(dp),intent(out) :: I1,I2m,I2f,I3
			
			!local variable
			real(dp) :: qu,qv
			
			!I1,I-2,I2,I3

			qu=sqrt((gu**2.d0)-epsls(1)*(ku**2.d0))
			qv=sqrt((gv**2.d0)-epsls(1)*(kv**2.d0))
			I1=1.d0/(qu+qv)

			qu=sqrt((gu**2.d0)-epsls(3)*(ku**2.d0))
			qv=sqrt((gv**2.d0)-epsls(3)*(kv**2.d0))
			I3=1.d0/(qu+qv)				
			
			
			qu=sqrt(epsls(2)*(ku**2.d0)-(gu**2.d0))
			qv=sqrt(epsls(2)*(kv**2.d0)-(gv**2.d0))

			
			I2f=sin(0.5d0*d*(qu+qv))/(0.5d0*(qu+qv))
			
			if (abs(qu-qv)<=tolval) then
				!lim of sinus cardinal function when x converge towards 0
				I2m=d
			else
				I2m=sin(0.5d0*d*(qu-qv))/(0.5d0*(qu-qv))
			end if
			
		end subroutine guided_integ
		
		!extract the fourier coefficient for an matrix element at coordinate (u(1)v(1),u(2)v(2))
		!Each coordinate is a pair of u and v
		subroutine gugv_coeff(u,v,dimfft,fftmap,val)
			
			implicit none
			integer,intent(in) ::u(2),v(2),dimfft
			complex(dp),intent(in):: fftmap(dimfft,dimfft)
			complex(dp),intent(out) ::val
			
			!local variable
			integer ::uv(2),xtmp,ytmp,redflag
			
			uv=u-v
			xtmp=uv(1)+1
			ytmp=uv(2)+1
			redflag=0
			if ((abs(xtmp)>dimfft) .or. (abs(ytmp)>dimfft)) then
				redflag=1
				val=cmplx(0.d0,0.d0,dp)
				print *, "index higher than matrix dimension.Process stopped"
			
			end if
			
			if (uv(2) <0) then
				ytmp=dimfft+uv(2)+1
			end if
			
			if (uv(1) <0) then
				xtmp=dimfft+uv(1)+1
			end if

			if (redflag==0) then
				val=fftmap(ytmp,xtmp)
			end if
			
			
		end subroutine gugv_coeff
		
		
		!find list of G vector inside cutoff value. Return Gls of G-vector
		!CAREFUL: Minv here is the inverse transformation matrix instead of Mtrans
		!Gls(2,Nk). flag =1 activate the writing of Gls to file
		!Gcut absolute value
		subroutine glist(g1,g2,Gcut,flag,Gls,output,kcount,flag1D)
		
		
			implicit none
			real(dp),intent(in) ::g1(2),g2(2),Gcut
			integer,intent(in) ::flag,flag1D
			character(len=120),intent(in) :: output
			real(dp),allocatable,intent(out) :: Gls(:,:)
			integer,intent(out) ::kcount
			
			!local variable
			real(dp) ::g1mod,g2mod,Gtol,grec(2),grecmod
			integer :: i,j,cutg1,cutg2
			integer,allocatable :: idx(:)
			real(dp),allocatable ::ktmp(:,:)
			logical :: exist
			
			
			g1mod=sqrt(g1(1)**2.d0+g1(2)**2.d0)
			g2mod=sqrt(g2(1)**2.d0+g2(2)**2.d0)
			
			cutg1=nint(Gcut/g1mod)+1
			if (flag1D==0) then
				cutg2=nint(Gcut/g2mod)+1
				allocate (ktmp(3,(2*cutg1+1)*(2*cutg2+1)))
			else
				cutg2=0
				allocate (ktmp(3,2*cutg1+1))
				
			end if
			kcount=0
			Gtol=1.d-3*Gcut
			!print *,"flag1D,flag1D",flag1D
			if (flag1D==0) then
				do i=-cutg1, cutg1
					do j=-cutg2,cutg2
						
						!separate G from G+k.
						!TODO: check carefully
						grec=1.d0*i*g1+1.d0*j*g2
						grecmod=sqrt(grec(1)**2.d0+grec(2)**2.d0)
						if (grecmod <= (Gcut+Gtol)) then
							kcount=kcount+1
							ktmp(:,kcount)=(/1.d0*i,1.d0*j,grecmod/)	
						end if
						
					end do
				end do
				
			else
					
					do i=-cutg1, cutg1
						!separate G from G+k.
						!TODO: check carefully
						grec=1.d0*i*g1
						grecmod=sqrt(grec(1)**2.d0+grec(2)**2.d0)
						if (grecmod <= (Gcut+Gtol)) then
							kcount=kcount+1
							ktmp(:,kcount)=(/1.d0*i,0.d0,grecmod/)	
						end if
				end do
				
			end if
			
			allocate (idx(kcount),Gls(3,kcount))
			
			call dheapsort(kcount,ktmp(3,1:kcount),idx) !increasing order
			
			Gls(:,1:kcount)=ktmp(:,idx(:))
			
			deallocate (ktmp,idx)
			
			if (flag==1) then
				
				inquire(file=output,exist=exist)
				if (exist) then
					open(1, file=output, status="replace", action="write",form="formatted")
				else
					open(1, file=output, status="new", action="write",form="formatted")

				end if
				
				do i=1,kcount
					write(1, '(f20.6,f20.6,f20.6)')Gls(1,i),Gls(2,i),Gls(3,i)
				end do
					close(1)
					
			end if

		end subroutine glist
		
	
		!splitted from glist to avoid recalculate k-vector below cutoff at each k-point
		subroutine gmodlist(kvect,Gls,kcount,g1,g2,flag,output,kmod)
		
			implicit none
			integer,intent(in) :: kcount,flag
			real(dp),intent(in) :: kvect(2),Gls(3,kcount),g1(2),g2(2)
			character(len=120),intent(in) :: output
			real(dp),allocatable,intent(out) :: kmod(:,:)
			
			!local variable
			real(dp) ::tmpg(2),modg,toto,gg(2)
			integer :: i
			logical ::exist
			
			allocate(kmod(3,kcount))
			do i =1,kcount
				
				tmpg=Gls(1,i)*g1+Gls(2,i)*g2 +kvect
				modg=sqrt(tmpg(1)**2.d0+tmpg(2)**2.d0)
				
				kmod(:,i)=(/(tmpg(1)/modg),(tmpg(2)/modg),modg/)
				
			end do
				
			if (flag==1) then
				
				inquire(file=output,exist=exist)
				if (exist) then
					open(1, file=output, status="replace", action="write",form="formatted")
				else
					open(1, file=output, status="new", action="write",form="formatted")

				end if
				
				do i=1,kcount
					gg=Gls(1,i)*g1+Gls(2,i)*g2
					toto=sqrt(gg(1)*gg(1)+gg(2)*gg(2))
					gg=gg+kvect
					write(1, '(f20.6,f20.6,f20.6,f20.6,f20.6,f20.6)')Gls(1,i),Gls(2,i),toto,gg(1),gg(2),kmod(3,i)
																																			
				end do
					close(1)
					
			end if
		
		end subroutine gmodlist
		
			
		! Calculate the inversion of dielectric matrix
		subroutine etamat(dimk,Gls,dimmap,Mtrans,epsmap,etamap)
		
		
			implicit none
			integer,intent(in) :: dimk,dimmap
			real(dp),intent(in) :: Gls(3,dimk),Mtrans(2,2)
			complex(dp),intent(in) :: epsmap(dimmap,dimmap)
			complex(dp),allocatable,intent(out) :: etamap(:,:)
			
			!local variable
			integer ::utmp(2),vtmp(2),i,j
			complex(dp) :: tmpval
			!complex(dp) :: token(dimk,dimk),tokcheck(dimk,dimk)
			
			!variable for inverse method
			
			integer ::INFO,LWORK
			integer ::IPIV(dimk)
			complex(dp),allocatable ::WORK(:)
			
			allocate (etamap(dimk,dimk))
			do i=1,dimk
				do j=1,dimk
					
					call iconvert_to_sc(Mtrans,Gls(1:2,i),utmp)
					call iconvert_to_sc(Mtrans,Gls(1:2,j),vtmp)
					
					call gugv_coeff(utmp,vtmp,dimmap,epsmap,tmpval)
					etamap(i,j)=tmpval
						
				end do
			end do
			

			call zgetrf(dimk,dimk,etamap,dimk,IPIV,INFO)
			!TODO: test+verbose part. check query LWORK=-1
			LWORK=-1
			allocate(WORK(1))
			call zgetri(dimk,etamap,dimk,IPIV,WORK,LWORK,INFO)
			LWORK=nint(real(WORK(1)))
			deallocate (WORK)
			allocate (WORK(LWORK))
			call zgetri(dimk,etamap,dimk,IPIV,WORK,LWORK,INFO)
			deallocate(WORK)

			
		end subroutine etamat
		
		
		!only store the non-zero part of g and e. kmod(3,:) store the module of vector g
		!inverse matrix will be added later for n(G,G*)
		!T=2 mark unfounded mode
		subroutine info_array(kcount,maxguide,kmod,d,refl,epsls,Tarr,earr,ABarr,karr,Hindex,dimH,verbose)
		
			
			implicit none
			integer,intent(in) :: refl,kcount,maxguide,verbose
			real(dp),intent(in) ::d,epsls(3),kmod(3,kcount)
			integer,allocatable,intent(out) :: Tarr(:,:),Hindex(:,:)
			complex(dp),allocatable,intent(out) :: ABarr(:,:)
			real(dp),allocatable,intent(out) :: earr(:,:),karr(:,:)
			integer,intent(inout) :: dimH
			
			
			!local variable 
			
			integer ::i
			
			
			!solver part
			integer :: tt,nbsol,Hcount
			real(dp),allocatable :: k(:)
			integer,allocatable :: pola(:)
			complex(dp) ::A2,A3,B1,B2
			real(dp) :: rawraw
			
			dimH=kcount*maxguide
			!print *,"ideal dimH: ",dimH
			!print *,"kcount",kcount
			allocate(Tarr(maxguide,kcount),ABarr(4*maxguide,kcount),karr(maxguide,kcount),earr(2,kcount))
			
			do i=1,kcount
				
				!check rank in other subroutines
				earr(:,i)=(/-kmod(2,i),kmod(1,i)/)
				!solver part
				call k_solver(maxguide,d,kmod(3,i),epsls,refl,k,pola,nbsol,verbose)
				do tt=1,nbsol
					

					karr(tt,i)=k(tt)
					Tarr(tt,i)=pola(tt)

					
					call coeff_abcd(k(tt),epsls,d,kmod(3,i),pola(tt),A2,A3,B1,B2,rawraw)
					ABarr(4*(tt-1)+1:4*(tt-1)+4,i)=(/A2,A3,B1,B2/)
				
				end do
				
				!filler (trans=2 signals guided mode cannot be found)
				! dimH recalculated each time there is a lack of guided mode (trans=2)
				if (nbsol <maxguide) then
				
					do tt=nbsol+1,maxguide
					
						karr(tt,i)=0.d0
						Tarr(tt,i)=2
						dimH=dimH-1
						ABarr(4*(tt-1)+1:4*(tt-1)+4,i)=(/cmplx(0.d0,0.d0,dp),&
									&cmplx(0.d0,0.d0,dp),cmplx(0.d0,0.d0,dp),cmplx(0.d0,0.d0,dp)/)
					
					end do
				
				end if

				deallocate (k,pola)
	
			end do
			
			!180411: dimensions of H calculated after all partial info array are formed
			!shape H index (containing the pw index and the mode order)
			
			allocate (Hindex(2,dimH))
			Hcount=0
			do i=1,kcount
				do tt=1,maxguide
					if (Tarr(tt,i) .ne. 2) then
						Hcount=Hcount+1
						Hindex(:,Hcount)=(/i,tt/)
					end if
				end do
			end do

			
		end subroutine info_array
			
			
		!DOS for loss calc
		subroutine dos(epsval,k0,gabs,dosval)
		
			implicit none
			real(dp),intent(in) ::epsval,k0,gabs
			real(dp),intent(out) :: dosval
			
			dosval=sqrt(epsval)/(4.d0*pi*sqrt(k0**2.d0-(gabs**2.d0)/epsval))
		
		end subroutine dos	
	

		!transfer matrix
		!only for symmetric slab atm. calculation based on and polarization(trans) 
		subroutine tmm_calc(q1,q2,q3,d,epsls,trans,M12,M23,M21,M32,verbose)
		
			implicit none
			real(dp),intent(in) :: q1,q2,q3,d,epsls(3)
			integer,intent(in) :: trans,verbose
			complex(dp),intent(out) :: M12(2,2),M23(2,2),M21(2,2),M32(2,2)
			
			!variable for inverse method
			
			integer ::INFO,LWORK
			integer ::IPIV(2)
			complex(dp),allocatable ::WORK(:)
			
			!TE polarization
			if (trans==0) then
				
				M12(1,:)=(0.5d0/q2)*(/(q1+q2)*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(-q1+q2)*exp(cmplx(0.d0,q2,dp)*d/2.d0)/)
				M12(2,:)=(0.5d0/q2)*(/(-q1+q2)*exp(-cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q1+q2)*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
									
				M23(1,:)=(0.5d0/q3)*(/(q3+q2)*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q3-q2)*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
				M23(2,:)=(0.5d0/q3)*(/(q3-q2)*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q3+q2)*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
									

			
			!TM polarization	
			else
				
				M12(1,:)=(0.5d0/(q2/epsls(2)))*(/(q1/epsls(1)+q2/epsls(2))*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(-q1/epsls(1)+q2/epsls(2))*exp(cmplx(0.d0,q2,dp)*d/2.d0)/)
				M12(2,:)=(0.5d0/(q2/epsls(2)))*(/(-q1/epsls(1)+q2/epsls(2))*exp(-cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q1/epsls(1)+q2/epsls(2))*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
									
				M23(1,:)=(0.5d0/(q3/epsls(3)))*(/(q3/epsls(3)+q2/epsls(2))*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q3/epsls(3)-q2/epsls(2))*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
				M23(2,:)=(0.5d0/(q3/epsls(3)))*(/(q3/epsls(3)-q2/epsls(2))*exp(cmplx(0.d0,q2,dp)*d/2.d0),&
									&(q3/epsls(3)+q2/epsls(2))*exp(-cmplx(0.d0,q2,dp)*d/2.d0)/)
				
			end if
			
			
			!inversing transfer matrix
			M21=M12
			M32=M23
			!inversing M12
			call zgetrf(2,2,M21,2,IPIV,INFO)
			!TODO: test+verbose part. check query LWORK=-1
			LWORK=-1
			allocate(WORK(1))
			call zgetri(2,M21,2,IPIV,WORK,LWORK,INFO)
			LWORK=nint(real(WORK(1)))
			deallocate (WORK)
			allocate (WORK(LWORK))
			call zgetri(2,M21,2,IPIV,WORK,LWORK,INFO)
			deallocate(WORK)
			
			!inversing M23
			call zgetrf(2,2,M32,2,IPIV,INFO)
			!TODO: test+verbose part. check query LWORK=-1
			LWORK=-1
			allocate(WORK(1))
			call zgetri(2,M32,2,IPIV,WORK,LWORK,INFO)
			LWORK=nint(real(WORK(1)))
			deallocate (WORK)
			allocate (WORK(LWORK))
			call zgetri(2,M32,2,IPIV,WORK,LWORK,INFO)
			deallocate(WORK)

			
		end subroutine tmm_calc
		
			
		subroutine zdot2(M,x,y)
			
			implicit none
			complex(dp),intent(in) :: M(2,2),x(2)
			complex(dp),intent(out) ::y(2)
			
			y(1)=M(1,1)*x(1)+M(1,2)*x(2)
			y(2)=M(2,1)*x(1)+M(2,2)*x(2)
		end subroutine zdot2

		! Construct the hamiltonian for the photonic crystal membrane
		subroutine h_guided(kcount,maxguide,dimH,kmod,d,refl,epsls,etamap,Hindex,Tarr,earr,ABarr,karr,H,flag)
		
			implicit none
			
			integer,intent(in) :: kcount,maxguide,refl,dimH,flag
			real(dp),intent(in) :: kmod(3,kcount),d,epsls(3)
			complex(dp),intent(in) :: etamap(kcount,kcount)
			!be careful !
			integer,intent(in) :: Hindex(2,dimH),Tarr(maxguide,kcount)
			complex(dp),intent(in) :: ABarr(4*maxguide,kcount)
			real(dp),intent(in) :: earr(2,kcount),karr(maxguide,kcount)
			complex(dp),allocatable,intent(out) :: H(:,:)
			
			!local variable
			integer ::i,j,uvi,uvj,ordi,ordj,tt
			complex(dp) ::eta2,eta1,eta3
			
			!for calculation
			real(dp) :: doteeb,dotggb,dotegb,dotgeb
			real(dp) ::I1,I2m,I2f,I3,ki,kj,giabs,gjabs
			real(dp) ::qi,qj,khi1i,khi1j,khi3i,khi3j
			complex(dp) ::A2(2),A3(2),B1(2),B2(2),tmpval
			
			
			

			
			allocate (H(dimH,dimH))

			do i=1,dimH

				do j=i,dimH
				
					uvi=Hindex(1,i)
					ordi=Hindex(2,i)
					uvj=Hindex(1,j)
					ordj=Hindex(2,j)
					ki=karr(Hindex(2,i),Hindex(1,i))
					kj=karr(Hindex(2,j),Hindex(1,j))
					giabs=kmod(3,Hindex(1,i))
					gjabs=kmod(3,Hindex(1,j))
					A2=(/ABarr(4*(ordi-1)+1,Hindex(1,i)),ABarr(4*(ordj-1)+1,Hindex(1,j))/)
					A3=(/ABarr(4*(ordi-1)+2,Hindex(1,i)),ABarr(4*(ordj-1)+2,Hindex(1,j))/)
					B1=(/ABarr(4*(ordi-1)+3,Hindex(1,i)),ABarr(4*(ordj-1)+3,Hindex(1,j))/)
					B2=(/ABarr(4*(ordi-1)+4,Hindex(1,i)),ABarr(4*(ordj-1)+4,Hindex(1,j))/)
					qi=sqrt(epsls(2)*(ki**2.d0)-(giabs**2.d0))
					qj=sqrt(epsls(2)*(kj**2.d0)-(gjabs**2.d0))
					khi1i=sqrt((giabs**2.d0)-epsls(1)*(ki**2.d0))
					khi1j=sqrt((gjabs**2.d0)-epsls(1)*(kj**2.d0))
					khi3i=sqrt((giabs**2.d0)-epsls(3)*(ki**2.d0))
					khi3j=sqrt((gjabs**2.d0)-epsls(3)*(kj**2.d0))
					eta2=etamap(uvi,uvj)
					
					if (uvi==uvj) then
						eta1=cmplx(1.d0,0.d0,dp)/epsls(1)
						!eta1=cmplx(13.d0,0.d0,dp)
						!eta3=cmplx(25.d0,0.d0,dp)
						eta3=cmplx(1.d0,0.d0,dp)/epsls(3)
					else
						eta1=cmplx(0.d0,0.d0,dp)
						eta3=cmplx(0.d0,0.d0,dp)
					end if
					
					
					
					call guided_integ(ki,kj,epsls,giabs,gjabs,d,I1,I2m,I2f,I3)
						
					!case TE-TE
					if ((Tarr(ordi,uvi)==0) .and.(Tarr(ordj,uvj)==0)) then
						doteeb=earr(1,uvi)*earr(1,uvj)+earr(2,uvi)*earr(2,uvj)
						tmpval=(ki**2.d0)*(kj**2.d0)*doteeb*(&
									&(epsls(1)**2.d0)*eta1*dconjg(B1(1))*B1(2)*I1+&
									&(epsls(3)**2.d0)*eta3*dconjg(A3(1))*A3(2)*I3+&
									&(epsls(2)**2.d0)*eta2*((dconjg(A2(1))*A2(2)+&
									&dconjg(B2(1))*B2(2))*I2m+(dconjg(A2(1))*B2(2)+&
									&dconjg(B2(1))*A2(2))*I2f))
						
						H(i,j)=tmpval
						

						if (i .ne. j) then
							H(j,i)=dconjg(tmpval)
						end if

					!case TM-TM
					else if ((Tarr(ordi,uvi)==1) .and.(Tarr(ordj,uvj)==1)) then
						dotggb=kmod(1,uvi)*kmod(1,uvj)+kmod(2,uvi)*kmod(2,uvj)
						tmpval=eta1*dconjg(B1(1))*B1(2)*(khi1i*khi1j*dotggb+giabs*gjabs)*I1+&
									&eta3*dconjg(A3(1))*A3(2)*(khi3i*khi3j*dotggb+giabs*gjabs)*I3+&
									&eta2*((dconjg(A2(1))*A2(2)+dconjg(B2(1))*B2(2))*(qi*qj*dotggb+&
									&giabs*gjabs)*I2m+(dconjg(A2(1))*B2(2)+dconjg(B2(1))*A2(2))*(&
									&-qi*qj*dotggb+giabs*gjabs)*I2f)
						H(i,j)=tmpval
						if (i .ne.j) then
							H(j,i)=dconjg(tmpval)
						end if

						
					!case TE-TM
					else if ((Tarr(ordi,uvi)==0) .and.(Tarr(ordj,uvj)==1)) then
						dotegb=earr(1,uvi)*kmod(1,uvj)+earr(2,uvi)*kmod(2,uvj)
						tmpval=(ki**2.d0)*dotegb*(-epsls(1)*eta1*dconjg(B1(1))*B1(2)*&
									&khi1j*I1+epsls(3)*eta3*dconjg(A3(1))*A3(2)*khi3j*I3+&
									&cmplx(0.d0,1.d0,dp)*epsls(2)*eta2*qj*((-dconjg(A2(1))*A2(2)+&
									&dconjg(B2(1))*B2(2))*I2m+(dconjg(A2(1))*B2(2)-&
									&dconjg(B2(1))*A2(2))*I2f))
									
						H(i,j)=tmpval
						H(j,i)=dconjg(tmpval)
						
					!case TM-TE
					else if ((Tarr(ordi,uvi)==1) .and.(Tarr(ordj,uvj)==0)) then
						dotgeb=kmod(1,uvi)*earr(1,uvj)+kmod(2,uvi)*earr(2,uvj)
						tmpval=(kj**2.d0)*dotgeb*(-epsls(1)*eta1*dconjg(B1(1))*B1(2)*&
									&khi1i*I1+epsls(3)*eta3*dconjg(A3(1))*A3(2)*khi3i*I3-&
									&cmplx(0.d0,1.d0,dp)*epsls(2)*eta2*qi*((-dconjg(A2(1))*A2(2)+&
									&dconjg(B2(1))*B2(2))*I2m+(dconjg(B2(1))*A2(2)-&
									&dconjg(A2(1))*B2(2))*I2f))
						
						
						H(i,j)=tmpval
						H(j,i)=dconjg(tmpval)

					end if
					


				end do
			end do		
		end subroutine h_guided
		


		!Calculation of loss in the photonic mode
		subroutine loss_calc(kcount,maxguide,dimH,d,epsls,kvect,&
						&eigrang,Hindex,Tarr,kmod,earr,ABarr,karr,Hmat,etamap,eigval,lossval,verbose,refl)
		
			implicit none
			integer,intent(in) :: eigrang(2),kcount,maxguide,dimH,verbose
			integer,intent(in) ::Hindex(2,dimH),Tarr(maxguide,kcount),refl
			real(dp),intent(in) :: kmod(3,kcount),earr(2,kcount)
			real(dp),intent(in) :: karr(maxguide,kcount),eigval(dimH),d,kvect(2),epsls(3)
			complex(dp),intent(in) :: Hmat(dimH,dimH),etamap(kcount,kcount),ABarr(4*maxguide,kcount)
			real(dp),intent(out) ::lossval(dimH)
			
			!local variable
			integer ::i,j,t,zz,Gid,Gjd,Ordid,pola
			real(dp) ::g00,llcrit,q1,q2,q3,ktmp,coeff
			real(dp) ::totalsum,locdos1,locdos3,locg
			complex(dp) ::WXTE1(2),WXTE2(2),WXTE3(2),YZTM1(2),YZTM2(2),YZTM3(2)
			complex(dp) :: M12TE(2,2),M23TE(2,2),M21TE(2,2),M32TE(2,2),MtotTE(2,2)
			complex(dp) :: M12TM(2,2),M23TM(2,2),M21TM(2,2),M32TM(2,2),MtotTM(2,2)
			real(dp) :: doteeb,dotggb,dotegb,dotgeb,ggb
			real(dp) :: qu,khi1u,khi3u,I2m,I2
			complex(dp) :: I1m,I1,I3m,I3,hval,A2,A3,B1,B2,eta2,eta1,eta3,locsumTE,locsumTM
			complex(dp) ::W1init,Y1init
			

			do i=1,dimH
				
				g00=sqrt(kvect(1)**2.d0+kvect(2)**2.d0)
				llcrit=g00/sqrt(epsls(1))
				
				!guided mode without any radiative coupling
				!180416: adding condition for eigrang limit of state
				if ((eigval(i) <= llcrit) .or. (i <eigrang(1)) .or. (i >eigrang(2)) ) then
					lossval(i)=9999.d0
				
				!radiative mode calculation
				else
					!initialize total loss sum
					totalsum=0.d0
					do j=1,kcount
						!only radiative modes above the lightlines are counted
						locg=kmod(3,j)
						llcrit=locg/sqrt(epsls(1))
						
						if ((eigval(i)> llcrit) .and. (isnan(kmod(1,j)).eqv. .false.) .and. (isnan(kmod(2,j)).eqv. .false.)) then
							
							!print *,"selected j",j
							q1=sqrt(epsls(1)*(eigval(i)**2.d0)-locg**2.d0)
							q2=sqrt(epsls(2)*(eigval(i)**2.d0)-locg**2.d0)
							q3=sqrt(epsls(3)*(eigval(i)**2.d0)-locg**2.d0)
							
							!contribution from each outgoing medium: 1 for lower slab, 2 for upper slab(3)
							call dos(epsls(1),eigval(i),locg,locdos1)
							call dos(epsls(3),eigval(i),locg,locdos3)
							!calculating transfer matrix for both TE and TM (1 frequency is capable to couple to both TE and TM)
							! Its implementation depends on which outgoing medium is being considered
							call tmm_calc(q1,q2,q3,d,epsls,0,M12TE,M23TE,M21TE,M32TE,verbose)
							call tmm_calc(q1,q2,q3,d,epsls,1,M12TM,M23TM,M21TM,M32TM,verbose)
							
							MtotTE=matmul(M23TE,M12TE)
							MtotTM=matmul(M23TM,M12TM)
							do t=1,1
								!initializing local sum for TE and TM radiative modes
								locsumTE=cmplx(0.d0,0.d0,dp)
								locsumTM=cmplx(0.d0,0.d0,dp)
								if (t==1) then
								
									!TEST for W and Y in TE and TM symmetric slab
									if (refl==1) then
										W1init=cmplx(1.d0/sqrt(2.d0*epsls(1)),0.d0,dp)*(1-MtotTE(1,2))/MtotTE(1,1)
										Y1init=-cmplx(1.d0/sqrt(2.d0),0.d0,dp)*(1+MtotTM(1,2))/MtotTM(1,1)
									else if (refl==-1) then
										W1init=-cmplx(1.d0/sqrt(2.d0*epsls(1)),0.d0,dp)*(1+MtotTE(1,2))/MtotTE(1,1)
										Y1init=cmplx(1.d0/sqrt(2.d0),0.d0,dp)*(1-MtotTM(1,2))/MtotTM(1,1)
									end if
									
									WXTE1=(/W1init,cmplx(1.d0/sqrt(2.d0*epsls(1)),0.d0,dp)/)
									call zdot2(M12TE,WXTE1,WXTE2)
									call zdot2(M23TE,WXTE2,WXTE3)
									YZTM1=(/Y1init,cmplx(1.d0/sqrt(2.d0),0.d0,dp)/)
									call zdot2(M12TM,YZTM1,YZTM2)
									call zdot2(M23TM,YZTM2,YZTM3)
									coeff=pi*locdos1

								else
									WXTE3=(/cmplx(1.d0/sqrt(epsls(3)),0.d0,dp),cmplx(0.d0,0.d0,dp)/)
									call zdot2(M32TE,WXTE3,WXTE2)
									call zdot2(M21TE,WXTE2,WXTE1)
									YZTM3=(/cmplx(1.d0,0.d0,dp),cmplx(0.d0,0.d0,dp)/)
									call zdot2(M32TM,YZTM3,YZTM2)
									call zdot2(M21TM,YZTM2,YZTM1)
									coeff=pi*locdos3
								end if
								
								do zz=1,dimH
									!calculating partial variable
									Gid=Hindex(1,zz)
									Ordid=Hindex(2,zz)
									Gjd=j
									ktmp=karr(Ordid,Gid)
									qu=sqrt(epsls(2)*(ktmp**2.d0)-kmod(3,Gid)**2.d0)
									khi1u=sqrt(-epsls(1)*(ktmp**2.d0)+kmod(3,Gid)**2.d0)
									khi3u=sqrt(-epsls(3)*(ktmp**2.d0)+kmod(3,Gid)**2.d0)
									I1m=cmplx(1.d0,0.d0,dp)/cmplx(khi1u,-q1,dp)
									I1=cmplx(1.d0,0.d0,dp)/cmplx(khi1u,q1,dp)
									I3m=cmplx(1.d0,0.d0,dp)/cmplx(khi3u,-q1,dp)
									I3=cmplx(1.d0,0.d0,dp)/cmplx(khi3u,q1,dp)
									
									if (abs(qu-q2)<=tolval) then
										I2m=d
									else	
										I2m=sin((qu-q2)*(d/2.d0))/((qu-q2)/2.d0)
									end if
									
									I2=sin((qu+q2)*(d/2.d0))/((qu+q2)/2.d0)
									eta2=etamap(Gid,Gjd)
									
									if (Gid==Gjd) then
										eta1=cmplx(1.d0,0.d0,dp)/epsls(1)
										eta3=cmplx(1.d0,0.d0,dp)/epsls(3)
									else
										eta1=cmplx(0.d0,0.d0,dp)/epsls(1)
										eta3=cmplx(0.d0,0.d0,dp)/epsls(3)
									end if
									
									!count TE and TM polarization
									do pola=0,1
										A2=ABarr(4*(ordid-1)+1,Gid)
										A3=ABarr(4*(ordid-1)+2,Gid)
										B1=ABarr(4*(ordid-1)+3,Gid)
										B2=ABarr(4*(ordid-1)+4,Gid)
										if ((Tarr(ordid,Gid)==0) .and. (pola==0)) then
											doteeb=earr(1,Gid)*earr(1,Gjd)+earr(2,Gid)*earr(2,Gjd)
											hval=(ktmp**2.d0)*eigval(i)*doteeb*(&
												&(epsls(1)**2.d0)*eta1*dconjg(B1)*(WXTE1(1)*I1+WXTE1(2)*I1m)+&
												&(epsls(3)**2.d0)*eta3*dconjg(A3)*(WXTE3(1)*I3m+WXTE3(2)*I3)+&
												&(epsls(2)**2.d0)*eta2*((dconjg(A2)*WXTE2(1)+dconjg(B2)*WXTE2(2))*I2m&
												&+(dconjg(A2)*WXTE2(2)+dconjg(B2)*WXTE2(1))*I2))
											
											locsumTE=locsumTE+dconjg(Hmat(zz,i))*hval
											
										else if ((Tarr(ordid,Gid)==1) .and. (pola==1)) then
											dotggb=kmod(1,Gid)*kmod(1,Gjd)+kmod(2,Gid)*kmod(2,Gjd)
											ggb=kmod(3,Gid)*kmod(3,Gjd)
											hval=eta1*dconjg(B1)*(cmplx(ggb,khi1u*q1*dotggb,dp)*YZTM1(1)*I1+&
												&cmplx(ggb,-khi1u*q1*dotggb,dp)*YZTM1(2)*I1m)+&
												&eta3*dconjg(A3)*(cmplx(ggb,-khi3u*q3*dotggb,dp)*YZTM3(1)*I3m+&
												&cmplx(ggb,khi3u*q3*dotggb,dp)*YZTM3(2)*I3)+eta2*(&
												&(dconjg(A2)*YZTM2(1)+dconjg(B2)*YZTM2(2))*&
												&(ggb+qu*q2*dotggb)*I2m +(dconjg(A2)*YZTM2(2)+dconjg(B2)*YZTM2(1))*&
												&(ggb-qu*q2*dotggb)*I2)
												
											locsumTM=locsumTM+dconjg(Hmat(zz,i))*hval
											
										
										else if ((Tarr(ordid,Gid)==0) .and. (pola==1)) then
											dotegb=earr(1,Gid)*kmod(1,Gjd)+earr(2,Gid)*kmod(2,Gjd)
											hval=cmplx(0.d0,ktmp**2.d0,dp)*dotegb*(epsls(1)*eta1*q1*&
												&dconjg(B1)*(-YZTM1(1)*I1+YZTM1(2)*I1m)+epsls(3)*eta3*q3*&
												&dconjg(A3)*(-YZTM3(1)*I3m+YZTM3(2)*I3)+epsls(2)*eta2*q2*&
												&((-dconjg(A2)*YZTM2(1)+dconjg(B2)*YZTM2(2))*I2m+&
												&(dconjg(A2)*YZTM2(2)-dconjg(B2)*YZTM2(1))*I2))

											locsumTM=locsumTM+dconjg(Hmat(zz,i))*hval
										
										else if ((Tarr(ordid,Gid)==1) .and. (pola==0)) then
											dotgeb=kmod(1,Gid)*earr(1,Gjd)+kmod(2,Gid)*earr(2,Gjd)
											hval=eigval(i)*dotgeb*(-epsls(1)*eta1*khi1u*dconjg(B1)*&
												&(WXTE1(1)*I1+WXTE1(2)*I1m)+epsls(3)*eta3*khi3u*dconjg(A3)*&
												&(WXTE3(1)*I3m+WXTE3(2)*I3)-cmplx(0.d0,epsls(2)*qu,dp)*eta2*&
												&((dconjg(B2)*WXTE2(2)-dconjg(A2)*WXTE2(1))*I2m+&
												&(dconjg(B2)*WXTE2(1)-dconjg(A2)*WXTE2(2))*I2))
											
											locsumTE=locsumTE+dconjg(Hmat(zz,i))*hval
										end if
										
										
									end do
									
								end do
								
								!add locsum contribution to imaginary part
								
								totalsum=totalsum+coeff*((abs(locsumTE)**2.d0)+(abs(locsumTM)**2.d0))

							end do
				
						end if

					end do
					!convert loss to Im (w/c)
					lossval(i)=totalsum/(2.d0*eigval(i))
					
				end if
			end do
			
		end subroutine loss_calc
		

		
		!list of kpts from fixed high sym pts
		subroutine listk(nbsympts,ksympts,nbdiv,outkls,nbpts)
			
			implicit none
			integer,intent(in) :: nbsympts
			real(dp),intent(in) :: ksympts(2,nbsympts)
			integer,intent(in) :: nbdiv
			real(dp),allocatable,intent(out) :: outkls(:,:)
			integer,intent(out) :: nbpts
			
			!local variable
			integer :: i,j,counter
			real(dp) ::dk(2),tmpk(2),startk(2)
			
			nbpts=nbdiv*(nbsympts-1)+1
			counter=0
			allocate(outkls(2,nbpts))
			
			do i=1,nbsympts-1
				startk=ksympts(:,i)
				dk=(ksympts(:,i+1)-ksympts(:,i))*(1.d0/(1.d0*nbdiv))
				
				do j=0,nbdiv-1
					counter=counter+1
					tmpk=startk+(1.d0*j)*dk
					outkls(:,counter)=tmpk
				end do
			end do
			
			counter=counter+1
			print *,"verifying counter and nbpts: ",counter,nbpts
			outkls(:,counter)=ksympts(:,nbsympts)
		
		end subroutine listk
		
		! Sampling the Brillouin zone
		subroutine bz_2d_sampling(ncap,per,g1,g2,nbdiv,outkls)
			
			implicit none
			real(dp),intent(in) ::g1(2),g2(2),ncap,per
			integer,intent(in) ::nbdiv(2)
			real(dp),allocatable,intent(out) :: outkls(:,:)
			
			!local variable
			integer ::ii,jj,counter
			real(dp) ::tmpk(2) !in cartesian coordinate
			!real(dp) ::normk !limfreq in norm freq unit
			
			allocate(outkls(2,nbdiv(1)*nbdiv(2)))
			
			counter=0
			do ii=1,nbdiv(1)
				do jj=1,nbdiv(2)
					
					counter=counter+1
					tmpk=(1.d0*(2*ii-nbdiv(1)-1)/(2.d0*nbdiv(1)))*g1+&
							&(1.d0*(2*jj-nbdiv(2)-1)/(2.d0*nbdiv(2)))*g2
					outkls(:,counter)=tmpk
				
				end do
			end do
			
		end subroutine bz_2d_sampling
		
end module mode_solver