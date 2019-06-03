
!This module allows user to extract the permittivity of a points in the photonic geometrical cell,
!based on the selected type of cell

module geometry_func
	
	
	use phys_cte
	
	implicit none
	
	!tolpos unit un um
	real(dp) :: tol=0.01
	real(dp) ::tolpos=0.002
	integer ::acellnb=8
	
	
	contains
	
	! This subroutine checks whether a points is located inside a polygone
	!INPUT:
	! pts(2): Coordinates (x,y) of the point to check
	! nbpoly: number of vertices of the polygon
	! polypts(2,nbpoly): List of polygon vertices's coordinates (x,y)
	!OUTPUT:
	! nbins: =1 if the points is inside the polygon, =0 otherwise
	subroutine polypts_check(pts,nbpoly,polypts,nbins)
	
		
		implicit none
		real(dp),intent(in) :: pts(2)
		integer,intent(in) :: nbpoly
		real(dp),intent(in) :: polypts(2,nbpoly)
		integer, intent(out) :: nbins
		
		!local var
		integer ::i,nbcutl,nbcutr,flagcheck,flagexit,remr,reml
		real(dp) ::polyb(2),polyf(2),cutpos(2),a,b
		real(dp) :: deltayb,deltayf,minybf,prod,deltax
		
		nbcutl=0
		nbcutr=0
		flagexit=0
		
		do i=1,nbpoly
			flagcheck=0
			polyb=polypts(:,i)
			if (i<nbpoly) then
				polyf=polypts(:,i+1)
			else
				polyf=polypts(:,1)
			end if
			
			cutpos(2)= pts(2)
			!normal case
			if ((abs(polyb(2)-polyf(2)) >tolpos) .and.(abs(polyb(1)-polyf(1)) >tolpos)) then
				
				a=(polyb(2)-polyf(2))/(polyb(1)-polyf(1))
				b=polyb(2)-a*polyb(1)
				cutpos(1)=(cutpos(2)-b)/a

				
			else if (abs(polyb(1)-polyf(1)) <= tolpos) then
				cutpos(1)= polyb(1)
			
			else if (abs(polyb(2)-polyf(2)) <= tolpos) then
				flagcheck=1
				
			end if
			
			!check condition
			deltax=cutpos(1)-pts(1)
			deltayb=cutpos(2)-polyb(2)
			deltayf=cutpos(2)-polyf(2)
			prod=deltayb*deltayf
			if (flagcheck==0) then
				if (abs(deltax)<=tolpos) then
					flagexit=1
					nbins=1
					exit
				else if (deltax >tolpos) then
					if ((abs(deltayb)<=tolpos) .and. (deltayf <deltayb)) then
						nbcutr=nbcutr+1
					else if ((abs(deltayf)<=tolpos) .and. (deltayb <deltayf)) then
						nbcutr=nbcutr+1
					else if ((abs(deltayf)>tolpos) .and. (abs(deltayb)>tolpos)&
									&.and. (prod <0.d0)) then
						nbcutr=nbcutr+1
					end if
					
						
						
				else if (deltax <-tolpos) then
					if ((abs(deltayb)<=tolpos) .and. (deltayf <deltayb)) then
						nbcutl=nbcutl+1
					else if ((abs(deltayf)<=tolpos) .and. (deltayb <deltayf)) then
						nbcutl=nbcutl+1
					else if ((abs(deltayf)>tolpos) .and. (abs(deltayb)>tolpos)&
									&.and. (prod <0.d0)) then
						nbcutl=nbcutl+1
					end if
				
				end if
			end if

		end do
		
		!check number of cutting points left and right
		if (flagexit .ne. 1) then
			remr=mod(nbcutr,2)
			reml=mod(nbcutl,2)
			if ((remr==0) .and. (reml==0)) then
				nbins =0
			else if ((remr == 1) .and. (reml == 1)) then
				nbins=1
			end if
		end if
		
	
	end subroutine polypts_check
	
	! 1D Bragg lattice (with finite membrane thickness):
	!INPUT
	! x,y: Point coordinates
	! inputls(period,wdiel,eps0,eps1,hdummy,d): Cell period, Filling factor of eps1,first permittivity, second permittivity, 
	!	cell height in the y direction, thickness of the cell in z-direction
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	
	subroutine bragg_latt(x,y,inputls,val)
	
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) :: val
		
		!local variable
		real(dp):: per,wdiel,d,eps0,eps1,hdummy
		real(dp):: ref(2) 
		
		per=inputls(1)
		wdiel=inputls(2)
		d=inputls(6)
		eps0=inputls(3)
		eps1=inputls(4)
		hdummy=inputls(5)
		
		val=eps0
		ref=(/(-per/2.d0),(-hdummy/2.d0)/)
		
		if ((x-ref(1)) <= wdiel) then
			val=eps1
		end if
		
	end subroutine bragg_latt
	
	
	! 2D photonic cell (with finite membrane thickness) with crosshatch pattern
	!INPUT
	! x,y: Point coordinates
	! inputls(period,h,angle,eps0,eps1,d,w): Cell period, length of each cross,angle between two cell basis vector,first permittivity, second permittivity, 
	!	thickness of the cell in z-direction, width of each cross
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	
	subroutine crosshatch_latt(x,y,inputls,val)
	
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) :: val
		
		!local variable
		real(dp) ::polyver(2,4),polyhor(2,4)
		real(dp) :: period,hval,angle,eps0,eps1,wval
		real(dp) :: centerver(2,5),centerhor(2,4)
		real(dp) :: dprime,outrad,tmpdist(2),dist,tmpx,tmpy
		integer :: i,nbins
		
		period=inputls(1)
		hval=inputls(2)
		angle=inputls(3)
		eps0=inputls(4)
		eps1=inputls(5)
		wval=inputls(7)
		
		dprime=period*s2/2.d0
		outrad=sqrt(hval*hval+wval*wval)/2.d0
		
		centerver(:,1)=(/0.d0,0.d0/)
		centerver(:,2)=(/period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		centerver(:,3)=(/period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		centerver(:,4)=(/-period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		centerver(:,5)=(/-period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		
		centerhor(:,1)=(/-dprime,0.d0/)
		centerhor(:,2)=(/0.d0,-dprime/)
		centerhor(:,3)=(/dprime,0.d0/)
		centerhor(:,4)=(/0.d0,dprime/)
		
		polyver(:,1)=(/-wval/2.d0,-hval/2.d0/)
		polyver(:,2)=(/wval/2.d0,-hval/2.d0/)
		polyver(:,3)=(/wval/2.d0,hval/2.d0/)
		polyver(:,4)=(/-wval/2.d0,hval/2.d0/)
		
		polyhor(:,1)=centerhor(:,1)+(/-hval/2.d0,-wval/2.d0/)
		polyhor(:,2)=centerhor(:,1)+(/hval/2.d0,-wval/2.d0/)
		polyhor(:,3)=centerhor(:,1)+(/hval/2.d0,wval/2.d0/)
		polyhor(:,4)=centerhor(:,1)+(/-hval/2.d0,wval/2.d0/)
		
		val=eps1
		
		!vertical box scan
		do i=1,5
			tmpdist=(/x-centerver(1,i),y-centerver(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			
			!check algorithm for points inside the outer circle
			if (dist <= outrad) then
			
				!translation
				tmpx= x-centerver(1,i)
				tmpy= y-centerver(2,i)
				call polypts_check((/tmpx,tmpy/),4,polyver,nbins)
				if (nbins==1) then 
					val=eps0
				end if
				exit
			end if
		end do
		
		!horizontal box scan
		do i=1,4
			tmpdist=(/x-centerhor(1,i),y-centerhor(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			
			!check algorithm for points inside the outer circle
			if (dist <= outrad) then
			
				!translation
				tmpx= x-(centerhor(1,i)-centerhor(1,1))
				tmpy= y-(centerhor(2,i)-centerhor(2,1))
				call polypts_check((/tmpx,tmpy/),4,polyhor,nbins)
				if (nbins==1) then 
					val=eps0
				end if
				exit
			end if
		end do
		

		
	end subroutine crosshatch_latt
	

	! 2D photonic cell (with finite membrane thickness) with equilateral triangle pattern
	!INPUT
	! x,y: Point coordinates
	! inputls(period,L,angle,eps0,eps1,d): Cell period, triangle edge width,angle between two cell basis vector,first permittivity, second permittivity, 
	!	thickness of the cell in z-direction
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	subroutine triang_gen_latt(x,y,inputls,val)
		
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) :: val
		
		!local variable
		
		real(dp) :: center(2,5),period,L,angle,eps0,eps1
		real(dp) :: polypts(2,3),tmpdist(2),dist,outrad,tmpx,tmpy
		integer :: nbins,i
		
		
		period=inputls(1)
		L=inputls(2)
		angle=inputls(3)
		eps0=inputls(4)
		eps1=inputls(5)
		
		center(:,1)=(/0.d0,0.d0/)
		center(:,2)=(/-period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		center(:,3)=(/period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		center(:,4)=(/period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		center(:,5)=(/-period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		

		
		polypts(:,1)=(/L*s3/3.d0,0.d0/)
		polypts(:,2)=(/(-L*s3/6.d0),(L/2.d0)/)
		polypts(:,3)=(/(-L*s3/6.d0),(-L/2.d0)/)
		!print *,"s3",s3
		outrad= L*s3/3.d0
		
		val=eps1
		
		do i=1,5
			tmpdist=(/x-center(1,i),y-center(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			
			!check algorithm for points inside the outer circle
			if (dist <= outrad) then
			
				!translation
				tmpx= x-center(1,i)
				tmpy= y-center(2,i)
				
				call polypts_check((/tmpx,tmpy/),3,polypts,nbins)
				if (nbins==1) then 
					val=eps0
				end if
				exit
			end if
		end do

	end subroutine triang_gen_latt

	! 2D photonic cell (with finite membrane thickness) with circle pattern
	!INPUT
	! x,y: Point coordinates
	! inputls(period,radius,angle,eps0,eps1,d): Cell period, circle radius,angle between two cell basis vector,first permittivity, second permittivity, 
	!	thickness of the cell in z-direction
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	subroutine gen_latt(x,y,inputls,val)
		
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) :: val
		
		!local variable
		
		real(dp) :: center(2,5),tmpdist(2),dist
		real(dp) :: period,radius,angle,eps0,eps1
		integer :: i
		
		period=inputls(1)
		radius=inputls(2)
		angle=inputls(3)
		eps0=inputls(4)
		eps1=inputls(5)
		
		
		
		center(:,1)=(/0.d0,0.d0/)
		center(:,2)=(/-period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		center(:,3)=(/period*cos(angle/2.d0),period*sin(angle/2.d0)/)
		center(:,4)=(/period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		center(:,5)=(/-period*cos(angle/2.d0),-period*sin(angle/2.d0)/)
		
		!print *,period,radius,angle,epsi0,epsi1
		val=eps1

		do i=1,5
			
			tmpdist=(/x-center(1,i),y-center(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			if (dist<=radius) then

				val=eps0

			exit
			end if
		end do
		
	end subroutine gen_latt
	
	
	
	! 2D photonic cell (with finite membrane thickness) with hybrid hexagonal-graphite pattern (2 holes per unit cell). angle=60 deg
	!INPUT
	! x,y: Point coordinates
	! inputls(period,rhex,rgraph,eps0,eps1,d): Cell period, hexagonal radius,graphite radius,angle between two cell basis vector,
	!first permittivity, second permittivity, thickness of the cell in z-direction
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	subroutine hexgraph_latt(x,y,inputls,val)
		
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) ::val
		
		!local variable
		real(dp) :: rhex,rgraph,period,eps0,eps1,centerh(2,5),centerg(2,6)
		real(dp) :: tmpdist(2),dist
		integer ::flag,i
		
		period=inputls(1)
		rhex=inputls(2)
		rgraph=inputls(3)
		eps0=inputls(4)
		eps1=inputls(5)
		
		centerh(:,1)=(/0.d0,0.d0/)
		centerh(:,2)=(/-period*cos(pi/6.d0),period*sin(pi/6.d0)/)
		centerh(:,3)=(/period*cos(pi/6.d0),period*sin(pi/6.d0)/)
		centerh(:,4)=(/period*cos(pi/6.d0),-period*sin(pi/6.d0)/)
		centerh(:,5)=(/-period*cos(pi/6.d0),-period*sin(pi/6.d0)/)
		
		centerg(:,1)=(/-period*s3/3.d0,0.d0/)
		centerg(:,2)=(/-period*s3/6.d0,period*sin(pi/6.d0)/)
		centerg(:,3)=(/period*s3/6.d0,period*sin(pi/6.d0)/)
		centerg(:,4)=(/period*s3/3.d0,0.d0/)
		centerg(:,5)=(/period*s3/6.d0,-period*sin(pi/6.d0)/)
		centerg(:,6)=(/-period*s3/6.d0,-period*sin(pi/6.d0)/)
		
		flag=1
		val=eps1
		
		do i=1,5
			
			tmpdist=(/x-centerh(1,i),y-centerh(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			if (dist<=rhex) then

				val=eps0
				flag=0
				exit
			end if
		end do
		
		if (flag==1) then
			do i=1,6
			
			tmpdist=(/x-centerg(1,i),y-centerg(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			if (dist<=rgraph) then

				val=eps0
				flag=0
				exit
			end if
			end do
		end if
	
	end subroutine hexgraph_latt
	
	! 2D photonic cell (with finite membrane thickness) with kagome pattern (3 holes per unit cell). angle=60 deg
	!INPUT
	! x,y: Point coordinates
	! inputls(period,a,rkag,eps0,eps1,d): Cell period, distance between each hole in the unit cell,hole radius,
	! first permittivity, second permittivity, thickness of the cell in z-direction
	!OUTPUT:
	! val: eps0 or eps1 depends on the position of the point
	subroutine kagome_latt(x,y,inputls,val)
		
		implicit none
		real(dp),intent(in) :: x,y,inputls(*)
		real(dp),intent(inout) ::val
		
		!local variable
		real(dp) :: period,a,rkag,eps0,eps1,centerk(2,9)
		real(dp) :: tmpdist(2),dist,dd
		integer ::flag,i
		
		period=inputls(1)
		a=period/2.d0
		rkag=inputls(3)
		eps0=inputls(4)
		eps1=inputls(5)
		
		centerk(:,1)=(/0.d0,0.d0/)
		centerk(:,2)=(/-period*cos(pi/6.d0),period*sin(pi/6.d0)/)
		centerk(:,3)=(/period*cos(pi/6.d0),period*sin(pi/6.d0)/)
		centerk(:,4)=(/period*cos(pi/6.d0),-period*sin(pi/6.d0)/)
		centerk(:,5)=(/-period*cos(pi/6.d0),-period*sin(pi/6.d0)/)
		centerk(:,6)=(/-a*cos(pi/6.d0),a*sin(pi/6.d0)/)
		centerk(:,7)=(/a*cos(pi/6.d0),a*sin(pi/6.d0)/)
		centerk(:,8)=(/a*cos(pi/6.d0),-a*sin(pi/6.d0)/)
		centerk(:,9)=(/-a*cos(pi/6.d0),-a*sin(pi/6.d0)/)

		
		val=eps1
		
		do i=1,9
			
			tmpdist=(/x-centerk(1,i),y-centerk(2,i)/)
			dist=sqrt((tmpdist(1)**2.d0)+(tmpdist(2)**2.d0))
			if (dist<=rkag) then

				val=eps0
				exit
			end if
		end do
		
	end subroutine kagome_latt
	

	
	! Supercell for H1 cavity
	! H2 cavity will be adressed later since shifted hole index is not regular
	! inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
	! full generation of supercell instead of (x,y) single coordinates
	! map will be allocated in this subroutine
	!INPUT:
	! dimmap(2): map size (x and y)
	! inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb):period,hole radius, angle between two basis vector,first permittivity, second permittivity,
	! membrane thickness, radius of the first neighbor ring of hole to the cavity, the shift of each hole in the first neighbor ring from ideal position 
	!OUTPUT:
	! map(dimmap(1),dimmap(2)): Output dielectric map
	! scper: Supercell period
	subroutine h1_latt(dimmap,inpls,map,scper)
		
		implicit none
		
		real(dp),intent(out) ::scper
		integer,intent(in) :: dimmap(2)
		real(dp),intent(in) :: inpls(*)
		real(dp),intent(inout) :: map(dimmap(1),dimmap(2))
		!local variable
		real(dp) :: per,rad,alpha,eps0,eps1,d,ord,radnb,shiftnb
		real(dp) :: tmprad,tmpshift(2),normnb(2),a1(2),a2(2)
		integer ::nbhole,i,j,t,counter,countnb,nbindex(2,6)
		real(dp) :: stepx,stepy,mapori(2),xwidth,ywidth,tmpx,tmpy,dxy,tmpcen(2)
		real(dp),allocatable ::holels(:,:)
		
		
		
		!load parameter
		per=inpls(1)
		rad=inpls(2)
		alpha=inpls(3)
		eps0=inpls(4)
		eps1=inpls(5)
		d=inpls(6)
		radnb=inpls(7)
		shiftnb=inpls(8)
		!1 hole removed
		nbhole=(2*acellnb+1)*(2*acellnb+1)-1
		xwidth=2.d0*per*acellnb*cos(alpha/2.d0)
		ywidth=2.d0*per*acellnb*sin(alpha/2.d0)
		!output scper
		scper=acellnb*per

		stepx=xwidth/((dimmap(2)-1)*1.d0)
		stepy=ywidth/((dimmap(1)-1)*1.d0)
		mapori=(/(-xwidth/2.d0),(-ywidth/2.d0)/)
		nbindex(1,:)=(/1,0,1,-1,-1,0/)
		nbindex(2,:)=(/0,1,-1,1,0,-1/)
		
		a1=(/per*cos(alpha/2.d0),per*sin(alpha/2.d0)/)
		a2=(/per*cos(alpha/2.d0),-per*sin(alpha/2.d0)/)
		
		allocate (holels(3,nbhole))
		counter=0
		countnb=0
		!supercell holes indexing
		do i=-acellnb,acellnb
			do j=-acellnb,acellnb
				
				tmprad=rad
				tmpshift=(/0.d0,0.d0/)
				if (countnb <6) then
				
					do t=1,6
						if ((i==nbindex(1,t)).and.(j==nbindex(2,t))) then
							
							normnb=i*1.d0*a1+j*1.d0*a2
							tmprad=radnb
							tmpshift=normnb*(shiftnb/sqrt(normnb(1)*normnb(1)+&
												&normnb(2)*normnb(2)))
							countnb=countnb+1
							exit
						end if
					end do
					
				end if
				if ((i .ne. 0).or.(j .ne. 0)) then
					counter=counter+1
					tmpcen=1.d0*i*a1+1.d0*j*a2+tmpshift
					holels(1,counter)=tmpcen(1)
					holels(2,counter)=tmpcen(2)
					holels(3,counter)=tmprad
				end if

			end do
		end do
		!end of supercell holes indexing
		
		
		do i=1,dimmap(2)
			do j=1,dimmap(1)
				
				map(j,i)=eps1
				tmpx=mapori(1)+1.d0*(i-1)*stepx
				tmpy=mapori(2)+1.d0*(j-1)*stepy
				do t=1,nbhole
					dxy=sqrt(((tmpx-holels(1,t))**2.d0)+((tmpy-holels(2,t))**2.d0))
					if (dxy <= holels(3,t)) then
						map(j,i)=eps0
						exit
					end if
				end do
				
			end do
			
		end do
		

	end subroutine h1_latt
	
	! Supercell for H2 cavity
	!INPUT:
	! dimmap(2): map size (x and y)
	! inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb):period,hole radius, angle between two basis vector,first permittivity, second permittivity,
	! membrane thickness, radius of the first neighbor ring of hole to the cavity, the shift of each hole in the first neighbor ring from ideal position 
	!OUTPUT:
	! map(dimmap(1),dimmap(2)): Output dielectric map
	! scper: Supercell period
	subroutine h2_latt(dimmap,inpls,map,scper)
		
		implicit none
		
		real(dp),intent(out) ::scper
		integer,intent(in) :: dimmap(2)
		real(dp),intent(in) :: inpls(*)
		real(dp),intent(inout) :: map(dimmap(1),dimmap(2))
		
		!local variable
		real(dp) :: per,rad,alpha,eps0,eps1,d,ord,radnb,shiftnb
		real(dp) :: tmprad,tmpshift(2),normnb(2),a1(2),a2(2)
		integer ::nbhole,i,j,t,qq,counter,countnb,nbindex(2,12),r1index(2,7),shiftid(6)
		real(dp) :: stepx,stepy,mapori(2),xwidth,ywidth,tmpx,tmpy,dxy
		real(dp),allocatable ::holels(:,:)
		integer :: flag
		
		
		!load parameter
		per=inpls(1)
		rad=inpls(2)
		alpha=inpls(3)
		eps0=inpls(4)
		eps1=inpls(5)
		d=inpls(6)
		radnb=inpls(7)
		shiftnb=inpls(8)
		!1 hole removed
		nbhole=(2*acellnb+1)*(2*acellnb+1)-7
		xwidth=2.d0*per*acellnb*cos(alpha/2.d0)
		ywidth=2.d0*per*acellnb*sin(alpha/2.d0)
		!output scper
		scper=acellnb*per

		stepx=xwidth/((dimmap(2)-1)*1.d0)
		stepy=ywidth/((dimmap(1)-1)*1.d0)
		mapori=(/(-xwidth/2.d0),(-ywidth/2.d0)/)
		nbindex(1,:)=(/2,0,1,-2,0,-1,2,2,1,-2,-1,-2/)
		nbindex(2,:)=(/0,2,1,0,-2,-1,-1,-2,-2,2,2,1/)
		r1index(1,:)=(/0,1,0,1,-1,-1,0/)
		r1index(2,:)=(/0,0,1,-1,1,0,-1/)
		shiftid=(/3,6,7,9,11,12/)
		a1=(/per*cos(alpha/2.d0),per*sin(alpha/2.d0)/)
		a2=(/per*cos(alpha/2.d0),-per*sin(alpha/2.d0)/)
		
		allocate (holels(3,nbhole))
		counter=0
		countnb=0
		!supercell holes indexing
		do i=-acellnb,acellnb
			do j=-acellnb,acellnb
				flag=0
				tmprad=rad
				tmpshift=(/0.d0,0.d0/)
				if (countnb <12) then
				
					do t=1,12
						if ((i==nbindex(1,t)).and.(j==nbindex(2,t))) then
							print *,""
							tmprad=radnb
							
							!check shift position
							do qq=1,6
								if (t==shiftid(qq)) then
									normnb=i*1.d0*a1+j*1.d0*a2
									tmpshift=normnb*(shiftnb/sqrt(normnb(1)*normnb(1)+&
														&normnb(2)*normnb(2)))
									exit
								end if
							end do
							countnb=countnb+1
							exit
						end if
					end do
					
				end if
				
				do qq=1,7
					if ((i == r1index(1,qq)).and.(j == r1index(2,qq))) then
						flag=1
						exit
					end if
				end do
				
				if (flag==0) then
						counter=counter+1
						holels(1:2,counter)=1.d0*i*a1+1.d0*j*a2+tmpshift
						holels(3,counter)=tmprad
				end if

			end do
		end do
		!end of supercell holes indexing
		
		
		do i=1,dimmap(2)
			do j=1,dimmap(1)
				
				map(j,i)=eps1
				tmpx=mapori(1)+1.d0*(i-1)*stepx
				tmpy=mapori(2)+1.d0*(j-1)*stepy
				do t=1,nbhole
					dxy=sqrt(((tmpx-holels(1,t))**2.d0)+((tmpy-holels(2,t))**2.d0))
					if (dxy <= holels(3,t)) then
						map(j,i)=eps0
						exit
					end if
				end do
				
			end do
		end do
		
		deallocate (holels)
	
	end subroutine h2_latt

	! Supercell for L3 cavity
	!INPUT:
	! dimmap(2): map size (x and y)
	! inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb):period,hole radius, angle between two basis vector,first permittivity, second permittivity,
	! membrane thickness, radius of the two holes at the begining and the end of the cavity, 
	! the shift of these two holes from ideal position 
	!OUTPUT:
	! map(dimmap(1),dimmap(2)): Output dielectric map
	! scper: Supercell period
	subroutine l3_latt(dimmap,inpls,map,scper)
	
		implicit none
		
		real(dp),intent(out) ::scper
		integer,intent(in) :: dimmap(2)
		real(dp),intent(in) :: inpls(*)
		real(dp),intent(inout) :: map(dimmap(1),dimmap(2))
		
		!local variable
		real(dp) :: per,rad,alpha,eps0,eps1,d,ord,radnb,shiftnb
		real(dp) :: tmprad,tmpshift(2),normnb(2),a1(2),a2(2)
		integer ::nbhole,i,j,t,qq,counter,countnb,nbindex(2,2),r1index(2,3),shiftid(2)
		real(dp) :: stepx,stepy,mapori(2),xwidth,ywidth,tmpx,tmpy,dxy
		real(dp),allocatable ::holels(:,:)
		integer :: flag
		
		
		!load parameter
		per=inpls(1)
		rad=inpls(2)
		alpha=inpls(3)
		eps0=inpls(4)
		eps1=inpls(5)
		d=inpls(6)
		radnb=inpls(7)
		shiftnb=inpls(8)
		!1 hole removed
		nbhole=(2*acellnb+1)*(2*acellnb+1)-3
		xwidth=2.d0*per*acellnb*cos(alpha/2.d0)
		ywidth=2.d0*per*acellnb*sin(alpha/2.d0)
		!output scper
		scper=acellnb*per

		stepx=xwidth/((dimmap(2)-1)*1.d0)
		stepy=ywidth/((dimmap(1)-1)*1.d0)
		mapori=(/(-xwidth/2.d0),(-ywidth/2.d0)/)
		nbindex(1,:)=(/2,-2/)
		nbindex(2,:)=(/0,0/)
		r1index(1,:)=(/0,1,-1/)
		r1index(2,:)=(/0,0,0/)
		shiftid=(/1,2/)
		a1=(/per*cos(alpha/2.d0),per*sin(alpha/2.d0)/)
		a2=(/per*cos(alpha/2.d0),-per*sin(alpha/2.d0)/)
		
		allocate (holels(3,nbhole))
		counter=0
		countnb=0
		!supercell holes indexing
		do i=-acellnb,acellnb
			do j=-acellnb,acellnb
				flag=0
				tmprad=rad
				tmpshift=(/0.d0,0.d0/)
				if (countnb <2) then
				
					do t=1,2
						if ((i==nbindex(1,t)).and.(j==nbindex(2,t))) then
							print *,""
							tmprad=radnb
							!check shift position
							do qq=1,2
								if (t==shiftid(qq)) then
									normnb=i*1.d0*a1+j*1.d0*a2
									tmpshift=normnb*(shiftnb/sqrt(normnb(1)*normnb(1)+&
														&normnb(2)*normnb(2)))
									exit
								end if
							end do
							countnb=countnb+1
							exit
						end if
					end do
					
				end if
				
				do qq=1,3
					if ((i == r1index(1,qq)).and.(j == r1index(2,qq))) then
						flag=1
						exit
					end if
				end do
				
				if (flag==0) then
						counter=counter+1
						holels(1:2,counter)=1.d0*i*a1+1.d0*j*a2+tmpshift
						holels(3,counter)=tmprad
				end if

			end do
		end do
		!end of supercell holes indexing
		
		
		!scanning map (x,y) (dimmap(2),dimmap(1))

		do i=1,dimmap(2)
			do j=1,dimmap(1)
				
				map(j,i)=eps1
				tmpx=mapori(1)+1.d0*(i-1)*stepx
				tmpy=mapori(2)+1.d0*(j-1)*stepy
				do t=1,nbhole
					dxy=sqrt(((tmpx-holels(1,t))**2.d0)+((tmpy-holels(2,t))**2.d0))
					if (dxy <= holels(3,t)) then
						map(j,i)=eps0
						exit
					end if
				end do
				
			end do
		end do
		
		deallocate (holels)
		
	end subroutine l3_latt		
			
end module geometry_func
