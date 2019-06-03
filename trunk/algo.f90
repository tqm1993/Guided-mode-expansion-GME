!This module contains heapsort algorithm

module algo
	
	use phys_cte
	
	contains 
	
	subroutine dsift_down(n,narr,idarr,arr)
		
		implicit none
		integer,intent(in) ::n,narr
		integer,intent(inout) ::idarr(narr)
		real(dp),intent(inout) :: arr(narr)
		
		!local variable
		integer ::nsift,ii,idleft,idright,idchild
		real(dp) ::dummy
		integer ::iddummy
		
		nsift=n/2
		
		do ii=nsift,1,-1
			
			idleft=ii*2
			idright=ii*2+1
			
			if (idright <(n+1)) then
				if (arr(idright)>arr(idleft)) then
					idchild=idright
				else
					idchild=idleft
				end if
			
			else
				idchild=idleft
			end if
			
			if (arr(idchild)>arr(ii)) then
				dummy=arr(ii)
				arr(ii)=arr(idchild)
				arr(idchild)=dummy
				
				!swap index
				iddummy=idarr(ii)
				idarr(ii)=idarr(idchild)
				idarr(idchild)=iddummy
				
			end if

		end do
		

		
	end subroutine dsift_down
	
	!heapsort in increasing direction for real array
	!narr: size of the array
	!arr(narr): array to sort
	!idarr(narr): index array after using heapsort
	subroutine dheapsort(narr,arr,idarr)
		
		implicit none
		integer,intent(in) ::narr
		real(dp),intent(inout)::arr(narr)
		integer,intent(out) ::idarr(narr)
		
		!local variable
		integer ::ii
		real(dp) ::dummy
		integer ::iddummy
		
		do ii=1,narr
			idarr(ii)=ii
		end do
		
		do ii=narr,2,-1
			call dsift_down(ii,narr,idarr,arr)
			!swap first and last elements
			dummy=arr(1)
			iddummy=idarr(1)
			arr(1)=arr(ii)
			arr(ii)=dummy
			idarr(1)=idarr(ii)
			idarr(ii)=iddummy
		end do
	
	end subroutine dheapsort
	

	subroutine isift_down(n,narr,idarr,arr)
		
		implicit none
		integer,intent(in) ::n,narr
		integer,intent(inout) ::idarr(narr)
		integer,intent(inout) :: arr(narr)
		
		!local variable
		integer ::nsift,ii,idleft,idright,idchild
		integer ::dummy
		integer ::iddummy
		
		nsift=n/2
		
		do ii=nsift,1,-1
			
			idleft=ii*2
			idright=ii*2+1
			
			if (idright <(n+1)) then
				if (arr(idright)>arr(idleft)) then
					idchild=idright
				else
					idchild=idleft
				end if
			
			else
				idchild=idleft
			end if
			
			if (arr(idchild)>arr(ii)) then
				dummy=arr(ii)
				arr(ii)=arr(idchild)
				arr(idchild)=dummy
				
				!swap index
				iddummy=idarr(ii)
				idarr(ii)=idarr(idchild)
				idarr(idchild)=iddummy
			end if
			
		end do
		
	end subroutine isift_down
	
	!heapsort in increasing direction for integer array
	!narr: size of the array
	!arr(narr): array to sort
	!idarr(narr): index array after using heapsort
	subroutine iheapsort(narr,arr,idarr)
		
		implicit none
		integer,intent(in) ::narr
		integer,intent(inout)::arr(narr)
		integer,intent(out) ::idarr(narr)
		
		!local variable
		integer ::ii
		integer ::dummy
		integer ::iddummy
		
		do ii=1,narr
			idarr(ii)=ii
		end do
		
		do ii=narr,2,-1
			call isift_down(ii,narr,idarr,arr)
			!swap first and last elements
			dummy=arr(1)
			iddummy=idarr(1)
			arr(1)=arr(ii)
			arr(ii)=dummy
			idarr(1)=idarr(ii)
			idarr(ii)=iddummy
		end do
	
	end subroutine iheapsort


end module algo
