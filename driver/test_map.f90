program test_map

	use phys_cte
	use geometry_func
	use geometry_defs
	
	implicit none
	
	real(dp) :: per,phi,inpls(6),center(2)
	integer :: divxy(2),mdim(2)
	real(dp),allocatable :: map(:,:)
	real(dp) :: alatt(2,2),glatt(2,2),asc(2,2),glattxy(2,2)
	complex(dp),allocatable :: fftmap(:,:)
	character(len=120) :: output,outfft
	
	!inpls(per,rad,angle,eps0,eps1,d,radnb,shiftnb)
	!inpls=(/1.d0,0.3d0,pi/3.d0,1.d0,12.d0,0.6d0,0.2d0,0.2d0/)
	!inputls(period,wdiel,eps0,eps1,hdummy,d)
	inpls=(/1.d0,0.3d0,1.d0,12.d0,2.d0,0.6d0/)
	per=inpls(1)
	divxy=(/127,127/)
	mdim=(/divxy(2)+1,divxy(1)+1/)
	output="test_Bragg"
	outfft="bragg_fft"
	center=(/0.d0,0.d0/)
	
	call unit_cell1D(per,center,inpls,bragg_latt,divxy,map,alatt,asc,glatt,glattxy)
	call map_fft(map,mdim,fftmap)
	call dwrite_map_tofile(mdim,map,output)
	call zwrite_map_tofile(mdim,fftmap,outfft)






end program test_map
