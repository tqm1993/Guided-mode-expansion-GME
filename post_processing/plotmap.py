

from numpy import *
from scipy import *
from extract_data import *
import numpy as np

import matplotlib.pyplot as plt


def plotmap (filename,period,phi,nbpts):

	x_ax=linspace(-period*cos(phi/2.0),period*cos(phi/2.0),num=nbpts[0])
	y_ax=linspace(-period*sin(phi/2.0),period*sin(phi/2.0),num=nbpts[1])
	qq=period*cos(phi/2.0)
	data=read_file(filename)
	
	#check this part
	x,y=meshgrid(x_ax,y_ax)
	#print y[0,:]
	#print data[0,:]
	fig=plt.figure()
	ax=fig.add_subplot(111)
	cp=ax.contourf(x,y,abs(data))
	plt.xlim([-qq,qq])
	plt.ylim([-qq,qq])
	plt.axes().set_aspect('equal')
	#plt.colorbar(cp)
	plt.xlabel(r'$\mu$m')
	plt.ylabel(r'$\mu$m')
	plt.show()
	
def plotmap1D (filename,period,hdummy,nbpts):

	x_ax=linspace(-period/2.0,period/2.0,num=nbpts[0])
	y_ax=linspace(-hdummy/2.0,hdummy/2.0,num=nbpts[1])
	qq=max(hdummy,period)
	data=read_file(filename)
	
	#check this part
	x,y=meshgrid(x_ax,y_ax)
	#print y[0,:]
	#print data[0,:]
	fig=plt.figure()
	ax=fig.add_subplot(111)
	cp=ax.contourf(x,y,abs(data))
	plt.xlim([-qq/2.0,qq/2.0])
	plt.ylim([-qq/2.0,qq/2.0])
	plt.axes().set_aspect('equal')
	#plt.colorbar(cp)
	plt.xlabel(r'$\mu$m')
	plt.ylabel(r'$\mu$m')
	plt.show()
	


#data=read_file(filename)

#fftdata=np.fft.fft2(data)
#databis=np.fft.ifft2(fftdata)

#print databis[1,:]
def plotfield (filename,flag=0):
	
	data=read_file(filename)
	nbrow=data.shape[0]
	nbcol=data.shape[1]
	
	deriv=np.empty((nbrow-1,nbcol))
	
	for i in range(0,nbrow-1):
		deriv[i,0]=data[i+1,0]
		deriv[i,1]=(data[i+1,1]-data[i,1])/(data[i+1,0]-data[i,0])
	
	fig=plt.figure()
	plt.hold(True)
	if flag==0:
		plt.plot(data[:,0],data[:,1],'b')
		plt.plot([data[0,0],data[-1,0]],[0.0,0.0],'r')
	else:
		plt.plot(deriv[:,0],deriv[:,1])
	
	plt.show()
	
def plotmulti (filename):
	
	fig=plt.figure()
	plt.hold(True)
	for i in range(0,len(filename)):
		data=read_file(filename[i])	
		plt.plot(data[:,0],data[:,1])

	
	plt.show()
	
def plot_gsphere(g1,g2,Gcut,filename):
	
	data=read_file(filename)
	
	nbrow=data.shape[0]
	nbcol=data.shape[1]
	
	cart=np.empty((nbrow,nbcol-1))
	
	for i in range (0,nbrow):
		gcart=data[i,0]*g1+data[i,1]*g2
		cart[i,:]=array([gcart[0],gcart[1]])
	
	
	
	fig=plt.figure()
	ax=fig.add_subplot(111)		
	plt.scatter (cart[:,0],cart[:,1])
	Gcutoff=plt.Circle((0.0, 0.0), Gcut, color='r', fill=False)
	ax.add_artist(Gcutoff)
	plt.xlim([-Gcut,Gcut])
	plt.ylim([-Gcut,Gcut])
	plt.axes().set_aspect('equal')
	plt.show()	
	
	
	
	

#***********************************************
#***********************************************
period=1.0
hdummy=2.0
phi=pi/3.0

#dimx,dimy
nbpts=[128,128]
filename="test_Bragg"
filenamer="fftmap_output_r"	
filenamei="fftmap_output_i"	
#plotmap(filename,period,phi,nbpts)	
plotmap1D (filename,period,hdummy,nbpts)
#************************************************
#************************************************
#test glist
g1=array([3.6276,-6.2832])
g2=array([3.6276,6.2832])
gabs=sqrt(g1[0]*g1[0]+g1[1]*g1[1])
Gcut=10.1*gabs
filename="outk"








