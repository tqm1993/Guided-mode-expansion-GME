from numpy import *
from scipy import *
from extract_data import *
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker

#add more color order later
base_color=array(['b','r','g','k','m','c','y'])
#for the file format, please check subroutine banddata in mode_solver.f90
def bandplot(filename):

	data=read_file(filename)
	
	lvlcap=data.shape[1]-4
	nbpts=data.shape[0]
	
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.hold(True)
	for i in range(0,lvlcap):
		ax.plot(data[:,0],data[:,i+4],'.r')
	
	#plotting the lightlines
	ax.plot(data[:,0],data[:,3],'-k')
	
	plt.xlim([0.0,1.0])
	plt.ylim([0.0,1.0])
	plt.xlabel('k-path')
	plt.ylabel(r'$\omega$a/2$\pi$c')
	plt.show()
	

def lossplot(filename):

	#read marker infos
	dummy=open(filename,'r')
	
	info=dummy.readline()
	
	markraw=info.split()
	marker=[]
	
	for i in range(1,len(markraw)):
		marker.append(int(markraw[i]))
		
	dummy.close()
	
	#read data
	nbloss=marker[0]
	data=read_file(filename,reject_el="#")
	
	figeig=plt.figure()
	axeig=figeig.add_subplot(111)
	figloss=plt.figure()
	axloss=figloss.add_subplot(111)
	figQ=plt.figure()
	axQ=figQ.add_subplot(111)
	
	axeig.hold(True)
	axloss.hold(True)
	axQ.hold(True)
	tmpQ=[]
	tmploss=[]
	axeig.set_ylim(0.0,1.0)
	nbpts=data.shape[0]
	nbrang=[marker[5],marker[6]]
	#print nbrang
	for i in range (0,nbloss):
		tmpk=[]
		tmpQ=[]
		tmploss=[]
		tmplabel="band "+str(i+nbrang[0])
		#print tmplabel
		axeig.plot(data[:,0],data[:,marker[1]+i],'*',c=base_color[i],label=tmplabel,markersize=6)
		for j in range (0,nbpts):
			if (data[j,marker[4]+i] >0.0):
				tmpk.append(data[j,0])
				tmploss.append(data[j,marker[2]+i])
				tmpQ.append(data[j,marker[3]+i])
		
		axloss.plot(tmpk,tmploss,'-*',c=base_color[i],label=tmplabel,markersize=6)
		axQ.plot(tmpk,tmpQ,'-*',c=base_color[i],label=tmplabel,markersize=6)
		
	
		
	axeig.plot(data[:,0],data[:,3],'-k')
	axeig.set_xlim(0.0,1.0)
	axloss.set_xlim(0.0,1.0)
	axQ.set_xlim(0.0,1.0)
	#axloss.set_yscale("log")
	axQ.set_yscale("log")
	
	axeig.legend()
	axQ.legend()
	axloss.legend()
	
	plt.show()

#plot all energy, Q factor and loss value
#filename refers to the prefix "prefix_e"
def scan_1D(filename):
	
	filenameE=filename+"_e"
	filenameQ=filename+"_q"
	filenameL=filename+"_l"
	#extract scan coordinates
	dummy=open(filenameE,'r')
	info=dummy.readline()
	
	lsraw=info.split()
	ls=[]
	
	for i in range(1,len(lsraw)):
		ls.append(float(lsraw[i]))
	
	#print ls
		
	dummy.close()
	
	#extract data

	
	dataE=read_file(filenameE,reject_el="#")
	dataQ=read_file(filenameQ,reject_el="#")
	dataL=read_file(filenameL,reject_el="#")
	

	figE=plt.figure()
	axEplot=figE.add_subplot(111)
	figQ=plt.figure()
	axQplot=figQ.add_subplot(111)
	figL=plt.figure()
	axLplot=figL.add_subplot(111)
	
	axEplot.plot(ls,dataE[:,0],'*b')
	axQplot.plot(ls,dataQ[:,0],'*b')
	axLplot.plot(ls,dataL[:,0],'*b')
	
	axQplot.set_yscale("log")
	axLplot.set_yscale("log")
	
	plt.show()
	
def scan_2D(filename):
	
	filenameE=filename+"_e"
	filenameQ=filename+"_q"
	filenameL=filename+"_l"
	#extract scan coordinates (scan2 is x axis; scan1 is y axis)
	dummy=open(filenameE,'r')
	
	info1=dummy.readline()
	scan1raw=info1.split()
	scan1=[]
	
	for i in range(1,len(scan1raw)):
		scan1.append(float(scan1raw[i]))
		
	info2=dummy.readline()
	scan2raw=info2.split()
	scan2=[]
	dummy.close()
	
	#extract_data
	
	for i in range(1,len(scan2raw)):
		scan2.append(float(scan2raw[i]))
	
	dataE=read_file(filenameE,reject_el="#")
	dataQ=read_file(filenameQ,reject_el="#")
	dataL=read_file(filenameL,reject_el="#")
	
	figE=plt.figure()
	axEplot=figE.add_subplot(111)
	figQ=plt.figure()
	axQplot=figQ.add_subplot(111)
	figL=plt.figure()
	axLplot=figL.add_subplot(111)

	cpE=axEplot.contourf(scan2,scan1,dataE)
	cpQ=axQplot.contourf(scan2,scan1,dataQ,locator=ticker.LogLocator())
	cpL=axLplot.contourf(scan2,scan1,dataL,locator=ticker.LogLocator())
	
	figE.colorbar(cpE)
	figQ.colorbar(cpQ)
	figL.colorbar(cpL)
	


#Activate to plot band structure or band structure + loss		

#bandplot("toto.dat") #to plot band structure only
#lossplot("toto.dat") #to plot band structure +loss 
