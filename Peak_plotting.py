#Code for plotting lightcurve of each peak in each quadrant. 
#Written by Priyanka Shahane 
#Date: 18-june-2018 Time 11:35 AM
###################################################################################################################
from astropy.table import Table
from numpy import arange
import numpy as np
import matplotlib.pyplot as plt
#----------------------------------User inputs------------------------------#
event_filename ='/home/cztipoc/AS1G06_029T01_9000000964_07040cztM0_level2_quad_clean.evt'
outputfits ='/home/cztipoc/9000000964_07040_vmean.fits'
pha_list = fits.open(event_filename, memmap=True)
pha_list1 = fits.open(outputfits, memmap=True)
#---------------------------------------------------------------------------#
#Get info of the fits file
pha_list.info()

# load the data in separate variable
pha_data1 = pha_list[1].data
pha_data2 = pha_list[2].data
pha_data3= pha_list[3].data
pha_data4 = pha_list[4].data
pha_data_op = pha_list1[1].data
#Peak Time
peak_time = pha_data_op['Ptime(s)'] 
#Vmean of Q1
Vm_Q1 = pha_data_op['VmQ1'] 
#Vmean of Q2
Vm_Q2 = pha_data_op['VmQ2']
#Vmean of Q3
Vm_Q3 = pha_data_op['VmQ3']
#Vmean of Q4
Vm_Q4 = pha_data_op['VmQ4']
#Time Q1
Time1 = pha_data1['TIME']
#Time Q2
Time2 = pha_data2['TIME']
#Time Q3
Time3 = pha_data3['TIME']
#Time Q4
Time4 = pha_data4['TIME']
#print Time1
  #------------------- Loop over  Peak_time --------------------#
for i in range(len(peak_time)):
    #Taking time 50s before peak time
    less=peak_time[i]-50.0
    #Taking time 50s after peak time
    greater=peak_time[i]+50.0
    #Array of time 50s before and after peak time for Q1. 
    op1= Time1[ (Time1 >= less) & (Time1 <= greater)]
    op1[:] = [x - peak_time[i] for x in op1]
    #Array of time 50s before and after peak time for Q2.
    op2= Time2[ (Time2 >= less) & (Time2 <= greater)]
    op2[:] = [x - peak_time[i] for x in op2]
    #Array of time 50s before and after peak time for Q3.
    op3= Time3[ (Time3 >= less) & (Time3 <= greater)]
    op3[:] = [x - peak_time[i] for x in op3]
    #Array of time 50s before and after peak time for Q4.
    op4= Time4[ (Time4 >= less) & (Time4 <= greater)]
    op4[:] = [x - peak_time[i] for x in op4]
    #-----------------------Plotting peak time in each quadrant----------------#
    x=np.linspace(-50,50,100)
    fig=plt.figure(1)
    plt.suptitle('Peak Time='+str(peak_time[i]))
    plt.subplot(221)
    y=np.ones(100)*float(Vm_Q1[i])
    plt.hist(op1,bins=(int)(op1[-1]-op1[0]) ,histtype='step',color="green")
    plt.plot(x,y,label=round(Vm_Q1[i]),color="red")
    plt.legend(loc='upper right',fontsize='small', framealpha=0.2)
    plt.xlabel('Time (s)')
    plt.ylabel('Counts/Bin')
    plt.title('Q0')
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.subplot(222)
    y=np.ones(100)*float(Vm_Q2[i])
    plt.hist(op2,bins=(int)(op2[-1]-op2[0]),histtype='step',color="green")
    plt.plot(x,y,label=round(Vm_Q2[i]),color="red")
    plt.legend(loc='upper right',fontsize='small', framealpha=0.2)
    plt.xlabel('Time (s)')
    plt.ylabel('Counts/Bin')
    plt.title('Q1')
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.subplot(223)
    y=np.ones(100)*float(Vm_Q3[i])
    plt.hist(op3,bins=(int)(op3[-1]-op3[0]),histtype='step',color="green")
    plt.plot(x,y,label=round(Vm_Q3[i]),color="red")
    plt.legend(loc='upper right',fontsize='small', framealpha=0.2)
    plt.xlabel('Time (s)')
    plt.ylabel('Counts/Bin')
    plt.title('Q2')
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.subplot(224)
    y=np.ones(100)*float(Vm_Q4[i])
    #plt.subplot(2, 2, 4)
    plt.hist(op4,bins=(int)(op4[-1]-op4[0]),histtype='step',color="green")
    plt.plot(x,y,label=round(Vm_Q4[i]),color="red")
    plt.legend(loc='upper right',fontsize='small', framealpha=0.2)
    plt.xlabel('Time (s)')
    plt.ylabel('Counts/Bin')
    plt.title('Q3')
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    fig.savefig('/home/cztipoc/peak_9000000964_07040'+str(i)+'.png',bbox_inches='tight') 
    plt.show()

