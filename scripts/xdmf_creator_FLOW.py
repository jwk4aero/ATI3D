#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os


######################################
#####        default values      #####
######################################

para_vec=False

# multiblock output? -Dsubspace_multiblock
#multiblock = False
multiblock = True

# precision in bytes
prec       = 4

# byte order: little endian = 'Little' , big endian = 'Big'
endian     = 'Little'

#output option 1:physical space / 2:spectral space
option     = 1 # (1,2)

#number of modes? (1 for zero mode)
modes      = 1 #(not necessary for option 1)

start_TS   = None
end_TS     = None
step_TS    = None


# get command line arguments
i=1
while i < (len(sys.argv)):
    if sys.argv[i]=='--prec' or sys.argv[i]=='-p':
        if sys.argv[i+1]=='4':
            prec=4
        elif sys.argv[i+1]=='8':
            prec=8
        else:
            print ('the next argument following the --prec or -p statement '
                    +'has to be either 4 or 8. Each command line argument has'
                    +' to be provided separately!')
            sys.exit("Quitting")
        i=i+2
        
    elif sys.argv[i]=='--multi' or sys.argv[i]== '-m':
        multiblock = True
        i=i+1
    elif sys.argv[i]=='--block' or sys.argv[i]== '-b':
        multiblock = False
        i=i+1
    elif sys.argv[i]=='--endian' or sys.argv[i]=='-e':
        if sys.argv[i+1]=='Little':
            endian='Little'
        elif sys.argv[i+1]=='Big':
            endian='Big'
        else:
            print ('the next argument following the --endian or -e statement '
                    +'has to be either "Little" or "Big. Each command line argument has'
                    +' to be provided separately!')
            sys.exit("Quitting")     
        i=i+2
    elif sys.argv[i]=='--option' or sys.argv[i]=='-o':
        if sys.argv[i+1]=='1':
            option=1
        elif sys.argv[i+1]=='2':
            option=2
        elif sys.argv[i+1]=='3':
            option=3
        elif sys.argv[i+1]=='11':
            option=1
            para_vec=True
        else:
            print ('the next argument following the --option or -o statement '
                    +'has to be either 1,2 or 3. Each command line argument has'
                    +' to be provided separately!')
            sys.exit("Quitting")
        i=i+2        
    elif sys.argv[i]=='--modes' or sys.argv[i]=='-k':    
        try:
            modes=int(sys.argv[i+1])
        except:
            print ('the next argument following the --modes or -k statement '
                    +'has to an integer. Each command line argument has'
                    +' to be provided separately!')
            sys.exit("Quitting")
        i=i+2              
    elif sys.argv[i]=='--time' or sys.argv[i]=='-t':    
        try:
            start_TS   = int(sys.argv[i+1])
            end_TS     = int(sys.argv[i+2])
            step_TS    = int(sys.argv[i+3])
        except:
            print ('the next argument following the --modes or -k statement '
                    +'has to an integer. Each command line argument has'
                    +' to be provided separately!')
            sys.exit("Quitting")
        i=i+4 
    elif sys.argv[i]=='--help' or sys.argv[i]=='-h':     
        print """
The following command line arguments can be specified:
    -b, --block:    write the xmf header for binary files written per block.
                    It is opposed by the --multi or -m option. No argument 
                    following this entry is expected
    -e, -- endian:  specify the endianess of the binary file. Can be either 
                    "Little" or "Big". It depends on the platform where the 
                    File was written. The default is "Little"
    -h,--help:      Shows this help
    -k,--modes:     The next argument specifies the number for modes in case 
                    spectral data is written. It only has an effect if option
                    2 is chosen.
    -m,--multi:     write the xmf header for multiblock binary files. It is 
                    opposed by the -b or --block flag. It does not expect any
                    following argument
    -o, --option:   The data content of the flow file. Can either be 1 
                    (physical data), 2 (spectral data) or 3 (test numerics 
                    output)    
    -p, --prec:     The next argument has to be an integer that specifies the
                    precission in which the binary files were written. 4 
                    (single) or 8 (double) can be specified
    -t,--time:      The next three (integer) arguments specify start time, 
                    end time and the time step 
    
    
    """
        sys.exit("Quitting")
    else:
        print 'Command line argument '+sys.argv[i]+' not know. It will be skipped' 
        print 'Use option -h for help'
        i=i+1
if option==3:
    start_str='test'            
else:
    start_str='FLOW'

#if time steps not given get them from directory
if start_TS==None or end_TS==None or step_TS==None:
    ls_dir=os.listdir('./')
    timesteps=[]
    for entry in ls_dir:
        if entry.endswith('.raw'):
            if entry.startswith(start_str) :
                temp=(int(entry.split('.')[-2].split('_')[-1]))
                if not any(temp==timesteps[i] for i in range(len(timesteps))):
                    timesteps.append(temp)
            else:
                timesteps.append(int(entry.split('.')[-2].split('_')[-1]))
    timesteps.sort()
else:
    timesteps=range(start_TS,end_TS+step_TS,step_TS)

print ' ' 
print 'The following options have been chosen:'
if multiblock:
    print '-)multiblock data set'
else:
    print '-)one file per block data set'
if prec==4:
    print '-)Single precision data set'
else:
    print '-)Double precision data set'
print '-)'+endian+' Endian'
print '-)option: %i'%option
if option ==2 or option ==3:
    print '-)number of modes: %i'%modes
print '-)number of timesteps %i' %len(timesteps)
print ' ' 

###################################
#####     END USER INPUT      #####
##################################
import struct,glob

if multiblock:
    size_grid_header = 4
    size_flow_header = 4 +4*prec
else:
    size_grid_header = 12
    size_flow_header = 3*4 + 4*prec


# set name of the files
if option == 1:
    # output
    XDMF_filename = 'Flow_phys.xmf'
    # input
    grid_file = 'FLOW_phys_GRID'
    #grid_file = 'test_numerics_q_GRID'
 
    # input - will be appended with <TS>.raw  => subspace_1_1_ leads to e.g. subspace_1_1_600010.raw
    data_file = 'FLOW_phys'
elif option == 2:
    # output
    XDMF_filename = 'Flow_spec.xmf'
    # input
    grid_file = 'FLOW_spec_GRID'
    # input - will be appended with <TS>.raw  => subspace_1_1_ leads to e.g. subspace_1_1_600010.raw
    data_file = 'FLOW_spec'
elif option == 3:
    # output
    XDMF_filename = 'test_numerics.xmf'
    # input
    grid_file = 'test_numerics_q_GRID'
    # input - will be appended with <TS>.raw  => subspace_1_1_ leads to e.g. subspace_1_1_600010.raw
    data_file = 'test_numerics_q'

else:
    import sys
    sys.exit("option not implemented")

# set tags of contained variables
# all variables contained by the Plot3D file MUST be specified!!!
var=[]
if option == 1:
    var.append('$\\rho$')
    var.append('$ u$')
    var.append('$ v$')
    var.append('$ w$')
    var.append('$ T$')
elif option == 2:
    for i in range(modes):
        var.append('$ |F( \\rho)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ |F( \\rho u)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ |F( \\rho v)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ |F( \\rho w)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ |F( \\rho E)|(%i)$'%(i))
elif option == 3:
    for i in range(modes):
        var.append('$ \Delta |F( \\rho)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ \Delta |F( \\rho u)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ \Delta |F( \\rho v)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ \Delta |F( \\rho w)|(%i)$'%(i))
    for i in range(modes):
        var.append('$ \Delta |F( \\rho E)|(%i)$'%(i))


# set the grid size
nxp=[]
nyp=[]
nzp=[]

if endian == 'Little':
    unpackformat_int = '<i'
    unpackformat_flo = '<f'
else:
    unpackformat_int = '>i'
    unpackformat_flo = '>f'


if multiblock:
    f = open(grid_file+'.xyz', 'r+')
    nblocks = struct.unpack(unpackformat_int,f.read(4))[0]
    size_grid_header += 12*nblocks
    size_flow_header += 12*nblocks
    print 'number of blocks: ',nblocks
    for i in range(nblocks):
        nxp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        nyp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        nzp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        print 'block ',i+1,': nxp ',nxp[i],' nyp ',nyp[i],' nzp ',nzp[i]
    f.close()
else:
    nblocks = len(glob.glob('*.xyz'))
    print 'number of blocks: ',nblocks
    for i in range(nblocks):
        f = open(grid_file+'_'+str(i+1)+'.xyz', 'r+')
        nxp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        nyp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        nzp.append(struct.unpack(unpackformat_int,f.read(4))[0])
        print 'block ',i+1,': nxp ',nxp[i],' nyp ',nyp[i],' nzp ',nzp[i]
        f.close()
        
num_var = len(var)

# set time
if multiblock:
    # physical time of first timestep
    f = open(data_file+'_'+str(timesteps[0])+'.raw')
    f.seek(size_flow_header-4)
    start_T = struct.unpack(unpackformat_flo,f.read(4))[0]
    f.close()
    # physical timestep
    if (len(timesteps)>1):
        f = open(data_file+'_'+str((timesteps[1]))+'.raw')
        f.seek(size_flow_header-4)
        dt = struct.unpack(unpackformat_flo,f.read(4))[0]-start_T
        dt=dt/float(timesteps[1]-timesteps[0])
        f.close()
    else:
        dt=0.0
else:
    # physical time of first timestep
    f = open(data_file+'_1_'+str(timesteps[0])+'.raw')
    f.seek(24)
    start_T = struct.unpack(unpackformat_flo,f.read(4))[0]
    f.close()
    # physical timestep
    if (len(timesteps)>1):
        f = open(data_file+'_1_'+str(timesteps[1])+'.raw')
        f.seek(24)
        dt = struct.unpack(unpackformat_flo,f.read(4))[0]-start_T
        
        dt=dt/float(timesteps[1]-timesteps[0])
        f.close()
    else:
        dt=0.0

print 'starting time: ',start_T
print 'timestep: ',dt


#start writing xdmf file
iofile = open(XDMF_filename, 'w')

iofile.write('<?xml version="1.0" ?>\n')
iofile.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
iofile.write('<Xdmf Version="2.0">\n')
iofile.write('<Domain>\n')
#start "grid" for whole domain
iofile.write('<Grid Name="Domain" GridType="Collection">\n')
entry=1
if multiblock:
    grid_offset = size_grid_header
    flow_offset = size_flow_header
    for iblock in range(nblocks):
        #start "grid" for block
        iofile.write('<Grid Name="Block %i" GridType="Collection" CollectionType="Temporal">\n'%(iblock+1))
        iofile.write('<Time TimeType="List">\n')
        iofile.write('<DataItem Format="XML" NumberType="Float" Precision="%i" Dimensions="%i">\n'%(prec,len(timesteps)))
        for istep in timesteps:
            iofile.write(' %f'%(start_T+(istep-timesteps[0])*dt))
        iofile.write('\n')
        iofile.write('</DataItem>\n')
        iofile.write('</Time>\n')
        for istep in timesteps:
            iofile.write('<Grid Name="Time %i">\n'%(istep))
            if istep == timesteps[0]:
                #setup grid for first timestep
                iofile.write('<Topology Name="MainTopology%i" TopologyType="3DSMesh" Dimensions="%i %i %i"/>\n'
                                        %(iblock+1,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('<Geometry Name="MainGeometry%i" GeometryType="X_Y_Z">\n'%(iblock+1))
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %(grid_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s.xyz\n'%(grid_file))
                iofile.write('</DataItem>\n')
                entry+=1
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %(nzp[iblock]*nyp[iblock]*nxp[iblock]*prec+grid_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s.xyz\n'%(grid_file))
                iofile.write('</DataItem>\n')
                entry+=1
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %(nzp[iblock]*nyp[iblock]*nxp[iblock]*prec*2+grid_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s.xyz\n'%(grid_file))
                iofile.write('</DataItem>\n')
                iofile.write('</Geometry>\n')
                entry+=1
            else:
                #setup grid for timesteps > 1 (just referencing to grid of first timestep)
                iofile.write('<Topology Reference="//Topology[@Name=\'MainTopology%i\']"/>\n'%(iblock+1))
                iofile.write('<Geometry Reference="//Geometry[@Name=\'MainGeometry%i\']"/>\n'%(iblock+1))
            if para_vec:
                iofile.write('<Attribute Name="$ \\rho$" Center="Node">\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*0+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s_%i.raw\n'%(data_file,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')
                iofile.write('<Attribute Name="$ \\vec{u}$" Center="Node" AttributeType="Vector"> \n')
                iofile.write('<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="%i %i %i 3"  >\n'%(nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*1+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s_%i.raw\n'%(data_file,istep))
                iofile.write('</DataItem>\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*2+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s_%i.raw\n'%(data_file,istep))
                iofile.write('</DataItem>\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*3+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s_%i.raw\n'%(data_file,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')
                iofile.write('<Attribute Name="$ T$" Center="Node">\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*4+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('%s_%i.raw\n'%(data_file,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')

            else:
                if len(var)>0:
                    #setup data
                    for ivar in range(len(var)):
                        iofile.write('<Attribute Name="%s" Center="Node">\n'%var[ivar])
                        iofile.write('<DataItem Format="Binary" Seek="%i" NumberType="Float" Precision="%i" Endian="%s" Dimensions="%i %i %i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*ivar+flow_offset,prec,endian,nzp[iblock],nyp[iblock],nxp[iblock]))
                        iofile.write('%s_%i.raw\n'%(data_file,istep))
                        iofile.write('</DataItem>\n')
                        iofile.write('</Attribute>\n')
                        entry+=1
            #end timestep grid
            iofile.write('</Grid>\n')
        #end block "grid"
        iofile.write('</Grid>\n')
        grid_offset += nzp[iblock]*nyp[iblock]*nxp[iblock]*prec *3
        flow_offset += nzp[iblock]*nyp[iblock]*nxp[iblock]*prec *len(var)
else:    
    for iblock in range(nblocks):
        #start "grid" for block
        iofile.write('<Grid Name="Block %i" GridType="Collection" CollectionType="Temporal">\n'%(iblock+1))
        iofile.write('<Time TimeType="List">\n')
        iofile.write('<DataItem Format="XML" NumberType="Float" Precision="%i" Dimensions="%i">\n'%(prec,len(timesteps)))
        for istep in timesteps:
            iofile.write(' %f'%(start_T+(istep-timesteps[0])*dt))
        iofile.write('\n')
        iofile.write('</DataItem>\n')
        iofile.write('</Time>\n')
        for istep in timesteps:
            iofile.write('<Grid Name="Time %i">\n'%(istep))
            if istep == timesteps[0]:
                #setup grid for first timestep
                iofile.write('<Topology Name="MainTopology%i" TopologyType="3DSMesh" Dimensions="%i %i %i"/>\n'
                                        %(iblock+1,nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('<Geometry Name="MainGeometry%i" GeometryType="X_Y_Z">\n'%(iblock+1))
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %(size_grid_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                iofile.write('%s_%i.xyz\n'%(grid_file,iblock+1))
                iofile.write('</DataItem>\n')
                entry+=1
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %(nzp[iblock]*nyp[iblock]*nxp[iblock]*prec+size_grid_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                iofile.write('%s_%i.xyz\n'%(grid_file,iblock+1))
                iofile.write('</DataItem>\n')
                entry+=1
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %(nzp[iblock]*nyp[iblock]*nxp[iblock]*prec*2+size_grid_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                iofile.write('%s_%i.xyz\n'%(grid_file,iblock+1))
                iofile.write('</DataItem>\n')
                entry+=1
                iofile.write('</Geometry>\n')
            else:
                #setup grid for timesteps > 1 (just referencing to grid of first timestep)
                iofile.write('<Topology Reference="//Topology[@Name=\'MainTopology%i\']"/>\n'%(iblock+1))
                iofile.write('<Geometry Reference="//Geometry[@Name=\'MainGeometry%i\']"/>\n'%(iblock+1))
            if option==para_vec:
                iofile.write('<Attribute Name="$ \\rho$" Center="Node">\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %(size_flow_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')
                iofile.write('<Attribute Name="$ \\vec{u}$" Center="Node" AttributeType="Vector"> \n')
                iofile.write('<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="%i %i %i 3"  >\n'%(nzp[iblock],nyp[iblock],nxp[iblock]))
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i " Precision="%i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*1+size_flow_header,endian,nxp[iblock],nyp[iblock],nzp[iblock],prec))
                iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                iofile.write('</DataItem>\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i " Precision="%i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*2+size_flow_header,endian,nxp[iblock],nyp[iblock],nzp[iblock],prec))
                iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                iofile.write('</DataItem>\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i " Precision="%i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*3+size_flow_header,endian,nxp[iblock],nyp[iblock],nzp[iblock],prec))
                iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')
                iofile.write('<Attribute Name="$ T$" Center="Node">\n')
                iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*4+size_flow_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                iofile.write('</DataItem>\n')
                iofile.write('</Attribute>\n')



            else:
                if len(var)>0:
                   #setup data
                    for ivar in range(len(var)):
                        iofile.write('<Attribute Name="%s" Center="Node">\n'%var[ivar])
                        iofile.write('<DataItem Format="Binary" Seek="%i" Endian="%s" Dimensions="%i %i %i" Precision="%i">\n'
                                        %((nzp[iblock]*nyp[iblock]*nxp[iblock]*prec)*ivar+size_flow_header,endian,nzp[iblock],nyp[iblock],nxp[iblock],prec))
                        iofile.write('%s_%i_%i.raw\n'%(data_file,iblock+1,istep))
                        iofile.write('</DataItem>\n')
                        iofile.write('</Attribute>\n')
                        entry+=1
                    #iofile.write('<Attribute Name="$ \\vec{u}$" Center="Node"> \n')
                    #iofile.write('<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="%i %i %i"  >\n'%(nzp[iblock],nyp[iblock],nxp[iblock]))
                    #iofile.write('<DataItem Reference="/Xdmf/DataItem[%i]" />\n'%(entry-4))
                    #iofile.write('<DataItem Reference="/Xdmf/DataItem[%i]" />\n'%(entry-3)) 
                    #iofile.write('<DataItem Reference="/Xdmf/DataItem[%i]" />\n'%(entry-2))
                    #iofile.write('</DataItem>\n') 
                    #iofile.write('</Attribute>\n') 
                    #entry+=1
            #end timestep grid
            iofile.write('</Grid>\n')
        #end block "grid"
        iofile.write('</Grid>\n')


#end global "grid"
iofile.write('</Grid>\n')
iofile.write('</Domain>\n')
iofile.write('</Xdmf>')

iofile.close()
print 'Finished succesfully'
print 'To open the data set load file "%s"'%XDMF_filename
