# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:11:21 2020
@author: Benjamin
"""
from pygimli.meshtools import readGmsh
import numpy as np
import pybert as pb
import os 
from pygimli.viewer.mpl import drawStreams
from scipy import interpolate
import matplotlib.pyplot as plt
import pygimli as pg
    
class fct_utils(): 

    def definePath(main,date,**kwargs):
        #%% define PATHS 
        # or set you working directory local path here 
        #main = os.getcwd()
        #os.chdir(main)
        
        geomPath= './geom/'
        meshPath= './mesh/'
        icsdPath= './icsd/'
        figPath= './fig/'
        processedPath= './processed_data/'
        
        rmvInvalid = False
        all_gates = False
        Nfix = False
        gateIP = False
        
        icsdPath += date 
        if rmvInvalid is False: 
            icsdPath += '_raw' 
        if all_gates: 
            icsdPath += '_AllG' 
        if Nfix: 
            icsdPath += '_Nfix' 
        if gateIP: 
            icsdPath += '_M' + str(gateIP)
        icsdPath += '/'
        
        figPath += date + '/'
        if not os.path.exists(figPath):
            os.makedirs(figPath)

        processedPath += date + '/'
        if not os.path.exists(processedPath):
            os.makedirs(processedPath)
            
    
        return geomPath, meshPath, icsdPath, figPath, processedPath


#%%  ------- MALM -----------------------------------------------------####

    def mk_seqMALM(a_injection = 72, b_return = 2, sensors=None, mesh=None, 
                   Syscal=None, Rec=None, 
                   PosXYZid=None,Ch_return=None, 
                   mfix=False, Bmove=(False,[]), mref=40,
                   savefig=False,SaveName=None):
        
        """
        write the rhizotron malm sequence
        parameters:
        ===
        * a_injection: number of the electrode used for injection (t: int) , default = 1
        * b_return: number of the electrode used as return (t: int), default = 2
        * b_return associated with Ch_return --> b_return_left = b_return[0], while b_return_right = b_return[1] 
        returns:
        ===
        * seq: malm sequence (t: np.ndarray)
        """
        print('*' * 20)

        # function    
        a = a_injection -1 
        b= b_return - 1 

        seq = []
        
        if mfix==True:
            m = mref -1
            electrodes = np.arange(1,len(sensors))
            for row in electrodes:
                print(row)
                n = row
                seq.append([a, b, m, n]) 

            
        else:
            # horizontal
            electrodes = np.arange(0,72)
            electrodes = electrodes.reshape((9,8))
            for row in electrodes:
                for i in range(len(row) -1):
                    m = row[i]
                    n = row[i + 1]
                    seq.append([a, b, m, n])
        
            # vertical
            rotated_electrodes = np.flipud(np.rot90(electrodes))
            for row in rotated_electrodes:
                for i in range(len(row) -1):
                    m = row[i]
                    n = row[i + 1]
                    seq.append([a, b, m, n])
        
            # diag \ -6 : 6
            for d in range(-6, 7):
                diagonal = np.diag(electrodes, k = d)
                for i in range(len(diagonal) - 1):
                    m = diagonal[i]
                    n = diagonal[i + 1]
                    seq.append([a, b, m, n])
        
            # diag / 
            flipped_electrodes = np.fliplr(electrodes)
            for d in range(-6, 7):
                diagonal = np.diag(flipped_electrodes, k = d)
                for i in range(len(diagonal) - 1):
                    m = diagonal[i]
                    n = diagonal[i + 1]
                    seq.append([a, b, m, n])
        
        # remove where the same electrode appears twice
        sequence = np.array(seq)

        indices_delete = []
        for i, row in enumerate(sequence):
            if len(np.unique(row)) != len(row):
                indices_delete.append(i)
        sequence_clean = np.delete(sequence, indices_delete, axis = 0)
        sequence_clean = np.vstack(sequence_clean)


        id2Remove=[]
        if PosXYZid is not None:
            print('-' * 5)
            print('Split case remove wrong quadripoles and adapt return electrode')
            PosXYZid=np.array(PosXYZid)
            Left=np.array(PosXYZid[0,:])   
            Right=np.array(PosXYZid[1,:]) 
            print('TEST1: IF M and N are on the same SIDE')
            print('TEST2: IF B return same side as M and N')
            for i in range(0,len(sequence_clean)):
                # check if line ismenber of left or right
                # ---------# 
                # 1st TEST #
                # TEST IF M and N are on the same SIDE
                # ---------# 
#                print(sequence_clean[i, 2:4])
                Left_Test =  any(elem in sequence_clean[i, 2:4]+1  for elem in Left)
                Right_Test =  any(elem in sequence_clean[i, 2:4]+1 for elem in Right) 
                #print('Left_Test=' + str(Left_Test) + '   Right_Test=' + str(Right_Test))
                if i in id2Remove:
                    continue
                try:
                    if all([Left_Test,Right_Test])==True:
                        id2Remove.append(i)
                except:
                    continue 

                if Ch_return is not None:
                    # ---------# 
                    # 2nd TEST #
                    # TEST IF B return same side as M and N 
                    # ---------# 
                    B_Test = b_return in Left
                    if (B_Test==False and Left_Test==True):
                        sequence_clean[i, 1]=Ch_return
                        
            print('sequence before removed split=' + str(len(sequence_clean)))
            sequence_clean= np.delete(sequence_clean,id2Remove, axis=0)
            print('sequence after removed split=' + str(len(sequence_clean)))
            
            # remove where the same electrode appears twice
            sequence = np.array(sequence_clean)
            indices_delete = []
            for i, row in enumerate(sequence):
                if len(np.unique(row)) != len(row):
                    indices_delete.append(i)
            print(indices_delete)
            
            sequence_clean = np.delete(sequence, indices_delete, axis = 0)
            sequence_clean = np.vstack(sequence_clean)

        if Bmove[0]==True:
            sequence_clean_base= sequence_clean
            for b in Bmove[1]:
                print('b=' + str(b))
                for i in range(0,len(sequence_clean_base)):
                    sequence_clean_b= sequence_clean_base
                    #print(i,b)
                    sequence_clean_b[i, 1]=b
                sequence_clean= np.vstack([sequence_clean,sequence_clean_b])

                # remove where the same electrode appears twice
                sequence = np.array(sequence_clean)
                indices_delete = []
                for i, row in enumerate(sequence):
                    if len(np.unique(row)) != len(row):
                        indices_delete.append(i)
                
                sequence_clean = np.delete(sequence, indices_delete, axis = 0)
                sequence_clean = np.vstack(sequence_clean)
        
        if Rec==True:
           print('add rec')
           sequence_cleanR=sequence_clean[:,[2,3,0,1]] 
           sequence_clean= np.vstack([sequence_clean,sequence_cleanR])
    
        schemeMALM = pb.DataContainerERT()
        schemeMALM.resize(len(sequence_clean))        
        for i, j in enumerate("abmn"):
            schemeMALM.set(j, sequence_clean[:, i])


        if mesh!=None:            
            sensors=[]
            posXYZ=[]
            for node in mesh.nodes():
                if node.marker() == -99:
                    sensors.append(node.pos())
                    posXYZ.append(np.array(node.pos()))

        if len(sensors)>0: 
            schemeMALM.setSensorPositions(sensors)
            schemeMALM.set("valid", np.ones(len(sequence_clean)))
            if PosXYZid is not None:
                schemeMALM.save('SequenceMALM_Rhizo_72Elecs_splitted.shm')
            else:
                schemeMALM.save('SequenceMALM_Rhizo_72Elecs.shm')
                
                
            schemeMALM.setSensorPositions(sensors)
            schemeMALM.set("valid", np.ones(len(sequence_clean)))
            if PosXYZid is not None:
                schemeMALM.save('SequenceMALM_Rhizo_72Elecs_splitted.shm')
            else:
                schemeMALM.save('SequenceMALM_Rhizo_72Elecs.shm')
        
        if Syscal==True:
            print('Write syscal sequence file')
            posXYZ= np.vstack(posXYZ)
            posXYZ_ElectrPro= [np.arange(1,len(posXYZ)+1),np.zeros([2,len(posXYZ[:,0])])]
            posXYZ_ElectrPro= np.vstack(posXYZ_ElectrPro).T
            nb = np.arange(1,len(posXYZ_ElectrPro)+1)
            posXYZ_nb= np.column_stack((nb,posXYZ_ElectrPro))
            measSyScal= sequence_clean+1
            nb = np.arange(1,len(measSyScal)+1)
            measSyScal_nb= np.column_stack((nb,measSyScal))
            if PosXYZid is not None:
                fileName='MALMR'+ str(np.amax(measSyScal)) +'_Splitted.txt'
            else:
                fileName='MALMR'+ str(np.amax(measSyScal)) +'.txt'
            f = open(fileName,'w')
            np.savetxt(f, posXYZ_nb, fmt='%d %1.2f %1.2f %1.2f', delimiter='\t',header='X Y Z')   # X is an array
            np.savetxt(f, measSyScal_nb, fmt='%d', delimiter='\t',header='A	 B	  M	   N')   # X is an array
            f.close()

        if savefig==True:
            plt.ioff()
            fig = plt.figure()
            sensorsArray = np.vstack(sensors)
            from celluloid import Camera
            camera = Camera(plt.figure())
            for qq in np.arange(0,len(sequence_clean),int(np.round(len(sequence_clean)/72))):
                plt.scatter(sensorsArray[:,0],sensorsArray[:,1],color='b',alpha=0.2)
                plt.scatter(sensors[sequence_clean[qq,1]][0],sensors[sequence_clean[qq,1]][1],color='r',alpha=0.9)
                plt.scatter(sensors[sequence_clean[qq,2]][0],sensors[sequence_clean[qq,2]][1],color='g',alpha=0.9)
                plt.scatter(sensors[sequence_clean[qq,3]][0],sensors[sequence_clean[qq,3]][1],color='g',alpha=0.9)
                camera.snap()
            anim = camera.animate(blit=True) #interval=200,
            if PosXYZid is not None:
                AnimName='SeqMALMGifDescription_Splitted.mp4'
            else:
                AnimName='SeqMALMGifDescription.mp4'
            anim.save(AnimName)
            plt.close(fig)

        return(sequence_clean)
  

    def PrepareMALMData(DataName, Rec=False, 
                        MinV=1e-6, MaxRc=1e4, DevErr=5, Kfact=1e3, MinMaxAppRes=1e3, Rscheck=50,
                        SwE=False, ExpRec=True, gIP=None, RmvE=None,
                        valid=None,**kwargs):
        """
        Prepare MALM data for ICSD analysis
        
        parameters:
        ===
        return:
        ===
        """
        Obs = pb.importer.importSyscalPro(DataName) 
        Obs.set('r', Obs('u')/Obs('i'))
        Obs.set('Rec',np.zeros(len(Obs('u'))))
        print("Data before filtering:", Obs)
        process_folder = os.getcwd()
        savefile=False
        for key, value in kwargs.items():
            if key == 'savefile':
                savefile = value
            if key == 'date':
                date = value
                process_folder = './processed_data/' + date +'/' 
                
                if not os.path.exists(process_folder):
                    os.makedirs(process_folder)
                
        if Rec==True:
            sC = Obs.sensorCount()
            # Build an unique Index list for one measurement direction:
            idxF = Obs('a')*sC**3 + Obs('b')*sC**2 + Obs('m')*sC**1 + Obs('n')
            # Build a second Index list for reciprocal measurement direction:
            idxR = Obs('m')*sC**3 + Obs('n')*sC**2 + Obs('a')*sC**1 + Obs('b')
            # for each measurement in idxF search the same  idx in idxR which is your reciprocal
            #threshold=5
            recip_err=[]
            count=0
            idc2Skip = [] # search for duplicates reciprocals 
            for i, dr in enumerate(idxF):
                if i in idc2Skip:
                    recip_err.append(0)
                else:
                    idc= np.where(dr==idxR)[0]
                    if idc.size>0:
                        count = count + 1
                        R_dif = np.abs(Obs('r')[i] - Obs('r')[idc[0]])
                        R_avg = (Obs('r')[i] + Obs('r')[idc[0]]) / 2
                        recip_err.append(np.abs(R_dif / R_avg) * 100)
                        idc2Skip.append(idc[0])
                    else:
                        recip_err.append(0)
            Dup_Vec_Bool= np.zeros([len(Obs('a'))])
            Dup_Vec_Bool[idc2Skip]=1


            idDouble = [] # search for doubles quadrupoles 
            for i, dr in enumerate(idxF):
                    idc= np.where(dr==idxF)[0]
                    if idc.size>1:
                        idDouble.append(1)
                    else:
                        idDouble.append(0)

            if ExpRec==True:
                print('Write reciprocal file')               


            Obs.set('RecBool',Dup_Vec_Bool)
            Obs.set('Double',idDouble)
            Obs.markInvalid(Obs('Double') == 1)
            Obs.removeInvalid()
            print("Data after filtering double", Obs)

            Obs.markInvalid(Obs('err') > DevErr)
            Obs.removeInvalid()
            print("Data after filtering contact resistance:", Obs)

            Obs.markInvalid(Obs('Rec') > MaxRc)
            Obs.removeInvalid()
            print("Data after filtering rec threshold:", Obs)
            Obs.markInvalid(Obs('RecBool') == 1)
            Obs.removeInvalid()
            print("Data after remove rec:", Obs)
            
        print("Data after filtering:", Obs)

        if RmvE is not None:
            Obs.set('RmvE_Vec_Bool',RmvE)
            Obs.markInvalid(Obs('RmvE_Vec_Bool') == 1)
            Obs.removeInvalid()
            print("Data after RmvE:", Obs)

            
        if SwE==True:
            print('Switch electrode pair')
            IdColElecSwE_1=[]
            IdColElecSwE_2=[]
            for iSwE in enumerate(SwE):
                print(iSwE)
                for i, j in enumerate("abmn"):
                    IdLineElecSwE_1= np.where(iSwE[1][0]-1==np.array(Obs(j)))[0]
                    IdLineElecSwE_2= np.where(iSwE[1][1]-1==np.array(Obs(j)))[0]
                    IdColElecSwE_1.append(IdLineElecSwE_1)
                    IdColElecSwE_2.append(IdLineElecSwE_2)
            for i, j in enumerate("abmn"):
                Obs(j)[IdColElecSwE_1[i][:]]=[iSwE[1][1]-1]*np.ones(len(IdColElecSwE_1[i][:]))
                Obs(j)[IdColElecSwE_2[i]]=[iSwE[1][0]-1]*np.ones(len(IdColElecSwE_2[i][:]))
                Obs.set('a', Obs('a'))
                Obs.set('b', Obs('b'))
                Obs.set('m', Obs('m'))
                Obs.set('n', Obs('n'))
                
       
        if valid is not None:
            Obs.set('valid',valid)
            Obs.removeInvalid()
            print("Data after valid:", Obs)

        
        if savefile==True: 
            if gIP is not None: 
                print('TDIP export')
                sep = '_'
                NameSave = process_folder + 'O'+ os.path.basename(DataName).split(sep, 1)[0] + 'M'+str(gIP) + '.txt'
                f = open(NameSave,'w')
                np.savetxt(f, np.array(Obs['M'+str(gIP)]), delimiter='\t',fmt='%1.6f')   # X is an array
                f.close()
            else:
                sep = '.'
                NameSave = process_folder + 'O'+ os.path.basename(DataName).split(sep, 1)[0] + '.txt'
                f = open(NameSave,'w')
                np.savetxt(f, np.array(Obs("r")), delimiter='\t',fmt='%1.6f')   # X is an array
                f.close()
                print(NameSave)
        
        
        dataABMN = [np.array(Obs('a'))+1, np.array(Obs('b'))+1,
                    np.array(Obs('m'))+1,np.array(Obs('n'))+1]
        dataABMN= np.vstack(dataABMN).T
        return [Obs, dataABMN]
    
    def VRTEpos(mesh=None,dim=3,MarkerVRTE=991):
        """
        Find position of VRTE into the mesh
        
        parameters:
        ===
        * mesh_VRTs: the mesh with all virtuals sources included

        return:
        ===
        * VRTeCoord.txt : adequate file to be inverted into ICSD
        """
        VRTE=[]
        Pos_vrtS=[]
        for node in mesh.nodes():
            if node.marker() == -MarkerVRTE:
                VRTE.append(node.pos())
                Pos_vrtS.append(np.array(node.pos()))
        print(str(len(VRTE)) + ' VRTEs found')
        Pos_vrtS= np.vstack(Pos_vrtS)
    
        # --- Writing the VRTe Coordinates into a file! ------------------------------
        
        #f = open('VRTeCoord.txt','w')
        #np.savetxt(f, Pos_vrtS[:,0:2], delimiter='\t',fmt='%1.3f')   # X is an array
        # f.close()
        #if dim==3:
        #     f = open('VRTeCoord.txt','w')
        #     np.savetxt(f, Pos_vrtS[:,0:3], delimiter='\t',fmt='%1.3f')   # X is an array
        #     f.close()
        return Pos_vrtS
    
    
    def mesh_import(fname):
        """
        Import a gmsh mesh and return mesh and sensors positions according to markers
        Returns
        -------
        None.
    
        """
        ## --------- Import mesh --------- ##
        mesh3d=readGmsh(fname, verbose=True)
        #mesh3d.exportVTK('mesh3d_Markers')
        #MR.PyMesh3d() # plot mesh3d_Markers marker in 3d
        sensors = []
        for node in mesh3d.nodes():
            if node.marker() == -99:
                sensors.append(node.pos())
            elif node.marker() == -999:
                print('-999')
                sensors.append(node.pos())
            elif node.marker() == -1000:
                print('-1000')
                sensors.append(node.pos())
        return mesh3d, np.vstack(sensors)


    def streamlines(coordObs, Obs, model, mesh=None, **kwargs):
        """
        Current streamlines
        -------
        - mesh: mesh to compute the quiver plot
    
        """
        mesh_inv = [] # medium conductivity
        for key, value in kwargs.items():
            if key == 'Xaxis':
                 Xaxis = value            
            if key == 'mesh_inv':
                 mesh_inv = value
            if key == 'date':
                 date = value
                 
        xn = 30
        yn = 30

        xx = np.linspace(min(coordObs[:,1]), max(coordObs[:,1]),xn)
        yy = np.linspace(min(coordObs[:,2]), max(coordObs[:,2]),yn)

        xx, yy = np.meshgrid(xx, yy)
        points = np.transpose(np.vstack((coordObs[:,1], coordObs[:,2])))


        if mesh is None:
            mesh = pg.createGrid(x=np.linspace(min(coordObs[:,1]), max(coordObs[:,1]),xn),
                     y=np.linspace(min(coordObs[:,2]), max(coordObs[:,2]),yn))

        u_interp = interpolate.griddata(points,
                                        Obs,
                                        (xx, yy), 
                                        method='cubic')
        uu = np.reshape(u_interp,[xn*yn])
        #uu = pg.interpolate(mesh, (xx, yy), Obs, method='spline')

        
        if isinstance(model, float):
            stream = -pg.solver.grad(mesh, uu)*(1/model)     
            jj = -uu*(1/model)
        else: 
            res = pg.interpolate(mesh, mesh_inv, model, method='spline')
            jj = -uu*(1/res).array()
            stream = -pg.solver.grad(mesh, uu)*(1/res).array()[:, None]
            

        if kwargs.get('vmin'):
           vmin = kwargs.get('vmin')
        else: 
           vmin = min(Obs)
           
        if kwargs.get('vmax'):
           vmax = kwargs.get('vmax')
        else: 
           vmax = max(Obs)

        if kwargs.get('ax'):
            ax = kwargs.get('ax')
        else:
            fig, ax = plt.subplots()
            
        sc=ax.scatter(coordObs[:,1], coordObs[:,2], c=Obs, 
                      cmap ='coolwarm',s=5e2, vmin=vmin, vmax=vmax) # norm=matplotlib.colors.Normalize()
        cbar = plt.colorbar(sc,ax=ax)
        cbar.set_label('V')   
        ax.set_ylabel('y [m]',fontsize=15)
        ax.set_xlabel('x [m]',fontsize=15)
        
        if len(kwargs.get('sensors'))>0:
            sensors = kwargs.get('sensors')
            ax.scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
            for i in range(len(sensors[:,0])):
                   ax.annotate(str(i+1), (sensors[i,0],sensors[i,1]))   
            if kwargs.get('A'):
                A = kwargs.get('A')
                ax.scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
            if kwargs.get('B'):
                B = kwargs.get('B')
                ax.scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
            if kwargs.get('Nfix'):
                Nfix = kwargs.get('Nfix')
                ax.scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')

        if kwargs.get('gridCoarse'):
            gridCoarse = pg.createGrid(x=np.linspace(min(sensors[:,0]), max(sensors[:,0]),xn/2),
                             y=np.linspace(min(sensors[:,1]), max(sensors[:,1]),yn/2))
        

        drawStreams(ax, mesh, stream,
                     color='green', quiver=True, linewidth=3.0)
        
        
        # save for TL analysis
        
        #j = -pg.solver.grad(mesh, uu) * (1/Res)
        #ax, _ = pg.show(mesh, hold=True, alpha=0.3)
        #drawStreams(ax, mesh, j)
        
        return mesh, uu, model

    def Ref_N_zeroP(ref_elec=71):
        """
        The potential at the reference electrode is used as a zero potential
        and is removed to the electrical potential collected at each electrode
        Returns
        -------
        None.

        """
        
    def plot_scatter_obs(x,y,V):
        """
        Plot scattered potential
        Returns
        -------
        None.

        """
        
    def load_geom(path):
        """ load the geometry of the acquisition (*geom file custum for Mise-à-la-masse data)
        
        Parameters
        ----------

        """
        geom_files = [f for f in os.listdir(path) if f.endswith('.geom')]
        if len(geom_files) != 1:
            raise ValueError('should be only one geom file in the current directory')
        
        fileNameElec = geom_files[0]  
        line_number = 0
        line_of_injection = []
        line_of_remotes = []
        # Open the file in read only mode
        with open(path + fileNameElec, 'r') as read_obj:
            # Read all lines in the file one by one
            for line in read_obj:
                # For each line, check if line contains the string
                # print(line)
                line_number += 1
                if ('#Remote') in line:
                    # If yes, then add the line number & line as a tuple in the list
                    line_of_remotes.append((line_number))
                if ('#Injection') in line:
                    line_of_injection.append((line_number))
        
        # print(line_of_injection)
        # print(line_of_remotes)
        RemLineNb= int(line_of_remotes[0])-1
        Injection= int(line_of_injection[0])-1
        
        coordE = np.loadtxt(path+ fileNameElec)
        pointsE= np.vstack(coordE[:RemLineNb,1:4])
        
        return RemLineNb, Injection, coordE, pointsE

    def filterTDIP(dataTDIP,id2rmv):
        valid = np.ones(len(dataTDIP.data('m')))
        valid[id2rmv] = 0

        
        #if remove == True:
            #dataTDIP.data.set('valid',valid)
            #dataTDIP.MA = dataTDIP.MA[:, dataTDIP.data['valid'].array()==1]
            #dataTDIP.data.removeInvalid()

        return dataTDIP, valid
    
#%%
# if Nfix is not None:
#     fig, ax = plt.subplots(nrows=1, ncols=4)
#     for i, g in enumerate(range(3,20,5)):
#         ax[i].set_ylabel('y [m]',fontsize=15)
#         ax[i].set_xlabel('x [m]',fontsize=15)
#         sc=ax[i].scatter(coordE_f[:,1], coordE_f[:,2], c=Obs('M'+str(g)).array(), 
#                       cmap ='coolwarm',s=5e2, vmin=-100)
#         cbar = plt.colorbar(sc,ax=ax[i]) 
#         ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')
#         for ie in range(len(sensors[:,0])):
#             ax[i].annotate(str(ie+1), (sensors[ie,0],sensors[ie,1]))


# #%%


#     if Nfix is not None:
#     fig, ax = plt.subplots(nrows=1, ncols=2)
#     ax[0].scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
#     ax[0].scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
#     ax[0].scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
#     ax[0].scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
#     for i in range(len(sensors[:,0])):
#         ax[0].annotate(str(i+1), (sensors[i,0],sensors[i,1]))
    
#     ax[0].legend(loc="upper right")    
#     ax[0].set_xlabel('x (m)')
#     ax[0].set_ylabel('y (m)')
    
#     # plot secondary voltage distribuion along a given profile
#     abmn=[]
#     for nn in range(len(IPcurves_f.data['a'])):
#         abmn.append([int(IPcurves_f.data(t)[nn]+1) for t in ['a', 'b', 'm', 'n']])
#     abmn = np.vstack(abmn)
    
#     # define the profile
#     pelecs = np.array([0,1,2,3,4,5,6,7]) +32
#     pelecs = np.array([0,1,2,3,4,5,6,7]) +16
#     pelecs = np.array([0,1,2,3,4,5,6,7]) +32
#     pelecs = np.array([5,13,21,29,37,45,53,61,69])
    
#     idp = np.zeros(len(pelecs))
#     Vg =  np.zeros(20)
#     Vi = np.zeros([20,len(pelecs)])
#     for i,p in enumerate(pelecs):
#         if(len(np.where(p+1==abmn)[0]))>0:
#             idp[i]= np.where(p+1==abmn)[0][0]
#         for g in range(20):
#             gatestr = 'M' + str(g+1)
#             Vg[g] = IPcurves_f.data[gatestr].array()[int(idp[i])]
#         Vi[:,i] = Vg
        
#     for g in range(3,20,5):
#         ax[1].plot(pelecs+1,Vi[g,:],'o-', alpha=0.5, 
#                  label='Gate t:' + str(IPcurves.t[g]) + 's')
#     plt.legend()    
#     plt.xlabel('# Electrode')
#     plt.ylabel('Secondary voltage (mV)')
#     plt.xticks(pelecs+1)

# if Nfix is not None:
#     fig, ax = plt.subplots(nrows=1, ncols=2)
#     sc=ax[0].scatter(coordE_f[:,1], coordE_f[:,2], c=abs(Obs('r').array()), 
#                   cmap ='coolwarm',s=5e2, vmin=0.1, vmax=200) # norm=matplotlib.colors.Normalize()
#     cbar = plt.colorbar(sc,ax=ax[0])
#     cbar.set_label('V')   
#     ax[0].scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
#     ax[0].scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
#     ax[0].scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
#     ax[0].scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
#     for i in range(len(sensors[:,0])):
#         ax[0].annotate(str(i+1), (sensors[i,0],sensors[i,1]))   
#     ax[0].set_ylabel('y [m]',fontsize=15)
#     ax[0].set_xlabel('x [m]',fontsize=15)
#     sc=ax[1].scatter(coordE_f[:,1], coordE_f[:,2], c=Obs('gm').array(), 
#                   cmap ='coolwarm',s=5e2, norm=matplotlib.colors.Normalize())
#     cbar = plt.colorbar(sc,ax=ax[1]) 
#     ax[1].set_title('Global chargeability')
#     