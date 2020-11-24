# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:11:21 2020
@author: Benjamin
"""
from pygimli.meshtools import readGmsh
import numpy as np
import pybert as pb
import os

class fct_utils(): 
    
    def PrepareMALMData(DataName, Rec=False, 
                        MinV=1e-6, MaxRc=1e4, DevErr=5, Kfact=1e3, MinMaxAppRes=1e3, Rscheck=50,
                        SwE=False, ExpRec=True, gIP=None, RmvE=None):
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
                
        Obs.save('ObservationData_Afterfilter.dat','a b m n r')        
        sep = '_'
        NameSave = 'O'+ os.path.basename(DataName).split(sep, 1)[0] + '.txt'
        f = open(NameSave,'w')
        np.savetxt(f, np.array(Obs("r")), delimiter='\t',fmt='%1.6f')   # X is an array
        f.close()
        
        if gIP: 
            print('TDIP export')
            sep = '_'
            NameSave = 'O'+ os.path.basename(DataName).split(sep, 1)[0] + 'M'+str(gIP) + '.txt'
            f = open(NameSave,'w')
            np.savetxt(f, np.array(Obs['M'+str(gIP)]), delimiter='\t',fmt='%1.6f')   # X is an array
            f.close()
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
        f = open('VRTeCoord.txt','w')
        np.savetxt(f, Pos_vrtS[:,0:2], delimiter='\t',fmt='%1.3f')   # X is an array
        f.close()
        if dim==3:
            f = open('VRTeCoord.txt','w')
            np.savetxt(f, Pos_vrtS[:,0:3], delimiter='\t',fmt='%1.3f')   # X is an array
            f.close()
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