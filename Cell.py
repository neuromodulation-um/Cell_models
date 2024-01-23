# Code to develop the myelinated sensory neurons using the provided hoc files. For details please refer to Graham et al. Neuromodulation: Technology at the Neural Interface 24 (4), 655-671

from __future__ import division
from neuron import h
import neuron as nrn
h.load_file("stdrun.hoc")
from numpy import pi,shape,array,ones
from find_node_coordinates import find_devor_node_coordinates, find_node_coordinates

import sys

class ABetaFiber(Cell):
    '''
    Hybrid model of an ABeta Sensory Neuron. MRG Myelination, with variable node spacing near 
    t junction seen in (Ito & Takahashi 1960, Amir and Devor 2003).
    '''
    def __init__(self,**kwargs):
        self.variables = kwargs
        if 'peripheral_trajectory' not in self.variables: raise TypeError('Need to specify peripheral axon trajectory!!!')
        if 'dorsal_trajectory' not in self.variables: raise TypeError('Need to specify dorsal axon trajectory!!!')
        if 'stem_trajectory' not in self.variables: raise TypeError('Need to specify stem axon trajectory!!!')
        if 'fiberD_central' not in self.variables:  raise TypeError('Need to specify central axon diameter!!!') #um, diameter of DR fiber
        if 'fiberD_peripheral' not in self.variables:  raise TypeError('Need to specify peripheral axon diameter!!!') #um, diameter of DR fiber
        if 'fiberD_stem' not in self.variables:  raise TypeError('Need to specify stem axon diameter!!!') #um, diameter of DR fiber
        if 'pain' not in self.variables: raise TypeError('Need to specify if this is a painful fiber!!!')
        if 'CELL_DIR' not in self.variables: raise TypeError('Need to specify path to cell file!!!')
        if 'CELL_FILE_NAME' not in self.variables: raise TypeError('Need to specify the .hoc file name!!!')

        self.CELL_DIR = self.variables['CELL_DIR']
        self.CELL_FILE_NAME = self.variables['CELL_FILE_NAME']

        super(ABetaFiber, self).__init__(**kwargs)

    def __str__(self):
        return "DRG ABeta Cell"

    def get_secs(self):
        secs = [sec for sec in h.nodeP]
        secs.extend([sec for sec in h.MYSAP])
        secs.extend([sec for sec in h.FLUTP])
        secs.extend([sec for sec in h.STINP])

        secs.extend([sec for sec in h.nodeC])
        secs.extend([sec for sec in h.MYSAC])
        secs.extend([sec for sec in h.FLUTC])
        secs.extend([sec for sec in h.STINC])

        secs.extend([sec for sec in h.nodeT])
        secs.extend([sec for sec in h.MYSAT])
        secs.extend([sec for sec in h.FLUTT])
        secs.extend([sec for sec in h.STINT])
        
        secs.append(h.iseg)
        secs.append(h.soma)

        return secs

    def get_secs_in_order(self):
        nstins = self.numberOfStinCompartmentsPerStretch
        secs = []
        for ii in range(self.axonnodesP-1):
            if ii == 0:
                secs.append(h.nodeP[ii])
            secs.append(h.MYSAP[2*ii])
            secs.append(h.FLUTP[2*ii])
            for jj in range(nstins):
                secs.append(h.STINP[nstins*ii+jj])
            secs.append(h.FLUTP[2*ii+1])
            secs.append(h.MYSAP[2*ii+1])
            secs.append(h.nodeP[ii+1])
        for ii in range(self.axonnodesT-1):
            if ii == 0:
                secs.append(h.nodeT[ii])
            secs.append(h.MYSAT[2*ii])
            secs.append(h.FLUTT[2*ii])
            for jj in range(nstins):
                secs.append(h.STINT[nstins*ii+jj])
            secs.append(h.FLUTT[2*ii+1])
            secs.append(h.MYSAT[2*ii+1])
            secs.append(h.nodeT[ii+1])
        secs.append(h.iseg)
        secs.append(h.soma)
        for ii in range(self.axonnodesC-1):
            if ii == 0:
                secs.append(h.nodeC[ii])
            secs.append(h.MYSAC[2*ii])
            secs.append(h.FLUTC[2*ii])
            for jj in range(nstins):
                secs.append(h.STINC[nstins*ii+jj])
            secs.append(h.FLUTC[2*ii+1])
            secs.append(h.MYSAC[2*ii+1])
            secs.append(h.nodeC[ii+1])

        return secs

    def _construct_cell(self):
        # self.central_trajectory = self.get_variable('central_trajectory')
        self.dorsal_trajectory = self.get_variable('dorsal_trajectory')
        self.peripheral_trajectory = self.get_variable('peripheral_trajectory')
        self.stem_trajectory = self.get_variable('stem_trajectory')
        self.fiberD_central = self.get_variable('fiberD_central')
        self.fiberD_peripheral = self.get_variable('fiberD_peripheral')
        self.fiberD_stem = self.get_variable('fiberD_stem')
        self.pain = self.get_variable('pain')

        self.NODE_COORDINATES_DR, dxDorsal = find_devor_node_coordinates(self.dorsal_trajectory, 'dorsal', axonType='hybrid', fiberD=self.fiberD_central)
        self.NODE_COORDINATES_PERIPHERAL, dxPeripheral = find_devor_node_coordinates(self.peripheral_trajectory, 'peripheral', axonType='hybrid', fiberD=self.fiberD_peripheral)
        self.NODE_COORDINATES_STEM, dxDontUse = find_devor_node_coordinates(self.stem_trajectory, 'stemMRG', axonType='hybrid', fiberD=self.fiberD_stem)

        # print(self.NODE_COORDINATES_DR)
        # print(self.NODE_COORDINATES_PERIPHERAL)
        # quit()

        self.numberOfStinCompartmentsPerStretch = 6

        self.axonnodesP = shape(self.NODE_COORDINATES_PERIPHERAL)[0]
        self.axonnodesC = shape(self.NODE_COORDINATES_DR)[0]
        self.axonnodesT = 5
        self.numStemCompartments = 7

        self.paranodes1P = 2*(self.axonnodesP-1)
        self.paranodes1C = 2*(self.axonnodesC-1)
        self.paranodes1T = 2*(self.axonnodesT-1)
        self.paranodes2P = 2*(self.axonnodesP-1)
        self.paranodes2C = 2*(self.axonnodesC-1)
        self.paranodes2T = 2*(self.axonnodesT-1)
        self.axoninterP = self.numberOfStinCompartmentsPerStretch*(self.axonnodesP-1)
        self.axoninterC = self.numberOfStinCompartmentsPerStretch*(self.axonnodesC-1)
        self.axoninterT = self.numberOfStinCompartmentsPerStretch*(self.axonnodesT-1)

        self.axonnodesAbC = 4
        self.axonnodesAbP = 4
        self.paranodes1AbC = 2*(self.axonnodesAbC-1)
        self.paranodes1AbP = 2*(self.axonnodesAbP-1)
        self.paranodes2AbC = 2*(self.axonnodesAbC-1)
        self.paranodes2AbP = 2*(self.axonnodesAbP-1)
        self.axoninterAbC = self.numberOfStinCompartmentsPerStretch*(self.axonnodesAbC-1)
        self.axoninterAbP = self.numberOfStinCompartmentsPerStretch*(self.axonnodesAbP-1)
        
        # define necessary parameters
        h('fiberD_central = %.1f' % self.fiberD_central)
        h('fiberD_peripheral = %.1f' % self.fiberD_peripheral)
        h('fiberD_stem = %.1f' % self.fiberD_stem)
        h('numberOfStinCompartmentsPerStretch = %i' % self.numberOfStinCompartmentsPerStretch)
        h('axonnodesP = %i' % self.axonnodesP)
        h('axonnodesC = %i' % self.axonnodesC)
        h('paranodes1P = %i' % self.paranodes1P)
        h('paranodes1C = %i' % self.paranodes1C)
        h('paranodes2P = %i' % self.paranodes2P)
        h('paranodes2C = %i' % self.paranodes2C)
        h('axoninterP = %i' % self.axoninterP)
        h('axoninterC = %i' % self.axoninterC)
        h('axonnodesT = %i' % self.axonnodesT)
        h('paranodes1T = %i' % self.paranodes1T)
        h('paranodes2T = %i' % self.paranodes2T)
        h('axoninterT = %i' % self.axoninterT)

        h('axonnodesAbC = %i' % self.axonnodesAbC)
        h('axonnodesAbP = %i' % self.axonnodesAbP)
        h('paranodes1AbC = %i' % self.paranodes1AbC)
        h('paranodes1AbP = %i' % self.paranodes1AbP)
        h('paranodes2AbC = %i' % self.paranodes2AbC)
        h('paranodes2AbP = %i' % self.paranodes2AbP)
        h('axoninterAbC = %i' % self.axoninterAbC)
        h('axoninterAbP = %i' % self.axoninterAbP)

        if (self.pain == True):
            h('pain = 1')
        else:
            h('pain = 0')

        # load node coordinates into hoc
        h('objref nxC')
        h('objref nxP')
        h('objref nxT')
        h('objref nyC')
        h('objref nyP')
        h('objref nyT')
        h('objref nzC')
        h('objref nzP')
        h('objref nzT')
        h('nxC = new Vector(%i)' % self.axonnodesC)
        h('nxP = new Vector(%i)' % self.axonnodesP)
        h('nxT = new Vector(%i)' % self.numStemCompartments)
        h('nyC = new Vector(%i)' % self.axonnodesC)
        h('nyP = new Vector(%i)' % self.axonnodesP)
        h('nyT = new Vector(%i)' % self.numStemCompartments)
        h('nzC = new Vector(%i)' % self.axonnodesC)
        h('nzP = new Vector(%i)' % self.axonnodesP)
        h('nzT = new Vector(%i)' % self.numStemCompartments)
        for i in range(self.axonnodesC):
            h.nxC.x[i] = self.NODE_COORDINATES_DR[i,0]
            h.nyC.x[i] = self.NODE_COORDINATES_DR[i,1]
            h.nzC.x[i] = self.NODE_COORDINATES_DR[i,2]
        for i in range(self.axonnodesP):
            h.nxP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,0]
            h.nyP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,1]
            h.nzP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,2]
        for i in range(self.numStemCompartments):
            h.nxT.x[i] = self.NODE_COORDINATES_STEM[i,0]
            h.nyT.x[i] = self.NODE_COORDINATES_STEM[i,1]
            h.nzT.x[i] = self.NODE_COORDINATES_STEM[i,2]

        # account for normalized variable internode lengths to put in for STIN lengths
        h('objref varLenP')
        h('objref varLenC')
        h('varLenP = new Vector(3)')
        h('varLenC = new Vector(3)')
        for ii in range(3):
            h.varLenP.x[ii] = dxPeripheral[ii]
            h.varLenC.x[ii] = dxDorsal[ii]

        h.load_file(1, self.CELL_DIR+self.CELL_FILE_NAME)
        self.endP = h.nodeP[self.axonnodesP-1] # end of peripheral axon
        self.endC = h.nodeC[self.axonnodesC-1] # end of central axon
        self.soma = h.soma # soma
        self.TjuncStem = h.nodeT[0] # t-junction at stem axon
        self.TjuncP = h.nodeP[0] # t-junction in peripheral axon
        self.TjuncC = h.nodeC[0] # t-junction in central axon
        self.midCentral = h.nodeC[int(self.axonnodesC/2)] # midway along central axon
        self.midPeripheral = h.nodeP[int(self.axonnodesP/2)] # midway along peripheral axon

class ADeltaFiber(Cell):

    def __init__(self,**kwargs):
        self.variables = kwargs
        if 'peripheral_trajectory' not in self.variables: raise TypeError('Need to specify peripheral axon trajectory!!!')
        if 'dorsal_trajectory' not in self.variables: raise TypeError('Need to specify dorsal axon trajectory!!!')
        if 'stem_trajectory' not in self.variables: raise TypeError('Need to specify stem axon trajectory!!!')
        if 'fiberD_central' not in self.variables:  raise TypeError('Need to specify central axon diameter!!!') #um, diameter of DR fiber
        if 'fiberD_peripheral' not in self.variables:  raise TypeError('Need to specify peripheral axon diameter!!!') #um, diameter of DR fiber
        if 'fiberD_stem' not in self.variables:  raise TypeError('Need to specify stem axon diameter!!!') #um, diameter of DR fiber
        if 'pain' not in self.variables: raise TypeError('Need to specify if this is a painful fiber!!!')
        if 'CELL_DIR' not in self.variables: raise TypeError('Need to specify path to cell file!!!')
        if 'CELL_FILE_NAME' not in self.variables: raise TypeError('Need to specify the .hoc file name!!!')
        if 'variable_STIN' not in self.variables: raise TypeError('Need to specify if this fiber has variable STIN compartments! Boolean True or 1 for affirmative.')

        self.CELL_DIR = self.variables['CELL_DIR']
        self.CELL_FILE_NAME = self.variables['CELL_FILE_NAME']

        super(ADeltaFiber, self).__init__(**kwargs)

    def __str__(self):
        return "DRG ADelta Cell"

    def get_secs(self):
        secs = [sec for sec in h.nodeP]
        secs.extend([sec for sec in h.MYSAP])
        secs.extend([sec for sec in h.FLUTP])
        secs.extend([sec for sec in h.STINP])
        if (self.variable_STIN == 1) or (self.variable_STIN == True):
            secs.extend([sec for sec in h.STINPvar])

        secs.extend([sec for sec in h.nodeC])
        secs.extend([sec for sec in h.MYSAC])
        secs.extend([sec for sec in h.FLUTC])
        secs.extend([sec for sec in h.STINC])
        if (self.variable_STIN == 1) or (self.variable_STIN == True):
            secs.extend([sec for sec in h.STINCvar])

        secs.extend([sec for sec in h.nodeT])
        secs.extend([sec for sec in h.MYSAT])
        secs.extend([sec for sec in h.FLUTT])
        secs.extend([sec for sec in h.STINT])
        secs.append(h.soma)

        return secs

    def _construct_cell(self):
        # self.central_trajectory = self.get_variable('central_trajectory')
        self.dorsal_trajectory = self.get_variable('dorsal_trajectory')
        self.peripheral_trajectory = self.get_variable('peripheral_trajectory')
        self.stem_trajectory = self.get_variable('stem_trajectory')
        self.fiberD_central = self.get_variable('fiberD_central')
        self.fiberD_peripheral = self.get_variable('fiberD_peripheral')
        self.fiberD_stem = self.get_variable('fiberD_stem')
        self.pain = self.get_variable('pain')
        self.variable_STIN = self.get_variable('variable_STIN')

        self.NODE_COORDINATES_DR, dxDorsal = find_devor_node_coordinates(self.dorsal_trajectory, 'dorsal', axonType='hybrid', fiberD=self.fiberD_central)
        self.NODE_COORDINATES_PERIPHERAL, dxPeripheral = find_devor_node_coordinates(self.peripheral_trajectory, 'peripheral', axonType='hybrid', fiberD=self.fiberD_peripheral)
        self.NODE_COORDINATES_STEM, dxDontUse = find_devor_node_coordinates(self.stem_trajectory, 'stemMRG_adelta', axonType='hybrid', fiberD=self.fiberD_stem)

        self.numberOfStinCompartmentsPerStretch = 6
        self.numberOfStinCompartmentsPerVariableStretch = 1

        self.axonnodesP = shape(self.NODE_COORDINATES_PERIPHERAL)[0]
        self.axonnodesC = shape(self.NODE_COORDINATES_DR)[0]
        self.axonnodesT = 5
        self.numStemCompartments = 6

        self.paranodes1P = 2*(self.axonnodesP-1)
        self.paranodes1C = 2*(self.axonnodesC-1)
        self.paranodes1T = 2*(self.axonnodesT-1)
        self.paranodes2P = 2*(self.axonnodesP-1)
        self.paranodes2C = 2*(self.axonnodesC-1)
        self.paranodes2T = 2*(self.axonnodesT-1)
        self.axoninterT = self.numberOfStinCompartmentsPerStretch*(self.axonnodesT-1)
        
        # define necessary parameters
        h('fiberD_central = %.1f' % self.fiberD_central)
        h('fiberD_peripheral = %.1f' % self.fiberD_peripheral)
        h('fiberD_stem = %.1f' % self.fiberD_stem)
        h('numberOfStinCompartmentsPerStretch = %i' % self.numberOfStinCompartmentsPerStretch)
        h('numberOfStinCompartmentsPerVariableStretch = %i' % self.numberOfStinCompartmentsPerVariableStretch)
        h('axonnodesP = %i' % self.axonnodesP)
        h('axonnodesC = %i' % self.axonnodesC)
        h('paranodes1P = %i' % self.paranodes1P)
        h('paranodes1C = %i' % self.paranodes1C)
        h('paranodes2P = %i' % self.paranodes2P)
        h('paranodes2C = %i' % self.paranodes2C)
        h('axonnodesT = %i' % self.axonnodesT)
        h('paranodes1T = %i' % self.paranodes1T)
        h('paranodes2T = %i' % self.paranodes2T)
        h('axoninterT = %i' % self.axoninterT)

        if (self.variable_STIN == 1) or (self.variable_STIN == True):
            h('variable_STIN = 1')
            dxC = dxDorsal[-1]
            dxP = dxPeripheral[-1]
            self.numNodes20mmCentral = int(20e-3/dxC)
            self.numNodes20mmPeripheral = int(20e-3/dxP)
            self.remainingNodesCentral = self.axonnodesC - self.numNodes20mmCentral
            self.remainingNodesPeripheral = self.axonnodesP - self.numNodes20mmPeripheral

            self.axoninterP = self.numberOfStinCompartmentsPerStretch*(self.numNodes20mmPeripheral)
            self.axoninterC = self.numberOfStinCompartmentsPerStretch*(self.numNodes20mmCentral)
            self.axonremaininginterP = self.numberOfStinCompartmentsPerVariableStretch*(self.remainingNodesPeripheral-1)
            self.axonremaininginterC = self.numberOfStinCompartmentsPerVariableStretch*(self.remainingNodesCentral-1)

            h('axoninterP = %i' % self.axoninterP)
            h('axoninterC = %i' % self.axoninterC)
            h('axonremaininginterP = %i' % self.axonremaininginterP)
            h('axonremaininginterC = %i' % self.axonremaininginterC)
            h('numNodes20mmPeripheral = %i' % self.numNodes20mmPeripheral)
            h('numNodes20mmCentral = %i' % self.numNodes20mmCentral)

        else:
            h('variable_STIN = 0')
            self.axoninterP = self.numberOfStinCompartmentsPerStretch*(self.axonnodesP-1)
            self.axoninterC = self.numberOfStinCompartmentsPerStretch*(self.axonnodesC-1)
            h('axoninterP = %i' % self.axoninterP)
            h('axoninterC = %i' % self.axoninterC)


        if (self.pain == True):
            h('pain = 1')
        else:
            h('pain = 0')

        # load node coordinates into hoc
        h('objref nxC')
        h('objref nxP')
        h('objref nxT')
        h('objref nyC')
        h('objref nyP')
        h('objref nyT')
        h('objref nzC')
        h('objref nzP')
        h('objref nzT')
        h('nxC = new Vector(%i)' % self.axonnodesC)
        h('nxP = new Vector(%i)' % self.axonnodesP)
        h('nxT = new Vector(%i)' % self.numStemCompartments)
        h('nyC = new Vector(%i)' % self.axonnodesC)
        h('nyP = new Vector(%i)' % self.axonnodesP)
        h('nyT = new Vector(%i)' % self.numStemCompartments)
        h('nzC = new Vector(%i)' % self.axonnodesC)
        h('nzP = new Vector(%i)' % self.axonnodesP)
        h('nzT = new Vector(%i)' % self.numStemCompartments)
        for i in range(self.axonnodesC):
            h.nxC.x[i] = self.NODE_COORDINATES_DR[i,0]
            h.nyC.x[i] = self.NODE_COORDINATES_DR[i,1]
            h.nzC.x[i] = self.NODE_COORDINATES_DR[i,2]
        for i in range(self.axonnodesP):
            h.nxP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,0]
            h.nyP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,1]
            h.nzP.x[i] = self.NODE_COORDINATES_PERIPHERAL[i,2]
        for i in range(self.numStemCompartments):
            h.nxT.x[i] = self.NODE_COORDINATES_STEM[i,0]
            h.nyT.x[i] = self.NODE_COORDINATES_STEM[i,1]
            h.nzT.x[i] = self.NODE_COORDINATES_STEM[i,2]

        # account for normalized variable internode lengths to put in for STIN lengths
        h('objref varLenP')
        h('objref varLenC')
        h('varLenP = new Vector(3)')
        h('varLenC = new Vector(3)')
        for ii in range(3):
            h.varLenP.x[ii] = dxPeripheral[ii]
            h.varLenC.x[ii] = dxDorsal[ii]

        h.load_file(1, self.CELL_DIR+self.CELL_FILE_NAME)
        self.endP = h.nodeP[self.axonnodesP-1] # end of peripheral axon
        self.endC = h.nodeC[self.axonnodesC-1] # end of central axon
        self.soma = h.soma # soma
        self.TjuncStem = h.nodeT[0] # t-junction at stem axon
        self.TjuncP = h.nodeP[0] # t-junction in peripheral axon
        self.TjuncC = h.nodeC[0] # t-junction in central axon
        self.midCentral = h.nodeC[int(self.axonnodesC/2)] # midway along central axon
        self.midPeripheral = h.nodeP[int(self.axonnodesP/2)] # midway along peripheral axon
