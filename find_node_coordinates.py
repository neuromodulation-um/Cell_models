# Function to create node coordinates, called from Cell.py


def find_devor_node_coordinates(axon_trajectory, axonName, axonType, fiberD):

    from numpy import array,sqrt,shape,sin,pi,ones,zeros,sqrt,multiply

    if axonType == 'hybrid':
        # specify deltax (i.e. node spacing) based on fiber diameter, try to work into axon class
        if (fiberD==2.0): deltax, nodeD, fiberD = 117, 1.4, 1.6
        if (fiberD==3.0): deltax, nodeD, fiberD = 309, 1.6, 2.3
        if (fiberD==5.7): deltax, nodeD, fiberD = 500, 1.9, 3.4
        if (fiberD==7.3): deltax, nodeD, fiberD = 750, 2.4, 4.6
        if (fiberD==8.7): deltax, nodeD, fiberD = 1000, 2.8, 5.8
        if (fiberD==10.0): deltax, nodeD, fiberD = 1150, 3.3, 6.9
        if (fiberD==11.5): deltax, nodeD, fiberD = 1250, 3.7, 8.1
        if (fiberD==12.8): deltax, nodeD, fiberD = 1350, 4.2, 9.2
        if (fiberD==14.0): deltax, nodeD, fiberD = 1400, 4.7, 10.4
        if (fiberD==15.0): deltax, nodeD, fiberD = 1450, 5.0, 11.5
        if (fiberD==16.0): deltax, nodeD, fiberD = 1500, 5.5, 12.7

        # ratios for normalizing variable node lengths based on "normal" internode lengths
        # from Amir and Devor 2003/Ito and Takahashi 1960
        ratio1stInternodeC = (358*1e-6)/(1450*1e-6)
        ratio2ndInternodeC = (780*1e-6)/(1450*1e-6)
        ratio3rdInternodeC = (1170*1e-6)/(1450*1e-6)
        ratio1stInternodeP = (461*1e-6)/(1567*1e-6)
        ratio2ndInternodeP = (670*1e-6)/(1567*1e-6)
        ratio3rdInternodeP = (1119*1e-6)/(1567*1e-6)

        peripheralInternodeRatios = array([ratio1stInternodeP,ratio2ndInternodeP,ratio3rdInternodeP])
        centralInternodeRatios = array([ratio1stInternodeC,ratio2ndInternodeC,ratio3rdInternodeC])

        dxDorsal = ((deltax) * 1e-6) * ones(len(axon_trajectory))
        dxPeripheral = ((deltax) * 1e-6) * ones(len(axon_trajectory))

        dxDorsal[0:3] = centralInternodeRatios * (deltax * 1e-6)
        dxPeripheral[0:3] = peripheralInternodeRatios * (deltax * 1e-6)

    elif axonType == 'devor':
        # 1e-6 scaling factors to convert microns to meters
        dxDorsal = ((1450+1.5) * 1e-6) * ones(len(axon_trajectory))
        dxPeripheral = ((1567+1.5) * 1e-6) * ones(len(axon_trajectory))
        dxStem = zeros(len(axon_trajectory))


        # 1e-6 scaling factors to convert microns to meters
        dxDorsal[0:3] = [(358+1.5)*1e-6, (780+1.5)*1e-6, (1170+1.5)*1e-6]
        dxPeripheral[0:3] = [(461+1.5)*1e-6, (670+1.5)*1e-6, (1119+1.5)*1e-6]

        # starting from junction node, moving towards soma
        dxStem = [(1.5+201)*1e-6, (1.5+168)*1e-6, (1.5+130)*1e-6, (1.5+85)*1e-6]


    ii, jj = 0, 0
    count = 0
    xx,yy,zz = list(), list(), list()
    P0 = axon_trajectory[ii,:]
    P1 = axon_trajectory[ii+1,:]
    xx.append(P0[0]); yy.append(P0[1]); zz.append(P0[2]) #define initial point of axon, i.e. first node
    zval = zz[0]
    if (axonName == 'dorsal'):
        while True:
            P = [xx[jj],yy[jj],zz[jj]] # node you're currently on
            PP1 = sqrt(sum((P1-P)**2)) # distance between current node and next point on trajectory

            # if you're on one of the first 3 nodes, your internode length is variable
            if jj <= 2:
                dx = dxDorsal[jj]
            else:
                dx = dxDorsal[3]

            if (PP1 > dx): # if the distance between current node and next point is greater than internode length
                P0P1 = sqrt(sum((P1-P0)**2)) # calculate distance between 
                tt = (P0P1-PP1+dx)/P0P1
                xx.append((1-tt)*P0[0]+tt*P1[0])
                yy.append((1-tt)*P0[1]+tt*P1[1])
                zz.append((1-tt)*P0[2]+tt*P1[2])
                jj+=1
            else:
                ii+=1                               #update P0 and P1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
        dxForReturn = dxDorsal

    elif (axonName == 'peripheral'):
        while True:
            P = [xx[jj],yy[jj],zz[jj]]
            PP1 = sqrt(sum((P1-P)**2))

            if jj <= 2:
                dx = dxPeripheral[jj]
            else:
                dx = dxPeripheral[3]

            if (PP1 > dx):
                P0P1 = sqrt(sum((P1-P0)**2))
                tt = (P0P1-PP1+dx)/P0P1
                xx.append((1-tt)*P0[0]+tt*P1[0])
                yy.append((1-tt)*P0[1]+tt*P1[1])
                zz.append((1-tt)*P0[2]+tt*P1[2])
                jj+=1
            else:
                ii+=1                               #update P0 and P1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
        dxForReturn = dxPeripheral
                
    elif (axonName == 'stem'): # note this is only for the cell model from Amir and Devor 2003
        '''
        Explicitly defined stem axon
        '''

        # second node, close to junction
        xx.append(0)
        yy.append((1.5+201)*1e-6)
        zz.append(zz[0])

        # third node
        xx.append(0)
        yy.append((1.5+201+1.5+168)*1e-6)
        zz.append(zz[0])

        # fourth node
        xx.append(0)
        yy.append((1.5+201+1.5+168+1.5+130)*1e-6)
        zz.append(zz[0])

        # fifth node, closest to soma
        xx.append(0)
        yy.append((1.5+201+1.5+168+1.5+130+1.5+85)*1e-6)
        zz.append(zz[0])

        ## hard code the initial segment and soma
        # iseg
        xx.append(0)
        yy.append((1.5+201+1.5+168+1.5+130+1.5+85+1.5+200)*1e-6)
        zz.append(zz[0])

        # soma
        xx.append(0)
        yy.append((1.5+201+1.5+168+1.5+130+1.5+85+1.5+200+80)*1e-6)
        zz.append(zz[0])

    elif (axonName == 'stemMRG'):
        from numpy import absolute

        xCoord = P0[0] + 1e-15
        yCoord = P0[1] + 1e-15
        zCoord = P0[2] + 1e-15

        # get unit vector that the middle of the arc is traveling in xy plane
        length = sqrt(xCoord**2 + yCoord**2)# + zCoord**2)
        xhat = xCoord/length
        yhat = absolute(yCoord/length)
        # xhat = 0
        # yhat = 1

        # second node
        inc = (1+201)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # third node
        inc += (1+168)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # fourth node
        inc += (1+130)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # fifth node, closest to soma
        inc += (1+85)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # iseg
        inc += (1+200)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # soma
        inc += (80)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        dxForReturn = None

    elif (axonName == 'stemMRG_adelta'):
        from numpy import absolute

        xCoord = P0[0] + 1e-15
        yCoord = P0[1] + 1e-15
        zCoord = P0[2] + 1e-15

        # get unit vector that the middle of the arc is traveling in xy plane
        length = sqrt(xCoord**2 + yCoord**2)# + zCoord**2)
        xhat = xCoord/length
        yhat = absolute(yCoord/length)
        # xhat = 0
        # yhat = 1

        # second node
        inc = (260)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # third node
        inc += (241)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # fourth node
        inc += (185)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # fifth node, closest to soma
        inc += (123)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        # soma
        inc += (29)*1e-6
        xInc = xhat*inc
        yInc = yhat*inc

        xx.append(xCoord + xInc)
        yy.append(yCoord + yInc)
        zz.append(zCoord)

        dxForReturn = None

    return array([xx,yy,zz]).T, dxForReturn
