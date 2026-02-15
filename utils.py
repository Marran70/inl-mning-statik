import numpy as np
import matplotlib.pyplot as plt

def coordxtr(edof, coords, dofs, nen=-1):
    """
    Create element coordinate matrices ex, ey, ez from edof
    coord and dofs matrices.
    
    Parameters:
    
        edof            [nel x (nen * nnd)], nnd = number of node dofs
        coords          [ncoords x ndims],   ndims = node dimensions
        dofs            [ncoords x nnd]
        
    Returns:
    
        ex              if ndims = 1
        ex, ey          if ndims = 2
        ex, ey, ez      if ndims = 3
    """

    # Create dictionary with dof indices

    dofDict = {}
    nDofs = np.size(dofs, 1)
    nElements = np.size(edof, 0)
    n_element_dofs = np.size(edof, 1)
    nDimensions = np.size(coords, 1)
    nElementDofs = np.size(edof, 1)

    if nen == -1:
        nElementNodes = int(nElementDofs/nDofs)
    else:
        nElementNodes = nen
    print(nElementNodes)
    if nElementNodes*nDofs != n_element_dofs:
        nDofs = nElementNodes*nDofs - n_element_dofs
        print(
            "dofs/edof mismatch. Using %d dofs per node when indexing." % nDofs)

    idx = 0
    for dof in dofs:
        #dofDict[dofHash(dof)] = idx
        dofDict[hash(tuple(dof[0:nDofs]))] = idx
        idx += 1

    # Loop over edof and extract element coords

    ex = np.zeros((nElements, nElementNodes))
    ey = np.zeros((nElements, nElementNodes))
    ez = np.zeros((nElements, nElementNodes))

    elementIdx = 0
    for etopo in edof:
        for i in range(nElementNodes):
            i0 = i*nDofs
            i1 = i*nDofs+nDofs-1
            dof = []
            if i0 == i1:
                dof = [etopo[i*nDofs]]
            else:
                dof = etopo[i*nDofs:(i*nDofs+nDofs)]

            nodeCoord = coords[dofDict[hash(tuple(dof[0:nDofs]))]]

            if nDimensions >= 1:
                ex[elementIdx, i] = nodeCoord[0]
            if nDimensions >= 2:
                ey[elementIdx, i] = nodeCoord[1]
            if nDimensions >= 3:
                ez[elementIdx, i] = nodeCoord[2]

        elementIdx += 1

    if nDimensions == 1:
        return ex

    if nDimensions == 2:
        return ex, ey

    if nDimensions == 3:
        return ex, ey, ez


def eldraw2(Ex, Ey, width = 1.0, color = "black"):
    """
    eldraw2(ex,ey,width=1.0,color="black")
    
     PURPOSE 
       Draw the undeformed 2D mesh for a number of bar elements

     INPUT  
        ex,ey:.......... nen:   number of element nodes
                         nel:   number of elements
        width:.......... width of the bars
        color:.......... color of the bars
    """

    nel = Ex.shape[0]
    for el in range(nel):
        plt.plot(Ex[el,:], Ey[el,:], linewidth=width, color=color)

def assem(edof, K, Ke, f=None, fe=None):
    """
    Assemble element matrices Ke ( and fe ) into the global
    stiffness matrix K ( and the global force vector f )
    according to the topology matrix edof.
    
    Parameters:
    
        edof        dof topology array
        K           the global stiffness matrix
        Ke          element stiffness matrix
        f           the global force vector
        fe          element force vector
        
    Output parameters:
    
        K           the new global stiffness matrix
        f           the new global force vector
        fe          element force vector
    
    """

    if edof.ndim == 1:
        idx = edof-1
        K[np.ix_(idx, idx)] = K[np.ix_(idx, idx)] + Ke
        if (not f is None) and (not fe is None):
            f[np.ix_(idx)] = f[np.ix_(idx)] + fe
    else:
        for row in edof:
            idx = row-1
            K[np.ix_(idx, idx)] = K[np.ix_(idx, idx)] + Ke
            if (not f is None) and (not fe is None):
                f[np.ix_(idx)] = f[np.ix_(idx)] + fe

    if f is None:
        return K
    else:
        return K, f

def solveq(K, f, bcdofs, bcvals):
    """
    Solve static FE-equations considering boundary conditions.
    
    Parameters:
    
        K           global stiffness matrix, dim(K)= nd x nd
        f           global load vector, dim(f)= nd x 1
    
        bcdofs      1-dim integer array containing prescribed dofs.
        bcvals      1-dim float array containing prescribed values.
                    If not given all prescribed dofs are assumed 0.
        
    Returns:
    
        a           solution including boundary values
    
    """

    ndofs = K.shape[0]
    print(ndofs)
    bcdofs = bcdofs - 1
    alldofs = np.arange(ndofs)
    freedofs = np.setdiff1d(alldofs, bcdofs)
    fsys = f[freedofs] - K[np.ix_(freedofs, bcdofs)] @ bcvals
    asys = np.linalg.solve(K[np.ix_(freedofs, freedofs)], fsys)
    a = np.zeros(ndofs)
    a[freedofs] = asys

    Q = K*a-f

    return (a, Q)

def extract_eldisp(edof, a):
    """
    Extract element displacements from the global displacement
    vector according to the topology matrix edof.
    
    Parameters:
    
        a           the global displacement vector
        edof        dof topology array
    
    Returns:
    
        ed:     element displacement array
    
    """

    ed = None

    if edof.ndim == 1:
        nDofs = len(edof)
        ed = np.zeros([nDofs])
        idx = edof-1
        ed[:] = a[np.ix_(idx)].T
    else:
        nElements = edof.shape[0]
        nDofs = edof.shape[1]
        ed = np.zeros([nElements, nDofs])
        i = 0
        for row in edof:
            idx = row-1
            ed[i, :] = a[np.ix_(idx)].T
            i += 1

    return ed

def eldisp2(Ex, Ey, Ed, sfac=1.0, width=1.0, color="g"):
    """
    eldisp2(Ex,Ey, Ed, sfac = 1.0, width=1.0, color="g")
    
     PURPOSE 
       Draw the deformed 2D mesh for a number of bar elements

     INPUT  
        Ex,Ey:.......... Element coordinetes
        Ed:   .......... Element displacements (obtained with e.g extract_eldisp)
        sfac: .......... Scale factor for deformation
        width:.......... width of the bars
        color:.......... color of the bars
    """
    nel = Ex.shape[0]
    for el in range(nel):
        ex =  Ex[el,:] + Ed[el,[0,2]] * sfac
        ey =  Ey[el,:] + Ed[el,[1,3]] * sfac
        plt.plot(ex, ey, linewidth=width, color=color)


def bar2s(ex, ey, ep, ed, eq=None, nep=None):
    """
    es = bar2s(ex, ey, ep, ed)
    -------------------------------------------------------------
    PURPOSE
    Compute normal force in two dimensional bar element.
    
    INPUT:  ex = [x1 x2]        element node coordinates

            ey = [y1 y2]        element node coordinates

            ep = [E A]          element properties,
                                E:  Young's modulus
                                A:  cross section area
 
            ed = [u1 ... u4]    element displacement vector 

            eq = [qX]           distributed load

            nep : number of evaluation points ( default=2 )

    OUTPUT: es = [N1 ;  section forces, local directions, in 
                  N2 ;  nep points along the beam, dim(es)= nep x 1
                  ...]  
           
            edi = [u1 ;    element displacements, local directions,
                   u2 ;    in n points along the bar, dim(edi)= nep x 1
                   ...]

            eci = [x1;     evaluation points on the local x-axis, 
                   x2;     (x1=0 and xn=L) 
                   ...] 
    -------------------------------------------------------------

    LAST MODIFIED: O Dahlblom  2015-12-04
                   O Dahlblom  2022-11-16 (Python version)
    Copyright (c)  Division of Structural Mechanics and
                   Division of Solid Mechanics.
                   Lund University
    -------------------------------------------------------------    
    """
    E, A = ep
    DEA = E*A
  
    qX = 0.
    if not eq is None:  
       qX = eq[0] 
     
    ne = 2
    if nep != None: 
       ne=nep

    x1, x2 = ex
    y1, y2 = ey
    dx = x2-x1
    dy = y2-y1
    L = np.sqrt(dx*dx+dy*dy)

    nxX = dx/L
    nyX = dy/L

    G = np.array([
        [nxX, nyX,   0,   0],  
        [  0,   0, nxX, nyX]
    ])
   
    a1 = G @ ed.reshape(4,1)

    C1 = np.array([
        [1.,      0.],
        [-1/L,   1/L]
    ]) 
    
    C1a = C1 @ a1

    X = np.arange(0., L+L/(ne-1), L/(ne-1)).reshape(ne,1) 
    zero = np.zeros(ne).reshape(ne,1)    
    one = np.ones(ne).reshape(ne,1)
  
    u = np.concatenate((one,  X), 1) @ C1a
    du = np.concatenate((zero,  one), 1) @ C1a
  
    if DEA != 0:
       u = u -(X**2-L*X)*qX/(2*DEA)
       du = du -(2*X-L)*qX/(2*DEA)
 
    N = DEA*du
    es = N
    edi = u
    eci = X

    if nep is None:
        return es[1]
    else:
        return es, edi, eci
