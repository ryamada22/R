def sinkhorn2(a,b, M, reg, numItermax = 1000, stopThr=1e-9, verbose=False, log=False):
    """
    Solve the entropic regularization optimal transport problem
    The function solves the following optimization problem:
    .. math::
        \gamma = arg\min_\gamma <\gamma,M>_F + reg\cdot\Omega(\gamma)
        s.t. \gamma 1 = a
             \gamma^T 1= b
             \gamma\geq 0
    where :
    - M is the (ns,nt) metric cost matrix
    - :math:`\Omega` is the entropic regularization term :math:`\Omega(\gamma)=\sum_{i,j} \gamma_{i,j}\log(\gamma_{i,j})`
    - a and b are source and target weights (sum to 1)
    The algorithm used for solving the problem is the Sinkhorn-Knopp matrix scaling algorithm as proposed in [2]_
    Parameters
    ----------
    a : np.ndarray (ns,)
        samples weights in the source domain
    b : np.ndarray (nt,)
        samples in the target domain
    M : np.ndarray (ns,nt)
        loss matrix
    reg : float
        Regularization term >0
    numItermax : int, optional
        Max number of iterations
    stopThr : float, optional
        Stop threshol on error (>0)
    verbose : bool, optional
        Print information along iterations
    log : bool, optional
        record log if True
    Returns
    -------
    gamma : (ns x nt) ndarray
        Optimal transportation matrix for the given parameters
    log : dict
        log dictionary return only if log==True in parameters
    Examples
    --------
    >>> import ot
    >>> a=[.5,.5]
    >>> b=[.5,.5]
    >>> M=[[0.,1.],[1.,0.]]
    >>> ot.sinkhorn(a,b,M,1)
    array([[ 0.36552929,  0.13447071],
           [ 0.13447071,  0.36552929]])
    References
    ----------
    .. [2] M. Cuturi, Sinkhorn Distances : Lightspeed Computation of Optimal Transport, Advances in Neural Information Processing Systems (NIPS) 26, 2013
    See Also
    --------
    ot.lp.emd : Unregularized OT
    ot.optim.cg : General regularized OT
    """

    a=np.asarray(a,dtype=np.float64)
    b=np.asarray(b,dtype=np.float64)
    M=np.asarray(M,dtype=np.float64)

    if len(a)==0:
        a=np.ones((M.shape[0],),dtype=np.float64)/M.shape[0]
    if len(b)==0:
        b=np.ones((M.shape[1],),dtype=np.float64)/M.shape[1]

    # init data
    Nini = len(a)
    Nfin = len(b)


    cpt = 0
    if log:
        log={'err':[]}

    # we assume that no distances are null except those of the diagonal of distances
    u = np.ones(Nini)/Nini
    v = np.ones(Nfin)/Nfin
    uprev=np.zeros(Nini)
    vprev=np.zeros(Nini)

    #print reg

    K = np.exp(-M/reg)
    #print np.min(K)

    Kp = np.dot(np.diag(1/a),K)
    transp = K
    cpt = 0
    err=1
    while (err>stopThr and cpt<numItermax):
        if np.any(np.dot(K.T,u)==0) or np.any(np.isnan(u)) or np.any(np.isnan(v)):
            # we have reached the machine precision
            # come back to previous solution and quit loop
            print('Warning: numerical errrors')
            if cpt!=0:
                u = uprev
                v = vprev
            break
        uprev = u
        vprev = v
        v = np.divide(b,np.dot(K.T,u))
        #u = 1./np.dot(Kp,v)
        u = np.divide(a,np.dot(K,v))
        if cpt%10==0:
            # we can speed up the process by checking for the error only all the 10th iterations
            transp = np.dot(np.diag(u),np.dot(K,np.diag(v)))
            err = np.linalg.norm((np.sum(transp,axis=0)-b))**2
            if log:
                log['err'].append(err)

            if verbose:
                if cpt%200 ==0:
                    print('{:5s}|{:12s}'.format('It.','Err')+'\n'+'-'*19)
                print('{:5d}|{:8e}|'.format(cpt,err))
        cpt = cpt +1
    if log:
        log['u']=u
        log['v']=v
          
    #print 'err=',err,' cpt=',cpt
    if log:
        return u.reshape((-1,1))*K*v.reshape((1,-1)),log
    else:
        return u.reshape((-1,1))*K*v.reshape((1,-1))