import numpy as np
import tensorflow as tf

def gaussian_kernel(x1, x2, beta = 1.0):
    r = tf.transpose(x1)
    r = tf.expand_dims(r, 2)
    return tf.reduce_sum(K.exp( -beta * K.square(r - x2)), axis=-1)
def MMD(x1, x2, beta):
    x1x1 = gaussian_kernel(x1, x1, beta)
    x1x2 = gaussian_kernel(x1, x2, beta)
    x2x2 = gaussian_kernel(x2, x2, beta)
    diff = tf.reduce_mean(x1x1) - 2 * tf.reduce_mean(x1x2) + tf.reduce_mean(x2x2)
    return diff


@tf.autograph.experimental.do_not_convert
def emd_loss(y_true, y_pred, cells_to_use):
    y_true = y_true[:cells_to_use, ]
    return tf.py_function(wasserstein_distance, [y_true, y_pred], tf.float32)
def compute_cosine_distances(a, b):
    normalize_a = tf.nn.l2_normalize(a,1)        
    normalize_b = tf.nn.l2_normalize(b,1)
    distance = 1 - tf.matmul(normalize_a, normalize_b, transpose_b=True)
    return distance

@tf.autograph.experimental.do_not_convert
def val_loss_mmd2(x1, x2):
    x1 = tf.transpose(x1)
    x2 = tf.transpose(x2)
    res= MMD(x2, x1, 1)
    return res
def KL_div(P, Q):
    P = np.array(P)
    Q = np.array(Q)
    mu_P = np.mean(P, axis=0)
    mu_Q = np.mean(Q, axis=0)    
    
    # Compute their covariance
    cov_P = np.cov(P, rowvar=False)
    cov_Q = np.cov(Q, rowvar=False)    
        
    cov_Q_inv = np.linalg.inv(cov_Q)
    
    # Compute KL divergence
    KL_div = np.log(np.linalg.det(cov_Q)/np.linalg.det(cov_P)) - mu_P.shape[0] + np.trace(cov_Q_inv@cov_P) + \
                (mu_P - mu_Q).T@cov_Q_inv@(mu_P - mu_Q)
    
    KL_div = 0.5 * KL_div
    
    return KL_div/P.shape[0]
@tf.autograph.experimental.do_not_convert
def kl_in_tf(A, B, cells_to_use):
    A = A[:cells_to_use, ]
    ans = tf.numpy_function(KL_div, [A, B], tf.float64)
    return ans

def mmd2(x1, x2):
    tf.shape(x1)
    x1 = tf.transpose(x1)
    x2 = tf.transpose(x2)
    res= MMD(x2, x1, 1)
    
@tf.autograph.experimental.do_not_convert
def mmd2(x1, x2):
    x1 = tf.transpose(x1)
    x2 = tf.transpose(x2)
    res= MMD(x2, x1, 1)
    return res
@tf.autograph.experimental.do_not_convert
def d_loss(x1, x2):
    ans = tf.reduce_mean(x2, axis = 0) - tf.reduce_mean(x1, axis = 0)
    ans = tf.reduce_mean(tf.square(ans))
    return tf.abs(ans)
@tf.autograph.experimental.do_not_convert
def mmd2_1(x1, x2, cells_to_use):
    x1 = x1[:cells_to_use, ]
    x1 = tf.transpose(x1)
    x2 = tf.transpose(x2)
    res= MMD(x2, x1, 1)
    return res
@tf.autograph.experimental.do_not_convert
def d_loss_1(x1, x2, cells_to_use):
    x1 = x1[:cells_to_use, ]
    ans = tf.reduce_mean(x2, axis = 0) - tf.reduce_mean(x1, axis = 0)
    ans = tf.reduce_mean(tf.square(ans))
    return tf.abs(ans)


@tf.autograph.experimental.do_not_convert
def distance_matrix(x1, x2, cells_to_use):
    x1 = x1[cells_to_use:, ]
    ac = tf.expand_dims(x1, axis =0)
    c = tf.expand_dims(x1, axis =1)
    res_x1 = tf.norm(ac-c, axis= 2)
    #res_x1 = compute_cosine_distances(ac, c)
    ac = tf.expand_dims(x2, axis =0)
    c = tf.expand_dims(x2, axis =1)
    res_x2 = tf.norm(ac-c, axis= 2)
    #res_x2 = compute_cosine_distances(ac, c)
    ans = tf.reduce_mean(tf.square(res_x2-res_x1))
    return ans

@tf.autograph.experimental.do_not_convert
def wasserstein(x1, x2, cells_to_use):
    x1 = x1[:cells_to_use, ]
    cdf_true = K.cumsum(x1, axis=0)
    cdf_pred = K.cumsum(x2, axis=0)
    emd = K.sqrt(K.mean(K.square(cdf_true - cdf_pred), axis=0))
    return K.mean(emd)

def maximum_mean_discrepancy(x, y):
    x_kernel = tf.reduce_mean(tf.exp(-tf.square(tf.math.subtract(x, tf.expand_dims(x, 1)))))
    y_kernel = tf.reduce_mean(tf.exp(-tf.square(tf.math.subtract(y, tf.expand_dims(y, 1)))))
    xy_kernel = tf.reduce_mean(tf.exp(-tf.square(tf.math.subtract(x, tf.expand_dims(y, 1)))))
    return x_kernel + y_kernel - 2 * xy_kernel

@tf.autograph.experimental.do_not_convert
def mmd2_1(x1, x2, cells_to_use):
    x1 = x1[:cells_to_use, ]
    x1 = tf.transpose(x1)
    x2 = tf.transpose(x2)
    #res= MMD(x2, x1, 1)
    res = maximum_mean_discrepancy(x1, x2)
    return res

def sample_integers(n, shape):
    sample = tf.random_uniform(shape, minval=0, maxval=tf.cast(n, 'float32'))
    sample = tf.cast(sample, 'int32')
    return sample

def resample_rows_per_column(x):
    """Permute all rows for each column independently."""
    n_batch = tf.shape(x)[0]
    n_dim = tf.shape(x)[1]
    row_indices = sample_integers(n_batch, (n_batch * n_dim,))
    col_indices = tf.tile(tf.range(n_dim), [n_batch])
    indices = tf.transpose(tf.stack([row_indices, col_indices]))
    x_perm = tf.gather_nd(x, indices)
    x_perm = tf.reshape(x_perm, (n_batch, n_dim))
    return x_perm

def z_score(x):
    """
    Z_scores each dimension of the data (across axis 0)
    """
    #mean_vals = tf.reduce_mean(x,axis=0,keep_dims=True)
    #std_vals = tf.sqrt(tf.reduce_var(x,axis=0,keep_dims=True))
    mean_vals,var_vals = tf.nn.moments(x,axes=[0],keep_dims=True)
    std_vals = tf.sqrt(var_vals)
    x_normalized = (x - mean_vals)/std_vals
    return x_normalized

def cost_matrix(x,y,p=2):
    "Returns the cost matrix C_{ij}=|x_i - y_j|^p"
    x_col = tf.expand_dims(x,1)
    y_lin = tf.expand_dims(y,0)
    c = tf.reduce_sum((tf.abs(x_col-y_lin))**p,axis=2)
    return c

def asymmetric_loss(alpha):
    def loss(y_true, y_pred):
        delta = y_pred - y_true
        return K.mean(K.square(delta) * 
                      K.square(K.sign(delta) + alpha), 
                      axis=-1)
    return loss

def outer_sinkhorn_loss(n,niter,epsilon, p=2):
    def sinkhorn_loss(x,y):
        """
        Given two emprical measures with n points each with locations x and y
        outputs an approximation of the OT cost with regularization parameter epsilon
        niter is the max. number of steps in sinkhorn loop

        Inputs:
            x,y:  The input sets representing the empirical measures.  Each are a tensor of shape (n,D)
            epsilon:  The entropy weighting factor in the sinkhorn distance, epsilon -> 0 gets closer to the true wasserstein distance
            n:  The number of support points in the empirical measures
            niter:  The number of iterations in the sinkhorn algorithm, more iterations yields a more accurate estimate
        Outputs:

        """
        # The Sinkhorn algorithm takes as input three variables :
        C = cost_matrix(x, y,p=p)  # Wasserstein cost function

        # both marginals are fixed with equal weights
        mu = tf.constant(1.0/n,shape=[n])
        nu = tf.constant(1.0/n,shape=[n])
        # Elementary operations
        def M(u,v):
            return (-C + tf.expand_dims(u,1) + tf.expand_dims(v,0) )/epsilon
        def lse(A):
            return tf.reduce_logsumexp(A,axis=1,keepdims=True)

        # Actual Sinkhorn loop
        u, v = 0. * mu, 0. * nu
        for i in range(niter):
            u = epsilon * (tf.math.log(mu) - tf.squeeze(lse(M(u, v)) )  ) + u
            v = epsilon * (tf.math.log(nu) - tf.squeeze( lse(tf.transpose(M(u, v))) ) ) + v

        u_final,v_final = u,v
        pi = tf.exp(M(u_final,v_final))
        cost = tf.reduce_sum(pi*C)
        return cost
    return sinkhorn_loss
    
def sinkhorn_from_product(x,epsilon,n,niter,z_score=False):
    y = resample_rows_per_column(x)
    if z_score:
        x = z_score(x)
        y = z_score(y)
    return sinkhorn_loss(x,y,epsilon,n,niter)




# Constraints
def positivity(f):
    '''
    Constraint 1: 
    Ensures flow moves from source to target
    '''
    return f 

def fromSrc(f, wp, i, shape):
    """
    Constraint 2: 
    Limits supply for source according to weight
    """
    fr = np.reshape(f, shape)
    f_sumColi = np.sum(fr[i,:])
    return wp[i] - f_sumColi

def toTgt(f, wq, j, shape):
    """
    Constraint 3: 
    Limits demand for target according to weight
    """
    fr = np.reshape(f, shape)
    f_sumRowj = np.sum(fr[:,j])
    return wq[j] - f_sumRowj

def maximiseTotalFlow(f, wp, wq): 
    """
    Constraint 4: 
    Forces maximum supply to move from source to target
    """
    return f.sum() - np.minimum(wp.sum(), wq.sum())

# Objective function
def flow(f, D):
    """
    The objective function
    The flow represents the amount of goods to be moved 
    from source to target
    """
    f = np.reshape(f, D.shape)
    return (f * D).sum()

# Distance
def groundDistance(x1, x2, norm = 2):
    """
    L-norm distance
    Default norm = 2
    """
    return np.linalg.norm(x1-x2, norm)

# Distance matrix
def getDistMatrix(s1, s2, norm = 2):
    """
    Computes the distance matrix between the source
    and target distributions.
    The ground distance is using the L-norm (default L2 norm)
    """
    # Slow method
    # rows = s1 feature length
    # cols = s2 feature length
    numFeats1 = s1.shape[0]
    numFeats2 = s2.shape[0]
    distMatrix = np.zeros((numFeats1, numFeats2))

    for i in range(0, numFeats1):
        for j in range(0, numFeats2):
            distMatrix[i,j] = groundDistance(s1[i], s2[j], norm)

    # Fast method (requires scipy.spatial)
    #import scipy.spatial
    #distMatrix = scipy.spatial.distance.cdist(s1, s2)

    return distMatrix

# Flow matrix
def getFlowMatrix(P, Q, D):
    """
    Computes the flow matrix between P and Q
    """
    numFeats1 = P[0].shape[0]
    numFeats2 = Q[0].shape[0]
    shape = (numFeats1, numFeats2)

    # Constraints  
    cons1 = [{'type':'ineq', 'fun' : positivity},
             {'type':'eq', 'fun' : maximiseTotalFlow, 'args': (P[1], Q[1],)}]

    cons2 = [{'type':'ineq', 'fun' : fromSrc, 'args': (P[1], i, shape,)} for i in range(numFeats1)]
    cons3 = [{'type':'ineq', 'fun' : toTgt, 'args': (Q[1], j, shape,)} for j in range(numFeats2)]

    cons = cons1 + cons2 + cons3

    # Solve for F (solve transportation problem) 
    F_guess = np.zeros(D.shape)
    F = scipy.optimize.minimize(flow, F_guess, args=(D,), constraints=cons)
    F = np.reshape(F.x, (numFeats1,numFeats2))

    return F

# Normalised EMD
def EMD(F, D):  
    """
    EMD formula, normalised by the flow
    """
    return (F * D).sum() / F.sum()

# Runs EMD program  
def getEMD(P,Q, norm = 2):
    """
    EMD computes the Earth Mover's Distance between
    the distributions P and Q

    P and Q are of shape (2,N)

    Where the first row are the set of N features
    The second row are the corresponding set of N weights

    The norm defines the L-norm for the ground distance
    Default is the Euclidean norm (norm = 2)
    """  

    D = getDistMatrix(P[0], Q[0], norm)
    F = getFlowMatrix(P, Q, D)

    return EMD(F, D)
@tf.autograph.experimental.do_not_convert
def get_loss(pointclouds1, pointclouds2):
    loss = tf.numpy_function(getEMD, [pointclouds1,pointclouds2], tf.float64)
    return loss