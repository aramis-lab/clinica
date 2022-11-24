def atlas_decomposition(dartel_input):
    """

    :param dartel_input: path to the dartel input
    :return: 3 atlases (gm, wm, csf)
    """

    import nibabel as nib

    dartel = nib.load(dartel_input)
    dartel = dartel.get_fdata(dtype="float32")
    atlas_1 = dartel[:, :, :, 0]
    atlas_2 = dartel[:, :, :, 1]
    atlas_3 = dartel[:, :, :, 2]
    atlas = [atlas_1, atlas_2, atlas_3]

    return atlas


def spm_read(fname):
    """
     Read the image and the header of fname
    :param fname: name of the image
    :return: it returns the image data as an array nibabel and the header
    """

    import nibabel as nib
    import numpy as np

    img = nib.load(fname)
    pico = img.get_fdata(dtype="float32")
    pico = np.array(pico, dtype="float32")
    mask = np.isnan(pico)
    pico[mask] = 0
    volu = img.header

    return [pico, volu]


def spm_write_vol(fname, regularized_features):
    """

    :param fname: name of the atlas image, necessary to take information about the affine matrix and the header
    :param regularized_features: new image data obtained with the regularization of the input image
    :return: new image data in NIFTI format
    """

    import nibabel as nib

    i = nib.load(fname)
    data = regularized_features
    img = nib.Nifti1Image(data, header=i.header, affine=i.affine)
    return img


def rescaleImage(image1, p):
    """
    Normalization of the histogram of intensity
    :param Image1: input image
    :param p: vector of minimum and maximum value for the normalization
    how the istogram is normalized:
    between [0 1] if there are no options
    between [1 p] if len(p) ==1
    between [p[0] p[1]] if len(p) == 2
    :return: image with the histogram normalized

    """
    import numpy as np

    eps = 2.2204e-16
    p = np.array(p)
    m = image1.min()
    M = (image1 - m).max()

    image2 = (image1 - m) / (M + eps)

    if len(p) == 1:
        image2 = image2 * (p - 1) + 1
    elif len(p) == 2:
        image2 = image2 * (p[1] - p[0]) + p[0]

    return image2


def tensor_scalar_product(sc, g1):
    """

    :param sc: scalar
    :param g1: 3 * 3 tensor
    :return: product between the tensor and the scalar
    """
    import numpy as np

    # we define the scalar and the tensor as complex

    sc = np.array(sc, dtype=np.complex128)
    g1 = np.array(g1, dtype=np.complex128)
    g = np.zeros(g1.shape, dtype=np.complex128)  # new vector

    for i in range(g1.shape[0]):
        for j in range(g1.shape[1]):
            for k in range(g1.shape[0]):
                g[i][j] = g1[i][j] * sc

    # g is the final tensor
    return g


def tensor_eye(atlas):
    """

    :param atlas: list of atlases
    :return: the identity matrix of a tensor
    """
    import numpy as np

    a = np.ones(atlas[0].shape)
    b = np.zeros(atlas[0].shape)
    # the new matrix is 1 on the diagonal and 0 in the other positions
    X = [a, b, b]
    Y = [b, a, b]
    Z = [b, b, a]

    matrix = [X, Y, Z]
    # it returns matrix, combination of array
    return matrix


def tensor_sum(g1, g2):
    """

    :param g1: first tensor
    :param g2: second tensor
    :return: sum of the two tensor
    """
    import numpy as np

    g = np.add(g1, g2)
    # numpy add to sum the tensors
    return g


def tensor_product(g1, g2):
    """

    :param g1: first tensor
    :param g2: second tensor
    :return: product between the two tensors
    """
    import numpy as np

    g1 = np.array(g1)
    g2 = np.array(g2)

    g = np.zeros(g1.shape)

    for i in range(g1.shape[0]):
        for j in range(g2.shape[0]):
            g[i][j] = 0
            for k in range(g1.shape[0]):
                g[i][j] = g[i][j] + np.multiply(g1[i][k], g2[k][j])

    # g = g1 * g2 (dim of the tensor: 3*3*xg*yg*zg)
    return g


def tensor_determinant(g):
    """

    :param g: tensor dim = 3*3*xg*yg*z
    :return: determinant of the tensor dim = xg*yg*zg
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    g = np.array(g)
    d = 0
    s = g.shape
    # recursive function we have a matrix of 3*3 and then we divide it in different blocks of 2*2 and we calculate the determinant
    # of them making a um of the different determinants. The resulting determinant it's a sum of the different blocks.
    if s[0] == 3:
        # if the tensor is 3*3

        for i in range(s[0]):

            if np.mod(i, 2) == 0:
                epsilon = 1
            else:
                epsilon = -1

            if i == 0:
                g1 = [[g[1][1], g[1][2]], [g[2][1], g[2][2]]]
            elif i == 1:
                g1 = [[g[0][1], g[0][2]], [g[2][1], g[2][2]]]
            else:
                g1 = [[g[0][1], g[0][2]], [g[1][1], g[1][2]]]
            # it's a recursive function
            prod = epsilon * g[i][0] * utils.tensor_determinant(g1)
            d = d + prod

    elif s[0] == 2:
        # if the tensor is 2*2
        for i in range(s[0]):
            if np.mod(i, 2) == 0:
                epsilon = 1
            else:
                epsilon = -1
            if i == 0:
                g1 = [g[1][1]]
            elif i == 1:
                g1 = [g[0][1]]
            prod = epsilon * g[i][0] * utils.tensor_determinant(g1)
            d = d + prod

    elif s[0] == 1:
        # if the tensor is 1 matrix
        d = [g[0]]
    return d


def tensor_trace(g):
    """

    :param g: tensor
    :return: trace of a tensor
    """
    import numpy as np

    trace = np.trace(g)
    # trace if the trace of the input tensor
    return trace


def roots_poly(C):
    """

    :param C: coefficients. If C has N+1 components, the polynomial is C(1)*X^N + ... + C(N) * X + C(N+1)
    :return: roots of the polynomial

    the functions find the polynomial roots. It computes the roots of the polynomail whose coefficients are the elements of the vector C?
    """
    import cmath
    import math

    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    C = np.array(C)
    if C.shape[0] < 2:
        rts = []

    elif C.shape[0] < 3:
        rts = -C[:, 1] * (1 / C[:, 0])

    elif C.shape[0] < 4:
        # implementation of delta
        delta = np.array(
            [
                cmath.sqrt((C[1, i] * C[1, i]) - (4 * C[0, i] * C[2, i]))
                for i in range(C.shape[1])
            ]
        )
        # two roots
        rts1 = (-C[1, :] + delta) * (1 / ((2 * C[0, :])))
        rts2 = (-C[1, :] - delta) * (1 / ((2 * C[0, :])))
        rts = np.array([rts1, rts2])

    elif C.shape[0] < 5:
        # implementation of the method of Cardan

        a = C[0, :]
        b = C[1, :]
        c = C[2, :]
        d = C[3, :]

        p = -b * b * (1 / (3 * a * a)) + c * (1 / a)
        q = b * (1 / (27.0 * a)) * (2 * b * b * (1 / (a * a)) - 9 * c * (1 / a)) + d * (
            1 / a
        )

        new_roots = np.array([np.ones((q.shape[0])), q, -(p * p * p) / 27])

        rts = utils.roots_poly(new_roots)

        u = rts[0, :]
        u_mod = abs(u) ** (1 / 3)
        u_angle = np.angle(u) * (1 / 3)

        v = rts[1, :]
        v_mod = abs(v) ** (1 / 3)
        v_angle = np.angle(v) * (1 / 3)

        rts = np.zeros([C.shape[1], 3], dtype="complex128")

        ind = np.zeros([u.shape[0]], dtype="int32")

        for k in [0, 1, 2]:
            u = u_mod * np.power(math.e, 1j * (u_angle + k * 2 * math.pi * (1 / 3)))
            u = np.array(u)

            for r in [0, 1, 2]:
                v = v_mod * np.power(math.e, 1j * (v_angle + r * 2 * math.pi * (1 / 3)))

                ind2 = abs(u * v + p * 1 / 3) < 1e-10

                rts[ind2, ind[ind2]] = u[ind2] + v[ind2] - b[ind2] * (1 / (3 * a[ind2]))

                ind[ind2] = np.minimum(2, ind[ind2] + 1)

    else:
        print("For degree > 3 use roots ")

    return rts


def tensor_eigenvalues(g):
    """

    :param g: tensor
    :return: eigenvalues of the tensor

    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    g = np.array(g)

    if g.shape[0] < 4:
        # condition if we have a tensor

        C1 = np.ones(len(np.ravel(g[0][0])))
        buff = -utils.tensor_trace(g)
        C2 = buff.flatten("F")
        buff = utils.tensor_trace(utils.tensor_product(g, g))
        buff = -0.5 * (buff.flatten("F") - np.multiply(C2, C2))
        C3 = buff.flatten("F")
        buff = -utils.tensor_determinant(g)
        C4 = buff.flatten("F")

        C = np.array([C1, C2, C3, C4])
        rts = utils.roots_poly(C)

    else:
        print("Degree too big : not still implemented")

    rts2 = rts.real.copy()
    rts2.sort()

    lamb = np.zeros(
        shape=(g.shape[0], g.shape[2], g.shape[3], g.shape[4]), dtype="complex128"
    )

    for i in range(g.shape[0]):
        lamb[i, :, :, :] = rts2[:, i].reshape(
            g.shape[2], g.shape[3], g.shape[4], order="F"
        )

    # lamb[0] is the smallest eigenvalues, lamb[2] is the biggest
    return lamb


def tensor_transpose(g):
    """

    :param g: tensor
    :return: tranpose of the tensor
    """
    import numpy as np

    g = np.array(g)
    tg = np.array(g)
    for i in range(g.shape[0]):
        for j in range(g.shape[0]):
            tg[i][j] = g[j][i]
    # tg is the transposed tensor
    return tg


def tensor_commatrix(g):
    """

    :param g: tensor
    :return: commatrix of the tensor
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    g = np.array(g)
    g_com = []

    for i in range(g.shape[0]):
        for j in range(g.shape[0]):
            if np.mod(i + j, 2) == 0:
                epsilon = 1
            else:
                epsilon = -1
            if i == 0:
                if j == 0:
                    n = [[g[1][1], g[1][2]], [g[2][1], g[2][2]]]
                    a0 = epsilon * utils.tensor_determinant(n)
                elif j == 1:
                    n = [[g[1][0], g[1][2]], [g[2][0], g[2][2]]]
                    a1 = epsilon * utils.tensor_determinant(n)
                else:
                    n = [[g[1][0], g[1][1]], [g[2][0], g[2][1]]]
                    a2 = epsilon * utils.tensor_determinant(n)
            elif i == 1:
                if j == 0:
                    n = [[g[0][1], g[0][2]], [g[2][1], g[2][2]]]
                    b0 = epsilon * utils.tensor_determinant(n)
                elif j == 1:
                    n = [[g[0][0], g[0][2]], [g[2][0], g[2][2]]]
                    b1 = epsilon * utils.tensor_determinant(n)
                else:
                    n = [[g[0][0], g[0][1]], [g[2][0], g[2][1]]]
                    b2 = epsilon * utils.tensor_determinant(n)
            else:
                if j == 0:
                    n = [[g[0][1], g[0][2]], [g[1][1], g[1][2]]]
                    c0 = epsilon * utils.tensor_determinant(n)
                elif j == 1:
                    n = [[g[0][0], g[0][2]], [g[1][0], g[1][2]]]
                    c1 = epsilon * utils.tensor_determinant(n)
                else:
                    n = [[g[0][0], g[0][1]], [g[1][0], g[1][1]]]
                    c2 = epsilon * utils.tensor_determinant(n)
    g_com = [[a0, a1, a2], [b0, b1, b2], [c0, c1, c2]]
    g_com = np.array(g_com)
    if len(g_com.shape) == 6:
        g_com = g_com[:, :, 0, :, :, :]

    return g_com


def create_fisher_tensor(atlas):
    """

    :param atlas: list of 3 atlases, the 3 probability maps from the template with 3 components
    :return: g = tensor
    """

    # create tensor for fisher metrics
    import numpy as np

    upper_bound = 0.999  # probabibilty limits to avoid log(0) and log(1)
    lower_bound = 0.001  # probability limits to avoid log(0) and log(1)
    epsilon = 1e-6  # regularization

    n = atlas[0].shape

    a = np.ones(atlas[0].shape)
    b = np.zeros(atlas[0].shape)
    X = [np.dot(a, epsilon), b, b]
    Y = [b, np.dot(a, epsilon), b]
    Z = [b, b, np.dot(a, epsilon)]

    g = [X, Y, Z]

    for i in range(3):  # for for each component of the tensor
        proba = 1
        proba = proba * atlas[i]
        proba = np.maximum(np.minimum(proba, upper_bound), lower_bound)
        gr = np.array(np.gradient(np.log(proba)))

        for x in range(3):
            for y in range(3):
                g[x][y] = g[x][y] + (proba * gr[x] * gr[y])

    return g


def tensor_helmholtz(x, h, detg, k):
    """

    :param x: 3D Array
    :param h: sqrt(det(g)) * inverse(g)^2 -> g is he metric tensor of the 3D manifold M
    :param detg: sqrt(det(g))
    :param k: constant( 0 - gives laplacian)
    :return:

    """

    import numpy as np

    detg = np.array(detg)
    if len(detg.shape) == 3:
        detg_ = detg
    else:
        detg_ = detg[0]

    weight = detg_[1:-1, 1:-1, 1:-1] * k
    if len(h.shape) == 6:
        h_ = h[
            :,
            :,
            0,
            :,
            :,
            :,
        ]
    else:
        h_ = h
    for i in range(len(h)):  # from 1 to 3
        mat = h_[i][i]
        weight = weight + mat[1:-1, 1:-1, 1:-1]

    mat_1 = h_[0][0]
    weight = weight + 0.5 * mat_1[:-2, 1:-1, 1:-1]
    weight = weight + 0.5 * mat_1[2:, 1:-1, 1:-1]

    mat_2 = h_[1][1]
    weight = weight + 0.5 * mat_2[1:-1, :-2, 1:-1]
    weight = weight + 0.5 * mat_2[1:-1, 2:, 1:-1]

    mat_3 = h_[2][2]
    weight = weight + 0.5 * mat_3[1:-1, 1:-1, :-2]
    weight = weight + 0.5 * mat_3[1:-1, 1:-1, 2:]

    y0 = weight * x[1:-1, 1:-1, 1:-1]

    mat1 = h_[0][1]
    mat2 = h_[1][0]
    mat3 = h_[1][2]
    mat4 = h_[2][1]
    mat5 = h_[0][2]
    mat6 = h_[2][0]

    y0 = (
        y0
        + (-1 * -1 * -0.25)
        * (mat1[:-2, 1:-1, 1:-1] + mat2[1:-1, :-2, 1:-1])
        * x[:-2, :-2, 1:-1]
    )
    y0 = (
        y0
        + (-1 * -1 * -0.25)
        * (mat3[1:-1, :-2, 1:-1] + mat4[1:-1, 1:-1, :-2])
        * x[1:-1, :-2, :-2]
    )
    y0 = (
        y0
        + (-1 * -1 * -0.25)
        * (mat5[:-2, 1:-1, 1:-1] + mat6[1:-1, 1:-1, :-2])
        * x[:-2, 1:-1, :-2]
    )
    y0 = (
        y0
        + (-1 * +1 * -0.25)
        * (mat1[:-2, 1:-1, 1:-1] + mat2[1:-1, 2:, 1:-1])
        * x[:-2, 2:, 1:-1]
    )
    y0 = (
        y0
        + (-1 * +1 * -0.25)
        * (mat3[1:-1, :-2, 1:-1] + mat4[1:-1, 1:-1, 2:])
        * x[1:-1, :-2, 2:]
    )
    y0 = (
        y0
        + (-1 * +1 * -0.25)
        * (mat5[:-2, 1:-1, 1:-1] + mat6[1:-1, 1:-1, 2:])
        * x[:-2, 1:-1, 2:]
    )
    y0 = (
        y0
        + (+1 * -1 * -0.25)
        * (mat1[2:, 1:-1, 1:-1] + mat2[1:-1, :-2, 1:-1])
        * x[2:, :-2, 1:-1]
    )
    y0 = (
        y0
        + (+1 * -1 * -0.25)
        * (mat3[1:-1, 2:, 1:-1] + mat4[1:-1, 1:-1, :-2])
        * x[1:-1, 2:, :-2]
    )
    y0 = (
        y0
        + (+1 * -1 * -0.25)
        * (mat5[2:, 1:-1, 1:-1] + mat6[1:-1, 1:-1, :-2])
        * x[2:, 1:-1, :-2]
    )
    y0 = (
        y0
        + (+1 * +1 * -0.25)
        * (mat1[2:, 1:-1, 1:-1] + mat2[1:-1, 2:, 1:-1])
        * x[2:, 2:, 1:-1]
    )
    y0 = (
        y0
        + (+1 * +1 * -0.25)
        * (mat3[1:-1, 2:, 1:-1] + mat4[1:-1, 1:-1, 2:])
        * x[1:-1, 2:, 2:]
    )
    y0 = (
        y0
        + (+1 * +1 * -0.25)
        * (mat5[2:, 1:-1, 1:-1] + mat6[1:-1, 1:-1, 2:])
        * x[2:, 1:-1, 2:]
    )

    y0 = (
        y0
        + (-0.5)
        * (mat_1[1:-1, 1:-1, 1:-1] + mat_1[:-2, 1:-1, 1:-1])
        * x[:-2, 1:-1, 1:-1]
    )
    y0 = (
        y0
        + (-0.5) * (mat_1[1:-1, 1:-1, 1:-1] + mat_1[2:, 1:-1, 1:-1]) * x[2:, 1:-1, 1:-1]
    )
    y0 = (
        y0
        + (-0.5)
        * (mat_2[1:-1, 1:-1, 1:-1] + mat_2[1:-1, :-2, 1:-1])
        * x[1:-1, :-2, 1:-1]
    )
    y0 = (
        y0
        + (-0.5) * (mat_2[1:-1, 1:-1, 1:-1] + mat_2[1:-1, 2:, 1:-1]) * x[1:-1, 2:, 1:-1]
    )
    y0 = (
        y0
        + (-0.5)
        * (mat_3[1:-1, 1:-1, 1:-1] + mat_3[1:-1, 1:-1, :-2])
        * x[1:-1, 1:-1, :-2]
    )
    y0 = (
        y0
        + (-0.5) * (mat_3[1:-1, 1:-1, 1:-1] + mat_3[1:-1, 1:-1, 2:]) * x[1:-1, 1:-1, 2:]
    )

    return y0


def tensor_inverse(g):
    """

    :param g: tensor
    :return: inverse of the tensor
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    h = utils.tensor_transpose(utils.tensor_commatrix(g))
    detg = utils.tensor_determinant(g)

    h = h * (1 / (detg))
    mask = np.isnan(h)
    h[mask] = 0
    return h


def operateur(x, ginv, detg):
    """

    :param x:
    :param ginv:
    :param detg:
    :return:
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    if len(x.shape) == 4:
        x = x[0, :, :, :]
    y = np.zeros([x.shape[0] + 2, x.shape[1] + 2, x.shape[2] + 2])
    y = np.array(y, dtype=np.complex_)
    y[1:-1, 1:-1, 1:-1] = x
    y = utils.tensor_helmholtz(y, ginv, detg, 0)

    return y


def largest_eigenvalue_heat_3D_tensor2(g, h, epsilon: float = 1e-6):
    """

    :param g: metric tensor
    :param h: space step
    :param epsilon: stop criterion (default: 1e-6)
    :return: lamba = the largest eigenvalues

    """
    import cmath

    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # parameters
    erreur = 1 + epsilon

    # tensors

    detg = utils.tensor_determinant(g)
    detg = np.array(detg, dtype=np.complex128)  # complex tensor
    detg = np.sqrt(detg)
    detg = detg[0]
    ginv = utils.tensor_inverse(g)

    if len(ginv.shape) == 6:
        ginv = ginv[:, :, 0, :, :, :]

    ginv = utils.tensor_scalar_product(detg, ginv)
    detg2 = detg[1:-1, 1:-1, 1:-1]  # 141*121*141
    detg2[np.isnan(detg2)] = 0
    detg[np.isnan(detg)] = 0
    ginv[np.isnan(ginv)] = 0

    # initialisation

    s = [g[0][0].shape[0] - 2, g[0][0].shape[1] - 2, g[0][0].shape[2] - 2]
    b1 = np.ones([s[0], s[1], s[2]])

    b1 = np.divide(
        b1,
        np.array(
            cmath.sqrt(np.dot(b1.flatten("F").transpose(), b1.flatten("F"))),
            dtype=np.complex128,
        ),
    )

    print("Computation of the largest eigenvalue ...")
    while erreur > epsilon:
        b0 = b1
        b2 = np.array(
            np.divide(np.array(utils.operateur(b1, ginv, detg)) * h, detg2) / h / h / h,
            dtype=np.complex128,
        )
        b1 = np.divide(
            b2,
            np.array(cmath.sqrt(np.dot(b2.flatten("F").transpose(), b2.flatten("F")))),
            dtype=np.complex128,
        )

        erreur = np.linalg.norm(b1.flatten("F") - b0.flatten("F"))

    print("done")

    lam = cmath.sqrt(np.dot(b2.flatten("F").transpose(), b2.flatten("F")))

    return lam


def heat_finite_elt_3D_tensor2(x0, t_final, t_step, h, g):
    """

    :param x0: vector x (at t = 0)
    :param t_final: time
    :param t_step: time step (must satisfy the CFL max(lambda) < 2)
    :param h:
    :param g: metric tensor
    :return: vector x (at t = t_final)

    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    if len(x0.shape) == 4:
        x0 = x0[0, :, :, :]

    # parameters
    nb_step = np.ceil(t_final / t_step)  # number of time step
    nb_step = nb_step.astype(int)
    t_step = t_final / nb_step

    # tensors
    detg = utils.tensor_determinant(g)
    detg = np.sqrt(detg)
    ginv = utils.tensor_inverse(g)
    ginv = utils.tensor_scalar_product(detg, ginv)
    detg2 = detg[:, 1:-1, 1:-1, 1:-1]
    if len(ginv.shape) == 6:
        ginv = ginv[:, :, 0, :, :, :]
    if len(detg.shape) == 4:
        detg = detg[0, :, :, :]
    ginv = np.array(ginv.real, dtype="float64")
    detg = np.array(detg.real, dtype="float64")
    detg2 = np.array(detg2.real, dtype="float64")

    # LOOP
    x = x0
    for i in range(nb_step):
        x = np.array(
            x
            - t_step
            * (np.divide(np.array(utils.operateur(x, ginv, detg)) * h, detg2))
            / h
            / h
            / h
        )

    return x


def heat_finite_elt_2D_tensor2(x0, t_final, t_step, h, g):
    """

    :param x0: vector x (at t = 0)
    :param t_final: time
    :param t_step: time step (must satisfy the CFL max(lambda) < 2)
    :param h:
    :param g: metric tensor
    :return: vector x (at t = t_final)

    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # parameters
    nb_step = np.ceil(t_final / t_step)  # number of time step
    t_step = t_final / nb_step

    # tensors
    detg = utils.tensor_determinant(g)
    detg = np.sqrt(detg)
    ginv = utils.tensor_inverse(g)
    ginv = utils.tensor_scalar_product(detg, ginv)
    detg2 = detg[1:-1, 1:-1]

    # LOOP
    x = x0
    for i in range(nb_step):
        m = t_step / h / h
        x = np.sum(
            x,
            -(
                utils.tensor_scalar_product(
                    m,
                    np.divide(
                        utils.tensor_scalar_product(h, utils.operateur(x, ginv, detg)),
                        detg2,
                        dtype=object,
                    ),
                )
            ),
        )

    return x


def heat_solver_tensor_3D_P1_grad_conj(
    f, g, t_final, h, t_step, CL_value=None, epsilon: float = 0.1
):
    """
    It solves the poisson's equation in 1D on the regular mesh (with mesh of size h)
    :param f: approximation of a function of L^(/Omega)
    :param g: tensor
    :param t_final:
    :param h:
    :param t_step:
    :param CL_value:
    :param epsilon:
    :return: u= solution of the poisson's equation
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # initialisation
    h = h or 1
    CL_value = CL_value or np.zeros(f.shape)

    # rigidity matrix
    b_h = f[1:-1, 1:-1, 1:-1] * (h * h * h)

    b_h[:, :, 0] = b_h[:, :, 0] + (
        CL_value[1:-1, 1:-1, 0] * h
    )  # not sure about b_h third value is 0 -> I need to avoid the column (HOW??)
    b_h[:, 0, :] = b_h[:, 0, :] + (CL_value[1:-1, 0, 1:-1] * h)
    b_h[0, :, :] = b_h[0, :, :] + (CL_value[0, 1:-1, 1:-1] * h)

    print("##########computation b_H#############@ ")

    # inversion of the linear system
    U_h = utils.heat_finite_elt_3D_tensor2(b_h, t_final, t_step, h, g)

    u = CL_value
    u[1:-1, 1:-1, 1:-1] = U_h

    return u


def heat_solver_tensor_2D_P1_grad_conj(
    f, g, t_final, h, t_step, CL_value=None, epsilon: float = 1e-4
):
    """
    It solves the poisson's equation in 1D on the regular mesh (with mesh of size h)
    :param f: approximation of a function of L^(/Omega)
    :param g: tensor
    :param t_final:
    :param h:
    :param t_step:
    :param CL_value:
    :param epsilon:
    :return: u= solution of the poisson's equation
    """
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # intiialisation
    h = h or 1
    CL_value = CL_value or np.zeros(f.shape)

    # rigidity matrix
    b_h = utils.tensor_scalar_product((h * h), f[1:-1, 1:-1])
    b_h[:, 0] = b_h[:, 0] + utils.tensor_scalar_product(h, CL_value[1:-1, 0])
    b_h[0, :] = b_h[0, :] + utils.tensor_scalar_product(h, CL_value[0, 1:-1])

    # inversion of the linear system
    U_h = utils.heat_finite_elt_2D_tensor2(b_h, t_final, t_step, h, g)

    u = CL_value
    u[1:-1, 1:-1] = U_h

    return u


def obtain_g_fisher_tensor(dartel_input, FWHM):
    """
    heat regularization based on the Fisher metric
    :param dartel_input: dartel template in MNI space
    :param sigma_loc: 10
    :param h: voxel size 1,5
    :param FWHM: mm of smoothing, parameters choosing by the user. default_value = 4
    :return: g: fisher tensor

    """

    import math
    import os

    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    #
    # PARAMETERS

    sigma_loc = 10
    error_tol = 0.001  # error for the estimation of the largest eigenvalue
    alpha_time = 0.9  # time_step = alpha_time * (time_step_max)
    max_proba = 0.999  # proba must be > 0 & < 1
    min_proba = 0.001
    h = 1.5  # voxel size

    # PARSE INPUTS/INIT
    sigma = FWHM / (2 * math.sqrt(2 * math.log(2)))  # sigma of voxels
    beta = sigma**2 / 2

    # SCALE MAPS
    xxx = []

    atlas = utils.atlas_decomposition(dartel_input)

    for i in atlas:

        image = utils.rescaleImage(i, [min_proba, max_proba])

        xxx.append(image)

    atlas = xxx
    si = atlas[0].shape

    # CREATE TENSOR
    g_atlas = utils.create_fisher_tensor(atlas)
    g_atlas = utils.tensor_scalar_product(h * h, g_atlas)
    g_pos = utils.tensor_eye(atlas)
    g_pos = utils.tensor_scalar_product(1 / float(sigma_loc**2), g_pos)

    g = utils.tensor_sum(g_atlas, g_pos)

    print("computing mean distance ... ")

    eigenv = utils.tensor_eigenvalues(g)

    print("done")

    dist_av = []

    for i in range(g.shape[0]):
        dist_av.append(np.sqrt(abs(eigenv[i])))
    dist_av = np.mean(dist_av)

    print("average distance ", dist_av)

    g = utils.tensor_scalar_product((1 / dist_av) / dist_av, g)

    np.save(os.path.abspath("./output_fisher_tensor.npy"), g)

    return g, os.path.abspath("./output_fisher_tensor.npy")


def obtain_time_step_estimation(dartel_input, FWHM, g):
    """

    :param h: 1,5 voxel size
    :param FWHM: mm of smoothing, defined by the user, default value = 4
    :param g: fisher tensor
    :return:
    """
    import json
    import math
    import os

    import nibabel as nib
    import numpy as np

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # obtain voxel size with dartel_input
    head = nib.load(dartel_input)
    head_ = head.header
    for i in range(len(head_["pixdim"])):
        if head_["pixdim"][i] > 0:
            h = head_["pixdim"][i]

    error_tol = 0.001  # error for the estimation of the largest eigenvalue
    alpha_time = 0.9  # time_step = alpha_time * (time_step_max)
    sigma = FWHM / (2 * math.sqrt(2 * math.log(2)))  # sigma of voxels
    beta = sigma**2 / 2

    lam = utils.largest_eigenvalue_heat_3D_tensor2(g, h, error_tol)
    print("lambda: ", lam)
    lam = np.array(lam.real, dtype="float64")
    t_step_max = 2 / lam

    t_step = alpha_time * t_step_max
    nbiter = np.ceil(beta / t_step)
    t_step = beta / nbiter

    # after t_step calculation: creation of json file
    data = {
        "MaxDeltaT": 0.0025,
        "Alpha": 0.9,
        "Epsilon": 1e-6,
        "BoundaryConditions": "TimeInvariant",
        "SigmaLoc": 10,
        "TimeStepMax": t_step_max,
        "SpatialPrior": "Tissues (GM, WM, CSF)",
        "RegularizationType": "Fisher",
        "FWHM": FWHM,
    }

    json_data = json.dumps(data)
    with open("./output_data.json", "w") as f:
        f.write(json_data)

    return t_step, os.path.abspath("./output_data.json")


def heat_solver_equation(input_image, g, FWHM, t_step, dartel_input):
    import math
    import os

    import nibabel as nib

    import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

    # obtain voxel size with dartel_input
    head = nib.load(dartel_input)
    head_ = head.header
    for i in range(len(head_["pixdim"])):
        if head_["pixdim"][i] > 0:
            h = head_["pixdim"][i]

    sigma = FWHM / (2 * math.sqrt(2 * math.log(2)))  # sigma of voxels
    beta = sigma**2 / 2

    input_image_read = nib.load(input_image)
    input_image_data = input_image_read.get_fdata(dtype="float32")

    u = utils.heat_solver_tensor_3D_P1_grad_conj(input_image_data, g, beta, h, t_step)

    img = utils.spm_write_vol(input_image, u)

    nib.save(img, "./regularized_" + os.path.basename(input_image))

    return os.path.abspath("./regularized_" + os.path.basename(input_image))
