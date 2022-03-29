import numpy as np



def forward(X, trans, emission, initial_distribution):
    alpha = np.zeros((X.shape[0], trans.shape[0]))
    c = np.zeros(X.shape[0])
    alpha[0, :] = initial_distribution * emission[0]
    c[0] = 1/alpha[0].sum(axis=0)
    alpha[0] *= c[0]
    for t in range(1, X.shape[0]):
        for j in range(trans.shape[0]):
            alpha[t, j] = alpha[t - 1].dot(trans[:, j]) * emission[t, j]
        c[t] = 1/alpha[t].sum(axis=0)
        alpha[t] *= c[t]

    return alpha , c
# **********************************************************************************************************************
def backward(X, trans, emission,c):
    beta = np.zeros((X.shape[0], trans.shape[0]))
    # setting beta(T) = 1
    beta[X.shape[0] - 1] = np.ones((trans.shape[0]))

    # Loop in backward way from T-1 to
    for t in range(X.shape[0] - 2, -1, -1):
        for j in range(trans.shape[0]):
            beta[t, j] = (beta[t + 1] * emission[t + 1, :]).dot(trans[j, :])
        beta[t] *= c[t]

    return beta
# **********************************************************************************************************************
def baum_welch(X, trans, emission, initial_distribution, n_iter=1):
    M = trans.shape[0]
    T = len(X)

    for n in range(n_iter):
        alpha,c  = forward(X, trans, emission, initial_distribution)
        beta = backward(X, trans, emission,c)

        gamma = np.zeros((M, T - 1))
        xi = np.zeros((M, M, T - 1))
        for t in range(T - 1):
            gamma[:, t] = (alpha[t, :] * beta[t, :]) / np.dot(alpha[t, :], beta[t, :])
            denominator = np.dot(np.dot(alpha[t, :].T, trans) * emission[t + 1].T, beta[t + 1, :])
            for i in range(M):
                numerator = alpha[t, i] * trans[i, :] * emission[t + 1].T * beta[t + 1, :].T
                xi[i, :, t] = numerator / denominator

        trans = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))

        # Add additional T'th element in gamma
        # gamma = np.hstack((gamma , (alpha[T-1, :] * beta[T-1, :] / np.dot(alpha[T-1, :] , beta[T-1, :])).reshape((-1, 1)) ))

        # denominator = np.sum(gamma, axis=1)
        # print(denominator)
        # print((np.sum(gamma, axis=0)))
        # for t in range(T):
        #     emission[t, :] = np.sum(gamma)

        # emission = np.divide(emission, denominator)

    return trans
# **********************************************************************************************************************
def viterbi(X, trans, emission, initial_distribution):
    T = X.shape[0]
    M = trans.shape[0]

    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * emission[:, X[0]])

    prev = np.zeros((T - 1, M))

    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(trans[:, j]) + np.log(emission[j, X[t]])

            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)

            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)

    # Path Array
    S = np.zeros(T)

    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])

    S[0] = last_state

    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1

    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)

    # Convert numeric values to actual hidden states
    result = []
    for s in S:
        if s == 0:
            result.append("A")
        else:
            result.append("B")

    return result
# **********************************************************************************************************************