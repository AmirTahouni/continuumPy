import sympy as sp

class tensor:
    def __init__(self, order, dimension=3):
        self.dimension_ = dimension
        self.order_ = order
        if self.order_==0:
            self.tensor_=0
        elif self.order_==1:
            self.tensor_=[0,0,0]
        elif self.order_==2:
            self.tensor_=[[0,0,0],[0,0,0],[0,0,0]]
        elif self.order_==3:
            self.tensor_=[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]]
        elif self.order_==4:
            self.tensor_=[[[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],
                     [[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]],
                     [[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]]]

    def flattenIndex(self, indices):
        flat_index = 0
        i=self.order_-1
        for element in indices:
            flat_index+=element*self.dimension_**(i)
            i=i-1
        return flat_index

    def print(self):
        print(self.tensor_)

    def __setitem__(self, key, value):
        self.tensor_[key]=value

    def __getitem__(self, index):
        return self.tensor_[index]

class coordinateSystem:
    def __init__(self, g, thetaVariables, initVariance = 'co'):
        self.thetaVariables_ = sp.Matrix(thetaVariables)
        self.initVariance = initVariance
        self.g_= {'co': None, 'contra': None, initVariance: g}
        self.metricTensor_ = {'co': None, 'contra': None, initVariance: sp.Matrix(3, 3, [0] * 9)}
        self.deriveMetricTensor()
        if self.initVariance == 'co':
            self.CC('co', 'contra')
        elif self.initVariance == 'contra':
            self.CC('contra', 'co')

        self.christoffel={'first':[[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)],
                            'second':[[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]}
        self.calculateChristoffelSymbols()
        # gMatCo = sp.Matrix(self.g_['co'])
        # gMatContra = sp.Matrix(self.g_['contra'])
        # for i in range(3):
        #     for j in range(3):
        #         for k in range(3):
        #             #self.christoffel['first'][i][j][k] = sp.simplify(sp.DotProduct(sp.diff(gMatCo.row(i), self.thetaVariables_[j]) , gMatCo.row(k)))
        #             #self.christoffel['second'][i][j][k] = sp.simplify(sp.DotProduct(sp.diff(gMatCo.row(i), self.thetaVariables_[j]) , gMatContra.row(k)))
        #             self.christoffel['first'][i][j][k] = 1/2*(sp.diff(self.metricTensor_['co'][j][k],self.thetaVariables_[i]) +
        #                                                       sp.diff(self.metricTensor_['co'][k][i],self.thetaVariables_[j]) -
        #                                                       sp.diff(self.metricTensor_['co'][i][j],self.thetaVariables_[k]))
        #             self.christoffel['second'][i][j][k] = sum(self.christoffel['first'][i][j][s] * self.metricTensor_['contra'][s][k] for s in range(3))

    def CC(self, base='co', target='contra'):
        self.g_[target] = self.metricTensor_[target] * self.g_[base]

    def deriveMetricTensor(self):
        for i in range(3):
            for j in range(3):
                # Manually compute dot product
                dot_product = sum(self.g_[self.initVariance][i,k] * self.g_[self.initVariance][j,k] for k in range(3))
                self.metricTensor_[self.initVariance][i, j] = sp.simplify(dot_product)

        if self.initVariance=='co':
            self.metricTensor_['contra'] = self.metricTensor_['co'].inv()
        else:
            self.metricTensor_['co'] = self.metricTensor_['contra'].inv()

    def calculateChristoffelSymbols(self):
        g = self.metricTensor_[self.initVariance]
        ginv = g.inv()
        christoffel_first = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]
        christoffel_second = [[[0 for _ in range(3)] for _ in range(3)] for _ in range(3)]

        # for i in range(3):
        #     for j in range(3):
        #         for k in range(3):
        #             # Calculating Christoffel symbols of the first kind
        #             term1 = sp.diff(g[j, k], self.thetaVariables_[i])
        #             term2 = sp.diff(g[i, k], self.thetaVariables_[j])
        #             term3 = sp.diff(g[i, j], self.thetaVariables_[k])
        #             christoffel_first[i][j][k] = sp.simplify((term1 + term2 - term3) / 2)

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        term1 = sp.diff(g[j, l], self.thetaVariables_[k])
                        term2 = sp.diff(g[k, l], self.thetaVariables_[j])
                        term3 = sp.diff(g[j, k], self.thetaVariables_[l])
                        christoffel_second[i][j][k] += ginv[i, l] * (term1 + term2 - term3) / 2

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    term1 = sp.diff(self.metricTensor_['co'][j,k], self.thetaVariables_[i])
                    term2 = sp.diff(self.metricTensor_['co'][k,i], self.thetaVariables_[j])
                    term3 = sp.diff(self.metricTensor_['co'][i,j], self.thetaVariables_[k])
                    christoffel_first[i][j][k] = (term1 + term2 - term3) / 2

        self.christoffel['second'] = christoffel_second
        self.christoffel['first'] = christoffel_first


    def getMetricTensor(self, variance = 'co'):
        return self.metricTensor_[variance]

    def getG(self,variance):
        return self.g_[variance]

class transformation:
    def __init__(self, R, targetParams):
        self.betaMatrix = self.deriveBetaMat(R, targetParams)

    def deriveBetaMat(self, R, targetCoordParams):
        bMat = sp.Matrix(3, 3, [0] * 9)
        for i in range(3):
            for j in range(3):
                bMat[i, j] = sp.diff(R[j], targetCoordParams[i])
        return bMat

    def transformVector(self, vector, inputOutputVariance = 'co'):
        vector_matrix = sp.Matrix(vector.vector_[inputOutputVariance])
        if inputOutputVariance == 'co':
            transformed_vector = self.betaMatrix * vector_matrix
        elif inputOutputVariance == 'contra':
            transformed_vector = self.betaMatrix.inv().T * vector_matrix

        # Convert the result back to a list
        return [transformed_vector[i] for i in range(transformed_vector.rows)]

    def transformBaseVectors(self, c1, variance='co'):
        gp = sp.Matrix(3, 3, [0] * 9)
        for i in range(3):
            # Extract the i-th row vector from c1.g_
            vector_i = c1.g_[variance].row(i)
            if variance == 'co':
                result = self.betaMatrix * vector_i.T  # Transpose to make it a column vector
            elif variance == 'contra':
                result = self.betaMatrix.inv().T * vector_i.T
            for j in range(3):
                gp[i, j] = result[j, 0]
        return gp.T

    def transformMatrix(self, matrix, outputVariance='contra-contra'):
        input_matrix = sp.Matrix(matrix.matrix_[matrix.initVariance])
        transformed_matrix = sp.Matrix(3,3, [0] * 9)
        if outputVariance == 'co-co':
            transformed_matrix = self.betaMatrix * input_matrix * self.betaMatrix.T
        elif outputVariance == 'contra-contra':
            transformed_matrix = self.betaMatrix.inv().T * input_matrix * self.betaMatrix.inv()
        elif outputVariance == 'co-contra':
            transformed_matrix = self.betaMatrix * input_matrix * self.betaMatrix.inv().T
        elif outputVariance == 'contra-co':
            transformed_matrix = self.betaMatrix.inv() * input_matrix * self.betaMatrix.T

        # Convert the result back to a list of lists
        #matrix.matrix_[outputVariance] = [[transformed_matrix[i, j] for j in range(transformed_matrix.cols)] for i in range(transformed_matrix.rows)]
        return [[transformed_matrix[i, j] for j in range(transformed_matrix.cols)] for i in range(transformed_matrix.rows)]

class tensor2ndOrder:
    def __init__(self, CS, inputMatrix, variance = 'co-co'):
        self.initVariance = variance
        self.matrix_ = {'co-co': [[None for _ in range(3)] for _ in range(3)],
                        'contra-contra': [[None for _ in range(3)] for _ in range(3)],
                        'co-contra': [[None for _ in range(3)] for _ in range(3)],
                        'contra-co': [[None for _ in range(3)] for _ in range(3)],
                        variance: inputMatrix}
        self.CS = CS
        variances = ['co-co', 'contra-contra', 'co-contra', 'contra-co']
        for variance in variances:
            if variance != self.initVariance:
                self.CC(self.initVariance, variance)

    def mat(self, variance = 'co-co'):
        return self.matrix_[variance]

    def coDerivative(self, variance = 'co-co'):
        derivative = [[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if variance == 'co-co':
                        derivative[i][j][k] = (sp.diff(self.matrix_['co-co'][i][j],self.CS.thetaVariables_[k]) -
                                               sum(self.matrix_['co-co'][l][j]*self.CS.christoffel['second'][i][k][l]
                                                   for l in range(3)) -
                                               sum(self.matrix_['co-co'][i][l]*self.CS.christoffel['second'][k][j][l]
                                                   for l in range(3)))
                    elif variance == 'contra-co':
                        derivative[i][j][k] = (sp.diff(self.matrix_['contra-co'][i][j], self.CS.thetaVariables_[k]) -
                                               sum(self.matrix_['contra-co'][l][j] * self.CS.christoffel['second'][i][k][l]
                                                   for l in range(3)) -
                                               sum(self.matrix_['contra-co'][i][l] * self.CS.christoffel['second'][k][l][j]
                                                   for l in range(3)))
                    elif variance == 'co-contra':
                        derivative[i][j][k] = (sp.diff(self.matrix_['co-contra'][i][j], self.CS.thetaVariables_[k]) -
                                               sum(self.matrix_['co-contra'][l][j] * self.CS.christoffel['second'][i][k][l]
                                                   for l in range(3)) -
                                               sum(self.matrix_['co-contra'][i][l] * self.CS.christoffel['second'][k][l][j]
                                                   for l in range(3)))
                    elif variance == 'contra-contra':
                        derivative[i][j][k] = (sp.diff(self.matrix_['contra-contra'][i][j], self.CS.thetaVariables_[k]) -
                                               sum(self.matrix_['contra-contra'][l][j] * self.CS.christoffel['second'][k][l][i]
                                                   for l in range(3)) -
                                               sum(self.matrix_['contra-contra'][i][l] * self.CS.christoffel['second'][k][l][j]
                                                   for l in range(3)))
        return derivative


    def CC(self, base='co-co', target='contra-contra'):
        targetSup = target.split('-')[0]
        targetSub = target.split('-')[1]
        self.matrix_[target] = self.CS.metricTensor_[targetSup].T * sp.Matrix(self.matrix_[base]) * self.CS.metricTensor_[targetSub]
        self.matrix_[target] = [[self.matrix_[target][i, j] for i in range(self.matrix_[target].rows)] for j in range(self.matrix_[target].cols)]

class vector:
    def __init__(self, CS, inputVector, variance='co'):
        self.initVariance = variance
        self.vector_ = {'co': [None for _ in range(3)], 'contra': [None for _ in range(3)], variance: inputVector}
        self.CS=CS
        if self.initVariance == 'co':
            self.CC('co', 'contra')
        elif self.initVariance == 'contra':
            self.CC('contra', 'co')

    def vec(self,variance = 'co'):
        return self.vector_[variance]

    def coDerivative(self, variance='contra'):
        derivative = [[0 for _ in range(3)] for _ in range(3)]  # Initialize a list of 3 zero vectors

        for i in range(3):
            for j in range(3):
                if variance == 'co':
                    derivative[i][j] = (sp.diff(self.vector_['co'][i], self.CS.thetaVariables_[j, 0]) -
                                        sum([self.vector_['co'][k] * self.CS.christoffel['second'][k][i][j]
                                             for k in range(3)]))
                elif variance == 'contra':
                    derivative[i][j] = (sp.diff(self.vector_['contra'][i], self.CS.thetaVariables_[j, 0]) +
                                        sum([self.vector_['contra'][k] * self.CS.christoffel['second'][j][k][i]
                                             for k in range(3)]))
        return derivative

    def getCo(self):
        return self.vector_['co']
    def getContra(self):
        return self.vector_['contra']

    def CC(self, base='co', target='contra'):
        self.vector_[target] = sp.Matrix(self.CS.metricTensor_[target] * sp.Matrix(self.vector_[base]))
        self.vector_[target] = [self.vector_[target][i, 0] for i in range(self.vector_[target].rows)]

