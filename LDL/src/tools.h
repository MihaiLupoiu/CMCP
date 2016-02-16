//
// Created by Mihaita Alexandru Lupoiu on 28/10/15.
//

#ifndef CMCP_LDL_FACTORIZATION_TOOLS_H
#define CMCP_LDL_FACTORIZATION_TOOLS_H

void printMatrix(double** matrix,size_t size);
void printMatrixCustom(double* matrix,size_t size_r,size_t size_c);
void printMatrixCustomCPU(double* matrix,size_t size_r,size_t size_c, int poces);

void printVector(double* vector, size_t size);

void initialize(double **matrix ,size_t size);
void initializeCust(double *matrix ,size_t size);

void copyMatrix(double **matrix, double **matrixToCopy ,size_t size);
void copyMatrixCustom(double *matrix, double *matrixToCopy ,size_t size);

double **make2dmatrix(size_t size);
double **make2dmatrixCustom(size_t size_r,size_t size_c);

void free2dmatrix(double **matrix, size_t size);

#endif //CMCP_LDL_FACTORIZATION_TOOLS_H
