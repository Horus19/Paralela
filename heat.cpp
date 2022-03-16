// heat.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include "mpi.h"

using namespace std;

void main(int argc, char** argv) {
    // Ejecuto la funcion de calor
    double Tmin;
    double Tmax;
    double tmax;
    double xmax;
    double dx;
    double dt;
    double k;
    int taskid, numtasks, numworkers, i, j, dest, Nt, Nx, xpertask, n;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    numworkers = numtasks-1;

    // Inicializar los valores
    Tmin = 20;
    Tmax = 100;
    tmax = 1;
    xmax = 2;
    dx = 0.2;
    dt = 0.01;
    k = 1;

    // Inicializar el vector T
    Nt = (int)(tmax / dt) + 1;
    Nx = (int)(xmax / dx) + 1;
    n = Nt * Nx;
    double* T = (double*)malloc(n * sizeof(double));

    // Master

    if (taskid == 0) {

        printf("MPI_mm ha empezado con %d procesos.\n",numtasks);

        for (i = 0; i < (Nt * Nx); i++) {
            T[i] = 0;
        }

        // Inicializar la primera fila

        for (j = 0; j < Nx; j++) {
            T[j] = Tmin;
        }
        // Inicializar la primera columna

        for (i = 0; i < Nt; i++) {
            T[i * Nx + 0] = Tmax;
        }
        // Inicializar la ultima columna

        for (i = 0; i < Nt; i++) {
            T[i * Nx + (Nx - 1)] = Tmax;
        }
        // Inicializar las esquinas
        T[0] = (T[1] + T[Nx]) / 2;
        T[Nx - 1] = (T[Nx - 2] + T[Nx + Nx - 1]) / 2;

        // Mandar el vector inicializado
        
        MPI_Send(T, n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        
        // Recibir el vector T calculado

        MPI_Recv(T, n , MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Imprimir el vector T

        ofstream heatsolution("heatsolution.txt", std::ofstream::out);

        if (heatsolution.is_open()) {
            double xt;
            double tx;
            tx = 0;
            xt = 0;
            int index = 0;
            for (i = 0; i < Nt; i++) {
                tx = (double)dt * i;
                for (j = 0; j < Nx; j++) {
                    xt = (double)dx * j;
                    heatsolution << xt << std::setw(10);
                    heatsolution << tx << std::setw(10);
                    heatsolution << T[index] << "\n";
                    index += 1;
                }
            }
            for (i = 0; i < Nt; i++) {
                for (j = 0; j < Nx; j++) {
                    cout << T[i * Nx + j] << " ";
                }
                cout << endl;

            heatsolution.close();
            }
        }

        // Libero la memoria asignada a T
        free(T);

    }

    // Worker
    
    else if (taskid == 1) {
        
        // Recibir el vector inicializado
        
        MPI_Status status;

        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);

        MPI_Get_count(&status, MPI_DOUBLE, &n);

        MPI_Recv(T, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Computar el resto del vector T

        for (i = 0; i < Nt - 1; i++) {

            for (j = 1; j < Nx - 1; j++) {
                T[((i + 1) * Nx) + j] = T[(i * Nx) + j] + (k * dt * ((T[(i * Nx) + j + 1] + T[(i * Nx) + j - 1] - (2 * T[(i * Nx) + j])) / pow(dx, 2)));
            }
        }
        
        // Mandar el vector calculado

        MPI_Send(T, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}

