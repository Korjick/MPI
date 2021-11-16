#include <iostream>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

namespace Task1 {
	int main(int argc, char** argv)
	{
		int rank;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::cout << "Hello World!\n";
		MPI_Finalize();

		return MPI_SUCCESS;
	}
}

namespace Task2 {
	int main(int argc, char** argv)
	{
		int arr[]{ 41, 18467, 6334, 26500, 19169, 15724, 11478, 29358, 26962, 24464 };
		int part = log(10) / log(2);

		int rank, size, div, tmp, dest, source;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Status status;

		if (size > part)
			div = part;
		else
			div = size;

		if (rank < div) {
			dest = rank + 1;
			source = rank - 1;

			int max = INT_MIN;
			if (source < 0) source = MPI_PROC_NULL;
			if (dest >= div) {
				dest = MPI_PROC_NULL;
				for (int i = rank * div; i < 10; i++)
				{
					if (arr[i] > max) max = arr[i];
				}

				MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
				cout << "max: " << (max > tmp ? max : tmp);
			}
			else {
				for (int i = rank * div; i < rank * div + div; i++)
				{
					if (arr[i] > max) max = arr[i];
				}

				if (source != MPI_PROC_NULL) {
					MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
					max = (max > tmp ? max : tmp);
					MPI_Send(&max, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
				}
				else {
					MPI_Send(&max, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
				}
			}
		}

		MPI_Finalize();

		return MPI_SUCCESS;
	}
}

namespace Task3 {
	int main(int argc, char** argv) {
		int i, id, np, N;
		double x, y, double_N, eTime, sTime, pTime;
		int lhit;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		MPI_Comm_size(MPI_COMM_WORLD, &np);

		sscanf_s(argv[1], "%lf", &double_N);
		N = lround(double_N);
		MPI_Barrier(MPI_COMM_WORLD);
		sTime = MPI_Wtime();
		lhit = 0;
		srand((unsigned)(time(0)));
		int lN = N / np;

		for (i = 0; i < lN; i++) {
			x = ((double)rand()) / ((double)RAND_MAX);
			y = ((double)rand()) / ((double)RAND_MAX);
			if (((x * x) + (y * y)) <= 1) lhit++;
		}

		int hit = 0;
		MPI_Allreduce(&lhit, &hit, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		double est;
		est = (hit * 4) / ((double)N);
		MPI_Barrier(MPI_COMM_WORLD);

		eTime = MPI_Wtime();
		pTime = fabs(eTime - sTime);

		if (id == 0) {
			printf("Number of Points Used:      %d\n", N);
			printf("Estimate of Pi:         %24.16f\n", est);
			printf("Elapsed Wall time:      %5.3e\n", pTime);
		}

		MPI_Finalize();

		return MPI_SUCCESS;
	}
}

namespace Task4 {
	int main(int argc, char** argv) {
		const int n = strtol(argv[1], NULL, 10);
		int size, * sendbuf;
		int rank, * rbuf, i, * displs, * scounts;
		int sum, count;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		sendbuf = (int*)malloc(n);
		rbuf = (int*)malloc(n);
		displs = (int*)malloc(size);
		scounts = (int*)malloc(size);

		for (i = 0; i < n; i++) {
			rbuf[i] = 0;
			sendbuf[i] = rand();
			if (rank == 0)
				printf("%d ", sendbuf[i]);
		}
		printf("\n");

		for (int i = 0; i < size; i++) {
			if (i == size - 1)
				scounts[i] = n - (i * (n / size));
			else
				scounts[i] = n / size;

			displs[i] = i * (n / size);
		}

		MPI_Scatterv(sendbuf, scounts, displs, MPI_INT, rbuf, n, MPI_INT, 0, MPI_COMM_WORLD);

		int lsum = 0, lcount = 0;
		for (int i = 0; i < n; i++) {
			if (rbuf[i] > 0) {
				lsum += rbuf[i];
				lcount++;
			}
		}

		MPI_Reduce(&lsum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&lcount, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("average: %f ", (1.0f * sum / count));
		}

		MPI_Finalize();
		return MPI_SUCCESS;
	}
}

namespace Task5 {
	int main(int argc, char** argv) {
		const int n = strtol(argv[1], NULL, 10);
		int size, * sendbuf1, * sendbuf2;
		int rank, * rbuf1, * rbuf2, * res, i, * displs, * scounts;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		rbuf1 = (int*)malloc(n);
		rbuf2 = (int*)malloc(n);
		res = (int*)malloc(n);
		sendbuf1 = (int*)malloc(n);
		sendbuf2 = (int*)malloc(n);
		displs = (int*)malloc(size);
		scounts = (int*)malloc(size);

		for (i = 0; i < n; i++) {
			rbuf1[i] = rbuf2[i] = res[i] = 0;

			sendbuf1[i] = 1 + rand() / ((RAND_MAX + 1u) / 6);
			sendbuf2[i] = 1 + rand() / ((RAND_MAX + 1u) / 6);
		}

		if (rank == 0) {
			for (i = 0; i < n; i++)
				printf("%d ", sendbuf1[i]);
			printf("\n");
			for (i = 0; i < n; i++)
				printf("%d ", sendbuf2[i]);
			printf("\n");
		}

		for (int i = 0; i < size; i++) {
			if (i == size - 1)
				scounts[i] = n - (i * (n / size));
			else
				scounts[i] = n / size;

			displs[i] = i * (n / size);
		}

		MPI_Scatterv(sendbuf1, scounts, displs, MPI_INT, rbuf1, n, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(sendbuf2, scounts, displs, MPI_INT, rbuf2, n, MPI_INT, 0, MPI_COMM_WORLD);

		for (i = 0; i < n; i++) {
			if (rbuf1[i] < 0) break;
			rbuf1[i] *= rbuf2[i];
		}

		MPI_Gatherv(rbuf1, scounts[rank], MPI_INT, res, scounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("res: ");
			for (i = 0; i < n; i++)
				printf("%d ", res[i]);
		}

		MPI_Finalize();
		return MPI_SUCCESS;
	}
}

namespace Task6 {
	int main(int argc, char** argv) {
		const int n = strtol(argv[1], NULL, 10);
		const int m = strtol(argv[2], NULL, 10);

		int size, * sendbuf;
		int rank, * rbuf, * displs, * scounts;
		int max, min, tmp;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Status status;

		sendbuf = (int*)malloc(n * m);
		rbuf = (int*)malloc(n * m);
		displs = (int*)malloc(size);
		scounts = (int*)malloc(size);

		for (int i = 0; i < n * m; i++) {
			rbuf[i] = 0;
			if (rank == 0) {
				sendbuf[i] = rand();
				printf_s("%d ", sendbuf[i]);
				if ((i + 1) % m == 0) printf("\n");
			}
		}

		if (rank == 0) {
			for (int i = 0; i < size; i++) {
				if (i == size - 1)
					scounts[i] = m * (n - (i * (n / size)));
				else
					scounts[i] = n / size * m;

				displs[i] = m * i * (n / size);
			}

			for (int i = 0; i < size; i++) {
				printf_s("%d ", scounts[i]);
			}
			printf_s("\n");
			for (int i = 0; i < size; i++) {
				printf_s("%d ", displs[i]);
			}
		}

		MPI_Scatterv(sendbuf, scounts, displs, MPI_INT, rbuf, n * m, MPI_INT, 0, MPI_COMM_WORLD);

		int lmax = rbuf[0], lmin = rbuf[0];
		for (int i = 0; i < n * m; i++) {
			if (rbuf[i] <= 0) break;
			if (rbuf[i] > lmax) lmax = rbuf[i];
			if (rbuf[i] < lmin) lmin = rbuf[i];
		}

		int dest = rank + 1;
		int source = rank - 1;
		if (dest >= size) {
			dest = MPI_PROC_NULL;
			MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
			max = (lmax > tmp ? lmax : tmp);
			MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
			min = (lmin < tmp ? lmin : tmp);
			printf_s("max: %d\n", max);
			printf_s("min: %d\n", min);
		}
		else {
			if (source < 0) {
				source = MPI_PROC_NULL;
				MPI_Send(&lmax, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
				MPI_Send(&lmin, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
			}
			else {
				MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
				max = (lmax > tmp ? lmax : tmp);
				MPI_Recv(&tmp, 1, MPI_INT, source, 42, MPI_COMM_WORLD, &status);
				min = (lmin < tmp ? lmin : tmp);

				MPI_Send(&max, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
				MPI_Send(&min, 1, MPI_INT, dest, 42, MPI_COMM_WORLD);
			}
		}

		MPI_Finalize();
		return MPI_SUCCESS;
	}
}

namespace Task7 {

	const int N_DIM = 4;

	void RowMatrixVectorMultiply(double* mat, double* vec, double* result);

	int main(int argc, char** argv) {
		int rank, size;

		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);


		if (N_DIM % size) {
			MPI_Finalize();
			return(0);
		}

		double matrix_data[N_DIM][N_DIM];
		double vector_data[N_DIM];
		double result[N_DIM] = { 0.0 };

		if (rank == 0) {
			for (int i = 0; i < N_DIM; i++)
			{
				vector_data[i] = 1 + rand() / ((RAND_MAX + 1u) / 6);
				for (int j = 0; j < N_DIM; j++)
				{
					matrix_data[i][j] = 1 + rand() / ((RAND_MAX + 1u) / 6);
				}
			}
		}
		RowMatrixVectorMultiply((double*)matrix_data, vector_data, result);

		/* Printing the Matrix*/
		if (rank == 0) {
			printf("Matrix  :\n");
			for (int i = 0; i < N_DIM; i++) {
				for (int j = 0; j < N_DIM; j++)
					printf("%.5f ", matrix_data[i][j]);
				printf("\n");
			}
			printf("Vector :\n");
			for (int i = 0; i < N_DIM; i++)
				printf("%.5f ", vector_data[i]);
			printf("\n\n");

			printf("Vector :\n");
			for (int i = 0; i < N_DIM; i++)
				printf("%.5f ", vector_data[i]);
			printf("\n\n");

			printf("Result :\n");
			for (int i = 0; i < N_DIM; i++)
				printf("%.5f ", result[i]);
			printf("\n\n");
		}

		MPI_Finalize();
		return(0);
	}

	void RowMatrixVectorMultiply(double* matrix_data, double* vector_data, double* result) {
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		double* localresult = new double[N_DIM / size]{};
		double matrix[N_DIM][N_DIM];   //local matrix

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Scatter(matrix_data, (N_DIM * N_DIM / size), MPI_DOUBLE, matrix, (N_DIM * N_DIM) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(vector_data, N_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);// Broadcast the Vector

		//Calculate the results
		for (int i = 0; i < (N_DIM / size); i++)
			for (int j = 0; j < N_DIM; j++)
				localresult[i] += vector_data[j] * matrix[i][j];

		MPI_Gather(localresult, (N_DIM) / size, MPI_DOUBLE, result, (N_DIM) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
}

namespace Task8 {
	int main(int argc, char** argv) {
		int rank, size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		const int n = 20;
		int part = n / size;
		int* recv = new int[part];
		if (rank == 0)
		{
			int array[n];
			for (int i = 0; i < n; i++)
			{
				array[i] = rand() % 1000;
				printf("%4d ", array[i]);
			}
			printf("\n");

			for (int i = 0; i < size; i++)
			{
				int* send_buf = new int[part];
				for (int j = 0; j < part; j++)
					send_buf[j] = array[part * i + j];

				if (i == 0)
					recv = send_buf;
				else
					MPI_Send(send_buf, part, MPI_INT, i, i, MPI_COMM_WORLD);

				delete[] send_buf;
			}
		}
		MPI_Status status;
		if (rank > 0)
		{
			MPI_Recv(recv, part, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
		}
		printf("Rank = %d: ", rank);
		for (int i = 0; i < part; i++)
		{
			printf("%4d ", recv[i]);
		}
		printf("\n");
		if (rank != size - 1)
		{
			MPI_Send(recv, part, MPI_INT, size - 1, size - 1, MPI_COMM_WORLD);
		}
		else
		{
			int result[n];
			for (int i = 0; i < part; i++)
				result[rank * part + i] = recv[i];
			delete[] recv;

			for (int i = 0; i < size - 1; i++)
				MPI_Recv(result + i * part, part, MPI_INT, i, rank, MPI_COMM_WORLD, &status);

			for (int i = 0; i < n; i++)
				printf("%4d ", result[i]);
		}
		MPI_Finalize();
		return MPI_SUCCESS;
	}
}

namespace Task9 {
	int main(int argc, char** argv) {
		int rank, size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		const int n = 20;
		int part_size = n / size;
		int array[n];
		int result[n];

		int* recv = new int[part_size];
		int* sendcounts = new int[size];
		int* displs = new int[size];
		int* rev = new int[size];

		if (rank == 0)
		{
			for (int i = 0; i < n; i++)
			{
				array[i] = rand() % 1000;
				printf("%4d ", array[i]);
			}
			printf("\n");
			for (int i = 0; i < size; i++)
			{
				sendcounts[i] = part_size;
				displs[i] = i * part_size;
				rev[size - i - 1] = i * part_size;
			}
		}

		MPI_Scatterv(array, sendcounts, displs, MPI_INT, recv, part_size, MPI_INT, 0, MPI_COMM_WORLD);
		for (int i = 0; i < part_size / 2; i++)
			swap(recv[i], recv[part_size - i - 1]);
		
		delete[] displs;
		MPI_Gatherv(recv, part_size, MPI_INT, result, sendcounts, rev, MPI_INT, 0, MPI_COMM_WORLD);

		if (rank == 0)
		{
			for (int i = 0; i < n; i++)
				printf("%4d ", result[i]);
		}

		delete[] recv;
		delete[] sendcounts;
		delete[] rev;
		MPI_Finalize();

		return MPI_SUCCESS;
	}
}

namespace Task10 {
	int main(int argc, char** argv) {
		int rank, size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Status status;

		const int n = 100000000;
		int* array = new int[n];
		double from, to;

		if (rank == 0)
			for (int i = 0; i < n; i++)
				array[i] = rand();

		// Send
		if (rank == 0)
		{
			from = MPI_Wtime();
			MPI_Send(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
			to = MPI_Wtime();
			printf("Send time: %f\n", to - from);
		}
		else if (rank == 1)
		{
			MPI_Recv(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Send(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//SSend
		if (rank == 0)
		{
			from = MPI_Wtime();
			MPI_Ssend(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
			to = MPI_Wtime();
			printf("SSend time: %f\n", to - from);
		}
		else if (rank == 1)
		{
			MPI_Recv(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Ssend(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//BSend
		if (rank == 0)
		{
			int buffer_attached_size = MPI_BSEND_OVERHEAD + sizeof(int) * n;
			int* buffer_attached = new int[buffer_attached_size];
			MPI_Buffer_attach(buffer_attached, buffer_attached_size);

			from = MPI_Wtime();
			MPI_Bsend(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
			to = MPI_Wtime();

			MPI_Buffer_detach(&buffer_attached, &buffer_attached_size);

			delete[] buffer_attached;
			printf("Bsend time: %f\n", to - from);
		}
		else if (rank == 1)
		{
			MPI_Recv(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			int buffer_attached_size = MPI_BSEND_OVERHEAD + sizeof(int) * n;
			int* buffer_attached = new int[buffer_attached_size];
			MPI_Buffer_attach(buffer_attached, buffer_attached_size);

			MPI_Bsend(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Buffer_detach(&buffer_attached, &buffer_attached_size);

			delete[] buffer_attached;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//Rsend
		if (rank == 0)
		{
			from = MPI_Wtime();
			MPI_Rsend(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(array, n, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
			to = MPI_Wtime();
			printf("Rsend time: %f\n", to - from);
		}
		else if (rank == 1)
		{
			MPI_Recv(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Rsend(array, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		delete[] array;
		MPI_Finalize();

		return MPI_SUCCESS;
	}
}

namespace Task11 {
	int main(int argc, char** argv) {
		int rank, size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Status status;
		int n = 1;
		if (rank == 0)
		{
			MPI_Send(&n, 1, MPI_INT, rank + 1, rank + 1, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Recv(&n, 1, MPI_INT, rank - 1, rank, MPI_COMM_WORLD, &status);
			printf("Rank = %d. N = %d\n", rank, n);
			n = n * 10;
			MPI_Send(&n, 1, MPI_INT, (rank + 1) % size, (rank + 1) % size, MPI_COMM_WORLD);
		}
		if (rank == 0)
		{
			MPI_Recv(&n, 1, MPI_INT, size - 1, rank, MPI_COMM_WORLD, &status);
			printf("Rank = %d. N = %d\n", rank, n);
		}
		MPI_Finalize();
		return MPI_SUCCESS;
	}
}

int main(int argc, char** argv) {
	Task8::main(argc, argv);
	return EXIT_SUCCESS;
}