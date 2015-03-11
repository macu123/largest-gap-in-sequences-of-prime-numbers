#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpi.h>

double MPI_Wtime();
long getLargestPrimegap(long start_num, long end_num, long upper_bound);

int main(int argc, char* argvp[]) {
	int my_rank;
	int p;
	long buffer_sd;
	long buffer_rv[2];
	MPI_Request request;
	MPI_Status status;

	double starttime;
   	double endtime;

    long lower_bound = strtol(argvp[1], NULL, 10);
	long upper_bound = strtol(argvp[2], NULL, 10);
	
	MPI_Init(&argc, &argvp);
	starttime = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if(my_rank == 0) {
		long range = upper_bound - lower_bound + 1;
		long segmentSize = 0;
	   	int i;
	   	//printf("Process %d start distributing segments\n", my_rank);
		if(range > 10000000*(p-1)) {
			segmentSize = 1000000;
	  	}
		else if(range > 1000000*(p-1)) {
			segmentSize = 100000;
		}
		else if(range > 100000*(p-1)) {
			segmentSize = 10000;
		}
		else if(range > 10000*(p-1)) {
			segmentSize = 1000;
		}
		else if(range > 1000*(p-1)) {
			segmentSize = 100;
		}
		else if(range > 100*(p-1)) {
			segmentSize = 10;
		}
	  	else {
			segmentSize = range/(p-1);
		}

		for(i=1; i<p; i++) {
			long buffer[2];
			buffer[0] = lower_bound;
			buffer[1] = lower_bound + segmentSize - 1;
	  		MPI_Isend(buffer, 2, MPI_LONG, i, 0, MPI_COMM_WORLD, &request);
	  		lower_bound += segmentSize;
		}
		//printf("Process %d finish distributing segments\n", my_rank);

		long buffers_rv[p];
		MPI_Request requests[p];

		for(i=1; i<p; i++) {
			MPI_Irecv(&buffers_rv[i-1], 1, MPI_LONG, i, 1, MPI_COMM_WORLD, &requests[i-1]);
		}

		//printf("Process %d is ready to receive data\n", my_rank);

		int indx;
		long largestPrimegap = 0;
		int nums_process_tonotified = p - 1;
		while(nums_process_tonotified > 0) {
			MPI_Waitany(p-1, requests, &indx, &status);
			//printf("Process %d received data from Process %d largest_primegap: %ld\n", my_rank, indx+1, buffers_rv[indx]);
			if(buffers_rv[indx] > largestPrimegap) {
				largestPrimegap = buffers_rv[indx];
				//printf("Process %d The largestPrimegap is %ld\n",my_rank, largestPrimegap);
			}
			MPI_Irecv(&buffers_rv[indx], 1, MPI_LONG, indx+1, 1, MPI_COMM_WORLD, &requests[indx]);
			long buffer[2];
			buffer[0] = lower_bound;
			if(upper_bound - lower_bound >= segmentSize) {
				buffer[1] = lower_bound + segmentSize - 1;
				lower_bound += segmentSize;
			}
			else if(upper_bound - lower_bound >= 0) {
				buffer[1] = upper_bound;
				lower_bound = upper_bound + 1;
			}
			else {
				nums_process_tonotified -= 1;
			}
			MPI_Isend(buffer, 2, MPI_LONG, indx+1, 0, MPI_COMM_WORLD, &request);
			//printf("Process %d has sent data to Process %d data: %ld %ld\n", my_rank, indx+1, buffer[0], buffer[1]);
		}
		printf("The largestPrimegap is %ld\n", largestPrimegap);
	}
	
	else {
		//printf("Process %d is ready to receive data from 0\n", my_rank);
		MPI_Recv(buffer_rv, 2, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
		//printf("Process %d has received data from Process 0: %ld %ld\n", my_rank, buffer_rv[0], buffer_rv[1]);
		while(buffer_rv[0] <= upper_bound) {
			buffer_sd = getLargestPrimegap(buffer_rv[0], buffer_rv[1], upper_bound);
			//printf("Process %d largestPrimegap: %ld\n", my_rank, buffer_sd);
			MPI_Isend(&buffer_sd, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD, &request);
			//printf("Process %d sent data to Process 0: %ld\n", my_rank, buffer_sd);
			MPI_Recv(buffer_rv, 2, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
			//printf("Process %d has received data from 0: %ld %ld\n", my_rank, buffer_rv[0], buffer_rv[1]);
			}
		}
		
		//printf("Process %d finishes!\n", my_rank);
		endtime = MPI_Wtime();
		MPI_Finalize();

		printf("Process %d took %f seconds\n", my_rank, endtime - starttime);
		return 0;
	}
	
	
long getLargestPrimegap(long start_num, long end_num, long upper_bound) {
	long next_prime;
	long largestPrimegap = 0;

	mpz_t a, b;

	mpz_init(a);
	mpz_init(b);
	mpz_set_si(a, start_num);
	
	if(mpz_probab_prime_p(a, 25) != 2) {
		mpz_nextprime(b, a);
		mpz_set(a, b);
		start_num = mpz_get_si(a);
	}
	
	while(start_num <= end_num) {
		mpz_nextprime(b, a);
		next_prime = mpz_get_si(b);
		if(next_prime > upper_bound) {
			break;
		}

		if(next_prime - start_num > largestPrimegap) {
			largestPrimegap = next_prime - start_num;
		}
		mpz_set(a, b);
		start_num = mpz_get_si(a);
	}
	
	//printf("Process %d largest_gap %ld\n", my_rank, largest_gap);	
	
	/*if(my_rank == 0) {
		printf("Process %d Result %ld\n", my_rank, result);
	}*/
	

	mpz_clear(a);
	mpz_clear(b);

	return largestPrimegap;
}
