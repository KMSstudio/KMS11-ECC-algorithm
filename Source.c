#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "KMS11.h"

int main(void) {
	CURVE cur;
	POINT PuA, PuB;
	POINT SeA, SeB;
	KMS A, B;
	clock_t st;
	double tval;
	double mean[3] = { 0, };
	int start_seed, lng;

	double maxv[3] = { 0, 0, 0 };
	double minv[3] = { 200, 200, 200 };

	// ************************* base setup ************************* //
	HEADER_SETUP(); //setup base value, array in header
	VLABL_SETUP(); //allocation memory to pointer in header
	CURVE_SETUP(&cur); //allocation memory to pointer in curve structure

	// ************************* pakg setup ************************* //
	KEY_PACKAGE_SETUP(&SeA, &PuA, &A); //allocation memory to pointer in key package
	KEY_PACKAGE_SETUP(&SeB, &PuB, &B); //allocation memory to pointer in key package

	printf("main seed : ");
	scanf("%d", &start_seed);
	srand(start_seed);

	printf("key length : ");
	scanf("%d", &lng);

	// ************************* prepareing ************************* //
	kECC_setSecp256r1(&cur); //set curve to secp256r1
	// ************************* strt setup ************************* //
	kECC_registArr(&cur, lng); //make structure (structure's form is curve's form)
	for (int now_seed = start_seed, cnt = 1; ; now_seed++, cnt++) {
		srand(now_seed);
		printf("{ seed %d }\n", now_seed);

		st = clock();
		TEMP_RAND(&A, &cur.n, lng); //temporary (lng bit) random function
		TEMP_RAND(&B, &cur.n, lng); //temporary (lng bit) random function
		tval = (double)(clock() - st) / 2;
		mean[0] += tval;
		printf(" : random number claculated\t   %4.1fms\t[%7.3f ms]\n", tval, mean[0] / cnt);

		st = clock();
		kECC_PublicKey(&PuA, &A, &cur); //make userA's public key
		kECC_PublicKey(&PuB, &B, &cur); //make userB's spublic key
		tval = (double)(clock() - st) / 2;
		mean[1] += tval;
		maxv[1] = (maxv[1] > tval) ? maxv[1] : tval;
		minv[1] = (minv[1] < tval) ? minv[1] : tval;
		printf(" : public key claculated\t%s %4.1fms\t[%7.3f ms] %f %f \n", (tval <= (mean[1] / cnt)) ? "¡å" : "¡ã", tval, mean[1] / cnt, maxv[1], minv[1]);

		st = clock();
		kECC_SecretKey(&SeA, &PuB, &A, &cur); //make userA's secret key
		kECC_SecretKey(&SeB, &PuA, &B, &cur); //make userB's secret key
		tval = (double)(clock() - st) / 2;
		mean[2] += tval;
		maxv[2] = (maxv[2] > tval) ? maxv[2] : tval;
		minv[2] = (minv[2] < tval) ? minv[2] : tval;
		printf(" : secret key claculated\t%s %4.1fms\t[%7.3f ms] %f %f \n", (tval <= (mean[2] / cnt)) ? "¡å" : "¡ã", tval, mean[2] / cnt, maxv[2], minv[2]);

		if ((kECC_cmp(&SeA.x, &SeB.x)) || (kECC_cmp(&SeA.y, &SeB.y))) {
			system("cls");
			printf("[ Here is some error to calculate ECC key pair ]\n");
			printf("srand seed = %d\n", now_seed);

			printf("\n== Pri b state ==\n");
			kECC_print(&A);
			kECC_print(&B);

			printf("\n== Pub x state ==\n");
			kECC_print(&PuA.x);
			kECC_print(&PuB.x);
			printf("\n== Pub y state ==\n");
			kECC_print(&PuA.y);
			kECC_print(&PuB.y);

			printf("\n== Sec x state ==\n");
			kECC_print(&SeA.x);
			kECC_print(&SeB.x);
			printf("\n== Sec y state ==\n");
			kECC_print(&SeA.y);
			kECC_print(&SeB.y);
			break;
		}
		else {
			printf(" : normal state\n\n");
		}
	}

	// ************************* unlod base ************************* //
	VLABL_FREE();
	CURVE_FREE(&cur);
	// ************************* unlod pakg ************************* //
	KEY_PACKAGE_FREE(&SeA, &PuA, &A);
	KEY_PACKAGE_FREE(&SeB, &PuB, &B);
	// ************************* unlod strt ************************* //
	ARRAY_FREE();

	return 0;
}