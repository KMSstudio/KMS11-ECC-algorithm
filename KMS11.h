#include <stdio.h>
#include <stdlib.h> //malloc, calloc, free
#include <math.h>
#include <string.h> // strlen, memcpy

#ifndef __KMS_11__3_0_ELEVATOR__
#define __KMS_11__3_0_ELEVATOR__

#define SML_SIZE	6
#define BIG_SIZE	9

#define DIGIT		60 //code 0x1
#define SEMI_DIGIT	30
#define PNUM		1152921504606846976 //code 0x2

//#define MIN(x, y)	(y + ((x - y) & ((x - y) >> (sizeof(__int64) * CHAR_BIT - 1))))
//#define MAX(x, y)	(x - ((x - y) & ((x - y) >> (sizeof(__int64) * CHAR_BIT - 1))))

#define MIN(x, y)	((y) ^ (((x) ^ (y)) & -((x) < (y))))
#define MAX(x, y)	((x) ^ (((x) ^ (y)) & -((x) < (y))))

typedef __int64 numData;

numData _binary[DIGIT + 1];
numData sem_binary[DIGIT + 1];

typedef struct N {
	numData* arr;
	short S;
	short regS;
}KMS;

typedef struct KMS_ECC_DIVIDE_STRUCTURE {
	KMS div_val;
	KMS mod_val;
}DIVIDE;

typedef struct P {
	KMS x;
	KMS y;
}POINT;

typedef struct KMS_ECC_CURVE {
	KMS p;		// standard prime for curve
	KMS n;		// smallest value whict satisfy G × n ≡ (0, 0) (mod p)
	POINT G;	// gemerator G (key is Q = xG)
	KMS a;		// y²= x³+ ax + b
	KMS b;		// y²= x³+ ax + b
}CURVE;

const int tab64[64] = { 63,  0, 58,  1, 59, 47, 53,  2, 60, 39, 48, 27, 54, 33, 42,  3, 61, 51, 37, 40, 49, 18, 28, 20, 55, 30, 34, 11, 43, 14, 22,  4, 62, 57, 46, 52, 38, 26, 32, 41, 50, 36, 17, 19, 29, 10, 13, 21, 56, 45, 25, 31, 35, 16,  9, 12, 44, 24, 15,  8, 23,  7,  6,  5 };
POINT* sem_mtpValue = NULL;
size_t Sze = 0, BSze = 0;
short Reg = 0;

typedef struct N2 {
	KMS num;
	short sign;
}SKMS;

SKMS EEAsup_semi;
numData* EEAsup_bigArr;

numData* Modmtp_bigArr;
numData* Modeql_bigArr;

SKMS Inv_a1, Inv_a2;
KMS Inv_c1, Inv_c2;
DIVIDE Inv_Div;

KMS Pnt_ceta, Pnt_t1;
POINT Pnt_semi;

POINT Sec_semi;

POINT Pnt_change;

void HEADER_SETUP(void) {
	Sze = sizeof(numData) * SML_SIZE;
	BSze = sizeof(numData) * BIG_SIZE;
	Reg = SML_SIZE - 1;
	numData KRIT = 1;
	for (int i = 0; i <= DIGIT; i++) {
		_binary[i] = KRIT;
		sem_binary[i] = KRIT - 1;
		KRIT <<= 1;
	}
	return;
}

void CURVE_SETUP(CURVE* cur) {
	cur->a.arr = (numData*)malloc(Sze);
	cur->b.arr = (numData*)malloc(Sze);
	cur->a.regS = cur->b.regS = Reg;

	cur->G.x.arr = (numData*)malloc(Sze);
	cur->G.y.arr = (numData*)malloc(Sze);
	cur->G.x.regS = cur->G.y.regS = Reg;

	cur->n.arr = (numData*)malloc(Sze);
	cur->p.arr = (numData*)malloc(Sze);
	cur->n.regS = cur->p.regS = Reg;
}

void CURVE_FREE(CURVE* cur) {
	free(cur->a.arr);
	free(cur->b.arr);
	free(cur->G.x.arr);
	free(cur->G.y.arr);
	free(cur->n.arr);
	free(cur->p.arr);
}

void VLABL_SETUP(void) {
	EEAsup_semi.num.arr = (numData*)malloc(Sze);
	EEAsup_bigArr = (numData*)malloc(BSze);
	EEAsup_semi.num.regS = Reg;

	Modmtp_bigArr = (numData*)malloc(BSze);
	Modeql_bigArr = (numData*)malloc(BSze);

	Inv_Div.div_val.arr = (numData*)malloc(Sze);
	Inv_Div.mod_val.arr = (numData*)malloc(Sze);
	Inv_Div.div_val.regS = Inv_Div.mod_val.regS = Reg;

	Inv_a1.num.arr = (numData*)malloc(Sze);
	Inv_a2.num.arr = (numData*)malloc(Sze);
	Inv_a1.num.regS = Inv_a2.num.regS = Reg;

	Inv_c1.arr = (numData*)malloc(Sze);
	Inv_c2.arr = (numData*)malloc(Sze);
	Inv_c1.regS = Inv_c2.regS = Reg;

	Pnt_ceta.arr = (numData*)malloc(Sze);
	Pnt_t1.arr = (numData*)malloc(Sze);
	Pnt_ceta.regS = Pnt_t1.regS = Reg;

	Pnt_semi.x.arr = (numData*)malloc(Sze);
	Pnt_semi.y.arr = (numData*)malloc(Sze);
	Pnt_semi.x.regS = Pnt_semi.y.regS = Reg;

	Sec_semi.x.arr = (numData*)malloc(Sze);
	Sec_semi.y.arr = (numData*)malloc(Sze);
	Sec_semi.x.regS = Sec_semi.y.regS = Reg;
}

void VLABL_FREE(void) {
	free(EEAsup_semi.num.arr);
	free(EEAsup_bigArr);
	free(Modmtp_bigArr);
	free(Modeql_bigArr);
	free(Inv_Div.div_val.arr);
	free(Inv_Div.mod_val.arr);
	free(Inv_a1.num.arr);
	free(Inv_a2.num.arr);
	free(Inv_c1.arr);
	free(Inv_c2.arr);
	free(Pnt_ceta.arr);
	free(Pnt_t1.arr);
	free(Pnt_semi.x.arr);
	free(Pnt_semi.y.arr);
	free(Sec_semi.x.arr);
	free(Sec_semi.y.arr);
}

void KEY_PACKAGE_SETUP(POINT* Sec, POINT* Pub, KMS* Pri) {
	Sec->x.arr = (numData*)malloc(Sze);
	Sec->y.arr = (numData*)malloc(Sze);
	Sec->x.regS = Sec->y.regS = Reg;

	Pub->x.arr = (numData*)malloc(Sze);
	Pub->y.arr = (numData*)malloc(Sze);
	Pub->x.regS = Pub->y.regS = Reg;

	Pri->arr = (numData*)malloc(Sze);
	Pri->regS = Reg;
}

void KEY_PACKAGE_FREE(POINT* Sec, POINT* Pub, KMS* Pri) {
	free(Sec->x.arr);
	free(Sec->y.arr);

	free(Pub->x.arr);
	free(Pub->y.arr);

	free(Pri->arr);
}

void kECC_print(KMS* KMS) {
	printf("[ address=%p | size=%d ] : ", KMS, KMS->S);
	for (int i = KMS->S; i >= 0; i--) { printf("%llX ", KMS->arr[i]); }
	printf("\n");
	return;
}

#define PROCESS_KECC_BDGT1(x) x |= x >> 1; x |= x >> 2; x |= x >> 4; x |= x >> 8; x |= x >> 16; x |= x >> 32;
#define PROCESS_KECC_BDGT2(x) ((unsigned long long)((x - (x >> 1)) * 0x07EDD5E59A4E28C2)) >> 58

int BitSize(KMS* factor) {
	int value = ((factor->S) ? (factor->S - 1) : 0) * DIGIT;
	unsigned long long logVal = (unsigned long long)factor->arr[factor->S];
	PROCESS_KECC_BDGT1(logVal)
		value += tab64[PROCESS_KECC_BDGT2(logVal)];
	return value;
}

void set(const char str[], KMS* KMS) {
	int i, idx;
	int size;
	int semi = 1;

	memset(KMS->arr, 0, sizeof(numData) * KMS->regS);
	size = (int)strlen(str);
	idx = semi &= 0;
	for (i = size - 1; i >= 0; i--) {
		if (semi >= DIGIT) { semi &= 0; idx++; }
		if (!(str[i] ^ '1')) { KMS->arr[idx] |= _binary[semi]; }
		semi++;
	}
	KMS->S = idx;
	return;
}

static void kECC_privateOP_A(KMS* k, numData a, numData b, numData c, numData d, numData e) {
	memset(k->arr, 0, sizeof(numData) * k->regS);

	k->arr[0] = a;
	k->arr[1] = b;
	k->arr[2] = c;
	k->arr[3] = d;
	k->arr[4] = e;

	k->S = 4;
}
void kECC_setSecp256r1(CURVE* cur) {
	kECC_privateOP_A(&(cur->p), 0xFFFFFFFFFFFFFFF, 0xFFFFFFFFF, 0x0, 0xFFFF00000001000, 0xFFFF);

	kECC_privateOP_A(&(cur->a), 0xFFFFFFFFFFFFFFC, 0xFFFFFFFFF, 0x0, 0xFFFF00000001000, 0xFFFF);
	kECC_privateOP_A(&(cur->b), 0xBCE3C3E27D2604B, 0x1D06B0CC53B0F63, 0xBBD55769886BC65, 0x35D8AA3A93E7B3E, 0x5AC6);

	kECC_privateOP_A(&(cur->G.x), 0x4A13945D898C296, 0x37D812DEB33A0F, 0xCE6E563A440F277, 0xD1F2E12C4247F8B, 0x6B17);
	kECC_privateOP_A(&(cur->G.y), 0xBB6406837BF51F5, 0xCE33576B315ECEC, 0x7EB4A7C0F9E162B, 0x42E2FE1A7F9B8EE, 0xAFE3);

	kECC_privateOP_A(&(cur->n), 0x3B9CAC2FC632551, 0xE6FAADA7179E84F, 0x7FFFFFFFFFFFFFBC, 0xFFFF00000000FFF, 0xFFFF);
	return;
}

// left ( ) right | [ > : 1 ] [ = : 0 ] [ < : -1 ]
int kECC_cmp(KMS* left, KMS* right) {
	if (left->S != right->S) { return (left->S > right->S) ? 1 : -1; }
	for (int i = left->S; i >= 0; i--) {
		if (left->arr[i] > right->arr[i]) { return 1; }
		else if (left->arr[i] < right->arr[i]) { return -1; }
	}
	return 0;
}

void TEMP_RAND(KMS* source, KMS* factor, int bit) {
	int i;
	do {
		for (i = 0; i <= source->regS; i++) { source->arr[i] = 0; }
		for (i = 0; i < bit; i++) { source->arr[i / DIGIT] |= (numData)(rand() & 1) << (i % DIGIT); }
		for (i = bit / DIGIT; (i > 0) && (!source->arr[i]); i--) { ; }
		source->S = i;
	} while (((!source->S) && (!source->arr[0])) || (kECC_cmp(source, factor) != -1));
	return;
}

//source(S) = (1 / n(S)) % mod(S)
void kECC_modInv(KMS* source, KMS* n, KMS* mod) {
	numData* pArr;
	SKMS storage;
	KMS pKms;
	numData semiA, semiB, semiC, semiD;
	int idxp1, i, j, fe, front, semi;
	unsigned long long bit1, bit2;

	memset(Inv_a1.num.arr, 0, Sze);
	memset(Inv_a2.num.arr, 0, Sze);

	Inv_a1.num.S = Inv_a2.num.S = 0;
	Inv_a1.sign = Inv_a2.sign = 1;
	Inv_a1.num.arr[0] = 0;
	Inv_a2.num.arr[0] = 1;

	memcpy(Inv_c1.arr, mod->arr, Sze);
	Inv_c1.S = mod->S;
	memcpy(Inv_c2.arr, n->arr, Sze);
	Inv_c2.S = n->S;
	// 호출에서의 mod 인자가 소수이기 때문에, 항상 최대공약수가 1임.
	while (Inv_c2.S || (Inv_c2.arr[0] ^ (numData)1)) {
		//kECC_stdDivide(&Inv_Div, &Inv_c1, &Inv_c2);
		memset(Inv_Div.div_val.arr, 0, Sze);
		Inv_Div.div_val.regS = SML_SIZE - 1;
		Inv_Div.div_val.S = 0;

		memcpy(Inv_Div.mod_val.arr, Inv_c1.arr, Sze);
		Inv_Div.mod_val.regS = SML_SIZE - 1;
		Inv_Div.mod_val.S = Inv_c1.S;

		fe = MAX(Inv_c2.S, Inv_c1.S);
		for (i = fe; i >= 0; i--) {
			if (Inv_c2.arr[i] < Inv_c1.arr[i]) { i = -2; break; }
			else if (Inv_c2.arr[i] > Inv_c1.arr[i]) { break; }
		}
		if (i < 0) {
			if (!(Inv_c2.S || Inv_c1.S)) {
				Inv_Div.div_val.arr[0] = Inv_c1.arr[0] / Inv_c2.arr[0];
				Inv_Div.mod_val.arr[0] = Inv_c1.arr[0] % Inv_c2.arr[0];
			}
			else {
				front = MAX(Inv_c1.S - Inv_c2.S - 1, 0);
				if ((Inv_c2.S < Inv_c1.S) && (Inv_c2.arr[Inv_c2.S] < Inv_c1.arr[Inv_c1.S])) { front++; }
				if (front) {
					//for (i = Inv_c2.S + front; i >= front; i--) { Inv_c2.arr[i] = Inv_c2.arr[i - front]; }
					memcpy(Inv_c2.arr + front, Inv_c2.arr, sizeof(numData) * Inv_c2.S);
					//for (i = 0; i < front; i++) { Inv_c2.arr[i] = 0; }
					memset(Inv_c2.arr, 0, sizeof(numData) * front);
					Inv_c2.S += front;
				}
				front *= DIGIT;
				//semi = bitDigit(Inv_c1.arr[Inv_c1.S]) - bitDigit(Inv_c2.arr[Inv_c2.S]);
				bit1 = Inv_c1.arr[Inv_c1.S];
				bit2 = Inv_c2.arr[Inv_c2.S];
				PROCESS_KECC_BDGT1(bit1);
				PROCESS_KECC_BDGT1(bit2);
				semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
				if (Inv_c1.S != Inv_c2.S) { semi += DIGIT; }
				if (semi && (Inv_c1.S != Inv_c2.S) ? ((Inv_c2.arr[Inv_c2.S] >> (DIGIT - semi)) > Inv_c1.arr[Inv_c1.S]) : ((Inv_c2.arr[Inv_c2.S] << semi) > Inv_c1.arr[Inv_c1.S])) { semi--; }
				if (semi >= 0) {
					i = Inv_c2.S;
					if (Inv_c2.arr[i] >= _binary[DIGIT - semi]) {
						Inv_c2.arr[i + 1] |= Inv_c2.arr[i] >> (DIGIT - semi);
						Inv_c2.arr[i] = (Inv_c2.arr[i] & sem_binary[DIGIT - semi]) << semi;
						Inv_c2.S++;
					}
					else { Inv_c2.arr[i] <<= semi; }
					for (i--; i >= 0; i--) {
						Inv_c2.arr[i + 1] |= Inv_c2.arr[i] >> (DIGIT - semi);
						Inv_c2.arr[i] = (Inv_c2.arr[i] & sem_binary[DIGIT - semi]) << (semi);
					}
					front += semi;

					for (;;) {
						for (i = MAX(Inv_Div.mod_val.S, Inv_c2.S); i >= 0; i--) {
							if (Inv_Div.mod_val.arr[i] > Inv_c2.arr[i]) { i = -1; break; }
							else if (Inv_Div.mod_val.arr[i] < Inv_c2.arr[i]) { break; }
						}
						if (i < 0) {
							//kECC_sub(&(Inv_Div.mod_val), &(Inv_Div.mod_val), factor);
							for (i = 0; i <= Inv_Div.mod_val.S; i++) {
								Inv_Div.mod_val.arr[i] -= Inv_c2.arr[i];
								if (Inv_Div.mod_val.arr[i] < 0) {
									Inv_Div.mod_val.arr[i] += PNUM;
									Inv_Div.mod_val.arr[i + 1]--;
								}
							}
							for (i = Inv_Div.mod_val.S; (i > 0) && (!Inv_Div.mod_val.arr[i]); i--) { ; }
							Inv_Div.mod_val.S = i;
							Inv_Div.div_val.arr[front / DIGIT] |= ((numData)1 << (front % DIGIT));
						}
						if (!front) { break; }

						for (i = 0; i < Inv_c2.S; i++) {
							Inv_c2.arr[i] >>= 1;
							Inv_c2.arr[i] |= (Inv_c2.arr[i + 1] & 1) << (DIGIT - 1);
						}
						Inv_c2.arr[Inv_c2.S] >>= 1;
						if (!(Inv_c2.arr[Inv_c2.S])) { Inv_c2.S--; }
						front--;
					}
					for (i = Inv_c1.S - Inv_c2.S + 1; (i > 0) && (!Inv_Div.div_val.arr[i]); i--) { ; }
					Inv_Div.div_val.S = i;
				}
			}
		}

		//signed_EEDA_supporter(&Inv_a1, &Inv_a2, &(Inv_Div.div_val)); //| a2 = a1 - (a2 * (c1 / c2))
		EEAsup_semi.sign = Inv_a2.sign;
		if ((!Inv_Div.div_val.S) && (Inv_Div.div_val.arr[0] < _binary[SEMI_DIGIT])) {
			EEAsup_semi.num.arr[0] = idxp1 = 0;
			for (i = 0; i <= Inv_a2.num.S; i++) {
				idxp1 = MIN(i + 1, EEAsup_semi.num.regS);
				semiA = Inv_a2.num.arr[i] & sem_binary[SEMI_DIGIT];
				semiC = Inv_a2.num.arr[i] >> SEMI_DIGIT;
				EEAsup_semi.num.arr[i] += semiA * Inv_Div.div_val.arr[0];
				EEAsup_semi.num.arr[i] += (semiC * Inv_Div.div_val.arr[0] & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
				EEAsup_semi.num.arr[idxp1] = semiC * Inv_Div.div_val.arr[0] >> SEMI_DIGIT;
				EEAsup_semi.num.arr[idxp1] += EEAsup_semi.num.arr[i] >> DIGIT;
				EEAsup_semi.num.arr[i] &= sem_binary[DIGIT];
			}
			for (i = idxp1 + 1; i <= EEAsup_semi.num.regS; i++) { EEAsup_semi.num.arr[i] = 0; }
			for (i = idxp1; (i > 0) && (!EEAsup_semi.num.arr[i]); i--) { ; }
			EEAsup_semi.num.S = i;
		}
		else {
			// make EEAsup_semi.num = Inv_a2.mi, * mtp | [mod] {mtp, mod}
			memset(EEAsup_semi.num.arr, 0, Sze);
			for (i = 0; i <= Inv_a2.num.S; i++) {
				for (j = 0; j <= Inv_Div.div_val.S; j++) {
					idxp1 = MIN(i + j + 1, EEAsup_semi.num.regS);
					semiA = (Inv_Div.div_val.arr[j] & sem_binary[SEMI_DIGIT]);
					semiB = (Inv_a2.num.arr[i] & sem_binary[SEMI_DIGIT]);
					semiC = (Inv_Div.div_val.arr[j] >> SEMI_DIGIT);
					semiD = (Inv_a2.num.arr[i] >> SEMI_DIGIT);
					EEAsup_semi.num.arr[i + j] += semiA * semiB;
					EEAsup_semi.num.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
					EEAsup_semi.num.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
					EEAsup_semi.num.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
					EEAsup_semi.num.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
					EEAsup_semi.num.arr[idxp1] += EEAsup_semi.num.arr[i + j] >> DIGIT;
					EEAsup_semi.num.arr[i + j] &= sem_binary[DIGIT];
					EEAsup_semi.num.arr[idxp1] += semiC * semiD;
				}
			}
			for (i = idxp1; (i > 0) && (!EEAsup_semi.num.arr[i]); i--) { ; }
			EEAsup_semi.num.S = i;
		}

		storage = Inv_a1;
		Inv_a1 = Inv_a2;
		Inv_a2 = storage;

		if (EEAsup_semi.sign == Inv_a2.sign) {
			//if (kECC_sub(&Inv_a2.num, &Inv_a2.num, &EEAsup_semi.num)) { Inv_a2.sign = !Inv_a2.sign; } , EEAsup_semi loss its value in this process | [sub] {cmp}
			fe = MAX(Inv_a2.num.S, EEAsup_semi.num.S);
			for (i = fe; i >= 0; i--) {
				if (Inv_a2.num.arr[i] < EEAsup_semi.num.arr[i]) { i = -2; break; }
				else if (Inv_a2.num.arr[i] > EEAsup_semi.num.arr[i]) { break; }
			}

			if (i ^ -2) {
				for (i = 0; i <= EEAsup_semi.num.S; i++) {
					Inv_a2.num.arr[i] -= EEAsup_semi.num.arr[i];
					if (Inv_a2.num.arr[i] < 0) {
						Inv_a2.num.arr[i] += PNUM;
						Inv_a2.num.arr[i + 1]--;
					}
				}
				for (i = Inv_a2.num.S; (i > 0) && (!Inv_a2.num.arr[i]); i--) { ; }
				Inv_a2.num.S = i;
			}
			else {
				pKms = Inv_a2.num;
				Inv_a2.num = EEAsup_semi.num;
				EEAsup_semi.num = pKms;

				for (i = 0; i <= EEAsup_semi.num.S; i++) {
					Inv_a2.num.arr[i] -= EEAsup_semi.num.arr[i];
					if (Inv_a2.num.arr[i] < 0) {
						Inv_a2.num.arr[i] += PNUM;
						Inv_a2.num.arr[i + 1]--;
					}
				}
				for (i = Inv_a2.num.S; (i > 0) && (!Inv_a2.num.arr[i]); i--) { ; }
				Inv_a2.num.S = i;

				Inv_a2.sign = !Inv_a2.sign;
			}
		}
		else {
			for (i = 0; i <= EEAsup_semi.num.S; i++) {
				Inv_a2.num.arr[i] += EEAsup_semi.num.arr[i];
				if (Inv_a2.num.arr[i] >> DIGIT) {
					Inv_a2.num.arr[i] &= sem_binary[DIGIT];
					Inv_a2.num.arr[i + 1]++;
				}
			}
			idxp1 = MAX(EEAsup_semi.num.S, Inv_a2.num.S);
			Inv_a2.num.S = (Inv_a2.num.arr[idxp1 + 1]) ? idxp1 + 1 : idxp1;
		}
		//return;

		Inv_c1.S = Inv_c2.S;
		Inv_c2.S = Inv_Div.mod_val.S;

		pArr = Inv_c1.arr;
		Inv_c1.arr = Inv_c2.arr;
		Inv_c2.arr = Inv_Div.mod_val.arr;
		Inv_Div.mod_val.arr = pArr;
	}

	if (!Inv_a2.sign) {
		//kECC_sub(&Inv_a2.num, mod, &Inv_a2.num);
		source->arr[0] = 0;
		for (i = 0; i <= mod->S; i++) {
			source->arr[i] += mod->arr[i] - Inv_a2.num.arr[i];
			if (source->arr[i] < 0) {
				source->arr[i] += PNUM;
				source->arr[i + 1] = -1;
			}
			else { source->arr[i + 1] &= 0; }
		}
		for (i = mod->S; (i > 0) && (!source->arr[i]); i--) { ; }
		source->S = i;
	}
	else {
		memcpy(source->arr, Inv_a2.num.arr, Sze);
		source->S = Inv_a2.num.S;
	}
	return;
}

void kECC_pointAddEql(POINT* source, POINT* factor, KMS* p) {
	KMS pKMS;
	numData semiA, semiB, semiC, semiD;
	numData* arrStorage;
	int i, j, fe, idxp1 = 0, front, semi;
	unsigned long long bit1, bit2;
	pKMS.arr = Modmtp_bigArr;

	//𝜆 = 𝑦₂- 𝑦₁ | kECC_modSub(&Pnt_ceta, &(factor->y), &(source->y), p); | [sub] {cmp, add}
	fe = MAX(factor->y.S, source->y.S);
	for (i = fe; i >= 0; i--) {
		if (factor->y.arr[i] < source->y.arr[i]) { i = -2; break; }
		else if (factor->y.arr[i] > source->y.arr[i]) { break; }
	}
	if (i == -2) {
		Pnt_ceta.arr[0] = 0;
		for (i = 0; i < p->S; i++) {
			Pnt_ceta.arr[i] += p->arr[i] + factor->y.arr[i] - source->y.arr[i];
			if (Pnt_ceta.arr[i] >= PNUM) {
				Pnt_ceta.arr[i] &= sem_binary[DIGIT];
				Pnt_ceta.arr[i + 1] = 1;
			}
			else if (Pnt_ceta.arr[i] < 0) {
				Pnt_ceta.arr[i] += PNUM;
				Pnt_ceta.arr[i + 1] = -1;
			}
			else { Pnt_ceta.arr[i + 1] = 0; }
		}
		Pnt_ceta.arr[i] += p->arr[i] + factor->y.arr[i] - source->y.arr[i];
		for (i++; i <= Pnt_ceta.regS; i++) { Pnt_ceta.arr[i] = 0; }
		for (i = p->S; (i > 0) && (!Pnt_ceta.arr[i]); i--) { ; }
		Pnt_ceta.S = i;
	}
	else {
		Pnt_ceta.arr[0] = 0;
		for (i = 0; i < factor->y.S; i++) {
			Pnt_ceta.arr[i] += factor->y.arr[i] - source->y.arr[i];
			if (Pnt_ceta.arr[i] < 0) {
				Pnt_ceta.arr[i] += PNUM;
				Pnt_ceta.arr[i + 1] = -1;
			}
			else { Pnt_ceta.arr[i + 1] = 0; }
		}
		Pnt_ceta.arr[i] += factor->y.arr[i] - source->y.arr[i];
		for (i++; i <= Pnt_ceta.regS; i++) { Pnt_ceta.arr[i] = 0; }
		for (i = factor->y.S; (i > 0) && (!Pnt_ceta.arr[i]); i--) { ; }
		Pnt_ceta.S = i;
	}
	//t1 = 𝑥₂- 𝑥₁ | kECC_modSub(&Pnt_t1, &(factor->x), &(source->x), p); | [sub] {cmp, add}
	fe = MAX(factor->x.S, source->x.S);
	for (i = fe; i >= 0; i--) {
		if (factor->x.arr[i] < source->x.arr[i]) { i = -2; break; }
		else if (factor->x.arr[i] > source->x.arr[i]) { break; }
	}
	if (i == -2) {
		Pnt_t1.arr[0] = 0;
		for (i = 0; i < p->S; i++) {
			Pnt_t1.arr[i] += p->arr[i] + factor->x.arr[i] - source->x.arr[i];
			if (Pnt_t1.arr[i] >= PNUM) {
				Pnt_t1.arr[i] &= sem_binary[DIGIT];
				Pnt_t1.arr[i + 1] = 1;
			}
			else if (Pnt_t1.arr[i] < 0) {
				Pnt_t1.arr[i] += PNUM;
				Pnt_t1.arr[i + 1] = -1;
			}
			else { Pnt_t1.arr[i + 1] = 0; }
		}
		Pnt_t1.arr[i] += p->arr[i] + factor->x.arr[i] - source->x.arr[i];
		for (i++; i <= Pnt_t1.regS; i++) { Pnt_t1.arr[i] = 0; }
		for (i = p->S; (i > 0) && (!Pnt_t1.arr[i]); i--) { ; }
		Pnt_t1.S = i;
	}
	else {
		Pnt_t1.arr[0] = 0;
		for (i = 0; i < factor->x.S; i++) {
			Pnt_t1.arr[i] += factor->x.arr[i] - source->x.arr[i];
			if (Pnt_t1.arr[i] < 0) {
				Pnt_t1.arr[i] += PNUM;
				Pnt_t1.arr[i + 1] = -1;
			}
			else { Pnt_t1.arr[i + 1] = 0; }
		}
		Pnt_t1.arr[i] += factor->x.arr[i] - source->x.arr[i];
		for (i++; i <= Pnt_t1.regS; i++) { Pnt_t1.arr[i] = 0; }
		for (i = factor->x.S; (i > 0) && (!Pnt_t1.arr[i]); i--) { ; }
		Pnt_t1.S = i;
	}
	//t1 = 1 / 𝑥₂- 𝑥₁ //////////////////////////////////////////////////////////
	kECC_modInv(&Pnt_t1, &Pnt_t1, p);

	//𝜆 = 𝑦₂- 𝑦₁ / 𝑥₂- 𝑥₁ | kECC_modMultEql(&Pnt_ceta, &Pnt_t1, p); [mod] {mtp, mod}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_t1.S; i++) {
		for (j = 0; j <= Pnt_ceta.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_ceta.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_t1.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_ceta.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_t1.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_ceta.S = pKMS.S;
	memcpy(Pnt_ceta.arr, pKMS.arr, Sze);

	//𝑥₃= 𝜆²  | kECC_modSquare(&(Pnt_semi.x), &Pnt_ceta, p); | [mod] {mod, mtp}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_ceta.S; i++) {
		for (j = 0; j <= Pnt_ceta.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_ceta.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_ceta.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_ceta.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_ceta.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_semi.x.S = pKMS.S;
	memcpy(Pnt_semi.x.arr, pKMS.arr, Sze);

	//t1 = 𝑥₁+ 𝑥₂ | kECC_modAdd(&Pnt_t1, &(source->x), &(factor->x), p); | [add] {cmp, add}
	fe = MAX(source->x.S, factor->x.S);
	Pnt_t1.arr[0] = 0;
	for (i = 0; i < fe; i++) {
		Pnt_t1.arr[i] += source->x.arr[i] + factor->x.arr[i];
		if (Pnt_t1.arr[i] >= PNUM) {
			Pnt_t1.arr[i] &= sem_binary[DIGIT];
			Pnt_t1.arr[i + 1] = 1;
		}
		else { Pnt_t1.arr[i + 1] = 0; }
	}
	Pnt_t1.arr[i] += source->x.arr[i] + factor->x.arr[i];
	if (Pnt_t1.arr[i] >> DIGIT) {
		Pnt_t1.arr[i] &= sem_binary[DIGIT];
		Pnt_t1.arr[++i] = 1;
	}
	Pnt_t1.S = i;
	for (i++; i <= Pnt_t1.regS; i++) { Pnt_t1.arr[i] = 0; }
	fe = MAX(Pnt_t1.S, p->S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_t1.arr[i] < p->arr[i]) { break; }
		else if (Pnt_t1.arr[i] > p->arr[i]) { i = -2; break; }
	}
	if (i < 0) {
		for (i = 0; i <= Pnt_t1.S; i++) {
			Pnt_t1.arr[i] -= p->arr[i];
			if (Pnt_t1.arr[i] < 0) {
				Pnt_t1.arr[i] += PNUM;
				Pnt_t1.arr[i + 1]--;
			}
		}
		for (i = Pnt_t1.S; (i > 0) && (!Pnt_t1.arr[i]); i--) { ; }
		Pnt_t1.S = i;
	}
	//𝑥₃= 𝜆²- 𝑥₁- 𝑥₂ | kECC_modSub(&(Pnt_semi.x), &(Pnt_semi.x), &Pnt_t1, p); | [sub] {cmp, add}
	fe = MAX(Pnt_semi.x.S, Pnt_t1.S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_semi.x.arr[i] < Pnt_t1.arr[i]) { i = -2; break; }
		else if (Pnt_semi.x.arr[i] > Pnt_t1.arr[i]) { break; }
	}
	if (i == -2) {
		for (i = 0; i <= p->S; i++) {
			Pnt_semi.x.arr[i] += p->arr[i] - Pnt_t1.arr[i];
			if (Pnt_semi.x.arr[i] >= PNUM) {
				Pnt_semi.x.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.x.arr[i + 1]++;
			}
			else if (Pnt_semi.x.arr[i] < 0) {
				Pnt_semi.x.arr[i] += PNUM;
				Pnt_semi.x.arr[i + 1]--;
			}
		}
		for (i = p->S; (i > 0) && (!Pnt_semi.x.arr[i]); i--) { ; }
		Pnt_semi.x.S = i;
	}
	else {
		for (i = 0; i <= Pnt_semi.x.S; i++) {
			Pnt_semi.x.arr[i] -= Pnt_t1.arr[i];
			if (Pnt_semi.x.arr[i] < 0) {
				Pnt_semi.x.arr[i] += PNUM;
				Pnt_semi.x.arr[i + 1]--;
			}
		}
		for (i = Pnt_semi.x.S; (i > 0) && (!Pnt_semi.x.arr[i]); i--) { ; }
		Pnt_semi.x.S = i;
	}

	//𝑦₃= 𝑥₃- 𝑥₁ | kECC_modSub(&(Pnt_semi.y), &(source->x), &(Pnt_semi.x), p); | [sub] {cmp, add}
	fe = MAX(source->x.S, Pnt_semi.x.S);
	for (i = fe; i >= 0; i--) {
		if (source->x.arr[i] < Pnt_semi.x.arr[i]) { i = -2; break; }
		else if (source->x.arr[i] > Pnt_semi.x.arr[i]) { break; }
	}
	if (i == -2) {
		Pnt_semi.y.arr[0] = 0;
		for (i = 0; i < p->S; i++) {
			Pnt_semi.y.arr[i] += p->arr[i] + source->x.arr[i] - Pnt_semi.x.arr[i];
			if (Pnt_semi.y.arr[i] >= PNUM) {
				Pnt_semi.y.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.y.arr[i + 1] = 1;
			}
			else if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1] = -1;
			}
			else { Pnt_semi.y.arr[i + 1] = 0; }
		}
		Pnt_semi.y.arr[i] += p->arr[i] + source->x.arr[i] - Pnt_semi.x.arr[i];
		for (i++; i <= Pnt_semi.y.regS; i++) { Pnt_semi.y.arr[i] = 0; }
		for (i = p->S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}
	else {
		Pnt_semi.y.arr[0] = 0;
		for (i = 0; i < source->x.S; i++) {
			Pnt_semi.y.arr[i] += source->x.arr[i] - Pnt_semi.x.arr[i];
			if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1] = -1;
			}
			else { Pnt_semi.y.arr[i + 1] = 0; }
		}
		Pnt_semi.y.arr[i] += source->x.arr[i] - Pnt_semi.x.arr[i];
		for (i++; i <= Pnt_semi.y.regS; i++) { Pnt_semi.y.arr[i] = 0; }
		for (i = source->x.S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}

	//𝑦₃= ( 𝑥₃- 𝑥₁)𝜆 | kECC_modMultEql(&(Pnt_semi.y), &Pnt_ceta, p); [mod] {mtp, mod}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_ceta.S; i++) {
		for (j = 0; j <= Pnt_semi.y.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_semi.y.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_ceta.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_semi.y.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_ceta.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_semi.y.S = pKMS.S;
	memcpy(Pnt_semi.y.arr, pKMS.arr, Sze);

	//𝑦₃= ( 𝑥₃- 𝑥₁)𝜆 - 𝑦₁ | kECC_modSub(&(Pnt_semi.y), &(Pnt_semi.y), &(source->y), p); | [sub]
	fe = MAX(Pnt_semi.y.S, source->y.S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_semi.y.arr[i] < source->y.arr[i]) { i = -2; break; }
		else if (Pnt_semi.y.arr[i] > source->y.arr[i]) { break; }
	}
	if (i == -2) {
		for (i = 0; i <= p->S; i++) {
			Pnt_semi.y.arr[i] += p->arr[i] - source->y.arr[i];
			if (Pnt_semi.y.arr[i] >= PNUM) {
				Pnt_semi.y.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.y.arr[i + 1]++;
			}
			else if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1]--;
			}
		}
		for (i = p->S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}
	else {
		for (i = 0; i <= Pnt_semi.y.S; i++) {
			Pnt_semi.y.arr[i] -= source->y.arr[i];
			if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1]--;
			}
		}
		for (i = Pnt_semi.y.S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}

	Pnt_change = Pnt_semi;
	Pnt_semi = *source;
	*source = Pnt_change;
	return;
}

void kECC_pointDbl(POINT* result, POINT* source, KMS* a, KMS* p) {
	KMS pKMS;
	numData semiA, semiB, semiC, semiD;
	numData* arrStorage;
	int i, j, fe, idxp1 = 0, front, semi;
	unsigned long long bit1, bit2;
	pKMS.arr = Modmtp_bigArr;

	//𝜆 = 𝑥² | kECC_modSquare(&Pnt_ceta, &(source->x), p); | [mod] {mod, mtp}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= source->x.S; i++) {
		for (j = 0; j <= source->x.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (source->x.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (source->x.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (source->x.arr[j] >> SEMI_DIGIT);
			semiD = (source->x.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_ceta.S = pKMS.S;
	memcpy(Pnt_ceta.arr, pKMS.arr, Sze);

	//𝜆 = 3𝑥² | kECC_modDouble(&Pnt_t1, &Pnt_ceta, p); kECC_modAdd(&Pnt_ceta, &Pnt_t1, &Pnt_ceta, p); | [add] {sub}
	Pnt_ceta.arr[0] *= 3;
	for (i = 1; i <= Pnt_ceta.S; i++) {
		Pnt_ceta.arr[i] *= 3;
		Pnt_ceta.arr[i] += Pnt_ceta.arr[i - 1] >> DIGIT;
		Pnt_ceta.arr[i - 1] &= sem_binary[DIGIT];
	}
	for (; i <= Pnt_ceta.regS; i++) { Pnt_ceta.arr[i] = 0; }
	if (Pnt_ceta.arr[Pnt_ceta.S] >> DIGIT) {
		Pnt_ceta.arr[Pnt_ceta.S + 1] = Pnt_ceta.arr[Pnt_ceta.S] >> DIGIT;
		Pnt_ceta.arr[Pnt_ceta.S] &= sem_binary[DIGIT];
		Pnt_ceta.S++;
	}

	fe = MAX(Pnt_ceta.S, p->S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_ceta.arr[i] > p->arr[i]) { i = -2; break; }
		else if (Pnt_ceta.arr[i] < p->arr[i]) { break; }
	}
	if (i < 0) {
		for (i = 0; i <= Pnt_ceta.S; i++) {
			Pnt_ceta.arr[i] -= p->arr[i];
			if (Pnt_ceta.arr[i] < 0) {
				Pnt_ceta.arr[i] += PNUM;
				Pnt_ceta.arr[i + 1]--;
			}
		}
		if (!Pnt_ceta.arr[Pnt_ceta.S]) {
			for (i = Pnt_ceta.S; (i > 0) && (!Pnt_ceta.arr[i]); i--) { ; }
			Pnt_ceta.S = i;
		}
		// part 2;
		fe = MAX(Pnt_ceta.S, p->S);
		for (i = fe; i >= 0; i--) {
			if (Pnt_ceta.arr[i] > p->arr[i]) { i = -2; break; }
			else if (Pnt_ceta.arr[i] < p->arr[i]) { break; }
		}
		if (i < 0) {
			for (i = 0; i <= Pnt_ceta.S; i++) {
				Pnt_ceta.arr[i] -= p->arr[i];
				if (Pnt_ceta.arr[i] < 0) {
					Pnt_ceta.arr[i] += PNUM;
					Pnt_ceta.arr[i + 1]--;
				}
			}
			if (!Pnt_ceta.arr[Pnt_ceta.S]) {
				for (i = Pnt_ceta.S; (i > 0) && (!Pnt_ceta.arr[i]); i--) { ; }
				Pnt_ceta.S = i;
			}
		}
	}

	//𝜆 = 3𝑥²+ a | kECC_modAdd(&Pnt_ceta, &Pnt_ceta, a, p); | [add] {sub}
	for (i = 0; i <= a->S; i++) {
		Pnt_ceta.arr[i] += a->arr[i];
		if (Pnt_ceta.arr[i] >> DIGIT) {
			Pnt_ceta.arr[i] &= sem_binary[DIGIT];
			Pnt_ceta.arr[i + 1]++;
		}
	}
	if (Pnt_ceta.arr[Pnt_ceta.S + 1]) { Pnt_ceta.S++; }
	fe = MAX(Pnt_ceta.S, p->S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_ceta.arr[i] > p->arr[i]) { i = -2; break; }
		else if (Pnt_ceta.arr[i] < p->arr[i]) { break; }
	}
	if (i < 0) {
		for (i = 0; i <= Pnt_ceta.S; i++) {
			Pnt_ceta.arr[i] -= p->arr[i];
			if (Pnt_ceta.arr[i] < 0) {
				Pnt_ceta.arr[i] += PNUM;
				Pnt_ceta.arr[i + 1]--;
			}
		}
		if (!Pnt_ceta.arr[Pnt_ceta.S]) {
			for (i = Pnt_ceta.S; (i > 0) && (!Pnt_ceta.arr[i]); i--) { i; }
			Pnt_ceta.S = i;
		}
	}
	//t1 = 2𝑦 | kECC_modDouble(&Pnt_t1, &source->y, p); | [double] {sub}
	Pnt_t1.arr[0] = source->y.arr[0] << 1;
	Pnt_t1.S = source->y.S;
	for (i = 1; i <= source->y.S; i++) {
		Pnt_t1.arr[i] = source->y.arr[i] << 1;
		Pnt_t1.arr[i] |= Pnt_t1.arr[i - 1] >> DIGIT;
		Pnt_t1.arr[i - 1] &= sem_binary[DIGIT];
	}
	for (; i <= Pnt_t1.regS; i++) { Pnt_t1.arr[i] = 0; }
	if (Pnt_t1.arr[Pnt_t1.S] >= PNUM) {
		Pnt_t1.arr[Pnt_t1.S] &= sem_binary[DIGIT];
		Pnt_t1.S++;
		Pnt_t1.arr[Pnt_t1.S] = 1;
	}
	fe = MAX(p->S, Pnt_t1.S);
	for (i = fe; i >= 0; i--) {
		if (p->arr[i] < Pnt_t1.arr[i]) { i = -2; break; }
		else if (p->arr[i] > Pnt_t1.arr[i]) { break; }
	}
	if (i < 0) {
		for (i = 0; i <= Pnt_t1.S; i++) {
			Pnt_t1.arr[i] -= p->arr[i];
			if (Pnt_t1.arr[i] < 0) {
				Pnt_t1.arr[i] += PNUM;
				Pnt_t1.arr[i + 1]--;
			}
		}
		for (i = Pnt_t1.S; (i > 0) && (!Pnt_t1.arr[i]); i--) { ; }
		Pnt_t1.S = i;
	}

	//t1 = 1 / t1 ////////////////////////////////////////////////////////////////////////
	kECC_modInv(&Pnt_t1, &Pnt_t1, p);

	//𝜆 = 3𝑥²+ a / 2𝑦 | kECC_modMultEql(&Pnt_ceta, &Pnt_t1, p); [mod] {mtp, mod}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_t1.S; i++) {
		for (j = 0; j <= Pnt_ceta.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_ceta.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_t1.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_ceta.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_t1.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_ceta.S = pKMS.S;
	memcpy(Pnt_ceta.arr, pKMS.arr, Sze);

	//t1 = 𝜆²
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_ceta.S; i++) {
		for (j = 0; j <= Pnt_ceta.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_ceta.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_ceta.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_ceta.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_ceta.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_semi.x.S = pKMS.S;
	memcpy(Pnt_semi.x.arr, pKMS.arr, Sze);

	//𝑥₃= 2𝑥₁ | kECC_modDouble(&Pnt_t1, &source->x, p); | [double] {sub}
	Pnt_t1.arr[0] = source->x.arr[0] << 1;
	for (i = 1; i <= source->x.S; i++) {
		Pnt_t1.arr[i] = source->x.arr[i] << 1;
		Pnt_t1.arr[i] |= Pnt_t1.arr[i - 1] >> DIGIT;
		Pnt_t1.arr[i - 1] &= sem_binary[DIGIT];
	}
	for (; i <= Pnt_t1.regS; i++) { Pnt_t1.arr[i] = 0; }
	Pnt_t1.S = source->x.S;
	if (Pnt_t1.arr[Pnt_t1.S] >> DIGIT) {
		Pnt_t1.arr[Pnt_t1.S] &= sem_binary[DIGIT];
		Pnt_t1.S++;
		Pnt_t1.arr[Pnt_t1.S] = 1;
	}
	fe = MAX(Pnt_t1.S, p->S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_t1.arr[i] < p->arr[i]) { break; }
		else if (Pnt_t1.arr[i] > p->arr[i]) { i = -2; break; }
	}
	if (i < 0) {
		for (i = 0; i <= Pnt_t1.S; i++) {
			Pnt_t1.arr[i] -= p->arr[i];
			if (Pnt_t1.arr[i] < 0) {
				Pnt_t1.arr[i] += PNUM;
				Pnt_t1.arr[i + 1]--;
			}
		}
		if (!Pnt_ceta.arr[Pnt_ceta.S]) {
			for (i = Pnt_t1.S; (i > 0) && (!Pnt_t1.arr[i]); i--) { ; }
			Pnt_t1.S = i;
		}
	}

	//𝑥₃= 𝜆²- 2𝑥₁ | kECC_modSub(&(Pnt_semi.x), &Pnt_semi.x, &Pnt_t1, p); | [sub] {add}
	fe = MAX(Pnt_semi.x.S, Pnt_t1.S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_semi.x.arr[i] < Pnt_t1.arr[i]) { i = -2; break; }
		else if (Pnt_semi.x.arr[i] > Pnt_t1.arr[i]) { break; }
	}
	if (i == -2) {
		for (i = 0; i <= p->S; i++) {
			Pnt_semi.x.arr[i] += p->arr[i] - Pnt_t1.arr[i];
			if (Pnt_semi.x.arr[i] >= PNUM) {
				Pnt_semi.x.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.x.arr[i + 1]++;
			}
			else if (Pnt_semi.x.arr[i] < 0) {
				Pnt_semi.x.arr[i] += PNUM;
				Pnt_semi.x.arr[i + 1]--;
			}
		}
		for (i = p->S; (i > 0) && (!Pnt_semi.x.arr[i]); i--) { ; }
		Pnt_semi.x.S = i;
	}
	else {
		for (i = 0; i <= Pnt_semi.x.S; i++) {
			Pnt_semi.x.arr[i] -= Pnt_t1.arr[i];
			if (Pnt_semi.x.arr[i] < 0) {
				Pnt_semi.x.arr[i] += PNUM;
				Pnt_semi.x.arr[i + 1]--;
			}
		}
		if (!Pnt_semi.x.arr[Pnt_semi.x.S]) {
			for (i = Pnt_semi.x.S; (i > 0) && (!Pnt_semi.x.arr[i]); i--) { ; }
			Pnt_semi.x.S = i;
		}
	}

	//𝑦 = 𝑥₁- 𝑥₃ | kECC_modSub(&(Pnt_semi.y), &(source->x), &(Pnt_semi.x), p); | [sub] {cmp, add}
	fe = MAX(source->x.S, Pnt_semi.x.S);
	for (i = fe; i >= 0; i--) {
		if (source->x.arr[i] < Pnt_semi.x.arr[i]) { i = -2; break; }
		else if (source->x.arr[i] > Pnt_semi.x.arr[i]) { break; }
	}
	if (i == -2) {
		Pnt_semi.y.arr[0] = 0;
		for (i = 0; i < p->S; i++) {
			Pnt_semi.y.arr[i] += p->arr[i] + source->x.arr[i] - Pnt_semi.x.arr[i];
			if (Pnt_semi.y.arr[i] >= PNUM) {
				Pnt_semi.y.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.y.arr[i + 1] = 1;
			}
			else if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1] = -1;
			}
			else { Pnt_semi.y.arr[i + 1] = 0; }
		}
		Pnt_semi.y.arr[i] += p->arr[i] + source->x.arr[i] - Pnt_semi.x.arr[i];
		for (i++; i <= Pnt_semi.y.regS; i++) { Pnt_semi.y.arr[i] = 0; }
		for (i = p->S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}
	else {
		Pnt_semi.y.arr[0] = 0;
		for (i = 0; i < source->x.S; i++) {
			Pnt_semi.y.arr[i] += source->x.arr[i] - Pnt_semi.x.arr[i];
			if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1] = -1;
			}
			else { Pnt_semi.y.arr[i + 1] = 0; }
		}
		Pnt_semi.y.arr[i] += source->x.arr[i] - Pnt_semi.x.arr[i];
		for (i++; i <= Pnt_semi.y.regS; i++) { Pnt_semi.y.arr[i] = 0; }
		for (i = source->x.S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}

	//𝑦 = ( 𝑥₃- 𝑥₁)𝜆 | kECC_modMultEql(&(Pnt_semi.y), &Pnt_ceta, p); [mod] {mtp, mod}
	pKMS.regS = BIG_SIZE - 1;
	memset(pKMS.arr, 0, BSze);
	for (i = 0; i <= Pnt_ceta.S; i++) {
		for (j = 0; j <= Pnt_semi.y.S; j++) {
			idxp1 = MIN(i + j + 1, pKMS.regS);
			semiA = (Pnt_semi.y.arr[j] & sem_binary[SEMI_DIGIT]);
			semiB = (Pnt_ceta.arr[i] & sem_binary[SEMI_DIGIT]);
			semiC = (Pnt_semi.y.arr[j] >> SEMI_DIGIT);
			semiD = (Pnt_ceta.arr[i] >> SEMI_DIGIT);
			pKMS.arr[i + j] += semiA * semiB;
			pKMS.arr[i + j] += (semiC * semiB & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiC * semiB >> SEMI_DIGIT;
			pKMS.arr[i + j] += (semiA * semiD & sem_binary[SEMI_DIGIT]) << SEMI_DIGIT;
			pKMS.arr[idxp1] += semiA * semiD >> SEMI_DIGIT;
			pKMS.arr[idxp1] += pKMS.arr[i + j] >> DIGIT;
			pKMS.arr[i + j] &= sem_binary[DIGIT];
			pKMS.arr[idxp1] += semiC * semiD;
		}
	}
	for (i = idxp1; (i > 0) && (!pKMS.arr[i]); i--) { ; }
	pKMS.S = i;
	// part of " | kECC_mmodEql(&pKMS, p);
	arrStorage = p->arr;
	if (p->S != pKMS.S) { i = (p->S > pKMS.S) ? -2 : 0; }
	else {
		for (i = p->S; i >= 0; i--) {
			if (p->arr[i] < pKMS.arr[i]) { break; }
			else if (p->arr[i] > pKMS.arr[i]) { i = -2; break; }
		}
	}
	if (i != -2) {
		p->arr = Modeql_bigArr;
		memset(p->arr + SML_SIZE, 0, BSze - Sze);
		memcpy(p->arr, arrStorage, Sze);
		front = MAX(pKMS.S - p->S - 1, 0);
		if ((p->S < pKMS.S) && (p->arr[p->S] < pKMS.arr[pKMS.S])) { front++; }
		if (front) {
			for (i = p->S + front; i >= front; i--) { p->arr[i] = p->arr[i - front]; }
			for (i = 0; i < front; i++) { p->arr[i] = 0; }
			p->S += front;
		}
		front *= DIGIT;
		//semi = bitDigit(pKMS.arr[pKMS.S]) - bitDigit(p->arr[p->S]);
		bit1 = pKMS.arr[pKMS.S];
		bit2 = p->arr[p->S];
		PROCESS_KECC_BDGT1(bit1);
		PROCESS_KECC_BDGT1(bit2);
		semi = tab64[PROCESS_KECC_BDGT2(bit1)] - tab64[PROCESS_KECC_BDGT2(bit2)];
		if (pKMS.S != p->S) { semi += DIGIT; }
		if (semi && (pKMS.S != p->S) ? ((p->arr[p->S] >> (DIGIT - semi)) > pKMS.arr[pKMS.S]) : ((p->arr[p->S] << semi) > pKMS.arr[pKMS.S])) { semi--; }
		if (semi >= 0) {

			i = p->S;
			if (p->arr[i] >= _binary[DIGIT - semi]) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << semi;
				p->S++;
			}
			else { p->arr[i] <<= semi; }
			for (i--; i >= 0; i--) {
				p->arr[i + 1] |= p->arr[i] >> (DIGIT - semi);
				p->arr[i] = (p->arr[i] & sem_binary[DIGIT - semi]) << (semi);
			}
			front += semi;

			for (;;) {
				if (pKMS.S != p->S) { i = (pKMS.S > p->S) ? -1 : 1; }
				else {
					for (i = MAX(pKMS.S, p->S); i >= 0; i--) {
						if (pKMS.arr[i] > p->arr[i]) { i = -1; break; }
						else if (pKMS.arr[i] < p->arr[i]) { break; }
					}
				}
				if (i < 0) {
					for (i = 0; i <= pKMS.S; i++) {
						pKMS.arr[i] -= p->arr[i];
						if (pKMS.arr[i] < 0) {
							pKMS.arr[i] += PNUM;
							pKMS.arr[i + 1]--;
						}
					}
					for (i = pKMS.S; (i > 0) && (!pKMS.arr[i]); i--) { ; }
					pKMS.S = i;
				}
				if (!front) { break; }

				for (i = 0; i < p->S; i++) {
					p->arr[i] >>= 1;
					p->arr[i] |= (p->arr[i + 1] & 1) << (DIGIT - 1);
				}
				p->arr[p->S] >>= 1;
				if (!(p->arr[p->S])) { p->S--; }
				front--;
			}
		}
		p->arr = arrStorage;
	}
	Pnt_semi.y.S = pKMS.S;
	memcpy(Pnt_semi.y.arr, pKMS.arr, Sze);

	//𝑦 = ( 𝑥₃- 𝑥₁)𝜆 - 𝑦₁ | kECC_modSub(&(Pnt_semi.y), &(Pnt_semi.y), &(source->y), p); | [sub]
	fe = MAX(Pnt_semi.y.S, source->y.S);
	for (i = fe; i >= 0; i--) {
		if (Pnt_semi.y.arr[i] < source->y.arr[i]) { i = -2; break; }
		else if (Pnt_semi.y.arr[i] > source->y.arr[i]) { break; }
	}
	if (i == -2) {
		for (i = 0; i <= p->S; i++) {
			Pnt_semi.y.arr[i] += p->arr[i] - source->y.arr[i];
			if (Pnt_semi.y.arr[i] >= PNUM) {
				Pnt_semi.y.arr[i] &= sem_binary[DIGIT];
				Pnt_semi.y.arr[i + 1]++;
			}
			else if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1]--;
			}
		}
		for (i = p->S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}
	else {
		for (i = 0; i <= Pnt_semi.y.S; i++) {
			Pnt_semi.y.arr[i] -= source->y.arr[i];
			if (Pnt_semi.y.arr[i] < 0) {
				Pnt_semi.y.arr[i] += PNUM;
				Pnt_semi.y.arr[i + 1]--;
			}
		}
		for (i = Pnt_semi.y.S; (i > 0) && (!Pnt_semi.y.arr[i]); i--) { ; }
		Pnt_semi.y.S = i;
	}

	Pnt_change = Pnt_semi;
	Pnt_semi = *result;
	*result = Pnt_change;
	return;
}

void kECC_registArr(CURVE* curve, int digit) {
	if (sem_mtpValue != NULL) { free(sem_mtpValue); }
	sem_mtpValue = (POINT*)malloc(sizeof(POINT) * digit);
	sem_mtpValue[0].x.arr = (numData*)malloc(Sze);
	sem_mtpValue[0].y.arr = (numData*)malloc(Sze);
	sem_mtpValue[0].x.regS = sem_mtpValue[0].y.regS = Reg;
	sem_mtpValue[0].x.S = curve->G.x.S;
	sem_mtpValue[0].y.S = curve->G.y.S;
	memcpy(sem_mtpValue[0].x.arr, curve->G.x.arr, Sze);
	memcpy(sem_mtpValue[0].y.arr, curve->G.y.arr, Sze);
	for (int i = 1; i < digit; i++) {
		sem_mtpValue[i].x.arr = (numData*)malloc(Sze);
		sem_mtpValue[i].y.arr = (numData*)malloc(Sze);
		sem_mtpValue[i].x.regS = sem_mtpValue[i].y.regS = Reg;
		kECC_pointDbl(&sem_mtpValue[i], &sem_mtpValue[i - 1], &(curve->a), &(curve->p));
	}
	return;
}

void ARRAY_FREE(void) { free(sem_mtpValue); }

void kECC_PublicKey(POINT* Public, KMS* Private, CURVE* curve) {
	int Eql = 1;
	int BitM, i, j, idx = 0;
	unsigned long long Bt;
	for (i = 0; i < Private->S; i++) {
		for (j = 0; j < DIGIT; j++) {
			if ((Private->arr[i] >> j) & 1) {
				if (Eql) {
					memcpy(Public->x.arr, sem_mtpValue[idx].x.arr, Sze);
					memcpy(Public->y.arr, sem_mtpValue[idx].y.arr, Sze);
					Public->x.S = sem_mtpValue[idx].x.S;
					Public->y.S = sem_mtpValue[idx].y.S;
					Eql = 0;
				}
				else { kECC_pointAddEql(Public, sem_mtpValue + idx, &(curve->p)); }
			}
			idx++;
		}
	}
	Bt = Private->arr[Private->S];
	PROCESS_KECC_BDGT1(Bt);
	BitM = tab64[PROCESS_KECC_BDGT2(Bt)];
	for (j = 0; j < BitM; j++) {
		if ((Private->arr[i] >> j) & 1) {
			if (Eql) {
				memcpy(Public->x.arr, sem_mtpValue[idx].x.arr, Sze);
				memcpy(Public->y.arr, sem_mtpValue[idx].y.arr, Sze);
				Public->x.S = sem_mtpValue[idx].x.S;
				Public->y.S = sem_mtpValue[idx].y.S;
				Eql = 0;
			}
			else { kECC_pointAddEql(Public, sem_mtpValue + idx, &(curve->p)); }
		}
		idx++;
	}
	return;
}

void kECC_SecretKey(POINT* Secret, POINT* Public, KMS* Private, CURVE* curve) {
	int Eql = 1;
	int BitM, i, j;
	unsigned long long Bt;
	Sec_semi.x.S = Public->x.S;
	Sec_semi.y.S = Public->y.S;
	memcpy(Sec_semi.x.arr, Public->x.arr, Sze);
	memcpy(Sec_semi.y.arr, Public->y.arr, Sze);
	for (i = 0; i < Private->S; i++) {
		for (j = 0; j < DIGIT; j++) {
			if ((Private->arr[i] >> j) & 1) {
				if (Eql) {
					memcpy(Secret->x.arr, Sec_semi.x.arr, Sze);
					memcpy(Secret->y.arr, Sec_semi.y.arr, Sze);
					Secret->x.S = Sec_semi.x.S;
					Secret->y.S = Sec_semi.y.S;
					Eql = 0;
				}
				else { kECC_pointAddEql(Secret, &Sec_semi, &(curve->p)); }
			}
			kECC_pointDbl(&Sec_semi, &Sec_semi, &(curve->a), &(curve->p));
		}
	}
	Bt = Private->arr[Private->S];
	PROCESS_KECC_BDGT1(Bt);
	BitM = tab64[PROCESS_KECC_BDGT2(Bt)];
	for (j = 0; j < BitM; j++) {
		if ((Private->arr[i] >> j) & 1) {
			if (Eql) {
				memcpy(Secret->x.arr, Sec_semi.x.arr, Sze);
				memcpy(Secret->y.arr, Sec_semi.y.arr, Sze);
				Secret->x.S = Sec_semi.x.S;
				Secret->y.S = Sec_semi.y.S;
				Eql = 0;
			}
			else { kECC_pointAddEql(Secret, &Sec_semi, &(curve->p)); }
		}
		kECC_pointDbl(&Sec_semi, &Sec_semi, &(curve->a), &(curve->p));
	}
	return;
}

#endif