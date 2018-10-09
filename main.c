#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define ROTL8(x,shift) ((uint8_t) ((x) << (shift)) | ((x) >> (8 - (shift))))

#define XTIME(x) ((uint8_t) ((x<<1) ^ ((x>>7)*0x1B)))

uint8_t RCON[10] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1B,0x36};

void initialize_aes_sbox(uint8_t sbox[256]) {
    uint8_t p = 1, q = 1;

    /* loop invariant: p * q == 1 in the Galois field */
    do {
        /* multiply p by 3 */
        p = p ^ (p << 1) ^ (p & 0x80 ? 0x1B : 0);

        /* divide q by 3 (equals multiplication by 0xf6) */
        q ^= q << 1;
        q ^= q << 2;
        q ^= q << 4;
        q ^= q & 0x80 ? 0x09 : 0;

        /* compute the affine transformation */
        uint8_t xformed = q ^ ROTL8(q, 1) ^ ROTL8(q, 2) ^ ROTL8(q, 3) ^ ROTL8(q, 4);

        sbox[p] = xformed ^ 0x63;
    } while (p != 1);

    /* 0 is a special case since it has no inverse */
    sbox[0] = 0x63;
}

uint8_t** strToMat(char* str){
    uint8_t **mat = (uint8_t **) malloc(sizeof(uint8_t *)*4);
    char *substr = (char*) malloc(sizeof(char)*2);

    for(uint8_t i =0; i < 4; i++){
        mat[i] = (uint8_t *) malloc(sizeof(uint8_t)*4);
    }

    for(uint8_t i =0; i<4; i++){
        for(uint8_t j = 0; j < 4; j++){
            memcpy(substr, &str[j*8+i*2], 2);
            mat[i][j] = (uint8_t) strtol(substr, NULL, 16);
        }
    }

    free(substr);

    return mat;
}

void initialize_sub_keys(uint8_t*** matKeys, char* key, const uint8_t* sBox){
    for (int i = 0; i < 11; ++i) {
        matKeys[i] = (uint8_t**) malloc(sizeof(uint8_t*)*4);
        if ( i != 0){
            for (int j = 0; j < 4; ++j) {
                matKeys[i][j] = (uint8_t*) malloc(sizeof(uint8_t)*4);
            }
        }
    }

    matKeys[0] = strToMat(key);

    for (int i = 1; i < 11; ++i) {
        for (int j = 0; j < 4; ++j) {
            if(j == 0){
                for (int k = 0; k < 4; ++k) {
                    matKeys[i][k][j] = matKeys[i-1][(k+1)%4][3];
                    matKeys[i][k][j] = (uint8_t) sBox[matKeys[i][k][j]/16*16+matKeys[i][k][j]%16];
                    matKeys[i][k][j] =  matKeys[i][k][j] ^ matKeys[i-1][k][j];
                    if (k == 0){
                        matKeys[i][k][j] = matKeys[i][k][j] ^ RCON[i-1];
                    }
                }
            } else {
                for (int k = 0; k < 4; ++k) {
                    matKeys[i][k][j] = matKeys[i-1][k][j] ^ matKeys[i][k][j-1];
                }
            }
        }
    }


};

void displayMat(uint8_t **mat){
    for (uint8_t i = 0; i < 4; ++i) {
        for (uint8_t j = 0; j < 4; ++j) {
            printf("%X\t", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");

}

void addRoundKey(uint8_t **mat1, uint8_t **mat2){
    for (uint8_t i = 0; i < 4; ++i) {
        for (uint8_t j = 0; j < 4; ++j) {
            mat1[i][j] = mat1[i][j] ^ mat2[i][j];
        }
    }
}

void subBytes(uint8_t **mat, const uint8_t* sBox){
    for (uint8_t i = 0; i < 4; ++i) {
        for (uint8_t j = 0; j < 4; ++j) {
            mat[i][j] = (uint8_t) sBox[mat[i][j]/16*16+mat[i][j]%16];
        }
    }
}

void shiftRows(uint8_t **mat){
    uint8_t temp;
    for (uint8_t i = 0; i < 4; ++i) {
        for (uint8_t j = 0; j < i; ++j) {
            temp = mat[i][0];

            mat[i][0] = mat[i][1];
            mat[i][1] = mat[i][2];
            mat[i][2] = mat[i][3];
            mat[i][3] = temp;
        }
    }
}

void mixedColumns(uint8_t **mat){
    uint8_t* temp = (uint8_t*) malloc(sizeof(uint8_t));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            temp[j] = mat[j][i];
        }

        for (int j= 0; j < 4; ++j) {
            mat[j][i] = XTIME(temp[j])^temp[(j+1)%4]^XTIME(temp[(j+1)%4])^temp[(j+2)%4]^temp[(j+3)%4];
        }
    }

    free(temp);
}

int main() {
    char *plain = "3243f6a8885a308d313198a2e0370734";
    char *key = "2b7e151628aed2a6abf7158809cf4f3c";

    uint8_t ** inputMat = strToMat(plain);
    uint8_t *** subKeysMat = (uint8_t***) malloc(sizeof(uint8_t**)*11);
    uint8_t* sBox = (uint8_t*) malloc(sizeof(uint8_t)*256);

    initialize_aes_sbox(sBox);
    initialize_sub_keys(subKeysMat, key, sBox);

    //initial round
    addRoundKey(inputMat, subKeysMat[0]); //

    //9 rounds
    for (int i = 0; i < 9; ++i) {
        subBytes(inputMat, sBox);
        shiftRows(inputMat);
        mixedColumns(inputMat);
        addRoundKey(inputMat, subKeysMat[i+1]);
    }

    //final round
    subBytes(inputMat, sBox);
    shiftRows(inputMat);
    addRoundKey(inputMat, subKeysMat[10]);

    displayMat(inputMat);

    return 0;
}