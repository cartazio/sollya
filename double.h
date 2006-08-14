#ifndef DOUBLE_H
#define DOUBLE_H


int mpfr_round_to_double(mpfr_t rop, mpfr_t op);
int mpfr_round_to_doubledouble(mpfr_t rop, mpfr_t op);
int mpfr_round_to_tripledouble(mpfr_t rop, mpfr_t op);
int printDoubleInHexa(mpfr_t x);
int readHexa(mpfr_t res, char *c);

#endif /* ifdef DOUBLE_H*/
