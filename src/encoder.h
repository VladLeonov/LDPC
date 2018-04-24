#include "matrix.h"
#include "ldpc_generator.h"

#ifndef ENCODER_H_INCLUDED
#define ENCODER_H_INCLUDED

matrix encode_message(ldpc ldpc_object, matrix message);

#endif // ENCODER_H_INCLUDED
