/*

  Copyright 2012-2015 by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  and by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors Ch. Lauter, S. Chevillard

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularily targeted to the automatized implementation of
  mathematical floating-point libraries (libm). Amongst other features,
  it offers a certified infinity norm, an automatic polynomial
  implementer and a fast Remez algorithm.

  This software is governed by the CeCILL-C license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-C
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-C license and that you accept its terms.

  This program is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*/

#ifndef SOLLYA_MESSAGES_H
#define SOLLYA_MESSAGES_H

/* Symbolic constants for the different messages in Sollya

   SOLLYA DEVELOPPERS: WHEN YOU ADD A NEW MESSAGE NUMBER, DO NOT
   FORGET TO ADD A CORRESPONDING MESSAGE TEXT IN sollya-messaging.h !

*/
#define SOLLYA_MSG_NO_MSG                                            0       /* Void error message, no error at all */
#define SOLLYA_MSG_CONTINUATION                                      1
#define SOLLYA_MSG_UNDEFINED_ERROR                                   2
#define SOLLYA_MSG_ABS_NOT_TWICE_DIFFERENTIABLE                      3
#define SOLLYA_MSG_HALF_NOT_DIFFERENTIABLE                           4
#define SOLLYA_MSG_SINGLE_NOT_DIFFERENTIABLE                         5
#define SOLLYA_MSG_DOUBLE_NOT_DIFFERENTIABLE                         6
#define SOLLYA_MSG_DOUBLEEXTENDED_NOT_DIFFERENTIABLE                 7
#define SOLLYA_MSG_DOUBLE_DOUBLE_NOT_DIFFERENTIABLE                  8
#define SOLLYA_MSG_TRIPLE_DOUBLE_NOT_DIFFERENTIABLE                  9
#define SOLLYA_MSG_QUAD_NOT_DIFFERENTIABLE                          10
#define SOLLYA_MSG_CEIL_NOT_DIFFERENTIABLE                          11
#define SOLLYA_MSG_FLOOR_NOT_DIFFERENTIABLE                         12
#define SOLLYA_MSG_NEARESTINT_NOT_DIFFERENTIABLE                    13
#define SOLLYA_MSG_UNDESIRED_ROUNDING_IN_ROUND_TO_FORMAT            14
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_TO_DOUBLE               15
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_TO_PREC                 16
#define SOLLYA_MSG_ROUNDING_OCCURRED_WHILE_CONVERTING_FROM_DOUBLE   17
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_TO_SINGLE               18
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_TO_DOUBLE_DOUBLE        19
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_TO_TRIPLE_DOUBLE        20
#define SOLLYA_MSG_ROUNDING_DOWN_BEFORE_PRINTING_DOUBLE             22
#define SOLLYA_MSG_ROUNDING_UP_BEFORE_PRINTING_DOUBLE               23
#define SOLLYA_MSG_COULD_NOT_FIGURE_OUT_ENDIANESS                   24
#define SOLLYA_MSG_ROUNDING_DOWN_BEFORE_PRINTING_SINGLE             25
#define SOLLYA_MSG_ROUNDING_UP_BEFORE_PRINTING_SINGLE               26
#define SOLLYA_MSG_SNAN_MIGHT_HAVE_BECOME_QNAN                      27
#define SOLLYA_MSG_UNABLE_TO_CONVERT_FROM_HEXADECIMAL_CONSTANT      28
#define SOLLYA_MSG_GIVEN_FUNCTION_IS_NO_POLYNOMIAL                  29
#define SOLLYA_MSG_NUM_OF_FORMATS_DOES_NOT_CORRESPOND_TO_DEGREE     30
#define SOLLYA_MSG_ERROR_WHILE_EXTRACTING_COEFFICIENTS_OF_POLY      31
#define SOLLYA_MSG_ERROR_POLY_COEFF_GETS_ROUNDED                    32
#define SOLLYA_MSG_ERROR_EVALUATION_OF_POLY_COEFF_NOT_FAITHFUL      33
#define SOLLYA_MSG_ERROR_UNKNOWN_EXPANSION_FORMAT                   34
#define SOLLYA_MSG_DOUBLE_ROUNDING_ON_HANDLING_POLY_COEFF           35
#define SOLLYA_MSG_AT_LEAST_ONE_POLY_COEFF_HAS_BEEN_ROUNDED         36
#define SOLLYA_MSG_NONE_OF_THE_POLY_COEFFS_HAS_BEEN_ROUNDED         37
#define SOLLYA_MSG_NON_REAL_NUMBER_CANNOT_BE_DOUBLE_EXPANSION       38
#define SOLLYA_MSG_DOUBLE_ROUNDING_ON_CONVERSION                    39
#define SOLLYA_MSG_DOUBLE_EXPANSION_INCOMPLETE                      40
#define SOLLYA_MSG_ROUNDING_WHILE_SIMPLIFYING_TO_POLYNOMIAL         41
#define SOLLYA_MSG_POLY_COEFF_IS_NOT_CONSTANT                       42
#define SOLLYA_MSG_SOME_EVALUATION_IS_NOT_FAITHFUL                  43
#define SOLLYA_MSG_DOUBLE_ROUNDING_IN_ROUND_IEEE_754_2008_OPERATOR  44
#define SOLLYA_MSG_STRING_CANNOT_BE_PARSED_BY_MINIPARSER            45
#define SOLLYA_MSG_ERROR_ON_RUNNING_GUESSDEGREE                     46
#define SOLLYA_MSG_FILE_EXECUTION_ASKED_FOR_QUIT_NOT_QUITTING       47
#define SOLLYA_MSG_EXPR_SHOULD_BE_CONSTANT_AND_SEEMS_CONSTANT       48
#define SOLLYA_MSG_EXPR_SHOULD_BE_CONSTANT_AND_IS_CONSTANT_ON_FP    49
#define SOLLYA_MSG_EXPR_SHOULD_BE_CONSTANT_AND_IS_NOT_FAITHFUL      50
#define SOLLYA_MSG_EXPR_SHOULD_BE_CONSTANT_NO_FAITHFUL_PLAIN_FP     51
#define SOLLYA_MSG_FAITHFUL_ROUNDING_FOR_EXPR_THAT_SHOULD_BE_CONST  52
#define SOLLYA_MSG_CONSTANT_IS_NOT_MACHINE_INTEGER_WILL_ROUND       53
#define SOLLYA_MSG_PRECISION_OF_NUMBERS_MUST_BE_AT_LEAST_TWO_BITS   54
#define SOLLYA_MSG_IDENTIFIER_ALREADY_BOUND                         55
#define SOLLYA_MSG_IDENTIFIER_REASSIGNMENT                          56
#define SOLLYA_MSG_DISPLAYED_VALUE_IS_FAITHFULLY_ROUNDED            57
#define SOLLYA_MSG_EXPRESSION_UNDEFINED_OR_UNSTABLE                 58
#define SOLLYA_MSG_ROUNDING_MAY_HAVE_HAPPENED_AND_NOT_FAITHFUL      59
#define SOLLYA_MSG_EVALUATION_WITH_PLAIN_FP_ARITHMETIC              60
#define SOLLYA_MSG_ROUNDING_MAY_HAVE_HAPPENED_SOMEWHERE             61
#define SOLLYA_MSG_EXPRESSION_TOO_BIG_FOR_HORNER_FORM               62
#define SOLLYA_MSG_EXPRESSION_TOO_BIG_FOR_CANONICAL_FORM            63
#define SOLLYA_MSG_DOUBLE_SIMPLIFICATION_NECESSARY                  64
#define SOLLYA_MSG_COMMAND_NOT_EXECUTABLE                           65
#define SOLLYA_MSG_TIMER_UNUSABLE                                   66
#define SOLLYA_MSG_CAN_MODIFY_ONLY_ELEMENTS_OF_STRUCTURES           67
#define SOLLYA_MSG_CONTROL_STRUCTURE_NOT_EXECUTABLE_EXPR_NO_BOOLEAN 68
#define SOLLYA_MSG_AT_END_OF_FOR_CNTRL_VAR_NO_LONGER_ASSIGNABLE     69
#define SOLLYA_MSG_AT_END_OF_FOR_CNTRL_VAR_NO_LONGER_CONSTANT       70
#define SOLLYA_MSG_TOOL_HAS_BEEN_RESTARTED_INSIDE_LOOP              71
#define SOLLYA_MSG_CNTRL_VAR_OF_LOOP_CANNOT_BE_ASSIGNED             72
#define SOLLYA_MSG_ARGS_OF_FOR_LOOP_NOT_CONSTANT                    73
#define SOLLYA_MSG_FOR_IN_LOOP_OVER_EMPTY_LIST                      74
#define SOLLYA_MSG_FOR_IN_LOOP_OVER_END_ELLIPTIC_LIST_NOT_ALLOWED   75
#define SOLLYA_MSG_AT_LEAST_ONE_OPERATION_MUST_BE_EXECUTED          76
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_MACHINE_INTEGER        77
#define SOLLYA_MSG_RESTART_IN_FILE_READ_INTO_ANOTHER                78
#define SOLLYA_MSG_QUIT_IN_FILE_READ_INTO_ANOTHER                   79
#define SOLLYA_MSG_IDENTIFIER_IS_FREE_VAR_CANNOT_BE_LOCAL_VAR       80
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_FUNC_CANNOT_BE_LOCAL_VAR   81
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_CONST_CANNOT_BE_LOCAL_VAR  82
#define SOLLYA_MSG_IDENTIFIER_IS_EXTERNAL_PROC_CANNOT_BE_LOCAL_VAR  83
#define SOLLYA_MSG_FRAME_SYSTEM_CORRUPTED_LOCAL_VAR_NOT_DECLARED    84
#define SOLLYA_MSG_FILE_COULD_NOT_BE_OPENED_FOR_WRITING             85
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_STRING                 86
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_LIST_OF_FUNCTIONS      87
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_RANGE                  88
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_STRING_NOR_DEFAULT     89
#define SOLLYA_MSG_IMPLEMENTATION_HAS_NOT_SUCCEEDED                 90
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_CONSTANT               91
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_AN_EXPRESSION          92
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_A_FUNCTION             93
#define SOLLYA_MSG_ROUNDING_WHILE_PRINTING                          94
#define SOLLYA_MSG_BASH_RETURNS_A_CERTAIN_RETURN_VALUE              95
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_ABS_NOR_REL            96
#define SOLLYA_MSG_FILE_COULD_NOT_BE_OPENED_FOR_APPENDING           97
#define SOLLYA_MSG_FILE_COULD_NOT_BE_OPENED_FOR_READING             98
#define SOLLYA_MSG_FREE_VARIABLE_HAS_BEEN_NAMED_SOMEHOW             99
#define SOLLYA_MSG_FREE_VARIABLE_HAS_BEEN_RENAMED                  100
#define SOLLYA_MSG_CAN_RENAME_ONLY_FREE_VARIABLE                   101
#define SOLLYA_MSG_IDENTIFIER_IS_FREE_VAR_CANNOT_BE_EXTERNAL       102
#define SOLLYA_MSG_IDENTIFIER_IS_BOUND_TO_VAR_CANNOT_BE_EXTERNAL   103
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_CONST_CANNOT_BE_EXTERNAL  104
#define SOLLYA_MSG_IDENTIFIER_IS_EXTERNAL_PROC_CANNOT_BE_EXTERNAL  105
#define SOLLYA_MSG_ERROR_OCCURRED_COMMAND_NOT_EXECUTED             106
#define SOLLYA_MSG_ASSIGNMENT_WILL_HAVE_NO_EFFECT                  107
#define SOLLYA_MSG_ASSIGNMENT_OF_INDEXED_ELEMENTS_NOT_IN_RANGE     108
#define SOLLYA_MSG_ASSIGNMENT_OF_INDEXED_EMPTY_LIST_ONLY_ON_ZERO   109
#define SOLLYA_MSG_STRING_NOT_OF_LENGTH_ONE                        110
#define SOLLYA_MSG_IDENTIFIER_NOT_BOUND_TO_LIST_OR_STRING          111
#define SOLLYA_MSG_IDENTIFIER_NOT_ASSIGNED_TO                      112
#define SOLLYA_MSG_FIRST_ELMENT_OF_LEFT_SIDE_NOT_IDENTIFIER        113
#define SOLLYA_MSG_LEFT_HAND_SIDE_NOT_ELEMENT_OF_STRUCTURED_TYPE   114
#define SOLLYA_MSG_IDENTIFIER_IS_FREE_VAR_CANNOT_BE_LIBRARY        115
#define SOLLYA_MSG_IDENTIFIER_IS_BOUND_TO_VAR_CANNOT_BE_LIBRARY    116
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_CONST_CANNOT_BE_LIBRARY   117
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_FUNC_CANNOT_BE_LIBRARY    118
#define SOLLYA_MSG_IDENTIFIER_IS_EXTERNAL_PROC_CANNOT_BE_LIBRARY   119
#define SOLLYA_MSG_PREC_MUST_BE_AT_LEAST_TWELVE_BITS               120
#define SOLLYA_MSG_POINTS_MUST_BE_AT_LEAST_THREE_POINTS            121
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_DISPLAY_TYPE          122
#define SOLLYA_MSG_VERBOSITY_MUST_NOT_BE_NEGATIVE                  123
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_ON_OR_OFF             124
#define SOLLYA_MSG_TAYLOR_RECURSIONS_MUST_NOT_BE_NEGATIVE          125
#define SOLLYA_MSG_HOPITAL_RECURSIONS_MUST_NOT_BE_NEGATIVE         126
#define SOLLYA_MSG_EXPR_OR_COMMAND_COULD_NOT_BE_HANDLED            127
#define SOLLYA_MSG_EXPR_NOT_CORRECTLY_TYPED                        128
#define SOLLYA_MSG_EVALUATION_CREATES_ERROR_SPECIAL_SYMBOL         129
#define SOLLYA_MSG_EXPR_TOO_BIG_FOR_AUTOMATIC_SIMPLIFICATION       130
#define SOLLYA_MSG_FPMINIMAX_LESS_FORMATS_THAN_MONOMIALS           131
#define SOLLYA_MSG_FPMINIMAX_LESS_MONOMIALS_THAN_FORMATS           132
#define SOLLYA_MSG_FPMINIMAX_FORMAT_LIST_MALFORMED                 133
#define SOLLYA_MSG_FPMINIMAX_FORMAT_NEGATIVE_FOR_FP_COEFFS         134
#define SOLLYA_MSG_QUIT_OR_RESTART_MUST_NOT_BE_IN_MATCH            135
#define SOLLYA_MSG_FRAME_SYSTEM_CORRUPTED_MATCH_NOT_EXECUTED       136
#define SOLLYA_MSG_IDENTIFIER_IS_FREE_VAR_CANNOT_BE_MATCHED        137
#define SOLLYA_MSG_IDENTIFIER_IS_EXTERNAL_PROC_CANNOT_BE_MATCHED   138
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_CONST_CANNOT_BE_MATCHED   139
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_FUNC_CANNOT_BE_MATCHED    140
#define SOLLYA_MSG_NO_MATCHING_CASE_FOUND                          141
#define SOLLYA_MSG_QUIT_OR_RESTART_MUST_NOT_BE_IN_PROC             142
#define SOLLYA_MSG_IDENTIFIER_IS_FREE_VAR_CANNOT_BE_PARAMETER      143
#define SOLLYA_MSG_IDENTIFIER_IS_EXTERNAL_PROC_CANNOT_BE_PARAMETER 144
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_CONST_CANNOT_BE_PARAMETER 145
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_FUNC_CANNOT_BE_PARAMETER  146
#define SOLLYA_MSG_FRAME_SYSTEM_CORRUPTED_PROC_NOT_EXECUTED        147
#define SOLLYA_MSG_REMEZ_MONOMIAL_DEGREES_MUST_NOT_BE_NEGATIVE     148
#define SOLLYA_MSG_REMEZ_MONOMIAL_DEGREE_GIVEN_TWICE               149
#define SOLLYA_MSG_REMEZ_TOO_MANY_ARGUMENTS                        151
#define SOLLYA_MSG_REMEZ_SECND_ARG_MUST_BE_NONNEGATIVE_INT_OR_LIST 152
#define SOLLYA_MSG_FPMINIMAX_TOO_MANY_ARGUMENTS                    153
#define SOLLYA_MSG_INVALID_FIFTH_ARGUMENT                          154
#define SOLLYA_MSG_FPMINIMAX_SECND_ARG_MUST_BE_NONNEG_INT_OR_LIST  155
#define SOLLYA_MSG_FPMINIMAX_THIRD_ARG_MUST_BE_FORMAT_INDICATIONS  156
#define SOLLYA_MSG_FPMINIMAX_FOURTH_ARG_INTERVAL_OR_LIST_OF_POINTS 157
#define SOLLYA_MSG_TEST_COMPARES_ERROR_TO_SOMETHING                158
#define SOLLYA_MSG_TEST_RELIES_ON_FP_RESULT_THAT_IS_NOT_FAITHFUL   159
#define SOLLYA_MSG_TEST_RELIES_ON_FP_RESULT                        160
#define SOLLYA_MSG_TEST_RELIES_ON_FP_RESULT_FAITHFUL_BUT_UNDECIDED 161
#define SOLLYA_MSG_TEST_RELIES_ON_FP_RESULT_FAITHFUL_BUT_NOT_REAL  162
#define SOLLYA_MSG_MIN_RELIES_ON_FP_RESULT_THAT_IS_NOT_FAITHFUL    163
#define SOLLYA_MSG_MIN_RELIES_ON_FP_RESULT                         164
#define SOLLYA_MSG_MIN_RELIES_ON_FP_RESULT_FAITHFUL_BUT_UNDECIDED  165
#define SOLLYA_MSG_MAX_RELIES_ON_FP_RESULT_THAT_IS_NOT_FAITHFUL    166
#define SOLLYA_MSG_MAX_RELIES_ON_FP_RESULT                         167
#define SOLLYA_MSG_MAX_RELIES_ON_FP_RESULT_FAITHFUL_BUT_UNDECIDED  168
#define SOLLYA_MSG_ERROR_WHILE_EXECUTING_A_PROCEDURE               169
#define SOLLYA_MSG_DIFFERENTIATING_FOR_EVAL_AS_START_PREC_LOW      170
#define SOLLYA_MSG_FREE_VAR_INTERPRETED_AS_IDENTITY_FUNCTION       171
#define SOLLYA_MSG_EXTERNAL_PROCEDURE_SIGNALED_FAILURE             172
#define SOLLYA_MSG_PRECISION_OF_NUMS_MUST_BE_AT_LEAST_TWELVE_BITS  173
#define SOLLYA_MSG_ROUNDING_OCCURRED_WHILE_READING_A_CONSTANT      174
#define SOLLYA_MSG_RANGE_BOUNDS_IN_INVERSE_ORDER                   175
#define SOLLYA_MSG_LITERAL_STRUCTURE_CONTAINS_ENTRY_TWICE          176
#define SOLLYA_MSG_NOT_LEAST_POSSIBLE_INCLUSION_INTERVAL           177
#define SOLLYA_MSG_ONLY_ONE_ENDPOINT_OF_RANGE_IS_NAN               178
#define SOLLYA_MSG_TREE_IS_CONSTANT_BUT_CANNOT_DO_FAITHFUL_EVAL    179
#define SOLLYA_MSG_TIMED_COMMAND_HAS_QUIT_THE_TOOL                 180
#define SOLLYA_MSG_PARAM_OF_PROCEDURE_DOES_NOT_EXIST               181
#define SOLLYA_MSG_AUTODIFF_DEGREE_MUST_NOT_BE_NEGATIVE            182
#define SOLLYA_MSG_EXPR_IS_NO_FRACTION                             183
#define SOLLYA_MSG_NO_ROUNDING_HAS_HAPPENED                        184
#define SOLLYA_MSG_ROUND_UP_HAS_HAPPENED                           185
#define SOLLYA_MSG_ROUND_DOWN_HAS_HAPPENED                         186
#define SOLLYA_MSG_XML_FILE_CANNOT_BE_READ                         187
#define SOLLYA_MSG_FILE_COULD_NOT_BE_OPENED_FOR_WRITING_IGNORING   188
#define SOLLYA_MSG_SUPNORM_DID_NOT_WORK_OUT_WELL                   189
#define SOLLYA_MSG_GUESSDEGREE_FIFTH_ARGUMENT_MUST_BE_NUMBER       190
#define SOLLYA_MSG_ROUND_TO_NEAREST_IMPOSSIBLE_WITH_BOUNDING       191
#define SOLLYA_MSG_OUT_OF_CURRENT_EXPONENT_RANGE                   192
#define SOLLYA_MSG_INADVERTED_ROUNDING_WHILE_DISPLAYING            193
#define SOLLYA_MSG_EXPRESSION_HAS_BEEN_SIMPLIFIED                  194
#define SOLLYA_MSG_EXPRESSION_HAS_BEEN_SIMPLIFIED_TO_ANOTHER_ONE   195
#define SOLLYA_MSG_FORMALLY_DIFFERENTIATING_AN_EXPRESSION          196
#define SOLLYA_MSG_FORMALLY_DIFFERENTIATING_A_PARTICULAR_EXPR      197
#define SOLLYA_MSG_EXPR_TOO_BIG_FOR_SIMPLIFICATION_BEFORE_DIFF     198
#define SOLLYA_MSG_DEGREE_OF_POLYNOMIAL_DOESNT_HOLD_ON_MACHINE_INT 199
#define SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_DOESNT_HOLD_ON_MACHINE_INT  200
#define SOLLYA_MSG_ROUNDING_UPON_BINOMIAL_COEFFICIENT_COMPUTATION  201
#define SOLLYA_MSG_ROUNDING_UPON_POW_EXPONENT_COMPUTATION          202
#define SOLLYA_MSG_RECURSION_ON_POLY_COEFFICIENTS_EXTRACTION       203
#define SOLLYA_MSG_POLY_COEFF_EXTRACTION_SPECIAL_ALGO_FOR_HORNER   204
#define SOLLYA_MSG_POLY_COEFF_EXTRACT_SPECIAL_ALGO_FOR_CANONICAL   205
#define SOLLYA_MSG_TRIED_TO_EXTRACT_COEFFS_OF_STH_NOT_POLYNOMIAL   206
#define SOLLYA_MSG_EXPR_NOT_HORNERIZED_AS_ALREADY_HORNERIZED       207
#define SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_HORNER_POLY    208
#define SOLLYA_MSG_DIFFERENTIATION_USES_SPECIAL_ALGO_FOR_HORNER    209
#define SOLLYA_MSG_DIFFERENTIATION_USES_SPECIAL_ALGO_FOR_CANONICAL 210
#define SOLLYA_MSG_ROUNDING_UPON_DIFFERENTIATION_OF_POLYNOMIAL     211
#define SOLLYA_MSG_EXPR_NOT_CANONICALIZED_AS_ALREADY_CANONICAL     212
#define SOLLYA_MSG_NO_COMMAND_PROVIDED                             213
#define SOLLYA_MSG_ERROR_WHILE_CREATING_A_PIPE                     214
#define SOLLYA_MSG_ERROR_WHILE_FORKING                             215
#define SOLLYA_MSG_UNABLE_TO_WRITE_TO_BASH                         216
#define SOLLYA_MSG_THE_EXIT_CODE_OF_CHILD_PROCESS_IS               216
#define SOLLYA_MSG_SAMPLING_PREC_MUST_BE_LOWER_THAN_CURR_PREC      217
#define SOLLYA_MSG_EXTERNALPLOT_COULD_NOT_OPEN_A_LIBRARY           218
#define SOLLYA_MSG_EXTERNALPLOT_DID_NOT_FIND_FUNCTION_F            219
#define SOLLYA_MSG_COULD_NOT_OPEN_PLOT_FILE                        220
#define SOLLYA_MSG_A_FUNCTION_COULD_NOT_BE_PLOTTED_AT_A_POINT      221
#define SOLLYA_MSG_TOOL_DIES_ON_ERROR_AS_PER_DIE_ON_ERROR_MODE     222
#define SOLLYA_MSG_ERROR_ON_INITIAL_SETUP                          223
#define SOLLYA_MSG_FRAME_STACK_HAS_BEEN_CORRUPTED                  224
#define SOLLYA_MSG_TIMING_STACK_HAS_BEEN_CORRUPTED                 225
#define SOLLYA_MSG_LAST_COMMAND_INTERRUPTED                        226
#define SOLLYA_MSG_RELEASING_FRAME_STACK                           227
#define SOLLYA_MSG_COEFF_NOT_TWICE_GREATER_THAN_SUBPOLY            228
#define SOLLYA_MSG_PREC_OF_HORNER_STEP_GREATER_THAN_FOR_PREV_ONE   229
#define SOLLYA_MSG_NO_AUTO_ROUND_FOR_COEFF_W_PREC_HIGHER_THAN_TD   230
#define SOLLYA_MSG_INFERED_COEFF_PREC_HIGHER_THAN_REQUIRED         231
#define SOLLYA_MSG_COEFF_HAS_BEEN_ROUNDED_TO_TRIPLE_DOUBLE         232
#define SOLLYA_MSG_COEFF_HAS_BEEN_ROUNDED_TO_DOUBLE_DOUBLE         233
#define SOLLYA_MSG_COEFF_HAS_BEEN_ROUNDED_TO_DOUBLE                234
#define SOLLYA_MSG_ERROR_ON_HANDLING_A_COEFFICIENT                 235
#define SOLLYA_MSG_IMPLEMENTED_POLY_DIFFERS_FROM_ORIGINAL_ONE      236
#define SOLLYA_MSG_COEFF_DOES_NOT_EVEN_HOLD_ON_TRIPLE_DOUBLE       237
#define SOLLYA_MSG_ROUNDING_ON_INTERNAL_HANDLING_OF_A_COEFFICIENT  238
#define SOLLYA_MSG_A_COEFF_COULD_NOT_BE_STORED_IN_ANY_KNOWN_FORMAT 239
#define SOLLYA_MSG_IMPLEMENTPOLY_FREE_VAR_HAS_UNKNOWN_FORMAT       240
#define SOLLYA_MSG_ERROR_IN_PRECISION_MANAGEMENT                   241
#define SOLLYA_MSG_CURRENT_PREC_INSUFFICIENT_FOR_TD_CODE           242
#define SOLLYA_MSG_TARGET_ACCURACY_GREATER_OR_EQUAL_THAN_ONE       243
#define SOLLYA_MSG_TARGET_ACCURACY_LESS_THAN_140_BITS              244
#define SOLLYA_MSG_INFERED_OUTPUT_PREC_LESS_THAN_VARIABLE_PREC     245
#define SOLLYA_MSG_COEFF_DOES_NOT_HOLD_ON_TD_USING_FAITHFUL_EVAL   246
#define SOLLYA_MSG_ERROR_ON_DETERMINING_THE_REQUIRED_PRECISIONS    247
#define SOLLYA_MSG_ERROR_ON_DETERMINING_THE_REQUIRED_POWERS        248
#define SOLLYA_MSG_THE_POLY_THAT_GETS_IMPLEMENTED_IS               249
#define SOLLYA_MSG_ERROR_ON_CODE_GENERATION_FOR_COEFFICIENTS       250
#define SOLLYA_MSG_COULD_NOT_WRITE_TO_THE_IMPLEMENTATION_FILE      251
#define SOLLYA_MSG_ERROR_ON_CODE_GENERATION_FOR_POWERS_OF_FREE_VAR 252
#define SOLLYA_MSG_ERROR_ON_CODE_GENERATION_FOR_HORNER_SCHEME      253
#define SOLLYA_MSG_EXPR_SEEMS_TO_BE_ZERO_INCREASE_PREC             254
#define SOLLYA_MSG_A_BASE_FUNC_IS_NOT_SUPPORTED_BY_IMPLEMENTCONST  255
#define SOLLYA_MSG_XML_PARSER_CHANGE                               256
#define SOLLYA_MSG_XML_PARSER_INDEX_CHANGE                         257
#define SOLLYA_MSG_XML_PARSER_STATE_INFORMATION                    258
#define SOLLYA_MSG_ROUNDING_ON_READING_CONSTANT_IN_XML_FILE        259
#define SOLLYA_MSG_XML_PARSER_FAILURE                              260
#define SOLLYA_MSG_XML_SYNCHRONIZATION_LOST_TRYING_TO_RESYNCH      261
#define SOLLYA_MSG_XML_PARSER_UNABLE_TO_OPEN_A_CERTAIN_FILE        262
#define SOLLYA_MSG_XML_PARSE_FUNCTIONALITY_NOT_COMPILED_IN         263
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_INTEGER               264
#define SOLLYA_MSG_ROUNDING_ON_CONVERTING_A_MACHINE_INTEGER        265
#define SOLLYA_MSG_HARMLESS_ERROR_OCCURRED_COMMAND_NOT_EXECUTED    266
#define SOLLYA_MSG_INPUT_AND_OUTPUT_PRECISION_MUST_BE_GREATER_TEN  267
#define SOLLYA_MSG_INTERNAL_PREC_LESS_THAN_IN_AND_OUT_PREC         268
#define SOLLYA_MSG_EPS_SPECIFIED_GREATER_THAN_HALFULP_OF_OUT_PREC  269
#define SOLLYA_MSG_GIVEN_EPS_MUST_BE_POSITIVE_TAKING_ABS           270
#define SOLLYA_MSG_CERTAIN_AMOUNT_OF_CASES_HANDLED                 271
#define SOLLYA_MSG_FUNC_EVALUATED_TO_ZERO_TAKING_ABSOLUTE_ERROR    272
#define SOLLYA_MSG_CANNOT_PERFORM_MORE_THAN_63_STEPS               273
#define SOLLYA_MSG_SEARCH_PREC_HIGHER_THAN_TOOL_PREC               274
#define SOLLYA_MSG_NUMBERS_OF_FUNCS_AND_FORMATS_DIFFER             275
#define SOLLYA_MSG_START_POINT_TOO_PRECISE_FOR_GIVEN_INPUT_PREC    276
#define SOLLYA_MSG_DEGREE_OF_TAYLORFORM_MUST_BE_AT_LEAST_ZERO      271
#define SOLLYA_MSG_ERROR_IN_TAYLORFORM_MULTIPLYING_INCOMPAT_MODELS 272
#define SOLLYA_MSG_ERROR_IN_TAYLORFORM_UNKNOWN_FUNC_FOR_ZUMKELLER  273
#define SOLLYA_MSG_ERROR_IN_TAYLORFORM_TRYING_TO_INCREASE_DEGREE   274
#define SOLLYA_MSG_DEVELOPMENT_POINT_NOT_CONSTANT                  275
#define SOLLYA_MSG_ROUNDING_ON_COMPUTATION_OF_TAYLOR_COEFFICIENT   276
#define SOLLYA_MSG_ROUNDING_ON_COMPUTATION_OF_TAYLOR_POWER         277
#define SOLLYA_MSG_A_CHARACTER_COULD_NOT_BE_RECOGNIZED             278
#define SOLLYA_MSG_A_FILE_COULD_NOT_BE_OPENED_FOR_READING          279
#define SOLLYA_MSG_SPACES_REMOVED_FROM_CONSTANT_IN_SCIENTIF_NOTAT  280
#define SOLLYA_MSG_SYNTAX_ERROR_ENCOUNTERED_WHILE_PARSING          281
#define SOLLYA_MSG_SUPNORM_NO_TAYLOR                               282
#define SOLLYA_MSG_SUPNORM_NOT_ENOUGH_WORKING_PRECISION            283
#define SOLLYA_MSG_SUPNORM_SINGULARITY_NOT_REMOVED                 284
#define SOLLYA_MSG_SUPNORM_COULD_NOT_SHOW_POSITIVITY               285
#define SOLLYA_MSG_SUPNORM_SINGULARITY_NOT_DETECTED                286
#define SOLLYA_MSG_SUPNORM_ANOTHER_SINGULARITY_IN_DOM              287
#define SOLLYA_MSG_SUPNORM_CANNOT_COMPUTE_LOWER_BOUND              288
#define SOLLYA_MSG_SUPNORM_CANNOT_COMPUTE_ABSOLUTE_INF             289
#define SOLLYA_MSG_SUPNORM_CANNOT_DETERMINE_SIGN_OF_T              290
#define SOLLYA_MSG_SUPNORM_CANNOT_DETERMINE_ORDER_OF_SINGU         291
#define SOLLYA_MSG_SUPNORM_GENERIC_ERROR                           292
#define SOLLYA_MSG_SUPNORM_ACCURACY_TOO_HIGH                       292
#define SOLLYA_MSG_SUPNORM_COULD_NOT_FAITHFULLY_EVAL_ERROR_FUNC    293
#define SOLLYA_MSG_DOMAIN_IS_NO_CLOSED_INTERVAL_ON_THE_REALS       294
#define SOLLYA_MSG_DOMAIN_IS_EMPTY                                 295
#define SOLLYA_MSG_DOMAIN_IS_REDUCED_TO_A_POINT_WILL_SIMPLY_EVAL   296
#define SOLLYA_MSG_SUPNORM_COULD_NOT_EVALUATE_ERROR_FUNC           297
#define SOLLYA_MSG_ACCUARCY_INDICATION_IS_NOT_A_REAL_NUMBER        298
#define SOLLYA_MSG_ACCUARCY_INDICATION_IS_ZERO                     299
#define SOLLYA_MSG_POLYNOMIAL_HAS_NON_DYADIC_COEFFICIENTS          300
#define SOLLYA_MSG_SUPNORM_SAFE_ENCLOSURE_COULD_NOT_BE_COMPUTED    301
#define SOLLYA_MSG_STURM_INTERVAL_A_CERTAIN_PREC_HAS_BEEN_CHOSEN   302
#define SOLLYA_MSG_STURM_COEFF_EVALUATED_TO_RATIONAL_NUMBER        303
#define SOLLYA_MSG_STURM_COEFF_NOT_CONSTANT_NOR_RATIONAL_ROUNDING  304
#define SOLLYA_MSG_STURM_COEFF_ROUNDED_TO_ZERO                     305
#define SOLLYA_MSG_STURM_USING_SLOWER_ALGORITHM_ON_RATIONALS       306
#define SOLLYA_MSG_STURM_POLY_IS_ZERO_POLY                         307
#define SOLLYA_MSG_CONSTANT_EXPR_CANNOT_BE_EVALUATED_AT_ALL        308
#define SOLLYA_MSG_PLOT_OVERFLOW_OCCURRED_ON_CONVERSION_TO_DOUBLE  309
#define SOLLYA_MSG_PLOT_FUNC_PROVEN_LESS_THAN_2_TO_MINUS_PREC      310
#define SOLLYA_MSG_PLOT_FUNC_UNDEFINED_OR_UNSTABLE_AT_POINT        311
#define SOLLYA_MSG_PLOT_NOT_FAITHFULLY_EVALUATED_AT_SOME_POINT     312
#define SOLLYA_MSG_DOMAIN_IS_REDUCED_TO_A_POINT_TRIVIAL_RESULT     313
#define SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_ALREADY_BOUND             314
#define SOLLYA_MSG_COULD_NOT_OPEN_LIBRARY_WITH_EXTERN_FUNC_OR_PROC 315
#define SOLLYA_MSG_EXTERNAL_FUNC_OR_PROC_NOT_FOUND_IN_LIBRARY      316
#define SOLLYA_MSG_COULD_NOT_CLOSE_LIBRARY                         317
#define SOLLYA_MSG_ENTERING_NEWTONS_ALGORITHM                      318
#define SOLLYA_MSG_NEWTON_ZERO_IS_EXACT_ZERO                       319
#define SOLLYA_MSG_NEWTON_AN_EXACT_ZERO_HAS_BEEN_FOUND             320
#define SOLLYA_MSG_ERROR_IN_TAYLORFORM_COPYING_INCOMPAT_MODELS     321
#define SOLLYA_MSG_PLOT_COULD_NOT_OPEN_FILE                        322
#define SOLLYA_MSG_NEWTON_FUNC_APPEARS_TO_HAVE_MORE_THAN_ONE_ZERO  323
#define SOLLYA_MSG_NEWTON_ZERO_TOO_CLOSE_TO_ZERO_TO_BE_ACCURATE    324
#define SOLLYA_MSG_NEWTON_ZERO_SEEMS_TO_BE_ZERO_NO_PROOF           325
#define SOLLYA_MSG_NEWTON_ALGORITHM_FAILS_DUE_TO_NUMERICAL_ISSUES  326
#define SOLLYA_MSG_NEWTON_ALGORITHM_FAILS_TO_LOCATE_ZERO           327
#define SOLLYA_MSG_NEWTON_PERFORMING_BISECTION_STEP                328
#define SOLLYA_MSG_NEWTON_PERFORMING_TRISECTION_STEP               329
#define SOLLYA_MSG_NEWTON_FINISHED_AFTER_NUMBER_OF_STEPS           330
#define SOLLYA_MSG_REMEZ_EXCHANGE_TAKE_A_CERTAIN_MINIMUM           331
#define SOLLYA_MSG_REMEZ_EXCHANGE_TAKE_A_CERTAIN_MAXIMUM           332
#define SOLLYA_MSG_REMEZ_FUNCTION_OSCILLATES_TOO_MUCH              333
#define SOLLYA_MSG_REMEZ_PERFORMING_AN_EXCHANGE_STEP               334
#define SOLLYA_MSG_REMEZ_COMPUTED_INFNORM_IS_A_CERTAIN_VALUE       335
#define SOLLYA_MSG_REMEZ_FAILED_TO_FIND_PSEUDOALTERNATING_POINTS   336
#define SOLLYA_MSG_REMEZ_CONSTRUCTING_THE_ERROR_TREE               337
#define SOLLYA_MSG_REMEZ_CONSTRUCTING_THE_ERROR_PRIME_TREE         338
#define SOLLYA_MSG_REMEZ_CONSTRUCTING_THE_ERROR_SECOND_TREE        339
#define SOLLYA_MSG_REMEZ_COMPUTING_THE_YI                          340
#define SOLLYA_MSG_REMEZ_THE_COMPUTED_YI_ARE_CERTAIN_VALUES        341
#define SOLLYA_MSG_REMEZ_ALGORITHM_IS_IN_A_CERTAIN_CASE            342
#define SOLLYA_MSG_REMEZ_THE_COMPUTED_SIGNS_ARE_CERTAIN_VALUES     343
#define SOLLYA_MSG_REMEZ_SIGNS_COULD_NOT_BE_EVALUATED              344
#define SOLLYA_MSG_REMEZ_MAIN_HEURISTIC_FAILED_USING_SLOWER_ALGO   345
#define SOLLYA_MSG_REMEZ_SLOWER_ALGORITHM_USED_FOR_A_STEP          346
#define SOLLYA_MSG_REMEZ_THE_NEW_POINTS_ARE_CERTAIN_VALUES         347
#define SOLLYA_MSG_REMEZ_THE_CURRENT_NORM_TAKES_A_CERTAIN_VALUE    348
#define SOLLYA_MSG_ENTERING_REMEZ_FUNCTION                         349
#define SOLLYA_MSG_REMEZ_COMPUTING_MONOMIALS                       350
#define SOLLYA_MSG_REMEZ_COMPUTING_INITIAL_POINT_SET               351
#define SOLLYA_MSG_REMEZ_THE_COMPUTED_POINT_SET_IS_CERTAIN_VALUES  352
#define SOLLYA_MSG_REMEZ_COMPUTING_THE_MATRIX                      353
#define SOLLYA_MSG_REMEZ_COMPUTAT_OF_MATRIX_ENTRY_USES_SLOWER_ALGO 354
#define SOLLYA_MSG_REMEZ_DEGENERATED_SYSTEM_IN_NON_HAAR_CONTEXT    355
#define SOLLYA_MSG_REMEZ_SIGNS_FOR_PSEUDO_ALTERN_ARE_CERTAIN_VALS  356
#define SOLLYA_MSG_REMEZ_THE_COMPUTED_MATRIX_HAS_A_CERTAIN_VALUE   357
#define SOLLYA_MSG_REMEZ_SOLVING_THE_SYSTEM                        358
#define SOLLYA_MSG_REMEZ_THE_COMPUTED_POLY_HAS_A_CERTAIN_VALUE     359
#define SOLLYA_MSG_REMEZ_CURRENT_EPSILON_HAS_A_CERTAIN_VALUE       360
#define SOLLYA_MSG_REMEZ_DIFFERENTIATING_THE_COMPUTED_POLYNOMIAL   361
#define SOLLYA_MSG_REMEZ_SEARCHING_FOR_EXTREMA_OF_ERROR_FUNCTION   362
#define SOLLYA_MSG_REMEZ_THE_BEST_POLY_GIVES_A_CERTAIN_ERROR       363
#define SOLLYA_MSG_REMEZ_CURRENT_QUALITY_HAS_A_CERTAIN_VALUE       364
#define SOLLYA_MSG_REMEZ_FINISHES_AS_TARGET_ERROR_IS_NOT_REACHABLE 365
#define SOLLYA_MSG_REMEZ_FINISHES_AS_TARGET_ERROR_HAS_BEEN_REACHED 366
#define SOLLYA_MSG_REMEZ_FINISHES_AS_QUALITY_HAS_BEEN_REACHED      367
#define SOLLYA_MSG_REMEZ_FAILS_AND_LOOPS_AGAIN                     368
#define SOLLYA_MSG_REMEZ_DOES_NOT_CONVERGE                         369
#define SOLLYA_MSG_REMEZ_MAY_HAPPEN_NOT_TO_CONVRG_AS_DOM_IS_POINT  370
#define SOLLYA_MSG_GUESSDEGREE_TRYING_A_CERTAIN_DEGREE             371
#define SOLLYA_MSG_GUESSDEGREE_NONE_OF_LESSER_DEGS_SATISFIES_ERROR 372
#define SOLLYA_MSG_GUESSDEGREE_TRYING_A_CERTAIN_DEG_WITHIN_BOUNDS  373
#define SOLLYA_MSG_GUESSDEG_NONE_OF_LESS_DEGS_SEEMS_TO_SATISFY_ERR 374
#define SOLLYA_MSG_FPMINIMAX_SINGULAR_MATRIX                       375
#define SOLLYA_MSG_FPMINIMAX_A_CERTAIN_COEFF_IS_EXACT_ZERO         376
#define SOLLYA_MSG_FPMINIMAX_MINIMAX_DOES_NOT_GIVE_ENOUGH_POINTS   377
#define SOLLYA_MSG_FPMINIMAX_THE_POINTS_ARE_CERTAIN_VALUES         378
#define SOLLYA_MSG_FPMINIMAX_FAILED_TO_RECOVER_COEFFS_FROM_POLY    379
#define SOLLYA_MSG_FPMINIMAX_THE_EXPONENTS_ARE_CERTAIN_VALUES      380
#define SOLLYA_MSG_FPMINIMAX_DID_NOT_CONVERGE                      381
#define SOLLYA_MSG_FPMINIMAX_NOT_ENOUGH_POINTS                     382
#define SOLLYA_MSG_FPMINIMAX_NOT_ENOUGH_FORMATS                    383
#define SOLLYA_MSG_FPMINIMAX_COMP_OF_MATRIX_ENTRY_USES_SLOWER_ALGO 384
#define SOLLYA_MSG_DIFFERENTIATING_FOR_DECORRELATION               385
#define SOLLYA_MSG_DECORRELATION_INTERVAL_ADDITION_OR_SUBTRACTION  386
#define SOLLYA_MSG_DIFFERENTIATING_FOR_HOPITALS_RULE               387
#define SOLLYA_MSG_USING_HOPITALS_RULE_ON_POINT_DIVISION           388
#define SOLLYA_MSG_SIMPLIFYING_INTERVAL_DIV_WITH_ZERO_POINT_NUMERA 389
#define SOLLYA_MSG_USING_HOPITALS_RULE_IN_GENERAL_CASE             390
#define SOLLYA_MSG_RECURSION_ON_USE_OF_HOPITALS_RULE               391
#define SOLLYA_MSG_AVOIDING_TAYLOR_EVALUATION_ON_POINT_INTERVAL    392
#define SOLLYA_MSG_NO_TAYLOR_EVALUATION_AS_NO_DERIVATIVE_GIVEN     393
#define SOLLYA_MSG_USING_TAYLOR_EVALUATION                         394
#define SOLLYA_MSG_NO_TAYLOR_EVALUATION_AS_NO_DERIVATIVE_GETS_HUGE 395
#define SOLLYA_MSG_DERIVATIVE_DOES_NOT_CHANGE_SIGN_ON_TAYLOR_EVAL  396
#define SOLLYA_MSG_NAN_OR_INF_ON_DERIVATIVE                        397
#define SOLLYA_MSG_INVOKING_RECURSIVE_INTERVAL_ZERO_SEARCH         398
#define SOLLYA_MSG_RECURSIVE_INTERVAL_ZERO_SEARCH_HAS_FINISHED     399
#define SOLLYA_MSG_CERTAIN_NUM_OF_INTVALS_ENCLOSING_ZEROS_OF_DERIV 400
#define SOLLYA_MSG_EXPRESSION_IS_CONSTANT                          401
#define SOLLYA_MSG_EVALUATION_AT_POINT_GIVES_NAN_EXCLUDING_POINT   402
#define SOLLYA_MSG_THE_CURRENT_MAXIMUM_IS_A_CERTAIN_VALUE          403
#define SOLLYA_MSG_EVALUATION_OF_DERIVATIVE_GIVES_NAN_NO_NEWTON    404
#define SOLLYA_MSG_NO_PROOF_WILL_BE_GENERATED                      405
#define SOLLYA_MSG_INFNORM_RESULT_IS_TRIVIAL                       406
#define SOLLYA_MSG_DERIVATIVE_IS_QUOTIENT                          407
#define SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_SINGULARITY            408
#define SOLLYA_MSG_DERIVATIVE_SEEMS_TO_HAVE_EXTENSIBLE_SINGULARITY 409
#define SOLLYA_MSG_DERIVATIVE_SEEMS_NOT_TO_HAVE_ANY_POLE           410
#define SOLLYA_MSG_INVOKING_INFNORM_SUBFUNCTION                    411
#define SOLLYA_MSG_INFNORM_SUBFUNCTION_HAS_FINISHED                412
#define SOLLYA_MSG_STARTING_TO_WRITE_THE_PROOF                     413
#define SOLLYA_MSG_THE_PROOF_HAS_BEEN_WRITTEN                      414
#define SOLLYA_MSG_THE_EXPRESSION_IS_NOT_CONSTANT                  415
#define SOLLYA_MSG_COULD_NOT_CHECK_INFNORM_ON_A_CERTAIN_INTERVAL   416
#define SOLLYA_MSG_REMOVING_A_POSSIBLE_ZERO_AT_SOME_POINT          417
#define SOLLYA_MSG_ZERO_FILTER_HAS_REMOVED_AT_LEAST_ONE_ZERO       418
#define SOLLYA_MSG_FAITHFUL_EVALUATION_RETURNS_NAN                 419
#define SOLLYA_MSG_INTERMEDIATE_PRECISION_HAS_BEEN_INCREASED       420
#define SOLLYA_MSG_TAYLOR_RECURSION_TEMPORARILY_SET_TO_A_VALUE     421
#define SOLLYA_MSG_ABS_DIAM_AND_PREC_SET_TO_CERTAIN_VALUES         422
#define SOLLYA_MSG_IDENTIFIER_IS_LIBRARY_FUNC_CANNOT_BE_EXTERNAL   423
#define SOLLYA_MSG_HANDLED_SIGSEGV                                 424
#define SOLLYA_MSG_HANDLED_SIGBUS                                  425
#define SOLLYA_MSG_HANDLED_SIGFPE                                  426
#define SOLLYA_MSG_HANDLED_SIGPIPE                                 427
#define SOLLYA_MSG_CANNOT_SUPPRESS_OR_UNSUPPRESS_A_MESSAGE         428
#define SOLLYA_MSG_EXPR_DOES_NOT_EVALUATE_TO_INT_OR_LIST_OF_INT    429
#define SOLLYA_MSG_SUPPRESSION_NUMBER_OMITTED                      430
#define SOLLYA_MSG_ROUNDING_ON_CONSTANT_RETRIEVAL                  431
#define SOLLYA_MSG_MIN_RELIES_ON_FP_RESULT_FAITHFUL_BUT_NOT_REAL   432
#define SOLLYA_MSG_MAX_RELIES_ON_FP_RESULT_FAITHFUL_BUT_NOT_REAL   433
#define SOLLYA_MSG_EXPRESSION_EVALUATES_TO_INFINITY                434
#define SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NOT_CONSTANT             435
#define SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NOT_INTEGER              436
#define SOLLYA_MSG_DEG_OF_MAX_POLY_DIV_IS_NEGATIVE                 437
#define SOLLYA_MSG_RATIONALAPPROX_SECOND_ARG_MUST_BE_GREATER_THAN_ONE 438
#define SOLLYA_MSG_ROUND_PREC_MUST_BE_AT_LEAST_TWO_BITS            439
#define SOLLYA_MSG_NAN_CONVERTED_TO_NUMBER_ON_CONSTANT_RETRIEVAL   440
#define SOLLYA_MSG_CHEBYSHEVFORM_DEGREE_MUST_NOT_BE_NEGATIVE       441
#define SOLLYA_MSG_CHEBYSHEVFORM_DOMAIN_MUST_NOT_BE_POINT_INTERVAL 442
#define SOLLYA_MSG_CHEBYSHEVFORM_ERROR_IN_COMPUTATION              443
#define SOLLYA_MSG_ERROR_IN_CHEBYSHEVFORM_COPYING_INCOMPAT_MODELS  444
#define SOLLYA_MSG_ERROR_IN_CHEBYSHEVFORM_UNKNOWN_FUNC_FOR_ZUMKELLER  445
#define SOLLYA_MSG_ERROR_IN_CHEBYSHEVFORM_NOT_A_POLYNOMIAL         446
#define SOLLYA_MSG_SPECIAL_ALGORITHM_USED_FOR_COEFF                447
#define SOLLYA_MSG_NO_CORRECT_ROUNDING_FOR_ROUND_OPERATOR          448
#define SOLLYA_MSG_ROUNDING_OF_BOUNDARY_INSTEAD_OF_CORRECT_ROUNDING 449
#define SOLLYA_MSG_NO_CORRECT_TERNARY_VALUE_FOR_ROUND_BUT_CORRECT_ROUNDING 450
#define SOLLYA_MSG_NO_CORRECT_TERNARY_VALUE_FOR_ROUND              451
#define SOLLYA_MSG_LIBRARY_CLOSER_ERROR                            452
#define SOLLYA_MSG_INF_CONVERTED_TO_NUMBER_ON_CONSTANT_RETRIEVAL   453
#define SOLLYA_MSG_ANNOTATION_COULD_NOT_BE_SET_UP                  454
#define SOLLYA_MSG_GENERIC_SOLLYA_LIBRARY_MSG                      455
#define SOLLYA_MSG_GUESSDEGREE_POSSIBLE_SINGULAR_WEIGHT            456
#define SOLLYA_MSG_SAFE_ROUNDING_FOR_EXPR_THAT_SHOULD_BE_CONST     457
#define SOLLYA_MSG_DEGREE_OF_POLYNOMIAL_LARGER_THAN_MULTIPRECISION_INT 458
#define SOLLYA_MSG_ANNOTATION_INCOHERENT                           459

#endif /* ifdef SOLLYA_MESSAGES_H*/
