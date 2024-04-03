#include <assert.h>
#include <inttypes.h>
#include <linux/limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define MIN_STEP 1e-9
// TODO:
// debug
//	calculate error
//	test lemma for some test (separate program)
//

typedef enum error_code {
  NO_ERROR = 0,
  READ_FILE_ERROR = 1,
  DIV_BY_ZERO = 2,
  INVALID_ARGUMENT = 3,
  NO_DATA = 4,
  INVALID_DATA = 5,
  STEP_IS_TO_SMALL = 6,
  MALLOC_ERROR = 7,
  WRITE_FILE_ERROR = 8
} error_code;

error_code read_input_file(const char *const file_path, double_t **input_data,
                           uint64_t *size);

double_t predict_left_border_point(double_t prev_point, double_t next_point);
double_t predict_left_border_result(double_t prev_result, double_t next_result);
double_t predict_right_border_point(double_t prev_point, double_t next_point);
double_t predict_right_border_result(double_t prev_result,
                                     double_t next_result);
error_code calculate_differences_between_points(
    const double_t *const input_data, uint64_t input_data_size,
    double_t **differences, uint64_t *differences_size);

error_code
calculate_divided_differences(const double_t *const input_data,
                              uint64_t input_data_size,
                              const double_t *const differences_between_points,
                              uint64_t differences_between_points_size,
                              double_t **div_diffs, uint64_t *size);
error_code calculate_omegas(const double_t *const div_diffs,
                            uint64_t div_diffs_size, double_t **omegas,
                            uint64_t *size);

error_code calculate_di(const double_t *const input_data,
                        uint64_t input_data_size,
                        const double_t *const div_diffs,
                        uint64_t div_diffs_size, const double_t *const omegas,
                        const double_t *const dbp, double_t **diths,
                        uint64_t *diths_size);

error_code calculate_a_1(const double_t *const input_data,
                         const uint64_t input_data_size, double_t **A_1,
                         uint64_t *size);
double_t const *calculate_a_2(const double_t *const dith);
// dbp means differences between points
error_code calculate_a_3(const double_t *const divided_differences,
                         const double_t *const dbp, const double_t *const dith,
                         uint64_t dith_size, double_t **A_3,
                         uint64_t *A_3_size);
error_code calculate_a_4(const double_t *const divided_differences,
                         const double_t *const dbp, const double_t *const dith,
                         uint64_t dith_size, double_t **A_4,
                         uint64_t *A_4_size);
uint64_t search(double_t x, const double_t *const input_data,
                uint64_t input_data_size);
double_t approximate(double_t **A_iths, const double_t *const input_data,
                     uint64_t input_data_size, double_t x);

double_t
calculate_maximal_step(const double_t *const differences_between_points);

double_t max_of_function_on_a_segment(const double_t *const input_data);

double_t test_approximation_error(double_t x, double_t fx, double_t h);

error_code
print_approximated_data_to_file(const char *const file_path,
                                const double_t *const approximated_data,
                                uint64_t approximated_data_size);
void usage();

extern char *optarg;

int main(int argc, char **argv) {
  int opt = 0;
  double_t *A_iths[4] = {NULL};
  uint64_t A_1_size = 0;
  uint64_t A_2_size = 0;
  uint64_t A_3_size = 0;
  uint64_t A_4_size = 0;
  double_t *input_data = NULL;
  uint64_t input_data_size = 0;
  double_t *divided_differences = NULL;
  uint64_t divided_differences_size = 0;
  double_t *omegas = NULL;
  uint64_t omegas_size = 0;
  double_t *differences_between_points = NULL;
  uint64_t differences_between_points_size = 0;
  double_t *diths = NULL;
  uint64_t diths_size = 0;

  double_t *approximated_data = NULL;
  char *input_file_path = NULL;
  double_t step = 1e-3;
  while ((opt = getopt(argc, argv, "f:s:")) != -1) {
    switch (opt) {
    case 'f':
      assert(read_input_file(optarg, &input_data, &input_data_size) ==
             NO_ERROR);
      input_file_path = optarg;
      break;
    case 's':
      step = (double_t)atof(optarg);
      if (fabs(step) < MIN_STEP) {
        usage();
        exit(INVALID_ARGUMENT);
      }
      break;
    default:
      usage();
      exit(INVALID_ARGUMENT);
    }
  }
  if (input_data == NULL) {
    exit(NO_DATA);
  }
  input_data[0] = predict_left_border_point(input_data[2], input_data[4]);
  input_data[1] = predict_left_border_result(input_data[3], input_data[5]);
  input_data[input_data_size - 2] = predict_right_border_point(
      input_data[input_data_size - 6], input_data[input_data_size - 4]);
  input_data[input_data_size - 1] = predict_right_border_result(
      input_data[input_data_size - 5], input_data[input_data_size - 3]);

  assert(calculate_differences_between_points(
             input_data, input_data_size, &differences_between_points,
             &differences_between_points_size) == NO_ERROR);
  assert(calculate_divided_differences(
             input_data, input_data_size, differences_between_points,
             differences_between_points_size, &divided_differences,
             &divided_differences_size) == NO_ERROR);
  assert(calculate_omegas(divided_differences, divided_differences_size,
                          &omegas, &omegas_size) == NO_ERROR);

  assert(calculate_di(input_data, input_data_size, divided_differences,
                      divided_differences_size, omegas,
                      differences_between_points, &diths,
                      &diths_size) == NO_ERROR);
  assert((calculate_a_1(input_data, input_data_size, &A_iths[0], &A_1_size)) ==
         NO_ERROR);
  A_iths[1] = (double_t *)calculate_a_2(diths);
  assert(calculate_a_3(divided_differences, differences_between_points, diths,
                       diths_size, &A_iths[2], &A_2_size) == NO_ERROR);
  assert(calculate_a_4(divided_differences, differences_between_points, diths,
                       diths_size, &A_iths[3], &A_4_size) == NO_ERROR);

  double_t min_x = input_data[2];
  double_t max_x = input_data[input_data_size - 4];
  uint64_t approximated_data_size =
      2 * (uint64_t)ceil((max_x - min_x + 1) / step);
  approximated_data =
      (double_t *)malloc(approximated_data_size * sizeof(double_t));
  if (approximated_data == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0; i < approximated_data_size - 1; i += 2) {
    approximated_data[i] = min_x;
    approximated_data[i + 1] =
        approximate(A_iths, input_data, input_data_size, min_x);
    min_x += step;
  }
  char out_file_path[PATH_MAX];
  snprintf(out_file_path, PATH_MAX, "%s_out", input_file_path);
  assert(print_approximated_data_to_file(out_file_path, approximated_data,
                                         approximated_data_size) == NO_ERROR);

  return 0;
}
// OK
error_code read_input_file(const char *const file_path, double_t **input_data,
                           uint64_t *size) {
  FILE *fp = fopen(file_path, "r");
  if (fp == NULL) {
    return READ_FILE_ERROR;
  }

  char buffer[256]; // maximum number of digits in double is 15 => 32 bytes per
                    // line MAX (double + " " + double + "\n")
  uint64_t number_of_lines = 0;
  while (fgets(buffer, sizeof buffer, fp) && !feof(fp)) {
    if (buffer[0] == '\n') {
      fclose(fp);
      return INVALID_DATA;
    }
    number_of_lines++;
  }
  if (number_of_lines == 0) {
    return NO_DATA;
  }

  double_t *input_data_local =
      (double_t *)malloc(2 * sizeof(double_t) * (number_of_lines + 2));
  if (input_data_local == NULL) {
    return MALLOC_ERROR;
  }
  uint64_t size_local = 2 * (number_of_lines + 2);

  if (fseek(fp, 0, SEEK_SET) != 0) {
    return READ_FILE_ERROR;
  }

  for (uint64_t i = 2; i < size_local - 2; i += 2) {
    fgets(buffer, sizeof buffer, fp);
    if (sscanf(buffer, "%lf %lf", &input_data_local[i],
               &input_data_local[i + 1]) != 2) {
      return INVALID_DATA;
    }
  }

  *input_data = input_data_local;
  *size = size_local;

  fclose(fp);
  return NO_ERROR;
}

inline double_t predict_left_border_point(double_t prev_point,
                                          double_t next_point) {
#ifdef _DEBUG
  double_t res = prev_point - (next_point - prev_point);

  printf("X_0 = %f\n", res);
  return res;
#endif

  return prev_point - (next_point - prev_point);
}
inline double_t predict_left_border_result(double_t prev_result,
                                           double_t next_result) {
#ifdef _DEBUG
  double_t res = prev_result - (next_result - prev_result);

  printf("Y_0 = %f\n", res);
  return res;
#endif

  return prev_result - (next_result - prev_result);
}
inline double_t predict_right_border_point(double_t prev_point,
                                           double_t next_point) {
#ifdef _DEBUG
  double_t res = next_point + (next_point - prev_point);

  printf("X_N+1 = %f\n", res);
  return res;
#endif

  return next_point + (next_point - prev_point);
}
inline double_t predict_right_border_result(double_t prev_result,
                                            double_t next_result) {
#ifdef _DEBUG
  double_t res = next_result + (next_result - prev_result);

  printf("Y_N+1 = %f\n", res);
  return res;
#endif

  return next_result + (next_result - prev_result);
}

error_code calculate_differences_between_points(
    const double_t *const input_data, uint64_t input_data_size,
    double_t **differences, uint64_t *differences_size) {
  uint64_t differences_size_local = input_data_size / 2 - 1;
  double_t *differences_local =
      (double_t *)malloc(sizeof(double_t) * differences_size_local);
  if (differences_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0, j = 0; i < differences_size_local; ++i, j += 2) {
    differences_local[i] = input_data[j + 2] - input_data[j];
    if (fabs(differences_local[i]) < MIN_STEP) {
      return STEP_IS_TO_SMALL;
    }
  }
  *differences = differences_local;
  *differences_size = differences_size_local;
#ifdef _DEBUG
  for (uint64_t i = 0; i < differences_size_local; ++i) {
    printf("x%" PRIu64 " - x%" PRIu64 " = %f\n", i + 1, i,
           differences_local[i]);
  }
  return NO_ERROR;
#endif

  return NO_ERROR;
}

error_code
calculate_divided_differences(const double_t *const input_data,
                              uint64_t input_data_size,
                              const double_t *const differences_between_points,
                              uint64_t differences_between_points_size,
                              double_t **div_diffs, uint64_t *size) {
  uint64_t size_local = input_data_size / 2 - 1;
  double_t *div_diffs_local = (double_t *)malloc(sizeof(double_t) * size_local);
  if (div_diffs_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0, j = 0; i < input_data_size - 3; i += 2, j++) {
    div_diffs_local[j] =
        (input_data[i + 3] - input_data[i + 1]) / differences_between_points[j];
  }
  *div_diffs = div_diffs_local;
  *size = size_local;
#ifdef _DEBUG
  for (uint64_t i = 0; i < size_local; ++i) {
    printf("f(x%" PRIu64 "; x%" PRIu64 ") = %f\n", i, i + 1,
           div_diffs_local[i]);
  }
  return NO_ERROR;
#endif

  return NO_ERROR;
}
error_code calculate_omegas(const double_t *const div_diffs,
                            uint64_t div_diffs_size, double_t **omegas,
                            uint64_t *size) {
  uint64_t size_local = div_diffs_size - 1;
  double_t *omegas_local = (double_t *)malloc(sizeof(double_t) * size_local);
  if (omegas_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0; i < size_local; ++i) {
    omegas_local[i] = fabs(div_diffs[i] - div_diffs[i + 1]);
  }
  *omegas = omegas_local;
  *size = size_local;
#ifdef _DEBUG
  for (uint64_t i = 0; i < size_local; ++i) {
    printf("\u03A9[%" PRIu64 "] = %f\n", i, i + 1, omegas_local[i]);
  }
  return NO_ERROR;
#endif

  return NO_ERROR;
}

error_code calculate_di(const double_t *const input_data,
                        uint64_t input_data_size,
                        const double_t *const div_diffs,
                        uint64_t div_diffs_size, const double_t *const omegas,
                        const double_t *const dbp, double_t **diths,
                        uint64_t *diths_size) {
  uint64_t diths_size_local = input_data_size / 2 - 2;
  double_t *diths_local =
      (double_t *)malloc(sizeof(double_t) * diths_size_local);
  if (diths_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 1; i < diths_size_local - 1; ++i) {
    if (fabs(omegas[i + 1]) < MIN_STEP && fabs(omegas[i - 1]) < MIN_STEP) {
      diths_local[i] = (dbp[i] * div_diffs[i - 1] + dbp[i - 1] * div_diffs[i]) /
                       (input_data[2 * (i + 1)] - input_data[2 * (i - 1)]);
    } else {
      if (fabs(omegas[i + 1] + omegas[i - 1]) >= MIN_STEP)
        diths_local[i] =
            (omegas[i + 1] * div_diffs[i - 1] + omegas[i - 1] * div_diffs[i]) /
            (omegas[i + 1] + omegas[i - 1]);
      else
        return DIV_BY_ZERO;
    }
  }
  diths_local[0] = (3 * div_diffs[1] - diths_local[1]) / 2;
  diths_local[diths_size_local - 1] =
      (3 * div_diffs[div_diffs_size - 2] - diths_local[diths_size_local - 2]) /
      2;
  *diths = diths_local;
  *diths_size = diths_size_local;
#ifdef _DEBUG
  for (uint64_t i = 0; i < diths_size_local; ++i) {
    printf("d[%" PRIu64 "] = %f\n", i, diths_local[i]);
  }
  return NO_ERROR;
#endif
  return NO_ERROR;
}

error_code calculate_a_1(const double_t *const input_data,
                         const uint64_t input_data_size, double_t **A_1,
                         uint64_t *A_1_size) {
  uint64_t A_1_local_size = (input_data_size) / 2;
  double_t *A_1_local = (double_t *)malloc(A_1_local_size * sizeof(double_t));
  if (A_1_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 1, j = 0; i < input_data_size; i += 2, ++j) {
    A_1_local[j] = input_data[i];
  }
  *A_1 = A_1_local;
  *A_1_size = A_1_local_size;
#ifdef _DEBUG
  for (uint64_t i = 0; i < A_1_local_size; ++i) {
    printf("A_1[%" PRIu64 "] = %f\n", i, A_1_local[i]);
  }
  return NO_ERROR;
#endif
  return NO_ERROR;
}
inline const double_t *calculate_a_2(const double_t *const dith) {
  return dith;
}

error_code calculate_a_3(const double_t *const divided_differences,
                         const double_t *const dbp, const double_t *const dith,
                         uint64_t dith_size, double_t **A_3,
                         uint64_t *A_3_size) {
  uint64_t A_3_local_size = (dith_size - 1);
  double_t *A_3_local = (double_t *)malloc(A_3_local_size * sizeof(double_t));
  if (A_3_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0; i < dith_size - 1; ++i) {
    A_3_local[i] = (divided_differences[i] - dith[i]) / dbp[i];
  }
  *A_3 = A_3_local;
  *A_3_size = A_3_local_size;
#ifdef _DEBUG
  for (uint64_t i = 0; i < A_3_local_size; ++i) {
    printf("A_3[%" PRIu64 "] = %f\n", i, A_3_local[i]);
  }
  return NO_ERROR;
#endif
  return NO_ERROR;
}

error_code calculate_a_4(const double_t *const divided_differences,
                         const double_t *const dbp, const double_t *const dith,
                         uint64_t dith_size, double_t **A_4,
                         uint64_t *A_4_size) {
  uint64_t A_4_local_size = dith_size - 1;
  double_t *A_4_local = (double_t *)malloc(sizeof(double_t) * (A_4_local_size));
  if (A_4_local == NULL) {
    return MALLOC_ERROR;
  }
  for (uint64_t i = 0; i < dith_size - 1; ++i) {
    A_4_local[i] = (dith[i] - dith[i + 1] - 2 * divided_differences[i]) /
                   (dbp[i] * dbp[i]);
  }
#ifdef _DEBUG
  for (uint64_t i = 0; i < A_4_local_size; ++i) {
    printf("A_4[%" PRIu64 "] = %f\n", i, A_4_local[i]);
  }
  return NO_ERROR;
#endif
  *A_4 = A_4_local;
  *A_4_size = A_4_local_size;
  return NO_ERROR;
}

uint64_t search(double_t x, const double_t *const input_data,
                uint64_t input_data_size) {
  double_t min_delta = __DBL_MAX__;
  uint64_t idx_left_bdr = 0;
  for (uint64_t i = 0; i < input_data_size - 1; i += 2) {
    if (x > input_data[i] && min_delta > x - input_data[i]) {
      min_delta = x - input_data[i];
      idx_left_bdr = i;
    }
  }
  return idx_left_bdr;
}

double_t approximate(double_t **A_iths, const double_t *const input_data,
                     uint64_t input_data_size, double_t x) {
  uint64_t idx = search(x, input_data, input_data_size);
  double_t dd = (x - input_data[2 * idx]);
  double_t dd_2 = dd * dd;
  return A_iths[0][idx] + A_iths[1][idx] * dd + A_iths[2][idx] * dd_2 +
         A_iths[3][idx] * dd_2 * (x - input_data[2 * (idx + 1)]);
}

error_code
print_approximated_data_to_file(const char *const file_path,
                                const double_t *const approximated_data,
                                uint64_t approximated_data_size) {
  FILE *out = fopen(file_path, "w");
  if (out == NULL) {
    return WRITE_FILE_ERROR;
  }
  for (uint64_t i = 0; i < approximated_data_size - 1; i += 2) {
    fprintf(out, "%f, %f\n", approximated_data[i], approximated_data[i + 1]);
  }
  fclose(out);
  return NO_ERROR;
}
void usage() { return; }
