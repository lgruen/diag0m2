#ifndef GEN_DIAG_0M2_H
#define GEN_DIAG_0M2_H

#include <cassert>

// This class generates permutation-net points as described in
// Leonhard Gruenschloss, Alexander Keller:
// "(t,m,s)-Nets and Maximized Minimum Distance, Part II".
//
// available at http://gruenschloss.org/diag0m2/gendiag0m2.h
//
// usage example:
//
// #include "gendiag0m2.h"
// #include <iostream>
// 
// int main(int argc, char** argv) {
//   // odd case
//   {
//     const GenDiag0m2 gen(7);
//     const double scale = 1.0 / (1ULL << gen.m);
//     for (unsigned int i = 0; i < gen.n; ++i) {
//       unsigned int x, y;
//       gen.get(i, x, y);
//       const double dx = x * scale;
//       const double dy = y * scale;
//       std::cout << dx << " " << dy << std::endl;
//     }
//   }
// 
//   // even case, with modified tiling
//   {
//     const GenDiag0m2 gen(6);
//     const double scale = 1.0 / (1ULL << gen.m);
//     for (unsigned int py = 0; py < 4; ++py) {
//       for (unsigned int px = 0; px < 4; ++px) {
//         for (unsigned int i = 0; i < gen.n; ++i) {
//           unsigned int x, y;
//           gen.getShiftedTiling(px, py, i, x, y);
//           const double dx = x * scale;
//           const double dy = y * scale;
//           std::cout << dx << " " << dy << std::endl;
//         }
//       }
//     }
//   }
// 
//   return 0;
// }
//
class GenDiag0m2 {
  public:
    const unsigned int m, n; // n = 2^m points are generated

    // Construct a generator for 2^m points.
    GenDiag0m2(const unsigned int m)
    : m(m), n(1 << m), m2((m + 1) >> 1), mask(~-(1 << m2)) {
      const unsigned int sqrtN = 1 << (m >> 1);
      d = new unsigned int[sqrtN];

      if (m & 1) { // odd
        dx = dy = m >> 1;
        // precompute the diagonal shifts
        for (unsigned int k = 0; k < sqrtN; k++)
          d[k] = (vdc(k) >> (32 - m)) + k;
      }
      else { // even
        dx = (m >> 1) - 2;
        dy = m >> 1;
        // precompute the diagonal shifts
        for (unsigned int k = 0; k < sqrtN; k++)
          d[k] = (vdc(k) >> (32 - m)) + (k >> 2) + (1 << dx);
      }
    }

    // Destructor.
    ~GenDiag0m2() { delete [] d; }

    // Return the integer coordinates (x, y)
    // for the i-th point with i in {0, ..., n - 1}.
    inline void get(const unsigned int i, unsigned int &x, unsigned int &y) const {
      const unsigned int k = i >> m2; // determine the diagonal
      const unsigned int j = i & mask; // j-th point on the k-th diagonal

      // multiplication by shift, modulo by bitwise and
      x = (d[k] + (j << dx)) & (n - 1);
      y = k + (j << dy);
    }

    // Compute the i-th point with i in {0, ..., n - 1} with a shifted
    // tiling for the even m case, inside the pixel (px, py).
    // The resulting integer-scaled point will already have a shift
    // corresponding to the pixel, so after dividing it by 2^m
    // it will be in [px, px + 1) x [py, py + 1).
    inline void getShiftedTiling(const unsigned int px, const unsigned int py,
        const unsigned int i, unsigned int &x, unsigned int &y) const {
      assert(!(m & 1)); // only makes sense for the even m case
      get(i, x, y);

      // modulo-wrap depending on row
      x += (py & 3) * (1 << (m - 2));
      x &= ~-(1 << m);

      // shift for pixel
      x += px << m;
      y += py << m;
    }

  private:

    // 32 bit version of van der Corput radical inverse in base 2
    inline static unsigned int vdc(unsigned int bits) {
      bits = (bits << 16) | (bits >> 16);
      bits = ((bits & 0x00ff00ff) << 8) | ((bits & 0xff00ff00) >> 8);
      bits = ((bits & 0x0f0f0f0f) << 4) | ((bits & 0xf0f0f0f0) >> 4);
      bits = ((bits & 0x33333333) << 2) | ((bits & 0xcccccccc) >> 2);
      bits = ((bits & 0x55555555) << 1) | ((bits & 0xaaaaaaaa) >> 1);
      return bits;
    }

    const unsigned int m2, mask;
    unsigned int dx, dy, *d;
};

#endif

