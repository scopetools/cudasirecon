// Self-contained, dependency-free unit test for dvfile.h (the IVE-replacement
// DeltaVision/MRC reader/writer).  Exercises the round-trip read/write + pixel-type
// conversion + big-endian rejection paths without needing IVE, a GPU, or fixtures.
//
// Built and registered with CTest by src/CMakeLists.txt.  Returns non-zero on any
// failure.
#include "dvfile.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(cond, msg)                                                      \
  do {                                                                        \
    if (!(cond)) {                                                            \
      std::printf("FAIL: %s\n", msg);                                         \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

// Build a minimal valid single-wave DV header for an nx*ny*nz stack of `mode`.
static IW_MRC_HEADER makeHeader(int nx, int ny, int nz, int mode) {
  IW_MRC_HEADER h{};
  h.nx = nx;
  h.ny = ny;
  h.nz = nz;
  h.mode = mode;
  h.mx = nx;
  h.my = ny;
  h.mz = nz;
  h.xlen = 1.f;
  h.ylen = 1.f;
  h.zlen = 1.f;
  h.mapc = 1;
  h.mapr = 2;
  h.maps = 3;
  h.nDVID = (int16_t)0xC0A0;  // little-endian DV magic (bytes A0 C0 at offset 96)
  h.num_times = 1;
  h.interleaved = ZTW_SEQUENCE;
  h.num_waves = 1;
  h.nlab = 0;
  return h;
}

// Round-trip a float stack through a given stored `mode`, reading back as float.
static void roundtrip(const std::string& path, int mode, bool intType) {
  const int nx = 8, ny = 6, nz = 3;
  std::vector<float> orig(nx * ny * nz);
  for (size_t i = 0; i < orig.size(); ++i)
    orig[i] = intType ? float((i * 7) % 200)        // integer-valued
                      : (float(i) - 50.f) * 1.5f;    // signed fractional

  IMOpen(0, path.c_str(), "new");
  IW_MRC_HEADER h = makeHeader(nx, ny, nz, mode);
  IMPutHdr(0, &h);
  IMAlCon(0, 1);  // conversion ON (float <-> stored type)
  for (int z = 0; z < nz; ++z)
    IMWrSec(0, orig.data() + (size_t)z * nx * ny);
  IMWrHdr(0, h.label, 0, 0.f, 1.f, 0.5f);
  IMClose(0);

  IMOpen(1, path.c_str(), "ro");
  int ixyz[3], mxyz[3], rmode;
  float mn, mx, me;
  IMRdHdr(1, ixyz, mxyz, &rmode, &mn, &mx, &me);
  CHECK(rmode == mode, ("stored mode mismatch for " + path).c_str());
  CHECK(ixyz[0] == nx && ixyz[1] == ny && ixyz[2] == nz,
        ("dims mismatch for " + path).c_str());
  IMAlCon(1, 1);
  std::vector<float> back(orig.size());
  for (int z = 0; z < nz; ++z)
    IMRdSec(1, back.data() + (size_t)z * nx * ny);
  IMClose(1);

  double maxerr = 0;
  for (size_t i = 0; i < orig.size(); ++i)
    maxerr = std::max(maxerr, (double)std::fabs(back[i] - orig[i]));
  // float/complex must be exact; integer types exact for integer-valued data.
  CHECK(maxerr == 0.0, ("round-trip not exact for mode " + std::to_string(mode)).c_str());
  std::printf("  mode %d round-trip maxerr=%g %s\n", mode, maxerr,
              maxerr == 0.0 ? "OK" : "FAIL");
}

int main() {
  std::string tmp = "dvfile_test_tmp";

  CHECK(sizeof(IW_MRC_HEADER) == 1024, "IW_MRC_HEADER must be 1024 bytes");

  std::printf("round-trip read/write + conversion:\n");
  roundtrip(tmp + "_f32.dv", IW_FLOAT, false);
  roundtrip(tmp + "_u16.dv", IW_USHORT, true);
  roundtrip(tmp + "_i16.dv", IW_SHORT, true);
  roundtrip(tmp + "_u8.dv", IW_BYTE, true);
  roundtrip(tmp + "_i32.dv", IW_LONG, true);

  // Big-endian files must be rejected, not silently mis-read.
  {
    std::string be = tmp + "_be.dv";
    IMOpen(0, (tmp + "_f32.dv").c_str(), "ro");  // reuse a valid LE file's bytes
    IMClose(0);
    std::vector<char> bytes;
    {
      std::ifstream in(tmp + "_f32.dv", std::ios::binary);
      bytes.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    }
    bytes[96] = (char)0xC0;  // flip magic to big-endian
    bytes[97] = (char)0xA0;
    std::ofstream(be, std::ios::binary).write(bytes.data(), bytes.size());
    bool rejected = false;
    try {
      DVFile f(be);
    } catch (const std::exception&) {
      rejected = true;
    }
    CHECK(rejected, "big-endian file was not rejected");
    std::printf("big-endian rejection: %s\n", rejected ? "OK" : "FAIL");
  }

  std::printf(failures ? "\n%d CHECK(S) FAILED\n" : "\nALL PASSED\n", failures);
  return failures ? 1 : 0;
}
