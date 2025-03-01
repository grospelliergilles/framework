﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SHA1HashAlgorithm.cc                                        (C) 2000-2023 */
/*                                                                           */
/* Calcule de fonction de hashage SHA-1.                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/SHA1HashAlgorithm.h"

#include "arcane/utils/Array.h"
#include "arcane/utils/FatalErrorException.h"
#include "arcane/utils/Ref.h"

#include <cstring>
#include <array>

// L'algorithme est decrit ici;
// https://en.wikipedia.org/wiki/SHA-3

// L'implémentation est issue du dépot suivant:
// https://github.com/rhash/RHash

/* sha1.c - an implementation of Secure Hash Algorithm 1 (SHA1)
 * based on RFC 3174.
 *
 * Copyright (c) 2008, Aleksey Kravchenko <rhash.admin@gmail.com>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE  INCLUDING ALL IMPLIED WARRANTIES OF  MERCHANTABILITY
 * AND FITNESS.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT,  OR CONSEQUENTIAL DAMAGES  OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE,  DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT,  NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF  OR IN CONNECTION  WITH THE USE  OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::SHA1Algorithm
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace
{
  constexpr int sha1_block_size = 64;
  constexpr int sha1_hash_size = 20;

  /* constants for SHA1 rounds */
  constexpr uint32_t K0 = 0x5a827999;
  constexpr uint32_t K1 = 0x6ed9eba1;
  constexpr uint32_t K2 = 0x8f1bbcdc;
  constexpr uint32_t K3 = 0xca62c1d6;

  uint32_t bswap_32(uint32_t x)
  {
    x = ((x << 8) & 0xFF00FF00u) | ((x >> 8) & 0x00FF00FFu);
    return (x >> 16) | (x << 16);
  }
  // ROTL/ROTR rotate a 32-bit word left by n bits
  static uint32_t _rotateLeft(uint32_t dword, int n)
  {
    return ((dword) << (n) ^ ((dword) >> (32 - (n))));
  }
  /**
   * Copy a memory block with simultaneous exchanging byte order.
   * The byte order is changed from little-endian 32-bit integers
   * to big-endian (or vice-versa).
   *
   * @param to the pointer where to copy memory block
   * @param index the index to start writing from
   * @param from  the source block to copy
   * @param length length of the memory block
   */
  void _swap_copy_str_to_u32(void* to, int index, const void* from, size_t length)
  {
    /* if all pointers and length are 32-bits aligned */
    if (0 == (((uintptr_t)to | (uintptr_t)from | (uintptr_t)index | length) & 3)) {
      /* copy memory as 32-bit words */
      const uint32_t* src = (const uint32_t*)from;
      const uint32_t* end = (const uint32_t*)((const char*)src + length);
      uint32_t* dst = (uint32_t*)((char*)to + index);
      for (; src < end; dst++, src++)
        *dst = bswap_32(*src);
    }
    else {
      const char* src = (const char*)from;
      for (length += index; (size_t)index < length; index++)
        ((char*)to)[index ^ 3] = *(src++);
    }
  }
} // namespace

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class SHA1
: public IHashAlgorithmContext
{
 private:

  //! Algorithm context.
  struct sha1_ctx
  {
    /* 512-bit buffer for leftovers */
    unsigned char message[sha1_block_size] = {};
    /* number of processed bytes */
    uint64_t length = 0;
    /* 160-bit algorithm internal hashing state */
    unsigned hash[5] = {};
  };
  sha1_ctx m_context;

 public:

  SHA1() { reset(); }

 public:

  void reset() override { sha1_init(); }

  void updateHash(Span<const std::byte> input) override
  {
    sha1_update(input);
  }

  void computeHashValue(HashAlgorithmValue& value) override
  {
    sha1_final(value);
  }

 private:

  void sha1_init();
  void sha1_update(Span<const std::byte> bytes);
  void sha1_final(HashAlgorithmValue& value);

  static void sha1_process_block(unsigned* hash, const unsigned* block);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Valide pour little-endian
//@{
#define be2me_32(x) bswap_32(x)
#define be32_copy(to, index, from, length) _swap_copy_str_to_u32((to), (index), (from), (length))
//@}

#define IS_ALIGNED_32(p) (0 == (3 & (uintptr_t)(p)))

/**
 * Initialize context before calculating hash.
 *
 * @param ctx context to initialize
 */
void SHA1::
sha1_init()
{
  sha1_ctx* ctx = &m_context;
  ctx->length = 0;

  /* initialize algorithm state */
  ctx->hash[0] = 0x67452301;
  ctx->hash[1] = 0xefcdab89;
  ctx->hash[2] = 0x98badcfe;
  ctx->hash[3] = 0x10325476;
  ctx->hash[4] = 0xc3d2e1f0;
}

/* round functions for SHA1 */
#define CHO(X, Y, Z) (((X) & (Y)) | ((~(X)) & (Z)))
#define PAR(X, Y, Z) ((X) ^ (Y) ^ (Z))
#define MAJ(X, Y, Z) (((X) & (Y)) | ((X) & (Z)) | ((Y) & (Z)))

#define ROTL32(a, b) _rotateLeft(a, b)
#define ROUND_0(a, b, c, d, e, FF, k, w) e += FF(b, c, d) + ROTL32(a, 5) + k + w
#define ROUND_1(a, b, c, d, e, FF, k, w) e += FF(b, ROTL32(c, 30), d) + ROTL32(a, 5) + k + w
#define ROUND_2(a, b, c, d, e, FF, k, w) e += FF(b, ROTL32(c, 30), ROTL32(d, 30)) + ROTL32(a, 5) + k + w
#define ROUND(a, b, c, d, e, FF, k, w) e = ROTL32(e, 30) + FF(b, ROTL32(c, 30), ROTL32(d, 30)) + ROTL32(a, 5) + k + w

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * The core transformation. Process a 512-bit block.
 * The function has been taken from RFC 3174 with little changes.
 *
 * @param hash algorithm state
 * @param block the message block to process
 */
void SHA1::
sha1_process_block(unsigned* hash, const unsigned* block)
{
  uint32_t W[80]; /* Word sequence */
  uint32_t A, B, C, D, E; /* Word buffers */

  A = hash[0];
  B = hash[1];
  C = hash[2];
  D = hash[3];
  E = hash[4];

  /* 0..19 */
  W[0] = be2me_32(block[0]);
  ROUND_0(A, B, C, D, E, CHO, K0, W[0]);
  W[1] = be2me_32(block[1]);
  ROUND_1(E, A, B, C, D, CHO, K0, W[1]);
  W[2] = be2me_32(block[2]);
  ROUND_2(D, E, A, B, C, CHO, K0, W[2]);
  W[3] = be2me_32(block[3]);
  ROUND(C, D, E, A, B, CHO, K0, W[3]);
  W[4] = be2me_32(block[4]);
  ROUND(B, C, D, E, A, CHO, K0, W[4]);

  W[5] = be2me_32(block[5]);
  ROUND(A, B, C, D, E, CHO, K0, W[5]);
  W[6] = be2me_32(block[6]);
  ROUND(E, A, B, C, D, CHO, K0, W[6]);
  W[7] = be2me_32(block[7]);
  ROUND(D, E, A, B, C, CHO, K0, W[7]);
  W[8] = be2me_32(block[8]);
  ROUND(C, D, E, A, B, CHO, K0, W[8]);
  W[9] = be2me_32(block[9]);
  ROUND(B, C, D, E, A, CHO, K0, W[9]);

  W[10] = be2me_32(block[10]);
  ROUND(A, B, C, D, E, CHO, K0, W[10]);
  W[11] = be2me_32(block[11]);
  ROUND(E, A, B, C, D, CHO, K0, W[11]);
  W[12] = be2me_32(block[12]);
  ROUND(D, E, A, B, C, CHO, K0, W[12]);
  W[13] = be2me_32(block[13]);
  ROUND(C, D, E, A, B, CHO, K0, W[13]);
  W[14] = be2me_32(block[14]);
  ROUND(B, C, D, E, A, CHO, K0, W[14]);

  W[15] = be2me_32(block[15]);
  ROUND(A, B, C, D, E, CHO, K0, W[15]);
  W[16] = ROTL32(W[13] ^ W[8] ^ W[2] ^ W[0], 1);
  ROUND(E, A, B, C, D, CHO, K0, W[16]);
  W[17] = ROTL32(W[14] ^ W[9] ^ W[3] ^ W[1], 1);
  ROUND(D, E, A, B, C, CHO, K0, W[17]);
  W[18] = ROTL32(W[15] ^ W[10] ^ W[4] ^ W[2], 1);
  ROUND(C, D, E, A, B, CHO, K0, W[18]);
  W[19] = ROTL32(W[16] ^ W[11] ^ W[5] ^ W[3], 1);
  ROUND(B, C, D, E, A, CHO, K0, W[19]);
  /* 20..39 */
  W[20] = ROTL32(W[17] ^ W[12] ^ W[6] ^ W[4], 1);
  ROUND(A, B, C, D, E, PAR, K1, W[20]);
  W[21] = ROTL32(W[18] ^ W[13] ^ W[7] ^ W[5], 1);
  ROUND(E, A, B, C, D, PAR, K1, W[21]);
  W[22] = ROTL32(W[19] ^ W[14] ^ W[8] ^ W[6], 1);
  ROUND(D, E, A, B, C, PAR, K1, W[22]);
  W[23] = ROTL32(W[20] ^ W[15] ^ W[9] ^ W[7], 1);
  ROUND(C, D, E, A, B, PAR, K1, W[23]);
  W[24] = ROTL32(W[21] ^ W[16] ^ W[10] ^ W[8], 1);
  ROUND(B, C, D, E, A, PAR, K1, W[24]);

  W[25] = ROTL32(W[22] ^ W[17] ^ W[11] ^ W[9], 1);
  ROUND(A, B, C, D, E, PAR, K1, W[25]);
  W[26] = ROTL32(W[23] ^ W[18] ^ W[12] ^ W[10], 1);
  ROUND(E, A, B, C, D, PAR, K1, W[26]);
  W[27] = ROTL32(W[24] ^ W[19] ^ W[13] ^ W[11], 1);
  ROUND(D, E, A, B, C, PAR, K1, W[27]);
  W[28] = ROTL32(W[25] ^ W[20] ^ W[14] ^ W[12], 1);
  ROUND(C, D, E, A, B, PAR, K1, W[28]);
  W[29] = ROTL32(W[26] ^ W[21] ^ W[15] ^ W[13], 1);
  ROUND(B, C, D, E, A, PAR, K1, W[29]);

  W[30] = ROTL32(W[27] ^ W[22] ^ W[16] ^ W[14], 1);
  ROUND(A, B, C, D, E, PAR, K1, W[30]);
  W[31] = ROTL32(W[28] ^ W[23] ^ W[17] ^ W[15], 1);
  ROUND(E, A, B, C, D, PAR, K1, W[31]);
  W[32] = ROTL32(W[29] ^ W[24] ^ W[18] ^ W[16], 1);
  ROUND(D, E, A, B, C, PAR, K1, W[32]);
  W[33] = ROTL32(W[30] ^ W[25] ^ W[19] ^ W[17], 1);
  ROUND(C, D, E, A, B, PAR, K1, W[33]);
  W[34] = ROTL32(W[31] ^ W[26] ^ W[20] ^ W[18], 1);
  ROUND(B, C, D, E, A, PAR, K1, W[34]);

  W[35] = ROTL32(W[32] ^ W[27] ^ W[21] ^ W[19], 1);
  ROUND(A, B, C, D, E, PAR, K1, W[35]);
  W[36] = ROTL32(W[33] ^ W[28] ^ W[22] ^ W[20], 1);
  ROUND(E, A, B, C, D, PAR, K1, W[36]);
  W[37] = ROTL32(W[34] ^ W[29] ^ W[23] ^ W[21], 1);
  ROUND(D, E, A, B, C, PAR, K1, W[37]);
  W[38] = ROTL32(W[35] ^ W[30] ^ W[24] ^ W[22], 1);
  ROUND(C, D, E, A, B, PAR, K1, W[38]);
  W[39] = ROTL32(W[36] ^ W[31] ^ W[25] ^ W[23], 1);
  ROUND(B, C, D, E, A, PAR, K1, W[39]);
  /* 40..59 */
  W[40] = ROTL32(W[37] ^ W[32] ^ W[26] ^ W[24], 1);
  ROUND(A, B, C, D, E, MAJ, K2, W[40]);
  W[41] = ROTL32(W[38] ^ W[33] ^ W[27] ^ W[25], 1);
  ROUND(E, A, B, C, D, MAJ, K2, W[41]);
  W[42] = ROTL32(W[39] ^ W[34] ^ W[28] ^ W[26], 1);
  ROUND(D, E, A, B, C, MAJ, K2, W[42]);
  W[43] = ROTL32(W[40] ^ W[35] ^ W[29] ^ W[27], 1);
  ROUND(C, D, E, A, B, MAJ, K2, W[43]);
  W[44] = ROTL32(W[41] ^ W[36] ^ W[30] ^ W[28], 1);
  ROUND(B, C, D, E, A, MAJ, K2, W[44]);

  W[45] = ROTL32(W[42] ^ W[37] ^ W[31] ^ W[29], 1);
  ROUND(A, B, C, D, E, MAJ, K2, W[45]);
  W[46] = ROTL32(W[43] ^ W[38] ^ W[32] ^ W[30], 1);
  ROUND(E, A, B, C, D, MAJ, K2, W[46]);
  W[47] = ROTL32(W[44] ^ W[39] ^ W[33] ^ W[31], 1);
  ROUND(D, E, A, B, C, MAJ, K2, W[47]);
  W[48] = ROTL32(W[45] ^ W[40] ^ W[34] ^ W[32], 1);
  ROUND(C, D, E, A, B, MAJ, K2, W[48]);
  W[49] = ROTL32(W[46] ^ W[41] ^ W[35] ^ W[33], 1);
  ROUND(B, C, D, E, A, MAJ, K2, W[49]);

  W[50] = ROTL32(W[47] ^ W[42] ^ W[36] ^ W[34], 1);
  ROUND(A, B, C, D, E, MAJ, K2, W[50]);
  W[51] = ROTL32(W[48] ^ W[43] ^ W[37] ^ W[35], 1);
  ROUND(E, A, B, C, D, MAJ, K2, W[51]);
  W[52] = ROTL32(W[49] ^ W[44] ^ W[38] ^ W[36], 1);
  ROUND(D, E, A, B, C, MAJ, K2, W[52]);
  W[53] = ROTL32(W[50] ^ W[45] ^ W[39] ^ W[37], 1);
  ROUND(C, D, E, A, B, MAJ, K2, W[53]);
  W[54] = ROTL32(W[51] ^ W[46] ^ W[40] ^ W[38], 1);
  ROUND(B, C, D, E, A, MAJ, K2, W[54]);

  W[55] = ROTL32(W[52] ^ W[47] ^ W[41] ^ W[39], 1);
  ROUND(A, B, C, D, E, MAJ, K2, W[55]);
  W[56] = ROTL32(W[53] ^ W[48] ^ W[42] ^ W[40], 1);
  ROUND(E, A, B, C, D, MAJ, K2, W[56]);
  W[57] = ROTL32(W[54] ^ W[49] ^ W[43] ^ W[41], 1);
  ROUND(D, E, A, B, C, MAJ, K2, W[57]);
  W[58] = ROTL32(W[55] ^ W[50] ^ W[44] ^ W[42], 1);
  ROUND(C, D, E, A, B, MAJ, K2, W[58]);
  W[59] = ROTL32(W[56] ^ W[51] ^ W[45] ^ W[43], 1);
  ROUND(B, C, D, E, A, MAJ, K2, W[59]);
  /* 60..79 */
  W[60] = ROTL32(W[57] ^ W[52] ^ W[46] ^ W[44], 1);
  ROUND(A, B, C, D, E, PAR, K3, W[60]);
  W[61] = ROTL32(W[58] ^ W[53] ^ W[47] ^ W[45], 1);
  ROUND(E, A, B, C, D, PAR, K3, W[61]);
  W[62] = ROTL32(W[59] ^ W[54] ^ W[48] ^ W[46], 1);
  ROUND(D, E, A, B, C, PAR, K3, W[62]);
  W[63] = ROTL32(W[60] ^ W[55] ^ W[49] ^ W[47], 1);
  ROUND(C, D, E, A, B, PAR, K3, W[63]);
  W[64] = ROTL32(W[61] ^ W[56] ^ W[50] ^ W[48], 1);
  ROUND(B, C, D, E, A, PAR, K3, W[64]);

  W[65] = ROTL32(W[62] ^ W[57] ^ W[51] ^ W[49], 1);
  ROUND(A, B, C, D, E, PAR, K3, W[65]);
  W[66] = ROTL32(W[63] ^ W[58] ^ W[52] ^ W[50], 1);
  ROUND(E, A, B, C, D, PAR, K3, W[66]);
  W[67] = ROTL32(W[64] ^ W[59] ^ W[53] ^ W[51], 1);
  ROUND(D, E, A, B, C, PAR, K3, W[67]);
  W[68] = ROTL32(W[65] ^ W[60] ^ W[54] ^ W[52], 1);
  ROUND(C, D, E, A, B, PAR, K3, W[68]);
  W[69] = ROTL32(W[66] ^ W[61] ^ W[55] ^ W[53], 1);
  ROUND(B, C, D, E, A, PAR, K3, W[69]);

  W[70] = ROTL32(W[67] ^ W[62] ^ W[56] ^ W[54], 1);
  ROUND(A, B, C, D, E, PAR, K3, W[70]);
  W[71] = ROTL32(W[68] ^ W[63] ^ W[57] ^ W[55], 1);
  ROUND(E, A, B, C, D, PAR, K3, W[71]);
  W[72] = ROTL32(W[69] ^ W[64] ^ W[58] ^ W[56], 1);
  ROUND(D, E, A, B, C, PAR, K3, W[72]);
  W[73] = ROTL32(W[70] ^ W[65] ^ W[59] ^ W[57], 1);
  ROUND(C, D, E, A, B, PAR, K3, W[73]);
  W[74] = ROTL32(W[71] ^ W[66] ^ W[60] ^ W[58], 1);
  ROUND(B, C, D, E, A, PAR, K3, W[74]);

  W[75] = ROTL32(W[72] ^ W[67] ^ W[61] ^ W[59], 1);
  ROUND(A, B, C, D, E, PAR, K3, W[75]);
  W[76] = ROTL32(W[73] ^ W[68] ^ W[62] ^ W[60], 1);
  ROUND(E, A, B, C, D, PAR, K3, W[76]);
  W[77] = ROTL32(W[74] ^ W[69] ^ W[63] ^ W[61], 1);
  ROUND(D, E, A, B, C, PAR, K3, W[77]);
  W[78] = ROTL32(W[75] ^ W[70] ^ W[64] ^ W[62], 1);
  ROUND(C, D, E, A, B, PAR, K3, W[78]);
  W[79] = ROTL32(W[76] ^ W[71] ^ W[65] ^ W[63], 1);
  ROUND(B, C, D, E, A, PAR, K3, W[79]);

  hash[0] += A;
  hash[1] += B;
  hash[2] += ROTL32(C, 30);
  hash[3] += ROTL32(D, 30);
  hash[4] += ROTL32(E, 30);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
 * Calculate message hash.
 * Can be called repeatedly with chunks of the message to be hashed.
 *
 * @param bytes message chunk
 */
void SHA1::
sha1_update(Span<const std::byte> bytes)
{
  sha1_ctx* ctx = &m_context;
  const unsigned char* msg = reinterpret_cast<const unsigned char*>(bytes.data());
  size_t size = bytes.size();
  unsigned index = (unsigned)ctx->length & 63;
  ctx->length += size;

  /* fill partial block */
  if (index) {
    unsigned left = sha1_block_size - index;
    memcpy(ctx->message + index, msg, (size < left ? size : left));
    if (size < left)
      return;

    /* process partial block */
    sha1_process_block(ctx->hash, (unsigned*)ctx->message);
    msg += left;
    size -= left;
  }
  while (size >= sha1_block_size) {
    unsigned* aligned_message_block;
    if (IS_ALIGNED_32(msg)) {
      /* the most common case is processing of an already aligned message
			without copying it */
      aligned_message_block = (unsigned*)msg;
    }
    else {
      memcpy(ctx->message, msg, sha1_block_size);
      aligned_message_block = (unsigned*)ctx->message;
    }

    sha1_process_block(ctx->hash, aligned_message_block);
    msg += sha1_block_size;
    size -= sha1_block_size;
  }
  if (size) {
    /* save leftovers */
    memcpy(ctx->message, msg, size);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1::
sha1_final(HashAlgorithmValue& value)
{
  value.setSize(sha1_hash_size);

  unsigned char* result = reinterpret_cast<unsigned char*>(value.bytes().data());
  sha1_ctx* ctx = &m_context;
  unsigned index = (unsigned)ctx->length & 63;
  unsigned* msg32 = (unsigned*)ctx->message;

  /* pad message and run for last block */
  ctx->message[index++] = 0x80;
  while ((index & 3) != 0) {
    ctx->message[index++] = 0;
  }
  index >>= 2;

  /* if no room left in the message to store 64-bit message length */
  if (index > 14) {
    /* then fill the rest with zeros and process it */
    while (index < 16) {
      msg32[index++] = 0;
    }
    sha1_process_block(ctx->hash, msg32);
    index = 0;
  }
  while (index < 14) {
    msg32[index++] = 0;
  }
  msg32[14] = be2me_32((unsigned)(ctx->length >> 29));
  msg32[15] = be2me_32((unsigned)(ctx->length << 3));
  sha1_process_block(ctx->hash, msg32);

  if (result)
    be32_copy(result, 0, &ctx->hash, sha1_hash_size);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::SHA1Algorithm

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
_computeHash(Span<const std::byte> input, HashAlgorithmValue& value)
{
  SHA1Algorithm::SHA1 sha1;
  sha1.updateHash(input);
  sha1.computeHashValue(value);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
_computeHash64(Span<const std::byte> input, ByteArray& output)
{
  HashAlgorithmValue value;
  _computeHash(input, value);
  output.addRange(value.asLegacyBytes());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
computeHash(Span<const std::byte> input, HashAlgorithmValue& value)
{
  _computeHash(input, value);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
computeHash(ByteConstArrayView input, ByteArray& output)
{
  Span<const Byte> input64(input);
  _computeHash64(asBytes(input64), output);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
computeHash64(Span<const Byte> input, ByteArray& output)
{
  Span<const std::byte> bytes(asBytes(input));
  _computeHash64(bytes, output);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SHA1HashAlgorithm::
computeHash64(Span<const std::byte> input, ByteArray& output)
{
  _computeHash64(input, output);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Ref<IHashAlgorithmContext> SHA1HashAlgorithm::
createContext()
{
  return makeRef<IHashAlgorithmContext>(new SHA1Algorithm::SHA1());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
