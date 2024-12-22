#include <iostream>
#include <bitset>
#include <sstream>
#include <vector>
#include <math.h>
#include "sha256.h"
#define SHA256_BIT_SIZE 256
#define SHA256_BLOCK_SIZE 512
#define BYTE 8
#define CHAR_BITS sizeof(char) * BYTE
#define WORD uint32_t

using namespace std;

/*
    PREPROCESSING

    1. for each character of the input string convert it to the binary (every char is 8 bits long / 1 byte )
        lets say our input string has 5 chars: 5 x 8 bits = 40 bits in total for our string

    2. we need to padd the message so that the number of bits is a multiple of 512 (it can also be 1024, depends on the implementation)
        to do that we append the bit "1" to the end of the message, followed by k (k > 0) is the smallest solution to the equation:
        l + 1 + k = 448 mod 512 -> k = 448 - , where l = 40 in our case

        then append the 64 bit block that is equal to the number l (40) expressed in binary, so a lot of 0000000...101000


*/

class SHA256
{
public:
    string hash(const string &msg);

private:
    void prepareBlocks(const string &str);
    void blockPadding(bitset<SHA256_BLOCK_SIZE> &block);
    void processBlocks();
    void processBlock(const bitset<SHA256_BLOCK_SIZE> &block);

    void prepareWORDChunks(WORD *chunks, const bitset<SHA256_BLOCK_SIZE> &block);
    void prepareMsgSchedule(WORD *msgSchedules, const WORD *chunks);

    void computeHash(WORD *chunks, WORD *msgSchedules);

    WORD ROTR(WORD x, int n) const;
    WORD CH(WORD e, WORD f, WORD g) const;
    WORD MAJ(WORD a, WORD b, WORD c) const;

    string createOutputHash() const;

    vector<bitset<SHA256_BLOCK_SIZE>> blocks;
    bitset<SHA256_BIT_SIZE>
        hashBits;

    void printMessageBlock() const;

    // initial hash values (first 32 bits of the fractional parts of the square roots of the first 8 primes 2..19):
    WORD hVals[8] = {
        0x6a09e667,
        0xbb67ae85,
        0x3c6ef372,
        0xa54ff53a,
        0x510e527f,
        0x9b05688c,
        0x1f83d9ab,
        0x5be0cd19};

    // constants (first 32 bits of the fractional parts of the cube roots of the first 64 primes 2..311):
    WORD k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

    int strBitSize;
};

int main()
{

    SHA256 obj;

    cout << endl
         << obj.hash("Olaf Dalach - Ph4sm4Solutions 12.23.2024") << endl
         << endl;

    return 0;
}

string SHA256::hash(const string &msg)
{
    this->prepareBlocks(msg);
    this->processBlocks();

    // combine the hash values into a single 256 bit bitset
    for (int i = 0; i < 8; i++)
    {
        bitset<32> hWord(hVals[i]);
        for (int j = SHA256_BIT_SIZE - 1 - i * 32; j >= SHA256_BIT_SIZE - 1 - (i + 1) * 32; j--)
        {
            if (j < 0)
                break;
            hashBits[j] = hWord[j % 32];
        }
    }

    return this->createOutputHash();
}

void SHA256::prepareBlocks(const string &str)
{
    const int charsPerBlock = 64;
    int blocksNumBeforePadding = ceil(float(str.length() * CHAR_BITS) / SHA256_BLOCK_SIZE);
    strBitSize = str.length() * CHAR_BITS;

    for (int blockIndex = 0; blockIndex < blocksNumBeforePadding; blockIndex++)
    {
        bitset<SHA256_BLOCK_SIZE> block;
        for (int i = charsPerBlock * blockIndex; i < charsPerBlock * (blockIndex + 1); i++)
        {
            if (i >= str.length())
                break;

            const char c = str[i];
            bitset<CHAR_BITS> charBits(c);
            for (int j = 0; j < CHAR_BITS; j++)
            {
                // the i % charsPerBlock is needed as we dont want to overflow the 512 range for higher indexed blocks
                block[SHA256_BLOCK_SIZE - 1 - CHAR_BITS * (i % charsPerBlock) - j] = charBits[CHAR_BITS - 1 - j]; // indexing is reversed - the higher the index, the more significant bit we deal with
            }
        }

        this->blocks.push_back(block);
    }
    blockPadding(blocks.at(blocks.size() - 1)); // padd the last block (may result in creating an additional one in case of an overflow)
}

void SHA256::blockPadding(bitset<SHA256_BLOCK_SIZE> &block)
{
    // we round up to the nearest multiple of 512
    int properMessageBlockSize = ceil(float(strBitSize + 1 + 64) / 512.f) * SHA256_BLOCK_SIZE;
    // calculate the number of additional 0's required between the actual message bits and the length bits
    int k = properMessageBlockSize - strBitSize - 1 - 64;

    block[SHA256_BLOCK_SIZE - ((strBitSize + 1) % SHA256_BLOCK_SIZE)] = 1; // add a '1' bit

    // 64 bit big endian number representing the length of the message that needs to get added at the end of the
    // message block
    bitset<64> bigEndianStrSize(strBitSize);

    bitset<SHA256_BLOCK_SIZE> additionalBlock;

    for (int i = 0; i < 64; i++)
    {
        // copy the bits from the 64 bit bigEndianStrSize to the proper block (attach the number at the end)
        // thus we are not reversing the assign index - we are starting from i = 0
        if (k > 447)
            additionalBlock[i] = bigEndianStrSize[i];
        else
            block[i] = bigEndianStrSize[i];
    }
    if (k > 447)
        this->blocks.push_back(additionalBlock);
}

void SHA256::processBlocks()
{
    for (bitset<SHA256_BLOCK_SIZE> block : this->blocks)
    {
        this->processBlock(block);
    }
}

void SHA256::processBlock(const bitset<SHA256_BLOCK_SIZE> &block)
{
    WORD chunks[16];
    // initially chunks are 0
    for (int i = 0; i < 16; i++)
        chunks[i] = 0;

    this->prepareWORDChunks(chunks, block);

    WORD msgSchedules[64];
    for (int i = 0; i < 16; i++)
        msgSchedules[i] = 0;

    this->prepareMsgSchedule(msgSchedules, chunks);

    this->computeHash(chunks, msgSchedules);
}

void SHA256::prepareWORDChunks(WORD *chunks, const bitset<SHA256_BLOCK_SIZE> &block)
{
    int ind = 0;
    for (int i = SHA256_BLOCK_SIZE - 1; i >= 0; i--)
    {
        // here we copy bit by bit the values from the block into our 32-bit chunks
        // by OR'ing the current WORD value with block's bit shifted proper number of times to the left
        // at first iter it would look like 1 (31 zeros), then 1 (30 zeros) etc.
        chunks[ind] = chunks[ind] | block[i] << i % 32;
        // switch to the next WORD
        if (i % 32 == 0)
            ind++;
    }
}

void SHA256::prepareMsgSchedule(WORD *msgSchedules, const WORD *chunks)
{
    // addition is performed using MOD 2^32

    // first 16 W (W0, W1, ..., W15) are just equal to our computed chunks
    for (int i = 0; i < 16; i++)
    {
        msgSchedules[i] = chunks[i];
    }

    /* the rest of the schedules W (W16, W17, ..., W63) are computed using this formula:
        sigma1(Wt-2) + Wt-7 + sigma0(Wt-15) + Wt-16

        where:

        sigma0(x) = ROTR7(x) | ROTR18(x) | SHR3(x)
        sigma1(x) = ROTR17(x) | ROTR19(x) | SHR10(x)

        ROTR[IND] (rotate right) means that we take the 32 bit WORD and shift every single bit IND times to the right
        while also looping the IND most right bits around to the start

        SHR[IND] is just shift right: x >> IND
    */

    for (int i = 16; i < 64; i++)
    {
        WORD x2 = msgSchedules[i - 2];
        WORD x7 = msgSchedules[i - 7];
        WORD x15 = msgSchedules[i - 15];
        WORD x16 = msgSchedules[i - 16];

        WORD sigma0 = ROTR(x15, 7) ^ ROTR(x15, 18) ^ (x15 >> 3);
        WORD sigma1 = ROTR(x2, 17) ^ ROTR(x2, 19) ^ (x2 >> 10);

        msgSchedules[i] = sigma1 + x7 + sigma0 + x16;
    }
}

void SHA256::computeHash(WORD *chunks, WORD *msgSchedules)
{
    int a = hVals[0], b = hVals[1], c = hVals[2], d = hVals[3],
        e = hVals[4], f = hVals[5], g = hVals[6], h = hVals[7];

    for (int i = 0; i < 64; i++)
    {
        WORD epsilon0 = ROTR(a, 2) ^ ROTR(a, 13) ^ ROTR(a, 22);
        WORD epsilon1 = ROTR(e, 6) ^ ROTR(e, 11) ^ ROTR(e, 25);

        WORD T1 = h + epsilon1 + CH(e, f, g) + k[i] + msgSchedules[i];
        WORD T2 = epsilon0 + MAJ(a, b, c);

        h = g;
        g = f;
        f = e;
        e = d + T1;
        d = c;
        c = b;
        b = a;
        a = T1 + T2;
    }

    hVals[0] += a;
    hVals[1] += b;
    hVals[2] += c;
    hVals[3] += d;
    hVals[4] += e;
    hVals[5] += f;
    hVals[6] += g;
    hVals[7] += h;
}

WORD SHA256::ROTR(WORD x, int n) const
{
    return (x >> n) | (x << (32 - n));
}

// CH - choose
WORD SHA256::CH(WORD e, WORD f, WORD g) const
{
    // this equation can be obtained by creating a truth table
    return ((~e) & g) ^ (e & f);
}

// MAJ - majority
WORD SHA256::MAJ(WORD a, WORD b, WORD c) const
{
    return (a & b) ^ (a & c) ^ (b & c);
}

/*
    creates an output hash using hexadecimal string representation
*/
string SHA256::createOutputHash() const
{
    std::stringstream hex_stream;

    // iterate every 4 bits
    for (int i = SHA256_BIT_SIZE - 1; i >= 0; i -= 4)
    {
        // 4 bit number (0 - 15), if the i, ..., i+3 bits are 1101, then each expression is gonna generate a number with only one bit switched on or off, based on whats in hashBits
        // and then all bits from those 4 numbers are gonna get put in a column, and we gonna perform an OR operation on every column of bits
        /*
          1000  (from hashBits[i] << 3)
        | 0100 (from hashBits[i + 1] << 2)
        | 0000 (from hashBits[i + 2] << 1)
        | 0001 (from hashBits[i + 3])
        resulting in 1101 (dec. 13)

        */
        unsigned int hexValue = (hashBits[i] << 3) | (hashBits[i - 1] << 2) | (hashBits[i - 2] << 1) | (hashBits[i - 3]);

        hex_stream << hex << hexValue;
    }

    return hex_stream.str();
}

void SHA256::printMessageBlock() const
{
    for (bitset<SHA256_BLOCK_SIZE> block : this->blocks)
    {
        int ct = 0;
        int colCt = 0;
        for (int i = SHA256_BLOCK_SIZE - 1; i >= 0; i--)
        {
            if (ct == 8)
            {
                ct = 0;
                colCt++;
                cout << " ";
            }

            if (colCt == 4)
            {
                cout << endl;
                colCt = 0;
            }

            if (ct < 8)
            {
                cout << block[i];
                ct++;
            }
        }
        cout << endl
             << endl;
    }
}
