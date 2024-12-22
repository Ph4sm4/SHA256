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
    void processBlock(bitset<SHA256_BLOCK_SIZE> &block);

    string createOutputHash() const;

    vector<bitset<SHA256_BLOCK_SIZE>> blocks;
    bitset<SHA256_BIT_SIZE>
        hashBits;

    void printMessageBlock() const;

    int blocksNumBeforePadding;
    int strBitSize;
};

int main()
{

    SHA256 obj;

    obj.hash("RedBlockBlue");

    return 0;
}

string SHA256::hash(const string &msg)
{
    this->prepareBlocks(msg);
    this->processBlocks();

    return this->createOutputHash();
}

void SHA256::prepareBlocks(const string &str)
{
    const int charsPerBlock = 64;
    blocksNumBeforePadding = ceil(float(str.length() * CHAR_BITS) / SHA256_BLOCK_SIZE);
    strBitSize = str.length() * CHAR_BITS;

    cout
        << "num of blocks: " << blocksNumBeforePadding << endl;

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

    printMessageBlock();
}

void SHA256::blockPadding(bitset<SHA256_BLOCK_SIZE> &block)
{
    // we round up to the nearest multiple of 512
    int properMessageBlockSize = ceil(float(strBitSize + 1 + 64) / 512.f) * SHA256_BLOCK_SIZE;
    // calculate the number of additional 0's required between the actual message bits and the length bits
    int k = properMessageBlockSize - strBitSize - 1 - 64;
    cout << "obtained: " << k << " " << properMessageBlockSize << endl;

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

void SHA256::processBlock(bitset<SHA256_BLOCK_SIZE> &block)
{
    WORD chunks[16];
    // initially chunks are 0
    for (int i = 0; i < 16; i++)
        chunks[i] = 0;

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
    cout << endl;
    for (int i = 0; i < 16; i++)
    {
        cout << bitset<32>(chunks[i]) << endl;
    }
    cout << endl;
}

/*
    creates an output hash using hexadecimal string representation
*/
string SHA256::createOutputHash() const
{
    std::stringstream hex_stream;

    // iterate every 4 bits
    for (int i = 0; i < SHA256_BIT_SIZE; i += 4)
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
        unsigned int hexValue = (hashBits[i] << 3) | (hashBits[i + 1] << 2) | (hashBits[i + 2] << 1) | (hashBits[i + 3]);

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
