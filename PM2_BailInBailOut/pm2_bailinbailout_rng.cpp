#include "pm2_bailinbailout.h"

/*<===========================================================================>
   
    Winter 2026
    Author: Kyryll Kotyk
    Advisor: Munehiro Fukuda 

    This file includes implementations of all the methods related to
    randomization, seeding, and item selection

    All randomization is seeded based on controlled factors
    that do not take rank / partitioning into account, ensuring
    same result regardless of how many processes are used to run the program

    Style Note:
        One parameter per line for function signatures
        Parantheses are closed on a new line, followed by body brace opening
        This rule does not apply if there are no parameters

<===========================================================================>*/

//<===========================================================================>
//<========================!Xoshiro256 Random Methods!========================>
//<===========================================================================>

uint64_t BailInBailOut::Xoshiro256::splitmix64(
    uint64_t& x
) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

BailInBailOut::Xoshiro256::Xoshiro256(
    uint64_t seed
) {
    for (int i = 0; i < 4; i++) {
        s[i] = splitmix64(seed);
    }
}

uint64_t BailInBailOut::Xoshiro256::rotl(
    uint64_t x,
    int k
) {
    return (x << k) | (x >> (64 - k));
}

uint64_t BailInBailOut::Xoshiro256::next() {
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rotl(s[3], 45);

    return result;
}

uint64_t BailInBailOut::Xoshiro256::nextInRange(
    uint64_t min, 
    uint64_t max
) {
    uint64_t range = max - min + 1;
    uint64_t x;
    uint64_t limit = numeric_limits<uint64_t>::max() - (numeric_limits
        <uint64_t>::max() % range);

    do {
        x = next();
    } while (x >= limit);

    return min + (x % range);
}

double BailInBailOut::Xoshiro256::nextDouble01() {
    return (next() >> 11) * (1.0 / 9007199254740992.0); // 2^53
}

double BailInBailOut::Xoshiro256::nextInRangeDouble(
    double min, 
    double max
) {
    return min + (max - min) * nextDouble01();
}

//<===========================================================================>
//<=====================!Misc Item Selection and Seeding!=====================>
//<===========================================================================>

uint64_t BailInBailOut::splitmix64Hash(
    uint64_t x
) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

unsigned int BailInBailOut::selectForWorker(
    uint64_t baseSeed,
    unsigned int workerGlobalId,
    unsigned int countTotal,
    uint64_t streamTag
) {

    uint64_t seed = makeSeed(baseSeed, 0, 0, workerGlobalId, streamTag);
    BailInBailOut::Xoshiro256 rng(seed);

    // Bank IDs are global
    return (unsigned int)rng.nextInRange(0, (uint64_t)countTotal - 1);
}


uint64_t BailInBailOut::makeSeed(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int timestep,
    unsigned int globalId,
    uint64_t streamTag
) {

    uint64_t key = baseSeed;

    // Mix in run number
    key ^= 0xD6E8FEB86659FD93ULL * (uint64_t)run;

    // Mix in timestep (0 for static parameters like graph/structure)
    key ^= 0xBF58476D1CE4E5B9ULL * (uint64_t)timestep;

    // Mix in entity global ID (bank/worker/firm)
    key ^= 0x9E3779B97F4A7C15ULL * (uint64_t)globalId;

    // Mix in stream tag to separate RNG 
    key ^= 0x94D049BB133111EBULL * streamTag;

    return splitmix64Hash(key);
}

unsigned short BailInBailOut::percent0to100FromMixedSeed(
    uint64_t seed
) {
    // Assumes mixed seed 

    uint64_t x = seed;

    const uint64_t range = 101;
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Retry by adding an odd constant 
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return (unsigned short)(x % range);
}

uint64_t BailInBailOut::randomInRangeFromSeed(
    uint64_t seed,
    uint64_t minVal,
    uint64_t maxVal
) {
    // Assumes mixed seed

    if (minVal > maxVal) {
        swap(minVal, maxVal);
    }

    uint64_t range = maxVal - minVal + 1;
    uint64_t x = seed;

    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Keep going until retry
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return minVal + (x % range);
}
