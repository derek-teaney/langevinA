#ifndef UUID_2112b883_2d99_4aea_8d23_476240d382da
#define UUID_2112b883_2d99_4aea_8d23_476240d382da

/* A C++ wrapper for xoroshiro128+ (see http://xoroshiro.di.unimi.it/),
 * designed as an alternative to std::mt19937.
 *
 * Usage:
 *     // with some fixed seed
 *     xoroshiro128plus prng;
 *     uint64_t random_number = prng();
 *
 *     // with a random 32 bit seed
 *     std::random_device rd;
 *     xoroshiro128plus prng{rd()};
 *     uint64_t random_number = prng();
 *
 *     // with a random seed for the whole 128 bits
 *     struct random_seeder {
 *         // https://reddit.com/comments/6e10b8/_/di759pn/
 *         template <typename T> void generate(T begin, T end) const {
 *             for (std::random_device r; begin != end; ++begin) *begin = r();
 *         }
 *     };
 *     //...
 *     random_seeder seeder;
 *     xoroshiro128plus prng{seeder};
 *     uint64_t random_number = prng();
 *
 *
 * xoroshiro128+ was
 * Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 * To the extent possible under law, the author has dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 *
 * See <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <type_traits>

template <class ...>
using void_t = void;

template <class Seq, class = void>
struct is_seed_sequence : std::false_type {};

template <class Seq>
struct is_seed_sequence<Seq, void_t<
    //typename Seq::result_type,
    decltype(
        std::declval<Seq>().generate(
            std::declval<uint32_t*>(),
            std::declval<uint32_t*>()),
        (void) 0)
    >> : std::true_type {};

// satisfies RandomNumberEngine
class xoroshiro128plus {
public:
    typedef uint64_t result_type;
    static constexpr uint64_t min() {
        return 0;
    }
    static constexpr uint64_t max() {
        return std::numeric_limits<uint64_t>::max();
    }
    xoroshiro128plus() {
        seed();
    }
    xoroshiro128plus(uint64_t s0) : s{s0} {
    }
    template <typename SeedSequence>
    xoroshiro128plus(SeedSequence seq) {
        seed(seq);
    }
    void seed() {
        // some randomly chosen constants
        s[0] = 0xe6a52bbd554f6446;
        s[1] = 0x4917ac5d1056b6d7;
    }
    void seed(uint64_t seed) {
        // s[0] = s0;
        // s[1] = 0;
        s[0] = 0x01d353e5f3993bb0+seed;
        s[1] = 0x7b9c0df6cb193b20*(seed+1);
    }
    template <typename SeedSequence, typename = typename std::enable_if<is_seed_sequence<SeedSequence>::value>::type>
    void seed(SeedSequence seq) {
        auto state = reinterpret_cast<uint32_t *>(s);
        seq.generate(state, state + 4);
    }
    uint64_t operator()() {
        const uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        const uint64_t result = s0 + s1;

        s1 ^= s0;
        s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
        s[1] = rotl(s1, 36); // c

        return result;
    }
    void discard(unsigned long long z) {
#if ULONGLONG_MAX > UINT64_MAX
        static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };
        for (auto i = 0ull, limit = z >> 64; i < limit; ++i) {
            uint64_t s0 = 0;
            uint64_t s1 = 0;
            for (std::size_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
                for(int b = 0; b < 64; b++) {
                    if (JUMP[i] & UINT64_C(1) << b) {
                        s0 ^= s[0];
                        s1 ^= s[1];
                    }
                    operator()();
                }
            s[0] = s0;
            s[1] = s1;
        }
        z &= 0xFFFFFFFFFFFFFFFFull;
#endif
        while (z--)
            operator()();
    }

    friend bool operator==(const xoroshiro128plus& lhs, const xoroshiro128plus& rhs) {
        return lhs.s[0] == rhs.s[0] && lhs.s[1] == rhs.s[1];
    }
    friend bool operator!=(const xoroshiro128plus& lhs, const xoroshiro128plus& rhs) {
        return lhs.s[0] != rhs.s[0] || lhs.s[1] != rhs.s[1];
    }
    friend std::ostream& operator<<(std::ostream& out, const xoroshiro128plus& x) {
        return out << x.s[0] << ' ' << x.s[1];
    }
    friend std::istream& operator>>(std::istream& in, xoroshiro128plus& x) {
        return in >> x.s[0] >> x.s[1];
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    uint64_t s[2];
};

// vim:set sw=4 et:
#endif
